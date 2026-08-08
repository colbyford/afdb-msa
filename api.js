/*
 * api.js -- the two remote services this site talks to. Both send
 * `access-control-allow-origin: *`, which is the only reason a page with no
 * backend can do any of this.
 *
 *   EBI Job Dispatcher (ncbiblast)  -- sequence search, returns UniProt hits
 *                                      *and* the query/hit pairwise alignment
 *   AlphaFold DB                    -- the precomputed a3m for a UniProt entry
 */

const EBI = 'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast';
const AFDB_MSA = acc => `https://alphafold.ebi.ac.uk/files/msa/AF-${String(acc).trim().toUpperCase()}-F1-msa_v6.a3m`;

/* ---------- EBI BLAST ---------- */

/*
 * EBI requires a contact email on every job and asks that clients poll no more
 * than once every few seconds and run no more than ~30 jobs at a time.
 * See https://www.ebi.ac.uk/jdispatcher/docs/
 */
async function blastSubmit(seq, opts) {
  const body = new URLSearchParams({
    email: opts.email,
    program: 'blastp',
    stype: 'protein',
    matrix: opts.matrix || 'BLOSUM62',
    database: opts.database || 'uniprotkb',
    exp: opts.evalue || '1e-3',
    alignments: String(opts.nHits || 50),
    scores: String(opts.nHits || 50),
    filter: 'F',
    sequence: `>query\n${seq}`,
  });
  const r = await fetch(`${EBI}/run`, { method: 'POST', body });
  const text = (await r.text()).trim();
  if (!r.ok) throw new Error(`BLAST submit failed (${r.status}): ${text.slice(0, 300)}`);
  return text; // job id
}

async function blastStatus(jobId) {
  const r = await fetch(`${EBI}/status/${jobId}`);
  return (await r.text()).trim();
}

async function blastResult(jobId) {
  const r = await fetch(`${EBI}/result/${jobId}/json`);
  if (!r.ok) throw new Error(`Could not fetch BLAST result (${r.status})`);
  return r.json();
}

/*
 * A SwissProt job finishes in roughly 3-4 seconds, so a flat 3 s poll spent up
 * to 3 s of that just not noticing. Measured on the same query: 3000 ms polling
 * took 70 s wall clock, 500 ms polling took 6.1 s.
 *
 * Start tight and back off, so short jobs are caught promptly while a slow
 * UniProtKB job (2-5 minutes) does not hammer EBI for its whole duration.
 */
async function blastWait(jobId, onTick, signal, intervalMs) {
  const t0 = Date.now();
  let every = intervalMs || 400;
  const max = intervalMs || 3000;
  for (;;) {
    if (signal && signal.aborted) throw new Error('cancelled');
    const s = await blastStatus(jobId);
    const elapsed = Math.round((Date.now() - t0) / 1000);
    if (onTick) onTick(s, elapsed);
    if (s === 'FINISHED') return;
    if (s === 'FAILURE' || s === 'ERROR' || s === 'NOT_FOUND') {
      throw new Error(`BLAST job ${s.toLowerCase()}`);
    }
    await new Promise(res => setTimeout(res, every));
    every = Math.min(Math.round(every * 1.4), max);
  }
}

/*
 * Flatten the EBI JSON into what the a3m transfer needs: a UniProt accession,
 * the hit length, and the aligned query/hit strings from the top HSP.
 */
function parseBlastHits(json) {
  const hits = [];
  for (const h of json.hits || []) {
    const hsp = (h.hit_hsps || [])[0];
    if (!hsp) continue;
    hits.push({
      acc: h.hit_acc,
      db: h.hit_db,                 // 'SP' (SwissProt) or 'TR' (TrEMBL)
      desc: h.hit_uni_de || h.hit_desc || '',
      organism: h.hit_uni_os || h.hit_os || '',
      taxid: h.hit_uni_ox || '',
      gene: h.hit_uni_gn || '',
      hitLen: h.hit_len,
      identity: hsp.hsp_identity,   // percent
      evalue: hsp.hsp_expect,
      bits: hsp.hsp_bit_score,
      alignLen: hsp.hsp_align_len,
      // fraction of the QUERY covered by this HSP
      queryCoverage: (hsp.hsp_query_to - hsp.hsp_query_from + 1) / (json.query_len || 1),
      hsp: {
        qseq: hsp.hsp_qseq,
        hseq: hsp.hsp_hseq,
        qFrom: hsp.hsp_query_from,
        hFrom: hsp.hsp_hit_from,
      },
    });
  }
  return hits;
}

/* ---------- exact-sequence shortcut via UniParc ---------- */

/*
 * UniParc indexes every sequence ever seen by its CRC64 checksum, so an exact
 * match is a hash lookup rather than a search: ~0.4 s against BLAST's 8 s
 * (SwissProt) or 250 s (full UniProtKB).
 *
 * This matters because pasting a real UniProt sequence is the common case, and
 * when the query IS the entry, the query/hit alignment is the identity map --
 * there is nothing for BLAST to compute.
 *
 * CRC64-ISO 3309, the variant UniProt uses.
 */
const CRC64_TABLE = (() => {
  const POLY_HI = 0xd8000000;
  const tbl = new Array(256);
  for (let i = 0; i < 256; i++) {
    let lo = i, hi = 0;
    for (let b = 0; b < 8; b++) {
      const carry = lo & 1;
      lo = (lo >>> 1) | ((hi & 1) << 31);
      hi = hi >>> 1;
      if (carry) hi = (hi ^ POLY_HI) >>> 0;
    }
    tbl[i] = [lo >>> 0, hi >>> 0];
  }
  return tbl;
})();

function crc64(str) {
  let lo = 0, hi = 0;
  for (let i = 0; i < str.length; i++) {
    const idx = (lo ^ str.charCodeAt(i)) & 0xff;
    const [tlo, thi] = CRC64_TABLE[idx];
    const nlo = ((lo >>> 8) | ((hi & 0xff) << 24)) >>> 0;
    const nhi = hi >>> 8;
    lo = (nlo ^ tlo) >>> 0;
    hi = (nhi ^ thi) >>> 0;
  }
  const h = n => n.toString(16).toUpperCase().padStart(8, '0');
  return h(hi) + h(lo);
}

/*
 * Every UniProt accession whose sequence is exactly `seq`, or [] if there is no
 * exact match. Ordered SwissProt-first, since curated entries are likelier to
 * have a well-populated AFDB MSA.
 */
async function uniparcAccessions(seq, signal) {
  const url = `https://rest.uniprot.org/uniparc/search?query=checksum:${crc64(seq.toUpperCase())}`
    + '&fields=upi,accession&format=json&size=1';
  const r = await fetch(url, { signal });
  if (!r.ok) return [];
  const j = await r.json();
  const hit = (j.results || [])[0];
  if (!hit) return [];
  // Entries look like "P69905" or "Q3MIF5.1" (versioned, and often obsolete).
  const accs = (hit.uniProtKBAccessions || [])
    .filter(a => !a.includes('.'))
    .map(a => a.trim())
    .filter(Boolean);
  // SwissProt accessions are 6 or 10 chars and never start with "A0A"; this is a
  // heuristic ordering only, so being wrong just costs a 404 and a retry.
  return accs.sort((a, b) => (a.startsWith('A0A') ? 1 : 0) - (b.startsWith('A0A') ? 1 : 0));
}

// Name/organism/taxid for one accession, so the exact-match path can label its
// row as informatively as the BLAST path does. Cheap (~0.2 s) and non-fatal.
async function uniprotSummary(acc, signal) {
  try {
    const r = await fetch(
      `https://rest.uniprot.org/uniprotkb/${acc}.json?fields=protein_name,organism_name,organism_id`,
      { signal });
    if (!r.ok) return {};
    const j = await r.json();
    return {
      desc: j.proteinDescription?.recommendedName?.fullName?.value
        || j.proteinDescription?.submissionNames?.[0]?.fullName?.value || '',
      organism: j.organism?.scientificName || '',
      taxid: j.organism?.taxonId ? String(j.organism.taxonId) : '',
    };
  } catch { return {}; }
}

/*
 * Sequences for a batch of accessions. The index returns accessions but the a3m
 * transfer needs an alignment, and an alignment needs the candidate's residues.
 * UniProt takes up to 100 accessions per request and answers in well under a
 * second, so a whole candidate list costs one round trip rather than one each.
 */
async function uniprotSequences(accs, signal) {
  const out = new Map();
  for (let i = 0; i < accs.length; i += 100) {
    const batch = accs.slice(i, i + 100);
    const url = `https://rest.uniprot.org/uniprotkb/accessions?accessions=${batch.join(',')}&format=fasta`;
    const r = await fetch(url, { signal });
    if (!r.ok) continue;                       // partial results beat none
    const text = await r.text();
    let acc = null, buf = [];
    const flush = () => { if (acc && buf.length) out.set(acc, buf.join('')); buf = []; };
    for (const line of text.split('\n')) {
      if (!line) continue;
      if (line[0] === '>') {
        flush();
        const m = /^>[a-z]{2}\|([^|]+)\|/i.exec(line);
        acc = m ? m[1] : null;
      } else if (acc) buf.push(line.trim());
    }
    flush();
  }
  return out;
}

/*
 * Sequence for one AFDB accession, straight from AlphaFold DB.
 *
 * UniProt is the wrong authority here. AFDB v6 was built on UniProt 2025_03 and
 * holds 91M more entries than UniProtKB does now; TrEMBL was pruned hard in
 * between, so a large share of AFDB accessions have since been deleted from
 * UniProtKB and return no sequence -- measured at 3 of 6 on a real candidate
 * list. They are all still in AFDB, which is the database we are borrowing from.
 * One request per accession, but only for the ones UniProt could not supply.
 */
// AFDB's /api/ endpoint answers 403 to requests without a browser User-Agent,
// unlike /files/, which is open. Browsers always send one and silently drop this
// header (it is forbidden to set from fetch), so this line is a no-op in the
// page and the thing that makes the node tooling work.
const BROWSER_UA = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36';

async function afdbSequence(acc, signal) {
  try {
    const r = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${acc}`,
      { signal, headers: { 'User-Agent': BROWSER_UA } });
    if (!r.ok) return null;
    const j = await r.json();
    const e = Array.isArray(j) ? j[0] : j;
    return (e && e.uniprotSequence) || null;
  } catch { return null; }
}

/*
 * Sequences for a candidate list: one batched UniProt request, then AFDB for
 * whatever that missed. Cheap in the common case, correct in all of them.
 */
async function candidateSequences(accs, signal, onNote) {
  const out = await uniprotSequences(accs, signal);
  const missing = accs.filter(a => !out.has(a));
  if (!missing.length) return out;

  const found = await Promise.all(missing.map(a => afdbSequence(a, signal)));
  let recovered = 0;
  missing.forEach((a, i) => { if (found[i]) { out.set(a, found[i]); recovered++; } });
  if (onNote && recovered) {
    onNote(`${recovered} of ${missing.length} candidate(s) missing from UniProtKB, recovered from AlphaFold DB`);
  }
  return out;
}

/* ---------- AlphaFold DB ---------- */

/*
 * These files run from tens of KB to tens of MB (EGFR is ~22 MB), so stream
 * with progress rather than blocking on one opaque await.
 *
 * There is deliberately no separate "does this exist" probe. AFDB ignores Range
 * requests and answers HEAD without CORS headers, so any probe would download
 * the whole file only for us to download it again -- doubling the bytes for
 * every hit. Callers fetch and handle NoMsaError instead.
 */
class NoMsaError extends Error {
  constructor(acc, status) {
    super(`No AlphaFold DB MSA for ${acc} (HTTP ${status})`);
    this.name = 'NoMsaError';
    this.acc = acc;
  }
}

async function afdbFetchMsa(acc, onProgress, signal) {
  const r = await fetch(AFDB_MSA(acc), { signal });
  if (r.status === 404) throw new NoMsaError(acc, 404);
  if (!r.ok) throw new NoMsaError(acc, r.status);
  const total = Number(r.headers.get('content-length')) || 0;

  if (!r.body || !r.body.getReader) return r.text();

  const reader = r.body.getReader();
  const dec = new TextDecoder();
  const parts = [];
  let got = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    got += value.length;
    parts.push(dec.decode(value, { stream: true }));
    if (onProgress) onProgress(got, total);
  }
  parts.push(dec.decode());
  return parts.join('');
}

/*
 * Try accessions in order, skipping those with no MSA. Used both by the exact
 * -match path (several accessions share one sequence) and by the BLAST path
 * (descend the diverse ranking when an entry has no MSA).
 */
async function afdbFetchFirstAvailable(accs, onProgress, signal, onSkip) {
  let last = null;
  for (const acc of accs) {
    try {
      const text = await afdbFetchMsa(acc, onProgress, signal);
      return { acc, text };
    } catch (e) {
      if (e.name !== 'NoMsaError') throw e;
      last = e;
      if (onSkip) onSkip(acc);
    }
  }
  throw last || new Error('No accessions to try');
}

/*
 * MSAs are immutable per release and often tens of MB, so a re-run with the
 * same query should not re-download them. Cache API persists across reloads;
 * the Map is the fallback when it is unavailable (file://, private mode).
 */
const memCache = new Map();
const CACHE_NAME = 'afdb-msa-v1';

async function afdbFetchMsaCached(acc, onProgress, signal) {
  if (memCache.has(acc)) return memCache.get(acc);
  const url = AFDB_MSA(acc);

  let cache = null;
  try { cache = typeof caches !== 'undefined' ? await caches.open(CACHE_NAME) : null; } catch { /* unavailable */ }
  if (cache) {
    const hit = await cache.match(url);
    if (hit) {
      const text = await hit.text();
      memCache.set(acc, text);
      return text;
    }
  }

  const text = await afdbFetchMsa(acc, onProgress, signal);
  memCache.set(acc, text);
  if (cache) {
    try { await cache.put(url, new Response(text)); } catch { /* quota; not fatal */ }
  }
  return text;
}

const Api = {
  EBI, AFDB_MSA,
  blastSubmit, blastStatus, blastResult, blastWait, parseBlastHits,
  crc64, uniparcAccessions, uniprotSummary, uniprotSequences, afdbSequence, candidateSequences,
  afdbFetchMsa, afdbFetchMsaCached, afdbFetchFirstAvailable, NoMsaError,
};

if (typeof module !== 'undefined' && module.exports) module.exports = Api;
if (typeof window !== 'undefined') window.Api = Api;
