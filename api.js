/*
 * api.js -- everything this page fetches from AlphaFold DB.
 *
 * There is one source of truth. An earlier version also asked UniProt for
 * candidate sequences and fell back to EBI BLAST when the index came up short;
 * both are gone. UniProt is the wrong authority here anyway -- AFDB v6 was built
 * on UniProt 2025_03 and holds 91M more entries than UniProtKB does now, so
 * measured on a real candidate list, 63% of AFDB accessions no longer resolve
 * there at all.
 */

const AFDB_MSA = acc =>
  `https://alphafold.ebi.ac.uk/files/msa/AF-${String(acc).trim().toUpperCase()}-F1-msa_v6.a3m`;

// AFDB's /api/ endpoint answers 403 without a browser User-Agent, unlike
// /files/. Browsers always send one and silently drop this header (fetch
// forbids setting it), so this line is a no-op in the page and the thing that
// makes the node tests work.
const BROWSER_UA = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36';

/* ---------- sequences ---------- */

/*
 * The residues for one AFDB accession. Needed because the index returns
 * accessions, and an alignment needs the actual sequence.
 */
async function sequenceOf(acc, signal) {
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
 * Sequences for a candidate list, fetched together. One request each, but they
 * are small and concurrent, so the whole list costs about one round trip.
 * Also returns the description and organism, which are worth showing and come
 * back in the same response.
 */
async function candidateInfo(accs, signal) {
  const out = new Map();
  await Promise.all(accs.map(async acc => {
    try {
      const r = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${acc}`,
        { signal, headers: { 'User-Agent': BROWSER_UA } });
      if (!r.ok) return;
      const j = await r.json();
      const e = Array.isArray(j) ? j[0] : j;
      if (e && e.uniprotSequence) {
        out.set(acc, {
          seq: e.uniprotSequence,
          desc: e.uniprotDescription || '',
          organism: e.organismScientificName || '',
          taxid: e.taxId ? String(e.taxId) : '',
        });
      }
    } catch { /* one missing candidate is not fatal */ }
  }));
  return out;
}

/* ---------- the MSA ---------- */

class NoMsaError extends Error {
  constructor(acc, status) {
    super(`No AlphaFold DB MSA for ${acc} (HTTP ${status})`);
    this.name = 'NoMsaError';
    this.acc = acc;
  }
}

/*
 * These run from tens of KB to tens of MB (EGFR is ~22 MB), so stream with
 * progress rather than blocking on one opaque await.
 *
 * There is deliberately no "does this exist" probe first: AFDB ignores Range
 * requests and answers HEAD without CORS headers, so probing would download the
 * whole file only to download it again.
 */
async function fetchMsa(acc, onProgress, signal) {
  const r = await fetch(AFDB_MSA(acc), { signal });
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

// MSAs are immutable per release and often tens of MB, so a re-run with the
// same query should not re-download them.
const memCache = new Map();
const CACHE_NAME = 'afdb-msa-v1';

async function fetchMsaCached(acc, onProgress, signal) {
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

  const text = await fetchMsa(acc, onProgress, signal);
  memCache.set(acc, text);
  if (cache) { try { await cache.put(url, new Response(text)); } catch { /* quota */ } }
  return text;
}

const Api = { AFDB_MSA, sequenceOf, candidateInfo, fetchMsa, fetchMsaCached, NoMsaError };

if (typeof module !== 'undefined' && module.exports) module.exports = Api;
if (typeof window !== 'undefined') window.Api = Api;
