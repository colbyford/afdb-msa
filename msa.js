/*
 * msa.js -- pure a3m logic. No DOM, no network, so it can be unit-tested in node.
 *
 * The trick this whole site rests on: AlphaFold DB publishes a precomputed a3m
 * for every UniProt entry. If a user's query is close to some AFDB entry, we can
 * borrow that entry's MSA instead of running a fresh homology search -- provided
 * we re-index it from the hit's residue frame into the query's. BLAST already
 * gives us the pairwise alignment needed to do that.
 */

/* ---------- a3m parsing ---------- */

// In a3m, lowercase letters are insertions relative to the query and do not
// occupy a column. Everything else (uppercase or '-') is a match state.
function countMatchStates(seq) {
  let n = 0;
  for (let i = 0; i < seq.length; i++) {
    const c = seq.charCodeAt(i);
    if (c < 97 || c > 122) n++;
  }
  return n;
}

function parseA3M(text) {
  const rows = [];
  let name = null;
  let buf = [];
  const flush = () => {
    if (name !== null) rows.push({ name, seq: buf.join('') });
    buf = [];
  };
  for (const line of text.split('\n')) {
    if (!line) continue;
    if (line[0] === '>') {
      flush();
      name = line.slice(1).trim();
    } else if (name !== null) {
      buf.push(line.trim());
    }
  }
  flush();
  return rows;
}

function formatA3M(rows, wrap) {
  const out = [];
  for (const r of rows) {
    out.push('>' + r.name);
    if (wrap) {
      for (let i = 0; i < r.seq.length; i += wrap) out.push(r.seq.slice(i, i + wrap));
    } else {
      out.push(r.seq);
    }
  }
  return out.join('\n') + '\n';
}

/* ---------- header metadata ---------- */

function taxidOf(header) {
  const m = /TaxID[=:](\d+)/i.exec(header);
  return m ? m[1] : null;
}

function organismOf(header) {
  const m = /Tax=([^=]*?)\s+TaxID=/i.exec(header);
  return m ? m[1].trim() : null;
}

// Fraction of the query's aligned residues that this row matches.
function identityToQuery(querySeq, rowSeq) {
  let same = 0, n = 0, qi = 0;
  for (let i = 0; i < rowSeq.length; i++) {
    const c = rowSeq[i];
    if (c >= 'a' && c <= 'z') continue; // insertion: no query column
    const q = querySeq[qi++];
    if (c === '-' || q === undefined) continue;
    n++;
    if (c === q) same++;
  }
  return n ? same / n : 0;
}

/* ---------- the core: re-index an a3m from hit frame to query frame ---------- */

/*
 * Build, from a BLAST HSP, the map hitResidueIndex -> queryResidueIndex (0-based),
 * where `null` means the hit residue has no query counterpart.
 *
 * hsp = { qseq, hseq, qFrom, hFrom }  (qFrom/hFrom are 1-based, as BLAST reports)
 */
function hspToMap(hsp, queryLen, hitLen) {
  const hitToQuery = new Array(hitLen).fill(null);
  const queryCovered = new Array(queryLen).fill(false);
  let hi = hsp.hFrom - 1;
  let qi = hsp.qFrom - 1;
  const { qseq, hseq } = hsp;
  for (let i = 0; i < qseq.length; i++) {
    const qc = qseq[i], hc = hseq[i];
    const qGap = qc === '-';
    const hGap = hc === '-';
    if (!qGap && !hGap) {
      if (hi < hitLen && qi < queryLen) {
        hitToQuery[hi] = qi;
        queryCovered[qi] = true;
      }
      hi++; qi++;
    } else if (qGap && !hGap) {
      hi++;          // hit residue falls in a query gap -> becomes an insertion
    } else if (!qGap && hGap) {
      qi++;          // query residue with no hit column -> all-gap column
    }
  }
  return { hitToQuery, queryCovered };
}

/*
 * Rewrite every row of a hit's a3m so its match states are the QUERY's residues.
 *
 * srcRows   : parseA3M output for the AFDB a3m (row 0 is the hit's own sequence)
 * queryLen  : length of the user's query
 * hitToQuery: from hspToMap
 *
 * Hit columns that map to no query position are demoted to lowercase insertions;
 * query positions the hit does not cover become '-' in every borrowed row.
 */
function transferRows(srcRows, queryLen, hitToQuery) {
  const out = [];
  for (const row of srcRows) {
    const cols = new Array(queryLen).fill('-');
    // inserts[i] is the lowercase run emitted immediately BEFORE query column i;
    // index queryLen holds any trailing run.
    const inserts = new Array(queryLen + 1).fill('');
    let pending = '';
    let matchIdx = 0;
    let any = false;
    const s = row.seq;

    for (let i = 0; i < s.length; i++) {
      const c = s[i];
      if (c >= 'a' && c <= 'z') { pending += c; continue; }

      const qpos = matchIdx < hitToQuery.length ? hitToQuery[matchIdx] : null;
      matchIdx++;
      if (qpos === null) {
        // dropped column: keep the residue as an insertion so nothing is lost
        if (c !== '-' && c !== '.') pending += c.toLowerCase();
        continue;
      }
      if (pending) { inserts[qpos] += pending; pending = ''; }
      cols[qpos] = c;
      if (c !== '-') any = true;
    }
    if (pending) inserts[queryLen] += pending;
    if (!any) continue; // row contributed nothing inside the query's span

    let seq = '';
    for (let i = 0; i < queryLen; i++) seq += inserts[i] + cols[i];
    seq += inserts[queryLen];
    out.push({ name: row.name, seq });
  }
  return out;
}

// The UniRef member id, e.g. ">UniRef90_A0A4U1FNC5/2-139 [subseq from] ..." -> "UniRef90_A0A4U1FNC5/2-139".
// Different AFDB entries search the same UniRef90 database, so the same member
// recurs verbatim across hits; this is the cheap half of merge deduplication.
function memberIdOf(header) {
  const m = /^(\S+)/.exec(header);
  return m ? m[1] : header;
}

/*
 * Full per-chain build: take one query + its BLAST hits (each with a fetched a3m)
 * and produce a single deduplicated, query-framed a3m.
 *
 * Merging several hits' MSAs is heavily redundant -- neighbouring UniProt
 * entries pull from the same UniRef90 clusters, so the second hit typically
 * contributes far less than its raw depth suggests. Two rounds handle it:
 * exact deduplication here (by UniRef member id and by transferred sequence),
 * and identity-threshold clustering afterwards in filter.js. `stats` reports
 * the per-hit yield so the UI can show what each hit was actually worth.
 *
 * chain = { name, seq }
 * hits  = [{ acc, hitLen, hsp:{qseq,hseq,qFrom,hFrom}, a3mText, identity, ... }]
 */
function buildChainA3M(chain, hits, opts) {
  opts = opts || {};
  const maxDepth = opts.maxDepth || 0;
  const queryLen = chain.seq.length;
  const rows = [{ name: chain.name, seq: chain.seq.toUpperCase() }];
  const seenId = new Set();
  // Seeded with the query: users routinely paste a sequence that IS a database
  // entry, and without this the top hit's own row would repeat it verbatim.
  const seenSeq = new Set([rows[0].seq]);
  const stats = [];
  let capped = false;

  for (const hit of hits) {
    const src = parseA3M(hit.a3mText);
    if (!src.length) { stats.push({ acc: hit.acc, added: 0, total: 0, dupId: 0, dupSeq: 0 }); continue; }
    // AFDB names row 0 ">query" -- it is the hit's own sequence. Label it usefully.
    const hitLen = hit.hitLen || countMatchStates(src[0].seq);
    src[0] = {
      name: `${hit.acc} ${hit.desc || ''} Tax=${hit.organism || '?'} TaxID=${hit.taxid || '0'}`.trim(),
      seq: src[0].seq,
    };
    const { hitToQuery } = hspToMap(hit.hsp, queryLen, hitLen);
    const moved = transferRows(src, queryLen, hitToQuery);

    let added = 0, dupId = 0, dupSeq = 0;
    for (const r of moved) {
      const id = memberIdOf(r.name);
      if (seenId.has(id)) { dupId++; continue; }
      // Transfer can map two different members onto the same query-frame string.
      const key = r.seq;
      if (seenSeq.has(key)) { dupSeq++; continue; }
      seenId.add(id);
      seenSeq.add(key);
      rows.push(r);
      added++;
      if (maxDepth && rows.length >= maxDepth) { capped = true; break; }
    }
    stats.push({ acc: hit.acc, added, total: moved.length, dupId, dupSeq });
    if (capped) break;
  }
  return { rows, stats, queryLen, capped };
}

/* ---------- choosing which hits to actually download ---------- */

/*
 * Project a hit's aligned sequence into query columns, so two hits of the same
 * query become directly comparable without any further alignment.
 */
function hitInQueryFrame(hit, queryLen) {
  const cols = new Array(queryLen).fill('-');
  let qi = hit.hsp.qFrom - 1;
  const { qseq, hseq } = hit.hsp;
  for (let i = 0; i < qseq.length; i++) {
    const qc = qseq[i], hc = hseq[i];
    if (qc === '-') continue;            // insertion relative to the query
    if (qi < queryLen && hc !== '-') cols[qi] = hc;
    qi++;
  }
  return cols.join('');
}

/*
 * AlphaFold DB builds one MSA per UniRef cluster, so UniProt entries that are
 * near-identical resolve to the SAME a3m -- BLAST's top hits for a given query
 * are exactly such a set. Empirically the top five hits for haemoglobin alpha
 * return five byte-identical 2783-sequence files, and merging them adds nothing.
 *
 * So pick greedily for diversity: accept a hit only if it is less than
 * `maxPairId` identical to every hit already accepted. The comparison uses the
 * hit sequences BLAST already returned, so redundant downloads never happen.
 *
 * The whole diverse ranking is returned, not just the first `maxHits` of it.
 * Callers must still drop hits that turn out to have no AFDB MSA, and when they
 * do they need the NEXT diverse hit to fall back on -- falling back to the raw
 * BLAST order would hand them exactly the duplicates this just rejected.
 */
function selectDiverseHits(hits, queryLen, opts) {
  opts = opts || {};
  const maxPairId = opts.maxPairId === undefined ? 0.9 : opts.maxPairId;
  const limit = opts.maxHits || 0; // 0 = rank everything

  const chosen = [];
  const frames = [];
  const skipped = [];

  for (const h of hits) {
    if (limit && chosen.length >= limit) break;
    const f = hitInQueryFrame(h, queryLen);
    let dupOf = null;
    for (let i = 0; i < frames.length; i++) {
      if (identityToQuery(frames[i], f) >= maxPairId) { dupOf = chosen[i].acc; break; }
    }
    if (dupOf) { skipped.push({ acc: h.acc, dupOf }); continue; }
    chosen.push(h);
    frames.push(f);
  }
  return { chosen, skipped };
}

/* ---------- pairing across chains (ColabFold-style) ---------- */

/*
 * chains = [{ queryLen, rows }] as returned by buildChainA3M.
 *
 * Rows are grouped by TaxID; the best-scoring row per taxid per chain is kept.
 * Taxa present in EVERY chain form the paired block (concatenated left-to-right);
 * everything else is stacked block-diagonally below, padded with gaps.
 */
function pairChains(chains, opts) {
  opts = opts || {};
  const maxPaired = opts.maxPaired || 0;
  const lens = chains.map(c => c.queryLen);
  const querySeq = chains.map(c => c.rows[0].seq).join('');

  const perChain = chains.map(c => {
    const q = c.rows[0].seq;
    const best = new Map();
    for (let i = 1; i < c.rows.length; i++) {
      const r = c.rows[i];
      const tax = taxidOf(r.name);
      if (!tax) continue;
      const id = identityToQuery(q, r.seq);
      const cur = best.get(tax);
      if (!cur || id > cur.id) best.set(tax, { row: r, id });
    }
    return best;
  });

  // taxa shared by all chains
  let shared = [...perChain[0].keys()];
  for (let i = 1; i < perChain.length; i++) {
    const s = perChain[i];
    shared = shared.filter(t => s.has(t));
  }
  shared.sort((a, b) => {
    const sa = Math.min(...perChain.map(p => p.get(a).id));
    const sb = Math.min(...perChain.map(p => p.get(b).id));
    return sb - sa;
  });
  if (maxPaired) shared = shared.slice(0, maxPaired);

  const out = [{ name: '101', seq: querySeq }];
  const usedPerChain = chains.map(() => new Set());

  for (const tax of shared) {
    const parts = [];
    const names = [];
    perChain.forEach((p, ci) => {
      const { row } = p.get(tax);
      parts.push(row.seq);
      names.push(row.name.split(/\s+/)[0]);
      usedPerChain[ci].add(row.name);
    });
    out.push({ name: names.join('\t'), seq: parts.join('') });
  }
  const nPaired = out.length - 1;

  // unpaired block-diagonal remainder
  chains.forEach((c, ci) => {
    const pad = chains.map((o, oi) => (oi === ci ? null : '-'.repeat(lens[oi])));
    for (let i = 1; i < c.rows.length; i++) {
      const r = c.rows[i];
      if (usedPerChain[ci].has(r.name)) continue;
      const parts = pad.map((p, oi) => (oi === ci ? r.seq : p));
      out.push({ name: r.name, seq: parts.join('') });
    }
  });

  const header = `#${lens.join(',')}\t${lens.map(() => 1).join(',')}`;
  return { rows: out, header, nPaired, lens };
}

function formatPairedA3M(paired, wrap) {
  return paired.header + '\n' + formatA3M(paired.rows, wrap);
}

/* ---------- FASTA input parsing ---------- */

const AA = /^[ACDEFGHIKLMNPQRSTVWYXBZJUO]+$/i;

function parseQueryInput(text) {
  const chains = [];
  const t = text.trim();
  if (!t) return chains;

  if (t[0] === '>') {
    for (const rec of t.split('\n>')) {
      const lines = rec.replace(/^>/, '').split('\n');
      const name = lines[0].trim() || `chain_${chains.length + 1}`;
      const seq = lines.slice(1).join('').replace(/[\s*]/g, '').toUpperCase();
      if (seq) chains.push({ name, seq });
    }
  } else {
    // bare sequence(s), ':' separates chains (ColabFold convention)
    const parts = t.replace(/[\s*]/g, '').toUpperCase().split(':');
    parts.forEach((p, i) => { if (p) chains.push({ name: `chain_${i + 1}`, seq: p }); });
  }
  for (const c of chains) {
    if (!AA.test(c.seq)) throw new Error(`"${c.name}" contains non-amino-acid characters`);
  }
  return chains;
}

const MSAKit = {
  parseA3M, formatA3M, countMatchStates, taxidOf, organismOf, identityToQuery,
  hspToMap, transferRows, buildChainA3M, pairChains, formatPairedA3M, parseQueryInput,
  memberIdOf, hitInQueryFrame, selectDiverseHits,
};

if (typeof module !== 'undefined' && module.exports) module.exports = MSAKit;
if (typeof window !== 'undefined') window.MSAKit = MSAKit;
