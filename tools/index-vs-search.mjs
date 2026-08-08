/*
 * node tools/index-vs-search.mjs <indexDir> <afdbPart.fa>
 *
 * The index has been shown to return candidates fast. This asks the only
 * question that matters: are they the same answers a real search would give?
 *
 * Ground truth is BLAST against full UniProtKB (~150M sequences, ~250 s each) --
 * the closest available reference at AFDB scale. SwissProt would be unfair to
 * both sides: it covers 570k, and a >=90% relative exists there for only 6% of
 * random AFDB sequences.
 *
 * Queries are real AFDB sequences with a few percent of positions mutated, which
 * is the realistic case: a user's sequence is not in the database verbatim but
 * has a close relative in it. Comparison is on the identity of the best hit each
 * method finds, aligned the same way (Smith-Waterman) so neither side is
 * flattered by its own scoring.
 */
import { createRequire } from 'node:module';
import { readFileSync, openSync, readSync, closeSync, statSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const require = createRequire(import.meta.url);
const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const A = require(join(root, 'align.js'));
const Api = require(join(root, 'api.js'));

const [, , indexDir, partPath] = process.argv;
if (!indexDir || !partPath) {
  console.error('usage: node tools/index-vs-search.mjs <indexDir> <afdbPart.fa>');
  process.exit(2);
}
const N_QUERIES = Number(process.env.N_QUERIES || 16);
const MUTATE = Number(process.env.MUTATE || 0.05);
const EMAIL = process.env.EMAIL || 'anonymous@example.org';
const DB = process.env.DB || 'uniprotkb';
const CONCURRENCY = Number(process.env.CONCURRENCY || 8);

/* ---------- the index ---------- */
const meta = JSON.parse(readFileSync(join(indexDir, 'meta.json'), 'utf8'));
const bk = new Uint32Array(readFileSync(join(indexDir, 'buckets.u32')).buffer);
const files = meta.files.map(f => ({ ...f, fd: openSync(join(indexDir, f.name), 'r') }));
const accFd = openSync(join(indexDir, 'acc.bin'), 'r');
const accBuf = Buffer.allocUnsafe(12);
const accOf = id => { readSync(accFd, accBuf, 0, 12, id * 12); return accBuf.toString('latin1').replace(/\0+$/, ''); };

const AA = 'ACDEFGHIKLMNPQRSTVWY';
const IDX = new Int8Array(128).fill(-1);
for (let i = 0; i < AA.length; i++) IDX[AA.charCodeAt(i)] = i;
const h32 = x => { x = (x ^ 0x9e3779b9) >>> 0; x = Math.imul(x ^ (x >>> 16), 0x85ebca6b) >>> 0; x = Math.imul(x ^ (x >>> 13), 0xc2b2ae35) >>> 0; return (x ^ (x >>> 16)) >>> 0; };

function seedsOf(seq) {
  const n = seq.length, c = new Int8Array(n);
  for (let i = 0; i < n; i++) c[i] = IDX[seq.charCodeAt(i) & 0x7f];
  const hs = new Int32Array(n).fill(-1);
  for (let i = meta.k - 1; i < n; i++) {
    let x = 0, ok = true;
    for (let j = i - meta.k + 1; j <= i; j++) { const q = c[j]; if (q < 0) { ok = false; break; } x = (Math.imul(x, 31) + q) | 0; }
    if (ok) hs[i] = h32(x >>> 0) & 0x7fffffff;
  }
  const w = Math.max(4, Math.min(128, Math.floor(n / meta.targetSeeds) || 4));
  const s = new Set();
  for (let i = meta.k - 1; i < n; i++) {
    if (i % w) continue;
    let b = -1; const stop = Math.min(i + w, n);
    for (let j = i; j < stop; j++) if (hs[j] >= 0 && (b < 0 || hs[j] < b)) b = hs[j];
    if (b >= 0) s.add(b);
  }
  return s;
}

const postBuf = Buffer.allocUnsafe(1 << 23);
function indexSearch(seq, topK) {
  const t0 = Date.now();
  const cnt = new Map();
  let bytes = 0;
  for (const seed of seedsOf(seq)) {
    const b = seed >>> meta.shift, lo = bk[b], hi = bk[b + 1];
    if (hi <= lo) continue;
    const len = (hi - lo) * meta.postBytes;
    if (len > postBuf.length) continue;                 // pathological seed
    const f = files.find(x => b >= x.fromBucket && b < x.toBucket);
    readSync(f.fd, postBuf, 0, len, (lo - f.postFrom) * meta.postBytes);
    bytes += len;
    const res = seed & 0xff;
    for (let i = 0; i < hi - lo; i++) {
      const o = i * meta.postBytes;
      if (postBuf[o] !== res) continue;
      const id = postBuf[o + 1] | (postBuf[o + 2] << 8) | (postBuf[o + 3] << 16) | (postBuf[o + 4] << 24);
      cnt.set(id, (cnt.get(id) || 0) + 1);
    }
  }
  const top = [...cnt.entries()].sort((a, b) => b[1] - a[1]).slice(0, topK);
  return { top: top.map(([id, c]) => ({ acc: accOf(id), shared: c })), bytes, ms: Date.now() - t0, candidates: cnt.size };
}

/* ---------- queries: real AFDB sequences, lightly mutated ---------- */
let st = 424242;
const rnd = () => { st ^= st << 13; st >>>= 0; st ^= st >>> 17; st ^= st << 5; st >>>= 0; return st / 4294967296; };
const queries = [];
{
  const fd = openSync(partPath, 'r');
  const size = statSync(partPath).size;
  const buf = Buffer.allocUnsafe(1 << 20);
  while (queries.length < N_QUERIES) {
    readSync(fd, buf, 0, 1 << 20, Math.floor(rnd() * (size - (1 << 20))));
    const lines = buf.toString('latin1').split('\n');
    for (let i = 1; i < lines.length - 1 && queries.length < N_QUERIES; i++) {
      if (lines[i].charCodeAt(0) !== 62) continue;
      const s = lines[i + 1];
      if (!s || s.charCodeAt(0) === 62) continue;
      if (s.length < 120 || s.length > 500 || !/^[ACDEFGHIKLMNPQRSTVWY]+$/.test(s)) continue;
      const m = / UA=([A-Z0-9]+)/.exec(lines[i]);
      const mut = s.split('');
      for (let p = 0; p < mut.length; p++) if (rnd() < MUTATE) mut[p] = AA[(rnd() * 20) | 0];
      queries.push({ origin: m ? m[1] : '?', seq: mut.join(''), native: s });
    }
  }
  closeSync(fd);
}
console.log(`${queries.length} queries: real AFDB sequences with ${(MUTATE * 100).toFixed(0)}% of positions mutated\n`);

/* ---------- run both ---------- */
async function blastBest(seq) {
  const t0 = Date.now();
  const job = await Api.blastSubmit(seq, { email: EMAIL, database: DB, nHits: 20, evalue: '1e-3' });
  await Api.blastWait(job, null, null);
  const hits = Api.parseBlastHits(await Api.blastResult(job));
  return { hits, ms: Date.now() - t0 };
}

const results = [];
for (let i = 0; i < queries.length; i += CONCURRENCY) {
  const batch = queries.slice(i, i + CONCURRENCY);
  const got = await Promise.all(batch.map(async q => {
    const idx = indexSearch(q.seq, 20);
    // align what the index returned, using AFDB as the sequence source
    let best = null;
    const seqs = await Api.candidateSequences(idx.top.map(t => t.acc), null, null);
    for (const t of idx.top) {
      const s = seqs.get(t.acc);
      if (!s) continue;
      const aln = A.smithWaterman(q.seq, s);
      if (aln && (!best || aln.identity > best.identity)) best = { acc: t.acc, identity: aln.identity, cov: aln.queryCoverage };
    }
    let bl = null;
    try {
      const r = await blastBest(q.seq);
      const h = r.hits[0];
      if (h) bl = { acc: h.acc, identity: h.identity, cov: h.queryCoverage, ms: r.ms };
      else bl = { acc: null, identity: 0, cov: 0, ms: r.ms };
    } catch (e) { bl = { error: e.message }; }
    return { q, idx, best, bl };
  }));
  results.push(...got);
  for (const r of got) {
    const i1 = r.best ? r.best.identity : 0;
    const i2 = r.bl && !r.bl.error ? r.bl.identity : 0;
    const flag = i1 >= i2 - 2 ? ' ' : '!';
    console.log(`  ${r.q.origin.padEnd(12)} index ${String((r.best?.acc) || '-').padEnd(12)} ${i1.toFixed(0).padStart(3)}%  `
      + `${String(r.idx.ms + 'ms').padStart(6)} ${(r.idx.bytes / 1024).toFixed(0).padStart(4)}KB   `
      + `| blast ${String((r.bl?.acc) || '-').padEnd(12)} ${i2.toFixed(0).padStart(3)}% ${String(Math.round((r.bl?.ms || 0) / 1000) + 's').padStart(5)} ${flag}`);
  }
}

/* ---------- summary ---------- */
const ok = results.filter(r => r.bl && !r.bl.error);
const n = ok.length;
const pct = x => `${(100 * x / n).toFixed(0)}%`;
const med = a => { const s = [...a].sort((x, y) => x - y); return s.length ? s[s.length >> 1] : 0; };
console.log(`\n${'='.repeat(72)}`);
console.log(`${n} queries, index over ${meta.nSeq.toLocaleString()} AFDB entries vs BLAST over ${DB}\n`);
console.log(`  index found a hit at all              ${pct(ok.filter(r => r.best).length)}`);
console.log(`  index identity >= BLAST identity - 2  ${pct(ok.filter(r => (r.best?.identity || 0) >= r.bl.identity - 2).length)}`);
console.log(`  index identity >= 90%                 ${pct(ok.filter(r => (r.best?.identity || 0) >= 90).length)}`);
console.log(`  BLAST identity >= 90%                 ${pct(ok.filter(r => r.bl.identity >= 90).length)}`);
console.log(`\n  median identity, index                ${med(ok.map(r => r.best?.identity || 0)).toFixed(0)}%`);
console.log(`  median identity, BLAST                ${med(ok.map(r => r.bl.identity)).toFixed(0)}%`);
console.log(`\n  median index time                     ${med(ok.map(r => r.idx.ms))} ms`);
console.log(`  median BLAST time                     ${(med(ok.map(r => r.bl.ms)) / 1000).toFixed(0)} s`);
console.log(`  median bytes read by index            ${(med(ok.map(r => r.idx.bytes)) / 1024).toFixed(0)} KB`);
