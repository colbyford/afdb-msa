/*
 * node tools/afdb-coverage.mjs <indexPart.fa> <queryPart.fa>
 *
 * SwissProt has a >=90% relative for only 6% of random AFDB sequences, so the
 * corpus has to be AFDB itself. The open question is how much of AFDB you must
 * index: if the database is redundant enough, a few million representatives
 * cover most of it and 1 GB is plenty; if it is genuinely diverse, only the full
 * index works and the budget is ~14 GB.
 *
 * So: index N random AFDB entries for several N, query with AFDB sequences drawn
 * from a *different* part of the file, and measure how often a >=90% relative is
 * found. The shape of coverage-vs-N is the answer.
 *
 * Sequences are held in one concatenated buffer with an offset table rather than
 * as N JavaScript strings -- at 10M entries the latter is tens of gigabytes of
 * object overhead for 3.5 GB of residues.
 */
import { createRequire } from 'node:module';
import { openSync, readSync, closeSync, statSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const require = createRequire(import.meta.url);
const A = require(join(dirname(fileURLToPath(import.meta.url)), '..', 'align.js'));

const [, , indexPart, queryPart] = process.argv;
if (!indexPart || !queryPart) {
  console.error('usage: node tools/afdb-coverage.mjs <indexPart.fa> <queryPart.fa>');
  process.exit(2);
}
const W = Number(process.env.W || 32);           // sparse seeding: 11 seeds/entry
const K = 10;
const N_QUERIES = Number(process.env.N_QUERIES || 300);
const SIZES = (process.env.SIZES || '250000,1000000,4000000,10000000').split(',').map(Number);
const MAXN = Math.max(...SIZES);

const t0 = Date.now();
const el = () => `${((Date.now() - t0) / 1000).toFixed(0)}s`;

const AA = 'ACDEFGHIKLMNPQRSTVWY';
const IDX = new Int8Array(128).fill(-1);
for (let i = 0; i < AA.length; i++) IDX[AA.charCodeAt(i)] = i;
const h32 = x => { x = (x ^ 0x9e3779b9) >>> 0; x = Math.imul(x ^ (x >>> 16), 0x85ebca6b) >>> 0; x = Math.imul(x ^ (x >>> 13), 0xc2b2ae35) >>> 0; return (x ^ (x >>> 16)) >>> 0; };

function minimizers(seq) {
  const n = seq.length;
  const c = new Int8Array(n);
  for (let i = 0; i < n; i++) c[i] = IDX[seq.charCodeAt(i) & 0x7f];
  const hs = new Int32Array(n).fill(-1);
  for (let i = K - 1; i < n; i++) {
    let x = 0, ok = true;
    for (let j = i - K + 1; j <= i; j++) { const q = c[j]; if (q < 0) { ok = false; break; } x = (Math.imul(x, 31) + q) | 0; }
    if (ok) hs[i] = h32(x >>> 0) & 0x7fffffff;
  }
  const s = new Set();
  for (let i = K - 1; i < n; i++) {
    if (i % W) continue;
    let b = -1;
    for (let j = i; j < Math.min(i + W, n); j++) if (hs[j] >= 0 && (b < 0 || hs[j] < b)) b = hs[j];
    if (b >= 0) s.add(b);
  }
  return s;
}

/* ---------- read the index corpus into a packed buffer ---------- */
console.log(`reading up to ${MAXN} sequences from ${indexPart}…`);
const CHUNK = 1 << 26;
let store = Buffer.allocUnsafe(1 << 30);
let storeLen = 0;
const offs = new Int32Array(MAXN + 1);
let nSeq = 0;
{
  const fd = openSync(indexPart, 'r');
  const buf = Buffer.allocUnsafe(CHUNK);
  let carry = '';
  let pos = 0;
  const size = statSync(indexPart).size;
  outer:
  while (pos < size) {
    const got = readSync(fd, buf, 0, CHUNK, pos);
    if (!got) break;
    pos += got;
    const txt = carry + buf.toString('latin1', 0, got);
    const lines = txt.split('\n');
    carry = lines.pop();
    for (let i = 0; i < lines.length - 1; i++) {
      if (lines[i].charCodeAt(0) !== 62) continue;
      const s = lines[i + 1];
      if (!s || s.charCodeAt(0) === 62) continue;
      const L = s.length;
      if (L < 50 || L > 1500) continue;
      if (storeLen + L > store.length) {
        const bigger = Buffer.allocUnsafe(store.length * 2);
        store.copy(bigger, 0, 0, storeLen);
        store = bigger;
      }
      store.write(s, storeLen, 'latin1');
      offs[nSeq] = storeLen;
      storeLen += L;
      offs[++nSeq] = storeLen;
      if (nSeq >= MAXN) break outer;
    }
  }
  closeSync(fd);
}
console.log(`  ${nSeq} sequences, ${(storeLen / 1e9).toFixed(2)} GB of residues (${el()})`);
const seqAt = i => store.toString('latin1', offs[i], offs[i + 1]);

/* ---------- queries from a different part of the file ---------- */
const queries = [];
{
  const fd = openSync(queryPart, 'r');
  const size = statSync(queryPart).size;
  const buf = Buffer.allocUnsafe(1 << 20);
  let st = 20260807;
  const rnd = () => { st ^= st << 13; st >>>= 0; st ^= st >>> 17; st ^= st << 5; st >>>= 0; return st / 4294967296; };
  while (queries.length < N_QUERIES) {
    readSync(fd, buf, 0, 1 << 20, Math.floor(rnd() * (size - (1 << 20))));
    const lines = buf.toString('latin1').split('\n');
    for (let i = 1; i < lines.length - 1 && queries.length < N_QUERIES; i++) {
      if (lines[i].charCodeAt(0) !== 62) continue;
      const s = lines[i + 1];
      if (!s || s.charCodeAt(0) === 62) continue;
      if (s.length >= 80 && s.length <= 1200 && /^[ACDEFGHIKLMNPQRSTVWY]+$/.test(s)) queries.push(s);
    }
  }
  closeSync(fd);
}
console.log(`  ${queries.length} query sequences from ${queryPart}\n`);

function radix(keys, vals, n) {
  let ki = keys, vi = vals, ko = new Uint32Array(n), vo = new Uint32Array(n);
  const c = new Uint32Array(257);
  for (let sh = 0; sh < 32; sh += 8) {
    c.fill(0);
    for (let i = 0; i < n; i++) c[((ki[i] >>> sh) & 255) + 1]++;
    for (let i = 0; i < 256; i++) c[i + 1] += c[i];
    for (let i = 0; i < n; i++) { const b = (ki[i] >>> sh) & 255, p = c[b]++; ko[p] = ki[i]; vo[p] = vi[i]; }
    let t = ki; ki = ko; ko = t; t = vi; vi = vo; vo = t;
  }
  return { keys: ki, vals: vi };
}

console.log('  index size    postings   index MB   >=90%   >=70%   >=50%   median id');
console.log('  ' + '-'.repeat(72));

for (const size of SIZES) {
  if (size > nSeq) continue;
  let tot = 0;
  const per = new Array(size);
  for (let i = 0; i < size; i++) { const s = minimizers(seqAt(i)); per[i] = s; tot += s.size; }
  let keys = new Uint32Array(tot), vals = new Uint32Array(tot);
  { let p = 0; for (let i = 0; i < size; i++) { for (const h of per[i]) { keys[p] = h; vals[p] = i; p++; } per[i] = null; } }
  ({ keys, vals } = radix(keys, vals, tot));
  const lb = k => { let lo = 0, hi = tot; while (lo < hi) { const m = (lo + hi) >>> 1; if (keys[m] < k) lo = m + 1; else hi = m; } return lo; };

  const best = [];
  for (const q of queries) {
    const qs = minimizers(q);
    const cnt = new Map();
    for (const h of qs) {
      const lo = lb(h);
      if (lo >= tot || keys[lo] !== h) continue;
      let hi = lo;
      while (hi < tot && keys[hi] === h) hi++;
      if (hi - lo > 2000) continue;
      for (let x = lo; x < hi; x++) { const e = vals[x]; cnt.set(e, (cnt.get(e) || 0) + 1); }
    }
    const top = [...cnt.entries()].sort((a, b) => b[1] - a[1]).slice(0, 10);
    let b2 = 0;
    for (const [e] of top) { const al = A.smithWaterman(q, seqAt(e)); if (al && al.identity > b2) b2 = al.identity; }
    best.push(b2);
  }
  best.sort((a, b) => a - b);
  const pct = x => (100 * best.filter(v => v >= x).length / best.length).toFixed(0) + '%';
  console.log(`  ${String(size).padStart(10)}  ${String((tot / 1e6).toFixed(0) + 'M').padStart(10)}  `
    + `${String((tot * 5 / 1e6).toFixed(0)).padStart(8)}   ${pct(90).padStart(5)}   ${pct(70).padStart(5)}   ${pct(50).padStart(5)}   `
    + `${String(best[best.length >> 1].toFixed(0) + '%').padStart(9)}   (${el()})`);
}
