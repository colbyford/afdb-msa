/*
 * node tools/minimizer-broad.mjs <swissprot.fasta> <famDir> [k] [w]
 *
 * The 100% recall figure for minimizer seeding came from 311 queries across six
 * families -- too narrow to trust, and measured only at >=90% identity, so it
 * said nothing about where recall falls off.
 *
 * This widens both: every AFDB MSA in famDir contributes queries (each family's
 * own SwissProt accession is the reference), and recall is reported per identity
 * bin. The bins below 90% are not a target -- they are what decides where the
 * fallback to BLAST has to kick in, which is currently set by assumption.
 */
import { createRequire } from 'node:module';
import { readFileSync, readdirSync } from 'node:fs';
import { dirname, join, basename } from 'node:path';
import { fileURLToPath } from 'node:url';

const require = createRequire(import.meta.url);
const here = dirname(fileURLToPath(import.meta.url));
const A = require(join(here, '..', 'align.js'));

const [, , fastaPath, famDir, kArg, wArg] = process.argv;
if (!fastaPath || !famDir) {
  console.error('usage: node tools/minimizer-broad.mjs <swissprot.fasta> <famDir> [k] [w]');
  process.exit(2);
}
const K = Number(kArg || 10);
const W = Number(wArg || 8);
const PER_FAM = Number(process.env.PER_FAM || 60);
const PER_BAND = Number(process.env.PER_BAND || 10);
const SCAN_LIMIT = Number(process.env.SCAN_LIMIT || 1500);
const BAND_EDGES = [0, 30, 50, 70, 80, 90, 95];
const TOPK = Number(process.env.TOPK || 10);
const CAP = Number(process.env.CAP || 2000);

const t0 = Date.now();
const el = () => `${((Date.now() - t0) / 1000).toFixed(0)}s`;

/* ---------- corpus ---------- */
console.log('loading corpus…');
const accs = [], seqs = [];
{
  const text = readFileSync(fastaPath, 'utf8');
  let acc = null, buf = [];
  const flush = () => { if (acc && buf.length) { const s = buf.join(''); if (s.length >= 30) { accs.push(acc); seqs.push(s); } } buf = []; };
  for (const line of text.split('\n')) {
    if (!line) continue;
    if (line.charCodeAt(0) === 62) { flush(); const m = /^>[a-z]{2}\|([^|]+)\|/i.exec(line); acc = m ? m[1] : null; }
    else if (acc) buf.push(line.trim());
  }
  flush();
}
const N = accs.length;
const accIndex = new Map(accs.map((a, i) => [a, i]));
console.log(`  ${N} sequences (${el()})`);

/* ---------- minimizers ---------- */
const AA = 'ACDEFGHIKLMNPQRSTVWY';
const IDX = new Int8Array(128).fill(-1);
for (let i = 0; i < AA.length; i++) IDX[AA.charCodeAt(i)] = i;

function h32(x) {
  x = (x ^ 0x9e3779b9) >>> 0;
  x = Math.imul(x ^ (x >>> 16), 0x85ebca6b) >>> 0;
  x = Math.imul(x ^ (x >>> 13), 0xc2b2ae35) >>> 0;
  return (x ^ (x >>> 16)) >>> 0;
}

function minimizers(seq, k, w) {
  const n = seq.length;
  const code = new Int8Array(n);
  for (let i = 0; i < n; i++) code[i] = IDX[seq.charCodeAt(i) & 0x7f];
  const hs = new Int32Array(n).fill(-1);
  for (let i = k - 1; i < n; i++) {
    let x = 0, ok = true;
    for (let j = i - k + 1; j <= i; j++) {
      const c = code[j];
      if (c < 0) { ok = false; break; }
      x = (Math.imul(x, 31) + c) | 0;
    }
    if (ok) hs[i] = h32(x >>> 0) & 0x7fffffff;
  }
  const set = new Set();
  for (let i = k - 1; i < n; i++) {
    if (i % w) continue;
    let best = -1;
    for (let j = i; j < Math.min(i + w, n); j++) if (hs[j] >= 0 && (best < 0 || hs[j] < best)) best = hs[j];
    if (best >= 0) set.add(best);
  }
  return set;
}

/* ---------- queries from every family ---------- */
console.log(`\nbuilding queries from ${famDir}…`);
let st = 20260807;
const rnd = () => { st ^= st << 13; st >>>= 0; st ^= st >>> 17; st ^= st << 5; st >>>= 0; return st / 4294967296; };

const queries = [];
let famUsed = 0;
for (const f of readdirSync(famDir)) {
  if (!f.endsWith('.a3m')) continue;
  const acc = basename(f, '.a3m');
  if (!accIndex.has(acc)) continue;
  const target = seqs[accIndex.get(acc)];
  const rows = readFileSync(join(famDir, f), 'utf8').split('\n');
  const mem = [];
  for (let i = 1; i < rows.length; i += 2) {
    if (!rows[i] || rows[i][0] === '>') continue;
    const s = rows[i].replace(/[-.]/g, '').toUpperCase();
    if (s.length > 60) mem.push(s);
  }
  if (mem.length < 10) continue;
  for (let i = mem.length - 1; i > 0; i--) { const j = (rnd() * (i + 1)) | 0; const t = mem[i]; mem[i] = mem[j]; mem[j] = t; }
  /*
   * MSA members are overwhelmingly divergent -- a flat sample gave 2224 queries
   * at 30-50% identity but only 24 above 90%, which is the band that matters.
   * So fill a quota per identity band instead, scanning further into the MSA
   * until the high bands are populated or the members run out.
   */
  const quota = new Map(BAND_EDGES.map(b => [b, 0]));
  const bandOf = id => BAND_EDGES.filter(b => id >= b).pop() ?? BAND_EDGES[0];
  let n = 0, scanned = 0;
  for (const s of mem) {
    if (n >= PER_FAM || scanned > SCAN_LIMIT) break;
    scanned++;
    const aln = A.smithWaterman(s, target);
    if (!aln || aln.queryCoverage < 0.5) continue;
    const b = bandOf(aln.identity);
    if (quota.get(b) >= PER_BAND) continue;
    quota.set(b, quota.get(b) + 1);
    queries.push({ seq: s, identity: aln.identity, fam: acc });
    n++;
  }
  famUsed++;
}
console.log(`  ${queries.length} queries from ${famUsed} families (${el()})`);

/* ---------- build the index ---------- */
console.log(`\nbuilding minimizer index k=${K} w=${W}…`);
function radixSortPairs(keys, vals, n) {
  let kIn = keys, vIn = vals;
  let kOut = new Uint32Array(n), vOut = new Uint32Array(n);
  const count = new Uint32Array(257);
  for (let shift = 0; shift < 32; shift += 8) {
    count.fill(0);
    for (let i = 0; i < n; i++) count[((kIn[i] >>> shift) & 0xff) + 1]++;
    for (let i = 0; i < 256; i++) count[i + 1] += count[i];
    for (let i = 0; i < n; i++) { const b = (kIn[i] >>> shift) & 0xff; const p = count[b]++; kOut[p] = kIn[i]; vOut[p] = vIn[i]; }
    let t = kIn; kIn = kOut; kOut = t; t = vIn; vIn = vOut; vOut = t;
  }
  return { keys: kIn, vals: vIn };
}

let total = 0;
const per = new Array(N);
for (let i = 0; i < N; i++) { const s = minimizers(seqs[i], K, W); per[i] = s; total += s.size; }
let keys = new Uint32Array(total), vals = new Uint32Array(total);
{ let p = 0; for (let i = 0; i < N; i++) { for (const s of per[i]) { keys[p] = s; vals[p] = i; p++; } per[i] = null; } }
({ keys, vals } = radixSortPairs(keys, vals, total));
console.log(`  ${(total / 1e6).toFixed(1)}M postings, ${(total * 4 / 1e6).toFixed(0)} MB (${el()})`);

const lowerBound = key => {
  let lo = 0, hi = total;
  while (lo < hi) { const mid = (lo + hi) >>> 1; if (keys[mid] < key) lo = mid + 1; else hi = mid; }
  return lo;
};

/* ---------- evaluate by identity ---------- */
const BINS = [[0, 30], [30, 50], [50, 70], [70, 80], [80, 90], [90, 95], [95, 101]];
const LAB = ['<30%', '30-50%', '50-70%', '70-80%', '80-90%', '90-95%', '95-100%'];
const S = BINS.map(() => ({ n: 0, hit: 0, cands: 0, lookups: 0, fams: new Set() }));

console.log('\nevaluating…');
for (let qi = 0; qi < queries.length; qi++) {
  const q = queries[qi];
  const b = BINS.findIndex(([lo, hi]) => q.identity >= lo && q.identity < hi);
  if (b < 0) continue;
  const st2 = S[b];
  st2.n++; st2.fams.add(q.fam);

  const qs = minimizers(q.seq, K, W);
  const counts = new Map();
  for (const s of qs) {
    const lo = lowerBound(s);
    if (lo >= total || keys[lo] !== s) continue;
    let hi = lo;
    while (hi < total && keys[hi] === s) hi++;
    if (hi - lo > CAP) continue;
    st2.lookups++;
    for (let x = lo; x < hi; x++) { const i = vals[x]; counts.set(i, (counts.get(i) || 0) + 1); }
  }
  st2.cands += counts.size;
  const top = [...counts.entries()].sort((a, c) => c[1] - a[1]).slice(0, TOPK);
  for (const [i] of top) {
    const aln = A.smithWaterman(q.seq, seqs[i]);
    if (aln && aln.identity >= q.identity - 2) { st2.hit++; break; }   // found something at least as good
  }
  if ((qi + 1) % 500 === 0) console.log(`  ${qi + 1}/${queries.length} (${el()})`);
}

console.log(`\nminimizer seeding, k=${K} w=${W}, top${TOPK} aligned, ${N} sequence corpus`);
console.log('"recall" = retrieved something at least as similar as the query\'s own reference\n');
console.log('  identity      n   fams   recall   cands/query   lookups');
console.log('  ' + '-'.repeat(60));
BINS.forEach((_, i) => {
  const x = S[i];
  if (!x.n) return;
  console.log(`  ${LAB[i].padEnd(9)} ${String(x.n).padStart(5)}  ${String(x.fams.size).padStart(5)}   `
    + `${String((100 * x.hit / x.n).toFixed(0) + '%').padStart(6)}   ${String((x.cands / x.n).toFixed(0)).padStart(11)}   ${String((x.lookups / x.n).toFixed(0)).padStart(7)}`);
});
console.log(`\ndone in ${el()}`);
