/*
 * node tools/afdb-shard.mjs <part.fa> <outDir> <shardIndex>
 *
 * One shard of the AFDB minimizer index: read a part of sequences.fasta, emit
 * (seed, localEntryId) postings sorted by seed, plus the accession list.
 * tools/afdb-merge.mjs then merges the shards and renumbers ids globally.
 *
 * Two things this has to get right:
 *
 * 1. Streaming. A part is ~4.9 GB and V8 caps strings near 512 MB, so
 *    readFileSync(path,'utf8') throws. Reads go through fixed Buffers with a
 *    carry for the line split across chunk boundaries.
 *
 * 2. Length-adaptive w. A fixed window undersamples short sequences -- at k=12
 *    w=16 a 142-residue globin got 8 seeds and shared 0 of them with a real 90%
 *    relative on two of three tries. Recall depends on the *number* of seeds,
 *    not the window, since a 10-mer survives 90% identity with p=0.9^10=0.35 and
 *    12 seeds give 1-0.65^12 = 99.2%. So w is chosen per sequence to land near
 *    TARGET_SEEDS regardless of length, which also stops long proteins from
 *    dominating the posting lists.
 */
import { openSync, readSync, closeSync, statSync, writeFileSync, createWriteStream } from 'node:fs';
import { join } from 'node:path';

const [, , partPath, outDir, shardArg] = process.argv;
if (!partPath || !outDir || shardArg === undefined) {
  console.error('usage: node tools/afdb-shard.mjs <part.fa> <outDir> <shardIndex>');
  process.exit(2);
}
const SHARD = Number(shardArg);
const K = Number(process.env.K || 10);
const TARGET_SEEDS = Number(process.env.TARGET_SEEDS || 12);
const MIN_W = 4, MAX_W = 128;
const MIN_LEN = Number(process.env.MIN_LEN || 40);
const MAX_LEN = Number(process.env.MAX_LEN || 2000);

const t0 = Date.now();
const el = () => `${((Date.now() - t0) / 1000).toFixed(0)}s`;

const AA = 'ACDEFGHIKLMNPQRSTVWY';
const IDX = new Int8Array(128).fill(-1);
for (let i = 0; i < AA.length; i++) IDX[AA.charCodeAt(i)] = i;

function h32(x) {
  x = (x ^ 0x9e3779b9) >>> 0;
  x = Math.imul(x ^ (x >>> 16), 0x85ebca6b) >>> 0;
  x = Math.imul(x ^ (x >>> 13), 0xc2b2ae35) >>> 0;
  return (x ^ (x >>> 16)) >>> 0;
}

// scratch reused across sequences so the hot loop allocates nothing
let codeBuf = new Int8Array(MAX_LEN);
let hashBuf = new Int32Array(MAX_LEN);

function seedsOf(seq, out) {
  const n = seq.length;
  if (codeBuf.length < n) { codeBuf = new Int8Array(n); hashBuf = new Int32Array(n); }
  for (let i = 0; i < n; i++) codeBuf[i] = IDX[seq.charCodeAt(i) & 0x7f];
  hashBuf.fill(-1, 0, n);
  for (let i = K - 1; i < n; i++) {
    let x = 0, ok = true;
    for (let j = i - K + 1; j <= i; j++) {
      const c = codeBuf[j];
      if (c < 0) { ok = false; break; }
      x = (Math.imul(x, 31) + c) | 0;
    }
    if (ok) hashBuf[i] = h32(x >>> 0) & 0x7fffffff;
  }
  const w = Math.max(MIN_W, Math.min(MAX_W, Math.floor(n / TARGET_SEEDS) || MIN_W));
  out.clear();
  for (let i = K - 1; i < n; i++) {
    if (i % w) continue;
    let best = -1;
    const stop = Math.min(i + w, n);
    for (let j = i; j < stop; j++) if (hashBuf[j] >= 0 && (best < 0 || hashBuf[j] < best)) best = hashBuf[j];
    if (best >= 0) out.add(best);
  }
  return out;
}

/* ---------- pass over the part ---------- */

const accOut = createWriteStream(join(outDir, `acc-${String(SHARD).padStart(2, '0')}.txt`));
let keys = new Uint32Array(1 << 24);
let vals = new Uint32Array(1 << 24);
let nPost = 0;
let nSeq = 0, tooShort = 0, tooLong = 0, badChars = 0;

function push(seed, id) {
  if (nPost === keys.length) {
    const k2 = new Uint32Array(keys.length * 2); k2.set(keys); keys = k2;
    const v2 = new Uint32Array(vals.length * 2); v2.set(vals); vals = v2;
  }
  keys[nPost] = seed; vals[nPost] = id; nPost++;
}

const accRe = / UA=([A-Z0-9]+)/;
{
  const fd = openSync(partPath, 'r');
  const size = statSync(partPath).size;
  const CH = 1 << 26;
  const buf = Buffer.allocUnsafe(CH);
  const seeds = new Set();
  let carry = '';
  let pos = 0;
  let pendingAcc = null;
  const accBatch = [];

  while (pos < size) {
    const got = readSync(fd, buf, 0, CH, pos);
    if (!got) break;
    pos += got;
    const txt = carry + buf.toString('latin1', 0, got);
    const lines = txt.split('\n');
    carry = lines.pop();

    for (const line of lines) {
      if (!line) continue;
      if (line.charCodeAt(0) === 62) {
        const m = accRe.exec(line);
        pendingAcc = m ? m[1] : null;
        continue;
      }
      if (pendingAcc === null) continue;
      const seq = line;
      const acc = pendingAcc;
      pendingAcc = null;
      if (seq.length < MIN_LEN) { tooShort++; continue; }
      if (seq.length > MAX_LEN) { tooLong++; continue; }
      seedsOf(seq, seeds);
      if (!seeds.size) { badChars++; continue; }
      const id = nSeq++;
      for (const s of seeds) push(s, id);
      accBatch.push(acc);
      if (accBatch.length >= 100000) { accOut.write(accBatch.join('\n') + '\n'); accBatch.length = 0; }
    }
    if (nSeq && nSeq % 2_000_000 < 5000) {
      console.error(`  [${SHARD}] ${(pos / 1e9).toFixed(1)} GB, ${nSeq} seqs, ${(nPost / 1e6).toFixed(0)}M postings (${el()})`);
    }
  }
  if (accBatch.length) accOut.write(accBatch.join('\n') + '\n');
  closeSync(fd);
}
accOut.end();

/* ---------- sort by seed ---------- */

function radix(k, v, n) {
  let ki = k, vi = v, ko = new Uint32Array(n), vo = new Uint32Array(n);
  const c = new Uint32Array(257);
  for (let sh = 0; sh < 32; sh += 8) {
    c.fill(0);
    for (let i = 0; i < n; i++) c[((ki[i] >>> sh) & 255) + 1]++;
    for (let i = 0; i < 256; i++) c[i + 1] += c[i];
    for (let i = 0; i < n; i++) { const b = (ki[i] >>> sh) & 255; const p = c[b]++; ko[p] = ki[i]; vo[p] = vi[i]; }
    let t = ki; ki = ko; ko = t; t = vi; vi = vo; vo = t;
  }
  return { keys: ki, vals: vi };
}

const sorted = radix(keys.subarray(0, nPost), vals.subarray(0, nPost), nPost);
const tag = String(SHARD).padStart(2, '0');
writeFileSync(join(outDir, `keys-${tag}.u32`), Buffer.from(sorted.keys.buffer, 0, nPost * 4));
writeFileSync(join(outDir, `vals-${tag}.u32`), Buffer.from(sorted.vals.buffer, 0, nPost * 4));
writeFileSync(join(outDir, `meta-${tag}.json`), JSON.stringify({
  shard: SHARD, part: partPath, nSeq, nPost, K, targetSeeds: TARGET_SEEDS,
  seedsPerSeq: nSeq ? nPost / nSeq : 0, tooShort, tooLong, badChars,
}, null, 2));

console.error(`  [${SHARD}] done: ${nSeq} seqs, ${(nPost / 1e6).toFixed(1)}M postings, `
  + `${nSeq ? (nPost / nSeq).toFixed(1) : 0} seeds/seq, skipped ${tooShort}/${tooLong} short/long (${el()})`);
