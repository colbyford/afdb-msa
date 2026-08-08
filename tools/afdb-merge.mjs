/*
 * node tools/afdb-merge.mjs <shardDir> <outDir>
 *
 * Merges the per-shard postings into the published index.
 *
 * Layout written, chosen so a browser can answer a query with one ranged GET
 * per seed and nothing else:
 *
 *   meta.json        entry count, k, target seeds, bucket bits, shard table
 *   buckets.u32      2^BITS+1 prefix offsets: seed's high bits -> posting range
 *   post-NNN.bin     5 bytes per posting: 1 byte key residual + 4 byte entry id
 *   acc.txt          accessions, one per line, indexed by entry id
 *
 * Seeds are 31-bit hashes, so the top BITS bits address a bucket directly and
 * the remaining bits are stored as a residual next to the entry id. A lookup is
 * therefore: read buckets[hi] and buckets[hi+1] (both already in memory from a
 * one-off 128 MB download), then one ranged GET of that byte range, then scan
 * for the matching residual. No binary search over the wire, which would have
 * cost ~25 round trips per seed.
 *
 * Shards are merged by bucket rather than by a k-way heap: every shard's
 * postings are already sorted by seed, so each contributes a contiguous run to
 * each bucket, and the runs concatenate. Order within a bucket does not matter
 * because the reader scans it.
 */
import { readFileSync, writeFileSync, openSync, readSync, closeSync, statSync, existsSync, mkdirSync, createWriteStream } from 'node:fs';
import { join } from 'node:path';

const [, , shardDir, outDir] = process.argv;
if (!shardDir || !outDir) {
  console.error('usage: node tools/afdb-merge.mjs <shardDir> <outDir>');
  process.exit(2);
}
const BITS = Number(process.env.BUCKET_BITS || 25);      // 33.5M buckets
const NBUCK = 1 << BITS;
const SHIFT = 31 - BITS;                                  // seeds are 31-bit
const SHARD_BYTES = Number(process.env.SHARD_BYTES || 1_500_000_000);

const t0 = Date.now();
const el = () => `${((Date.now() - t0) / 1000).toFixed(0)}s`;
if (!existsSync(outDir)) mkdirSync(outDir, { recursive: true });

/* ---------- shard inventory and global id offsets ---------- */
const metas = [];
for (let i = 0; i < 64; i++) {
  const tag = String(i).padStart(2, '0');
  const f = join(shardDir, `meta-${tag}.json`);
  if (!existsSync(f)) continue;
  const m = JSON.parse(readFileSync(f, 'utf8'));
  m.tag = tag;
  metas.push(m);
}
metas.sort((a, b) => a.shard - b.shard);
let base = 0;
for (const m of metas) { m.idBase = base; base += m.nSeq; }
const nSeq = base;
const nPost = metas.reduce((a, m) => a + m.nPost, 0);
console.log(`${metas.length} shards, ${nSeq.toLocaleString()} sequences, ${nPost.toLocaleString()} postings (${el()})`);

/* ---------- pass 1: bucket histogram ---------- */
console.log(`counting into ${NBUCK.toLocaleString()} buckets…`);
const counts = new Uint32Array(NBUCK + 1);
const CH = 1 << 24;
const kbuf = Buffer.allocUnsafe(CH * 4);
for (const m of metas) {
  const fd = openSync(join(shardDir, `keys-${m.tag}.u32`), 'r');
  let off = 0;
  const total = m.nPost * 4;
  while (off < total) {
    const got = readSync(fd, kbuf, 0, Math.min(kbuf.length, total - off), off);
    if (!got) break;
    const arr = new Uint32Array(kbuf.buffer, kbuf.byteOffset, got / 4);
    for (let i = 0; i < arr.length; i++) counts[(arr[i] >>> SHIFT) + 1]++;
    off += got;
  }
  closeSync(fd);
  console.log(`  shard ${m.tag} counted (${el()})`);
}
for (let i = 0; i < NBUCK; i++) counts[i + 1] += counts[i];
console.log(`  prefix sums done, total ${counts[NBUCK].toLocaleString()} (${el()})`);

/* ---------- plan output shards on bucket boundaries ---------- */
const POST_BYTES = 5;
const files = [];
{
  let startBucket = 0, startPost = 0;
  for (let b = 1; b <= NBUCK; b++) {
    const bytes = (counts[b] - startPost) * POST_BYTES;
    if (bytes >= SHARD_BYTES || b === NBUCK) {
      files.push({ from: startBucket, to: b, postFrom: startPost, postTo: counts[b] });
      startBucket = b; startPost = counts[b];
    }
  }
}
console.log(`  ${files.length} output files, <= ${(SHARD_BYTES / 1e9).toFixed(1)} GB each (${el()})`);

/* ---------- pass 2: scatter postings into place ---------- */
const cursor = Uint32Array.from(counts.subarray(0, NBUCK));
const vbuf = Buffer.allocUnsafe(CH * 4);

for (const [fi, f] of files.entries()) {
  const n = f.postTo - f.postFrom;
  const out = Buffer.allocUnsafe(n * POST_BYTES);
  for (const m of metas) {
    const kfd = openSync(join(shardDir, `keys-${m.tag}.u32`), 'r');
    const vfd = openSync(join(shardDir, `vals-${m.tag}.u32`), 'r');
    let off = 0;
    const total = m.nPost * 4;
    while (off < total) {
      const want = Math.min(kbuf.length, total - off);
      readSync(kfd, kbuf, 0, want, off);
      readSync(vfd, vbuf, 0, want, off);
      const ka = new Uint32Array(kbuf.buffer, kbuf.byteOffset, want / 4);
      const va = new Uint32Array(vbuf.buffer, vbuf.byteOffset, want / 4);
      for (let i = 0; i < ka.length; i++) {
        const b = ka[i] >>> SHIFT;
        if (b < f.from || b >= f.to) continue;
        const p = cursor[b]++ - f.postFrom;
        const o = p * POST_BYTES;
        out[o] = ka[i] & 0xff;                       // key residual
        const id = m.idBase + va[i];
        out[o + 1] = id & 0xff; out[o + 2] = (id >>> 8) & 0xff;
        out[o + 3] = (id >>> 16) & 0xff; out[o + 4] = (id >>> 24) & 0xff;
      }
      off += want;
    }
    closeSync(kfd); closeSync(vfd);
  }
  const name = `post-${String(fi).padStart(3, '0')}.bin`;
  writeFileSync(join(outDir, name), out);
  f.name = name;
  console.log(`  wrote ${name}: ${(out.length / 1e9).toFixed(2)} GB, buckets ${f.from}-${f.to} (${el()})`);
}

/* ---------- bucket table + accessions + meta ---------- */
writeFileSync(join(outDir, 'buckets.u32'), Buffer.from(counts.buffer, 0, (NBUCK + 1) * 4));

// Concatenate the per-shard accession lists. Writing multi-hundred-MB Buffers
// into a stream without awaiting drain queues them all in memory and dies with
// 'writev failed' -- and leaves a file padded to exactly 2^32 bytes. Copy in
// bounded chunks and respect backpressure instead.
{
  const out = createWriteStream(join(outDir, 'acc.txt'));
  const CHUNK = 1 << 24;
  const buf = Buffer.allocUnsafe(CHUNK);
  for (const m of metas) {
    const path = join(shardDir, `acc-${m.tag}.txt`);
    const fd = openSync(path, 'r');
    const size = statSync(path).size;
    let off = 0;
    while (off < size) {
      const got = readSync(fd, buf, 0, Math.min(CHUNK, size - off), off);
      if (!got) break;
      if (!out.write(Buffer.from(buf.subarray(0, got)))) {
        await new Promise(r => out.once('drain', r));
      }
      off += got;
    }
    closeSync(fd);
  }
  await new Promise(r => out.end(r));
}

writeFileSync(join(outDir, 'meta.json'), JSON.stringify({
  nSeq, nPost, k: metas[0].K, targetSeeds: metas[0].targetSeeds,
  bucketBits: BITS, shift: SHIFT, postBytes: POST_BYTES,
  files: files.map(f => ({ name: f.name, fromBucket: f.from, toBucket: f.to, postFrom: f.postFrom, postTo: f.postTo })),
}, null, 2));

console.log(`\nindex: ${nSeq.toLocaleString()} sequences, ${nPost.toLocaleString()} postings`);
console.log(`  postings ${(nPost * POST_BYTES / 1e9).toFixed(1)} GB + buckets ${((NBUCK + 1) * 4 / 1e6).toFixed(0)} MB`);
console.log(`done in ${el()}`);
