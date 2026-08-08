/*
 * node tools/select-entries.mjs <msa_depths.csv> <clusters.tsv.gz> <out.txt> [maxEntries]
 *
 * Chooses which AFDB entries to index. The structural-cluster representative
 * list turned out to be the wrong criterion -- measured MSA depth over a 3.2M
 * row sample:
 *
 *   cluster representatives            median depth   200   12.8% under 10 seqs
 *   everything else (mostly singleton) median depth 9,140    5.2% under 10 seqs
 *
 * Depth is what actually matters, because a donor with a 10-sequence MSA is
 * useless to borrow from no matter how well retrieval finds it. So rank by depth
 * instead, and use the cluster file only to suppress redundancy: from any one
 * structural cluster keep a single entry -- the deepest -- rather than the
 * arbitrary representative Foldseek picked, and let singletons through on their
 * own merit.
 *
 * Output: one accession per line, deepest first.
 */
import { readFileSync, writeFileSync, createReadStream } from 'node:fs';
import { createGunzip } from 'node:zlib';
import { createInterface } from 'node:readline';

const [, , depthsPath, clustersPath, outPath, maxArg] = process.argv;
if (!depthsPath || !clustersPath || !outPath) {
  console.error('usage: node tools/select-entries.mjs <msa_depths.csv> <clusters.tsv.gz> <out.txt> [maxEntries]');
  process.exit(2);
}
const MAX = Number(maxArg || 10_000_000);
const MIN_DEPTH = Number(process.env.MIN_DEPTH || 100);

const t0 = Date.now();
const el = () => `${((Date.now() - t0) / 1000).toFixed(0)}s`;

/*
 * Pass 1: entry -> cluster representative, but ONLY for entries in a
 * multi-member cluster.
 *
 * The cluster file has 199,373,567 rows -- one per AFDB entry, including
 * single-member clusters written as "X <tab> X". Structural clustering
 * therefore reduces almost nothing on its own: just 42M of those 199M entries
 * live in a cluster with more than one member, and the other ~157M are their
 * own. Tracking all of them needs a half-billion-slot table and buys nothing.
 *
 * So the reps file (2.6M genuine multi-member clusters) is loaded first and used
 * as a filter: only rows whose repId is one of those are recorded. That is 42M
 * entries, which fits comfortably, and every other entry is by definition a
 * singleton that competes on depth alone.
 */
console.log(`reading cluster membership from ${clustersPath}…`);

function fnv1a(s) {
  let h = 0x811c9dc5 >>> 0;
  for (let i = 0; i < s.length; i++) h = Math.imul(h ^ s.charCodeAt(i), 0x01000193) >>> 0;
  return h >>> 0;
}

// open-addressed table: entryHash -> repHash. 42M entries at a 0.33 load
// factor; linear probing degrades badly past ~0.7, and a table that fills
// completely makes put() spin forever rather than fail, so leave real headroom.
const CAP = 1 << 27;                       // 134M slots for ~42M entries
const MASK = CAP - 1;
const keyTab = new Uint32Array(CAP);
const repTab = new Uint32Array(CAP);
let nClustered = 0;

function put(entryHash, repHash) {
  let i = entryHash & MASK;
  for (;;) {
    const k = keyTab[i];
    if (k === 0) { keyTab[i] = entryHash || 1; repTab[i] = repHash; nClustered++; return; }
    if (k === (entryHash || 1)) return;
    i = (i + 1) & MASK;
  }
}
function get(entryHash) {
  let i = entryHash & MASK;
  for (;;) {
    const k = keyTab[i];
    if (k === 0) return 0;
    if (k === (entryHash || 1)) return repTab[i];
    i = (i + 1) & MASK;
  }
}

// the 2.6M genuine multi-member representatives, as a hash set
const repsPath = process.env.REPS_TSV;
const multiReps = new Set();
if (repsPath) {
  const { gunzipSync } = await import('node:zlib');
  const text = gunzipSync(readFileSync(repsPath)).toString('utf8');
  for (const line of text.split('\n')) {
    if (!line) continue;
    const t = line.indexOf('\t');
    multiReps.add(fnv1a(t < 0 ? line : line.slice(0, t)));
  }
  console.log(`  ${multiReps.size} multi-member cluster representatives loaded`);
}

{
  const rl = createInterface({ input: createReadStream(clustersPath).pipe(createGunzip()), crlfDelay: Infinity });
  let n = 0, skipped = 0;
  for await (const line of rl) {
    if (!line) continue;
    // repId \t entryId \t cluFlag \t taxId
    const t1 = line.indexOf('\t');
    if (t1 < 0) continue;
    const t2 = line.indexOf('\t', t1 + 1);
    const rep = line.slice(0, t1);
    const entry = t2 < 0 ? line.slice(t1 + 1) : line.slice(t1 + 1, t2);
    const rh = fnv1a(rep);
    // single-member clusters are written as "X <tab> X"; skip them and anything
    // whose representative is not a known multi-member cluster
    if (rep === entry || (multiReps.size && !multiReps.has(rh))) { skipped++; continue; }
    put(fnv1a(entry), rh);
    if (nClustered > CAP * 0.6) { console.error('cluster table over 60% full -- raise CAP'); process.exit(1); }
    if (++n % 25_000_000 === 0) console.log(`  ${n / 1e6}M rows, ${nClustered} clustered (${el()})`);
  }
  console.log(`  ${n} rows scanned, ${skipped} singleton/unknown, ${nClustered} entries in multi-member clusters (${el()})`);
}

/*
 * Pass 2 in two sweeps, because the obvious way does not fit.
 *
 * Collecting every candidate as a [depth, accession] pair and sorting would mean
 * ~157M JS arrays -- tens of gigabytes for what is ultimately a threshold. Depth
 * is a small integer, so sweep once to histogram it, solve for the depth cutoff
 * that yields MAX entries, then sweep again keeping only what clears it. Two
 * passes over a 3.8 GB file beats one pass over a 30 GB heap.
 */
const MAXD = 1 << 20;
const hist = new Float64Array(MAXD + 1);
const bestDepthPerRep = new Map();          // repHash -> depth (accession comes later)

console.log(`\nreading depths from ${depthsPath} (pass 1: histogram)…`);
{
  const rl = createInterface({ input: createReadStream(depthsPath), crlfDelay: Infinity });
  let n = 0, low = 0;
  for await (const line of rl) {
    if (!line) continue;
    const c = line.indexOf(',');
    if (c < 0) continue;
    const depth = +line.slice(c + 1);
    if (!(depth >= MIN_DEPTH)) { low++; continue; }
    const rep = get(fnv1a(line.slice(0, c)));
    if (rep) {
      const cur = bestDepthPerRep.get(rep);
      if (cur === undefined || depth > cur) bestDepthPerRep.set(rep, depth);
    } else {
      hist[Math.min(depth, MAXD)]++;        // singleton
    }
    if (++n % 50_000_000 === 0) console.log(`  ${n / 1e6}M rows (${el()})`);
  }
  console.log(`  ${n} rows above depth ${MIN_DEPTH}, ${low} below (${el()})`);
}

// cluster picks are always kept; singletons fill the rest of the budget
for (const d of bestDepthPerRep.values()) hist[Math.min(d, MAXD)]--;   // not singletons
const nClusters = bestDepthPerRep.size;
const budget = Math.max(0, MAX - nClusters);
let cutoff = MAXD, running = 0;
while (cutoff > 0 && running + hist[cutoff] <= budget) { running += hist[cutoff]; cutoff--; }
console.log(`\n  multi-member cluster picks : ${nClusters}`);
console.log(`  singleton budget           : ${budget}`);
console.log(`  singleton depth cutoff     : ${cutoff} (yields ~${running})`);

console.log(`\n(pass 2: collect)…`);
const out = [];
{
  const rl = createInterface({ input: createReadStream(depthsPath), crlfDelay: Infinity });
  for await (const line of rl) {
    if (!line) continue;
    const c = line.indexOf(',');
    if (c < 0) continue;
    const depth = +line.slice(c + 1);
    if (!(depth >= MIN_DEPTH)) continue;
    const acc = line.slice(0, c);
    const rep = get(fnv1a(acc));
    if (rep) {
      if (bestDepthPerRep.get(rep) === depth) { out.push([depth, acc]); bestDepthPerRep.delete(rep); }
    } else if (depth > cutoff) {
      out.push([depth, acc]);
    }
  }
}
out.sort((a, b) => b[0] - a[0]);
const kept = out.slice(0, MAX);
writeFileSync(outPath, kept.map(x => x[1]).join('\n') + '\n');

const depths = kept.map(x => x[0]);
console.log(`\n  selected ${kept.length} entries -> ${outPath}`);
console.log(`  depth: median ${depths[depths.length >> 1]}, min ${depths[depths.length - 1]}, max ${depths[0]}`);
console.log(`  (structural representatives had median depth 200 for comparison)`);
console.log(`done in ${el()}`);
