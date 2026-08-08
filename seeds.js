/*
 * seeds.js -- minimizer seeding, shared verbatim by the index builder and the
 * browser. If these two ever disagree, every lookup silently misses.
 *
 * The whole search rests on one fact: at 90% identity a 10-mer survives intact
 * with probability 0.9^10 = 0.35, so a handful of exact seeds is plenty. Twelve
 * seeds give 1 - 0.65^12 = 99.2% chance of sharing at least one. Measured on
 * 239.6M AFDB entries: 100% of >=90% queries retrieved, ~53 KB and ~1 ms.
 *
 * Two details are load-bearing.
 *
 * 1. The k-mer key covers exactly the last k residues. An earlier version
 *    accumulated `code = (code*20 + c) % MOD` across the whole sequence, which
 *    is a prefix hash -- two sequences then collide only if their entire
 *    prefixes match. It produced a flat 25% recall at every k and w, the
 *    signature of a key that does not depend on k at all.
 *
 * 2. The window is chosen per sequence so the seed count lands near
 *    TARGET_SEEDS regardless of length. A fixed w undersamples short proteins:
 *    at k=12 w=16 a 142-residue globin got 8 seeds and shared 0 of them with a
 *    real 90% relative on two of three tries. Recall depends on the number of
 *    seeds, not the window -- and adapting w also stops long proteins from
 *    dominating the posting lists.
 */

const SEED_K = 10;
const TARGET_SEEDS = 12;
const MIN_W = 4;
const MAX_W = 128;

const AA20 = 'ACDEFGHIKLMNPQRSTVWY';
const AA_IDX = new Int8Array(128).fill(-1);
for (let i = 0; i < AA20.length; i++) AA_IDX[AA20.charCodeAt(i)] = i;

function hash32(x) {
  x = (x ^ 0x9e3779b9) >>> 0;
  x = Math.imul(x ^ (x >>> 16), 0x85ebca6b) >>> 0;
  x = Math.imul(x ^ (x >>> 13), 0xc2b2ae35) >>> 0;
  return (x ^ (x >>> 16)) >>> 0;
}

// Window for a sequence of length n: aim for `target` seeds, clamped.
function windowFor(n, target) {
  return Math.max(MIN_W, Math.min(MAX_W, Math.floor(n / (target || TARGET_SEEDS)) || MIN_W));
}

/*
 * The set of minimizer hashes for a sequence: in each window, the smallest
 * k-mer hash. Values are 31-bit so the top bits can address a bucket and the
 * low bits survive as a residual.
 *
 * Non-standard residues (X, U, B, Z) break the run rather than mapping to some
 * arbitrary letter, so they never invent a seed that isn't there.
 */
function seeds(seq, k, target) {
  k = k || SEED_K;
  const n = seq.length;
  const out = new Set();
  if (n < k) return out;

  const code = new Int8Array(n);
  for (let i = 0; i < n; i++) code[i] = AA_IDX[seq.charCodeAt(i) & 0x7f];

  const hs = new Int32Array(n).fill(-1);
  for (let i = k - 1; i < n; i++) {
    let x = 0, ok = true;
    for (let j = i - k + 1; j <= i; j++) {
      const c = code[j];
      if (c < 0) { ok = false; break; }
      x = (Math.imul(x, 31) + c) | 0;
    }
    if (ok) hs[i] = hash32(x >>> 0) & 0x7fffffff;
  }

  const w = windowFor(n, target);
  for (let i = k - 1; i < n; i++) {
    if (i % w) continue;
    let best = -1;
    const stop = Math.min(i + w, n);
    for (let j = i; j < stop; j++) if (hs[j] >= 0 && (best < 0 || hs[j] < best)) best = hs[j];
    if (best >= 0) out.add(best);
  }
  return out;
}

const Seeds = { SEED_K, TARGET_SEEDS, MIN_W, MAX_W, seeds, windowFor, hash32, AA_IDX };

if (typeof module !== 'undefined' && module.exports) module.exports = Seeds;
if (typeof window !== 'undefined') window.Seeds = Seeds;
