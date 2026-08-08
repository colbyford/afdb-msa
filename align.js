/*
 * align.js -- Smith-Waterman local alignment with affine gaps, BLOSUM62.
 *
 * The index replaces BLAST's *search*, but the a3m transfer needs BLAST's other
 * output: the aligned query/hit strings that say which hit column maps to which
 * query residue. This produces exactly that, in the same shape as the HSP object
 * parsed out of the EBI JSON, so msa.js cannot tell the two apart.
 *
 * One query against a handful of candidates is a few hundred thousand cells --
 * microseconds. Nothing here needs to be clever.
 */

const B62_ORDER = 'ARNDCQEGHILKMFPSTWYVBZX*';

// BLOSUM62, row-major over B62_ORDER.
const B62_RAW = [
   4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4,
  -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4,
  -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4,
  -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4,
   0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4,
  -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4,
  -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4,
   0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4,
  -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4,
  -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4,
  -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4,
  -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4,
  -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4,
  -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4,
  -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4,
   1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4,
   0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4,
  -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4,
  -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4,
   0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4,
  -2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4,
  -1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4,
   0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4,
  -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1,
];

const NA = B62_ORDER.length;
const B62_IDX = new Int8Array(128).fill(B62_ORDER.indexOf('X'));
for (let i = 0; i < NA; i++) B62_IDX[B62_ORDER.charCodeAt(i)] = i;

const GAP_OPEN = 11;   // BLAST protein defaults
const GAP_EXT = 1;

/*
 * Local alignment. Returns null when nothing scores above zero, otherwise an
 * HSP shaped exactly like Api.parseBlastHits produces:
 *   { qseq, hseq, qFrom, hFrom, qTo, hTo, score, identity, alignLen }
 * qFrom/hFrom are 1-based, matching BLAST's convention (and hspToMap's).
 */
function smithWaterman(query, target, opts) {
  opts = opts || {};
  const go = opts.gapOpen === undefined ? GAP_OPEN : opts.gapOpen;
  const ge = opts.gapExtend === undefined ? GAP_EXT : opts.gapExtend;

  const n = query.length, m = target.length;
  if (!n || !m) return null;

  const qi = new Int8Array(n), ti = new Int8Array(m);
  for (let i = 0; i < n; i++) qi[i] = B62_IDX[query.charCodeAt(i) & 0x7f];
  for (let j = 0; j < m; j++) ti[j] = B62_IDX[target.charCodeAt(j) & 0x7f];

  // One row at a time; traceback comes from a compact direction matrix.
  // 0 = stop, 1 = diagonal, 2 = up (gap in target), 3 = left (gap in query)
  //
  // The two affine states recurse along different axes, which decides their
  // storage. Getting this backwards still produces plausible-looking alignments
  // with subtly wrong endpoints, so:
  //   E[i][j] = max(H[i][j-1] - open, E[i][j-1] - extend)   gap in query, along j
  //             -> a scalar carried across the row
  //   F[i][j] = max(H[i-1][j] - open, F[i-1][j] - extend)   gap in target, along i
  //             -> one slot per column, carried across rows
  const dir = new Uint8Array((n + 1) * (m + 1));
  const H = new Int32Array(m + 1);
  const F = new Int32Array(m + 1).fill(-1e8);
  let best = 0, bestI = 0, bestJ = 0;

  for (let i = 1; i <= n; i++) {
    let diagPrev = 0;      // H[i-1][j-1]
    let E = -1e8;          // E[i][j-1]
    let Hleft = 0;         // H[i][j-1]
    const rowOff = i * (m + 1);
    const qrow = qi[i - 1] * NA;
    for (let j = 1; j <= m; j++) {
      const up = H[j];     // H[i-1][j], before this row overwrites it

      const fj = (up - go) > (F[j] - ge) ? (up - go) : (F[j] - ge);
      F[j] = fj;
      E = (Hleft - go) > (E - ge) ? (Hleft - go) : (E - ge);

      const diag = diagPrev + B62_RAW[qrow + ti[j - 1]];
      let h = 0, d = 0;
      if (diag > h) { h = diag; d = 1; }
      if (fj > h) { h = fj; d = 2; }
      if (E > h) { h = E; d = 3; }

      dir[rowOff + j] = d;
      diagPrev = up;
      Hleft = h;
      H[j] = h;
      if (h > best) { best = h; bestI = i; bestJ = j; }
    }
  }

  if (!best) return null;

  let i = bestI, j = bestJ;
  const qs = [], hs = [];
  for (;;) {
    const d = dir[i * (m + 1) + j];
    if (!d) break;
    if (d === 1) { qs.push(query[i - 1]); hs.push(target[j - 1]); i--; j--; }
    else if (d === 2) { qs.push(query[i - 1]); hs.push('-'); i--; }
    else { qs.push('-'); hs.push(target[j - 1]); j--; }
    if (i === 0 || j === 0) break;
  }
  qs.reverse(); hs.reverse();
  const qseq = qs.join(''), hseq = hs.join('');

  let same = 0;
  for (let k = 0; k < qseq.length; k++) if (qseq[k] === hseq[k] && qseq[k] !== '-') same++;

  return {
    qseq, hseq,
    qFrom: i + 1, hFrom: j + 1,
    qTo: bestI, hTo: bestJ,
    score: best,
    alignLen: qseq.length,
    identity: qseq.length ? (100 * same) / qseq.length : 0,
    queryCoverage: (bestI - i) / n,
  };
}

const Align = { smithWaterman, B62_ORDER, B62_RAW, B62_IDX, GAP_OPEN, GAP_EXT };

if (typeof module !== 'undefined' && module.exports) module.exports = Align;
if (typeof window !== 'undefined') window.Align = Align;
