/*
 * filter.js -- MSA filtering, following the semantics in GREMLIN-GUI's msa.js
 * (sokrypton/GREMLIN-GUI), which in turn follow py2Dmol's thresholds.
 *
 * Definitions kept deliberately identical to that code so numbers are comparable:
 *   isGap      : '-', '.', ' ' and 'X' all count as gaps
 *   coverage   : fraction of columns holding a real residue
 *   identity   : fraction of ALL columns where both rows hold the same real
 *                residue (denominator is the full length, not just aligned pairs)
 *
 * Everything here operates on match-state strings; a3m insertions (lowercase)
 * are stripped for measurement and restored on output, so filtering never
 * changes the alignment itself.
 */

function isGap(ch) {
  return !ch || ch === '-' || ch === '.' || ch === ' ' || ch === 'X';
}

// Drop a3m insertion characters, leaving one character per match column.
function matchStates(seq) {
  let out = '';
  for (let i = 0; i < seq.length; i++) {
    const c = seq.charCodeAt(i);
    if (c < 97 || c > 122) out += seq[i]; // not a-z
  }
  return out;
}

function coverageOf(seq) {
  if (!seq || !seq.length) return 0;
  let n = 0;
  for (let i = 0; i < seq.length; i++) if (!isGap(seq[i])) n++;
  return n / seq.length;
}

function identityTo(seq, query) {
  if (!seq || !query || seq.length !== query.length) return 0;
  let m = 0;
  const t = seq.length;
  for (let i = 0; i < t; i++) {
    const a = seq[i], b = query[i];
    if (!isGap(a) && !isGap(b) && a === b) m++;
  }
  return t > 0 ? m / t : 0;
}

/*
 * Greedy redundancy clustering, as in GREMLIN-GUI's filterRedundancy: walk the
 * rows in order, keep each as a representative, and drop any later row more than
 * `maxId` identical to a representative already kept. Row 0 (the query) is
 * always a representative. Identity counts equal characters over all L columns,
 * gap-vs-gap included -- the same definition Meff reweighting uses.
 *
 * Plain O(N*R) with an early exit; the byte-packed version in GREMLIN-GUI is
 * worth it there because N reaches 100k, here N is bounded by maxDepth.
 */
function filterRedundancy(ms, maxId) {
  const L = ms[0].length;
  const need = Math.ceil(maxId * L);
  const keep = [0];
  const reps = [ms[0]];
  for (let n = 1; n < ms.length; n++) {
    const s = ms[n];
    let redundant = false;
    for (let r = 0; r < reps.length; r++) {
      const t = reps[r];
      let id = 0;
      let ok = true;
      for (let k = 0; k < L; k++) {
        if (s[k] === t[k]) id++;
        else if (id + (L - k - 1) < need) { ok = false; break; }
      }
      if (ok && id >= need) { redundant = true; break; }
    }
    if (!redundant) { keep.push(n); reps.push(s); }
  }
  return keep;
}

/*
 * rows  : [{name, seq}] in a3m form, row 0 = query
 * opts  : { minCoverage, minIdentity, maxIdentity, sortByIdentity, maxDepth }
 *
 * Returns { rows, warnings, stats } with the query always retained as row 0.
 */
function filterRows(rows, opts) {
  opts = opts || {};
  const minCoverage = opts.minCoverage === undefined ? 0 : opts.minCoverage;
  const minIdentity = opts.minIdentity === undefined ? 0 : opts.minIdentity;
  const maxIdentity = opts.maxIdentity === undefined ? 1 : opts.maxIdentity;
  const maxDepth = opts.maxDepth || 0;
  const warnings = [];

  if (!rows.length) return { rows: [], warnings, stats: {} };

  const ms = rows.map(r => matchStates(r.seq));
  const query = ms[0];

  const scored = [];
  let nDropCov = 0, nDropId = 0;
  for (let i = 1; i < rows.length; i++) {
    const cov = coverageOf(ms[i]);
    const idn = identityTo(ms[i], query);
    if (cov < minCoverage) { nDropCov++; continue; }
    if (idn < minIdentity) { nDropId++; continue; }
    scored.push({ i, cov, idn });
  }
  if (nDropCov) warnings.push(`Dropped ${nDropCov} sequence(s) below ${minCoverage.toFixed(2)} coverage.`);
  if (nDropId) warnings.push(`Dropped ${nDropId} sequence(s) below ${minIdentity.toFixed(2)} identity to the query.`);

  /*
   * Redundancy clustering is greedy and order-dependent: it keeps whichever
   * member of a cluster it meets first. So always run it over identity-sorted
   * rows, which makes the survivor the one closest to the query -- otherwise the
   * survivor is an accident of input order. Measured on a 600-row AFDB MSA, the
   * two orders keep different members 13 times out of ~500.
   *
   * `sortByIdentity` then governs only how the *output* is ordered, which is a
   * separate question from which rows survive.
   */
  const byIdentity = [...scored].sort((a, b) => b.idn - a.idn);
  let order = [0, ...byIdentity.map(s => s.i)];

  if (maxIdentity > 0 && maxIdentity < 1 && order.length > 1) {
    const sub = order.map(i => ms[i]);
    const keep = filterRedundancy(sub, maxIdentity);
    const before = order.length;
    order = keep.map(k => order[k]);
    if (before !== order.length) {
      warnings.push(`Removed ${before - order.length} sequence(s) above ${maxIdentity.toFixed(2)} identity to a kept sequence; ${order.length} remain.`);
    }
  }

  if (maxDepth && order.length > maxDepth) {
    warnings.push(`Truncated to the top ${maxDepth} of ${order.length} sequences.`);
    order = order.slice(0, maxDepth);
  }

  // Restore input order unless the caller asked for identity order. The query
  // stays first either way.
  if (!opts.sortByIdentity) {
    const keep = new Set(order);
    order = [0, ...scored.map(s => s.i).filter(i => keep.has(i))];
  }

  const byIdx = new Map(scored.map(s => [s.i, s]));
  const out = order.map(i => {
    const r = rows[i];
    const s = byIdx.get(i);
    return { name: r.name, seq: r.seq, cov: s ? s.cov : 1, idn: s ? s.idn : 1 };
  });

  return {
    rows: out,
    warnings,
    stats: { input: rows.length, kept: out.length, droppedCoverage: nDropCov, droppedIdentity: nDropId },
  };
}

/*
 * Meff / Neff: number of sequence clusters at `thresh` identity, the standard
 * measure of how much independent information an MSA actually carries.
 * Same clustering as filterRedundancy, reported rather than applied.
 */
function neff(rows, thresh) {
  thresh = thresh === undefined ? 0.8 : thresh;
  const ms = rows.map(r => matchStates(r.seq));
  if (!ms.length) return 0;
  return filterRedundancy(ms, thresh).length;
}

const Filter = { isGap, matchStates, coverageOf, identityTo, filterRedundancy, filterRows, neff };

if (typeof module !== 'undefined' && module.exports) module.exports = Filter;
if (typeof window !== 'undefined') window.Filter = Filter;
