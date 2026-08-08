/*
 * search.js -- client-side lookup over the static AFDB minimizer index.
 *
 * The index is ~17 GB and a query touches about 50 KB of it. Nothing is
 * downloaded up front, not even the bucket table:
 *
 *   1. compute the query's ~12 minimizer seeds
 *   2. one 8-byte ranged GET into buckets.u32 per seed gives its posting range
 *      (buckets is a prefix-sum array, so slots b and b+1 are adjacent and
 *      arrive in a single read)
 *   3. one ranged GET of that posting range per seed
 *   4. count how many seeds each entry shares; the top few are the candidates
 *
 * ~24 small ranged requests and ~50 KB, independent of index size. An earlier
 * draft had the browser download all 134 MB of buckets.u32 once; range-reading
 * it costs 96 bytes instead and removes the cold start entirely.
 *
 * Seeds are 31-bit. The top `bucketBits` address a bucket; the stored 1-byte
 * residual is the seed's low 8 bits. That overlaps the bucket by 2 bits -- a
 * harmless redundancy which still pins down the 6 low bits the bucket does not
 * carry. Residual collisions are filtered by the alignment that follows, since
 * shared-seed count is only a prefilter.
 */

const DEFAULT_TOPK = 20;
// Skip repeat / low-complexity seeds rather than pull megabytes for one lookup.
const MAX_POSTINGS = 20000;

class MinimizerIndex {
  constructor(baseUrl, opts) {
    this.base = baseUrl.replace(/\/$/, '') + '/';
    this.opts = opts || {};
    this.meta = null;
    this._loading = null;
  }

  async load() {
    if (this.meta) return this;
    if (this._loading) return this._loading;
    this._loading = (async () => {
      const r = await fetch(this.base + 'meta.json');
      if (!r.ok) throw new Error(`index meta.json: HTTP ${r.status}`);
      const m = await r.json();
      // A mismatch here means the page and the published index disagree about
      // how seeds are computed; every lookup would silently miss.
      if (m.k !== Seeds.SEED_K || m.targetSeeds !== Seeds.TARGET_SEEDS) {
        throw new Error(`index/seeder mismatch: index k=${m.k} targetSeeds=${m.targetSeeds}, `
          + `page k=${Seeds.SEED_K} targetSeeds=${Seeds.TARGET_SEEDS}`);
      }
      this.meta = m;
      return this;
    })();
    return this._loading;
  }

  _range(url, from, len, signal) {
    return fetch(this.base + url, { headers: { Range: `bytes=${from}-${from + len - 1}` }, signal })
      .then(r => {
        if (r.status !== 206 && r.status !== 200) throw new Error(`range read ${url}: HTTP ${r.status}`);
        return r.arrayBuffer();
      });
  }

  _fileFor(bucket) {
    for (const f of this.meta.files) {
      if (bucket >= f.fromBucket && bucket < f.toBucket) return f;
    }
    return null;
  }

  /* One seed -> the entry ids carrying it. Two ranged reads. */
  async _postingsFor(seed, signal) {
    const b = seed >>> this.meta.shift;
    const head = new Uint32Array(await this._range('buckets.u32', b * 4, 8, signal));
    const lo = head[0], hi = head[1];
    const count = hi - lo;
    if (!count) return { ids: [], bytes: 8 };
    if (count > MAX_POSTINGS) return { ids: [], bytes: 8, skipped: true };

    const f = this._fileFor(b);
    if (!f) return { ids: [], bytes: 8 };
    const PB = this.meta.postBytes;
    const buf = new Uint8Array(await this._range(f.name, (lo - f.postFrom) * PB, count * PB, signal));

    const residual = seed & 0xff;
    const ids = [];
    for (let i = 0; i < count; i++) {
      const o = i * PB;
      if (buf[o] !== residual) continue;
      ids.push((buf[o + 1] | (buf[o + 2] << 8) | (buf[o + 3] << 16) | (buf[o + 4] << 24)) >>> 0);
    }
    return { ids, bytes: 8 + count * PB };
  }

  /* Accessions are fixed-width records, so an entry id IS its byte offset. */
  async accessionsFor(ids, signal) {
    const w = this.meta.accBytes || 12;
    return Promise.all(ids.map(id =>
      this._range('acc.bin', id * w, w, signal)
        .then(b => {
          let s = '';
          const u = new Uint8Array(b);
          for (let i = 0; i < u.length && u[i]; i++) s += String.fromCharCode(u[i]);
          return s;
        })
        .catch(() => null)));
  }

  /*
   * Returns [{acc, id, shared}] best-first, where `shared` counts how many of
   * the query's seeds that entry carries. This ranking is a prefilter only --
   * the caller aligns the top few, and alignment decides the actual best hit.
   */
  async search(seq, opts) {
    await this.load();
    opts = opts || {};
    const topK = opts.topK || this.opts.topK || DEFAULT_TOPK;
    const signal = opts.signal;

    const seedList = [...Seeds.seeds(seq, this.meta.k, this.meta.targetSeeds)];
    if (!seedList.length) return Object.assign([], { bytesRead: 0, seeds: 0, candidates: 0 });

    // Seeds are independent; issuing them together makes the lookup ~2 RTT.
    const lists = await Promise.all(seedList.map(s => this._postingsFor(s, signal)));

    const counts = new Map();
    let bytes = 0, found = 0, skipped = 0;
    for (const l of lists) {
      bytes += l.bytes;
      if (l.skipped) { skipped++; continue; }
      if (l.ids.length) found++;
      for (const id of l.ids) counts.set(id, (counts.get(id) || 0) + 1);
    }

    const top = [...counts.entries()].sort((a, b) => b[1] - a[1]).slice(0, topK);
    const accs = await this.accessionsFor(top.map(t => t[0]), signal);

    const out = top.map(([id, shared], i) => ({ acc: accs[i], id, shared })).filter(x => x.acc);
    out.bytesRead = bytes;
    out.seeds = seedList.length;
    out.seedsFound = found;
    out.seedsSkipped = skipped;
    out.candidates = counts.size;
    return out;
  }
}

const Search = { MinimizerIndex, DEFAULT_TOPK, MAX_POSTINGS };

if (typeof module !== 'undefined' && module.exports) {
  if (typeof Seeds === 'undefined') global.Seeds = require('./seeds.js');
  module.exports = Search;
}
if (typeof window !== 'undefined') window.Search = Search;
