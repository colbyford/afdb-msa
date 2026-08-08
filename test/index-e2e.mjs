/*
 * node test/index-e2e.mjs <indexDir>
 *
 * The whole path against the real minimizer index, over HTTP with real Range
 * requests, so the client code under test is exactly what ships:
 *
 *   sequence -> seeds -> postings -> accessions -> Smith-Waterman
 *            -> AlphaFold DB a3m -> transfer -> filtered a3m
 *
 * A local server stands in for Hugging Face / R2 and answers Range the same way
 * (206 + Content-Range). AFDB is hit for real.
 */
import { createRequire } from 'node:module';
import { createServer } from 'node:http';
import { openSync, readSync, closeSync, statSync, existsSync } from 'node:fs';
import { join } from 'node:path';
const require = createRequire('file:///home/ubuntu/msa_afdb/x.js');
const R = '/home/ubuntu/msa_afdb';
globalThis.Seeds = require(R + '/seeds.js');
const A = require(R + '/align.js'), K = require(R + '/msa.js'), F = require(R + '/filter.js'), Api = require(R + '/api.js');
const { MinimizerIndex } = require(R + '/search.js');
const DIR = process.argv[2];

let reqs = 0, served = 0;
const fds = new Map();
const srv = createServer((q, res) => {
  const p = join(DIR, decodeURIComponent(q.url.replace(/^\//, '')));
  if (!existsSync(p)) { res.writeHead(404); res.end(); return; }
  const size = statSync(p).size;
  const m = /bytes=(\d+)-(\d+)?/.exec(q.headers.range || '');
  if (!fds.has(p)) fds.set(p, openSync(p, 'r'));
  if (m) {
    const from = +m[1], to = m[2] ? Math.min(+m[2], size - 1) : size - 1;
    const len = to - from + 1;
    const b = Buffer.allocUnsafe(len);
    readSync(fds.get(p), b, 0, len, from);
    reqs++; served += len;
    res.writeHead(206, { 'Content-Range': `bytes ${from}-${to}/${size}`, 'Content-Length': len, 'Accept-Ranges': 'bytes' });
    res.end(b);
  } else {
    const b = Buffer.allocUnsafe(size); readSync(fds.get(p), b, 0, size, 0);
    reqs++; served += size;
    res.writeHead(200, { 'Content-Length': size }); res.end(b);
  }
});
await new Promise(r => srv.listen(0, '127.0.0.1', r));
const base = `http://127.0.0.1:${srv.address().port}/`;

const QUERIES = {
  'HBA human (exact)': 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR',
  'HBA variant (~85%)': 'MVLSAADKSNVKATWDKIGSHAGDYGGEALDRTFQSFPTTKTYFPHFDLSHGSAQVKAHGKKVAAALVEAVNHIDDIAGALSKLSDLHAQKLRVDPVNFKLLGQCFLVVVAIHHPSALTPEVHASLDKFLCAVGNVLTAKYR',
  'EGFR human': null,
};
let pass = 0, fail = 0;
const ok = (c, m) => { c ? (pass++, console.log('  ok   ' + m)) : (fail++, console.log('  FAIL ' + m)); };
try {
  const idx = new MinimizerIndex(base);
  await idx.load();
  console.log(`index: ${idx.meta.nSeq.toLocaleString()} entries, ${idx.meta.files.length} posting files\n`);
  for (const [name, q] of Object.entries(QUERIES)) {
    if (!q) continue;
    console.log(`# ${name} (${q.length} aa)`);
    const r0 = reqs, s0 = served, t0 = Date.now();
    const hits = await idx.search(q, { topK: 20 });
    const dt = Date.now() - t0;
    ok(hits.length > 0, `${hits.length} candidates from ${hits.candidates} scored, ${hits.seedsFound}/${hits.seeds} seeds found`);
    console.log(`       ${reqs - r0} ranged requests, ${((served - s0) / 1024).toFixed(0)} KB, ${dt} ms`);
    console.log(`       top: ${hits.slice(0, 4).map(h => h.acc + '(' + h.shared + ')').join(' ')}`);

    const seqs = await Api.candidateSequences(hits.map(h => h.acc), null, n => console.log('       ' + n));
    ok(seqs.size > 0, `resolved ${seqs.size}/${hits.length} candidate sequences`);
    let best = null;
    for (const h of hits) {
      const s = seqs.get(h.acc); if (!s) continue;
      const aln = A.smithWaterman(q, s);
      if (aln && (!best || aln.score > best.aln.score)) best = { acc: h.acc, aln, len: s.length };
    }
    ok(!!best, best ? `best ${best.acc} at ${best.aln.identity.toFixed(1)}% identity, ${(best.aln.queryCoverage * 100).toFixed(0)}% coverage` : 'no alignment');
    // Only queries that actually have a >=90% relative can clear the bar. The
    // ~85% variant is built not to, and must still yield a usable alignment
    // rather than nothing -- that is the graceful-degradation case.
    const wantHigh = name.includes('exact');
    if (wantHigh) ok(best && best.aln.identity >= 90, `clears the 90% bar (${best ? best.aln.identity.toFixed(1) : 0}%)`);
    else ok(best && best.aln.identity >= 80, `no >=90% relative exists; still found ${best ? best.aln.identity.toFixed(1) : 0}%`);

    let got = null;
    for (const h of hits.slice(0, 6)) {
      try { got = { acc: h.acc, text: await Api.afdbFetchMsaCached(h.acc) }; break; }
      catch (e) { if (e.name !== 'NoMsaError') throw e; }
    }
    ok(!!got, got ? `AFDB MSA for ${got.acc}: ${(got.text.length / 1e6).toFixed(2)} MB` : 'no MSA');
    if (!got) continue;
    const dseq = seqs.get(got.acc);
    const aln = A.smithWaterman(q, dseq);
    const built = K.buildChainA3M({ name: 'query', seq: q }, [{
      acc: got.acc, a3mText: got.text, hitLen: dseq.length,
      hsp: { qseq: aln.qseq, hseq: aln.hseq, qFrom: aln.qFrom, hFrom: aln.hFrom },
      desc: '', organism: '', taxid: '',
    }]);
    const bad = built.rows.filter(x => K.countMatchStates(x.seq) !== q.length).length;
    ok(bad === 0, `transferred ${built.rows.length} sequences, all with ${q.length} match states`);
    const filt = F.filterRows(built.rows, { minCoverage: 0.5, maxIdentity: 0.9, sortByIdentity: true });
    ok(filt.rows.length > 50, `${filt.rows.length} after filtering, Neff@0.8 = ${F.neff(filt.rows, 0.8)}`);
    const a3m = K.formatA3M(filt.rows);
    ok(a3m.startsWith('>query\n') && a3m.split('\n')[1] === q, 'well-formed a3m, query first');
    console.log(`       -> ${(a3m.length / 1e6).toFixed(2)} MB a3m, ${filt.rows.length} sequences\n`);
  }
} finally { srv.close(); for (const fd of fds.values()) closeSync(fd); }
console.log(`total: ${reqs} ranged requests, ${(served / 1024).toFixed(0)} KB from a 17 GB index`);
console.log(`${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
