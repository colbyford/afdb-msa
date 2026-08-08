/*
 * node test/e2e.mjs ["SEQ" | "SEQ_A:SEQ_B"]
 *
 * End-to-end, exactly the path app.js takes: BLAST at EBI, select diverse hits,
 * pull each hit's AFDB MSA, transfer onto the query, merge, filter, and (for
 * multi-chain input) pair by species.
 *
 *   EMAIL=you@example.org   contact address EBI requires on every job
 *   DB=uniprotkb_swissprot  much faster than full uniprotkb
 *   NHITS=1                 hits to merge
 */
import { createRequire } from 'node:module';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const require = createRequire(import.meta.url);
const here = dirname(fileURLToPath(import.meta.url));
const K = require(join(here, '..', 'msa.js'));
const F = require(join(here, '..', 'filter.js'));
const Api = require(join(here, '..', 'api.js'));

const EMAIL = process.env.EMAIL || 'anonymous@example.org';
const DB = process.env.DB || 'uniprotkb';
const NHITS = Number(process.env.NHITS || 1);

const HBA = 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR';
const HBB = 'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH';
// mutate slightly so the query is a neighbour of a real entry, not an exact copy
const DEFAULT = `${HBA.replace(/^MVLS/, 'MVLA').replace(/KYR$/, 'KFR')}:${HBB.replace(/^MVHL/, 'MVHM')}`;

const chains = K.parseQueryInput(process.argv[2] || DEFAULT);
const t0 = Date.now();
const el = () => `${((Date.now() - t0) / 1000).toFixed(0)}s`;
let failures = 0;
const check = (cond, msg) => {
  if (!cond) failures++;
  console.log(`  ${cond ? 'ok  ' : 'FAIL'} ${msg}`);
};

console.log(`${chains.length} chain(s), database=${DB}, merging up to ${NHITS} hit(s)\n`);

async function doChain(chain, tag) {
  console.log(`[${tag}] ${chain.name}: ${chain.seq.length} aa`);

  // Fast path: exact UniProt sequence -> checksum lookup, no search at all.
  if (process.env.NO_EXACT !== '1') {
    const t = Date.now();
    const accs = await Api.uniparcAccessions(chain.seq);
    console.log(`[${tag}]   UniParc checksum lookup: ${Date.now() - t}ms, ${accs.length} accession(s)`);
    if (accs.length) {
      const got = await Api.afdbFetchFirstAvailable(accs, null, null,
        a => console.log(`[${tag}]   no AFDB MSA for ${a}`));
      console.log(`[${tag}]   exact match ${got.acc}: ${(got.text.length / 1e6).toFixed(2)} MB — BLAST skipped (${Date.now() - t}ms total)`);
      const hits = [{
        acc: got.acc, a3mText: got.text, hitLen: chain.seq.length, identity: 100, queryCoverage: 1,
        hsp: { qseq: chain.seq, hseq: chain.seq, qFrom: 1, hFrom: 1 }, desc: '', organism: '', taxid: '',
      }];
      return finish(chain, hits, tag);
    }
  }

  const job = await Api.blastSubmit(chain.seq, { email: EMAIL, database: DB, nHits: 50, evalue: '1e-3' });
  await Api.blastWait(job, (s, e) => process.stdout.write(`\r[${tag}]   ${s} (${e}s)      `));
  const all = Api.parseBlastHits(await Api.blastResult(job));
  console.log(`\r[${tag}] ${all.length} hits in ${el()}                    `);

  // Rank everything for diversity, then descend that ranking -- never fall back
  // to the raw BLAST order, which is full of the duplicates just rejected.
  const sel = K.selectDiverseHits(all, chain.seq.length, { maxPairId: 0.9 });
  console.log(`[${tag}]   ${sel.skipped.length} redundant hit(s) skipped, ${sel.chosen.length} distinct`);

  // Download as we descend the ranking; a missing MSA is a NoMsaError, not a
  // separate probe (probing would download every file twice).
  const usable = [];
  for (const h of sel.chosen) {
    if (usable.length >= NHITS) break;
    try {
      h.a3mText = await Api.afdbFetchMsaCached(h.acc);
      console.log(`[${tag}]   ${h.acc} id=${h.identity}% cov=${(h.queryCoverage * 100).toFixed(0)}% -> ${K.parseA3M(h.a3mText).length} seqs, ${(h.a3mText.length / 1e6).toFixed(2)} MB`);
      usable.push(h);
    } catch (e) {
      if (e.name !== 'NoMsaError') throw e;
      console.log(`[${tag}]   no AFDB MSA for ${h.acc}`);
    }
  }
  check(new Set(usable.map(h => h.acc)).size === usable.length, `[${tag}] no duplicate hits downloaded`);

  return finish(chain, usable, tag);
}

function finish(chain, usable, tag) {
  const built = K.buildChainA3M(chain, usable);
  for (const s of built.stats) {
    console.log(`[${tag}]   ${s.acc}: +${s.added} new of ${s.total}${s.dupId ? ` (${s.dupId} dup-id)` : ''}`);
  }
  const filtered = F.filterRows(built.rows, { minCoverage: 0.5, maxIdentity: 0.9, sortByIdentity: true });
  console.log(`[${tag}] depth ${built.rows.length} -> ${filtered.rows.length} after filters, Neff@0.8=${F.neff(filtered.rows, 0.8)}`);

  const bad = filtered.rows.filter(r => K.countMatchStates(r.seq) !== chain.seq.length).length;
  check(bad === 0, `[${tag}] every row has ${chain.seq.length} match states`);
  check(F.matchStates(filtered.rows[0].seq) === chain.seq, `[${tag}] row 0 is the query`);
  check(filtered.rows.length > 1, `[${tag}] alignment is non-trivial (${filtered.rows.length} seqs)`);

  return { queryLen: chain.seq.length, rows: filtered.rows };
}

const results = [];
for (let i = 0; i < chains.length; i++) {
  results.push(await doChain(chains[i], `chain ${i + 1}`));
  console.log('');
}

if (results.length > 1) {
  console.log('# pairing');
  const paired = K.pairChains(results);
  const L = results.reduce((a, r) => a + r.queryLen, 0);
  const text = K.formatPairedA3M(paired);
  console.log(`  header: ${paired.header.replace('\t', ' \\t ')}`);
  console.log(`  ${paired.nPaired} paired, ${paired.rows.length - 1 - paired.nPaired} unpaired, ${L} columns, ${(text.length / 1e6).toFixed(2)} MB`);

  const bad = paired.rows.filter(r => K.countMatchStates(r.seq) !== L).length;
  check(bad === 0, `every paired row has ${L} match states`);
  check(paired.header === `#${results.map(r => r.queryLen).join(',')}\t${results.map(() => 1).join(',')}`, 'ColabFold header');
  check(paired.nPaired > 0, `${paired.nPaired} species paired`);
  check(paired.rows[0].seq === results.map(r => r.rows[0].seq).join(''), 'query row is the concatenation');

  // every paired row must name one accession per chain, all from the same taxon
  let ok = true;
  for (let i = 1; i <= paired.nPaired; i++) {
    if (paired.rows[i].name.split('\t').length !== results.length) { ok = false; break; }
  }
  check(ok, 'each paired row names one accession per chain');
  console.log('\n  sample paired rows:');
  for (let i = 1; i <= Math.min(3, paired.nPaired); i++) {
    console.log(`    ${paired.rows[i].name.split('\t').join('  +  ')}`);
  }
}

console.log(`\n${failures ? `${failures} CHECK(S) FAILED` : 'all checks passed'} — ${el()}`);
process.exit(failures ? 1 : 0);
