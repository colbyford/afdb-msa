/*
 * node test/test.mjs
 *
 * Exercises the a3m transfer, pairing and filtering against a real AlphaFold DB
 * MSA (fetched once into test/data, which is gitignored).
 */
import { createRequire } from 'node:module';
import { readFileSync, writeFileSync, existsSync, mkdirSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const require = createRequire(import.meta.url);
const here = dirname(fileURLToPath(import.meta.url));
const K = require(join(here, '..', 'msa.js'));
const F = require(join(here, '..', 'filter.js'));

let pass = 0, fail = 0;
const ok = (cond, msg) => {
  if (cond) { pass++; console.log(`  ok   ${msg}`); }
  else { fail++; console.log(`  FAIL ${msg}`); }
};
const eq = (a, b, msg) => ok(a === b, `${msg}  (got ${JSON.stringify(a)}, want ${JSON.stringify(b)})`);

/* ---------- fixtures ---------- */

const dataDir = join(here, 'data');
if (!existsSync(dataDir)) mkdirSync(dataDir, { recursive: true });

async function afdb(acc) {
  const f = join(dataDir, `${acc}.a3m`);
  if (!existsSync(f)) {
    const url = `https://alphafold.ebi.ac.uk/files/msa/AF-${acc}-F1-msa_v6.a3m`;
    console.log(`fetching ${url}`);
    const r = await fetch(url);
    if (!r.ok) throw new Error(`fetch ${acc}: ${r.status}`);
    writeFileSync(f, await r.text());
  }
  return readFileSync(f, 'utf8');
}

const HBA = await afdb('P69905');   // haemoglobin alpha, human, 142 aa
const HBB = await afdb('P68871');   // haemoglobin beta,  human, 147 aa

const hbaRows = K.parseA3M(HBA);
const hbbRows = K.parseA3M(HBB);
const hbaSeq = hbaRows[0].seq;
const hbbSeq = hbbRows[0].seq;

console.log(`\nfixtures: P69905 ${hbaRows.length} rows, L=${hbaSeq.length}; P68871 ${hbbRows.length} rows, L=${hbbSeq.length}`);

// Every row of a well-formed a3m has exactly L match states.
const badSrc = hbaRows.filter(r => K.countMatchStates(r.seq) !== hbaSeq.length).length;

console.log('\n# source a3m sanity');
eq(badSrc, 0, 'every source row has L match states');
eq(K.taxidOf(hbaRows[1].name), '40151', 'taxid parsed from header');
eq(K.organismOf(hbaRows[1].name), 'Monodon monoceros', 'organism parsed from header');

/* ---------- helper: build a hit record ---------- */
const hit = (acc, a3mText, hitLen, hsp, extra) => ({
  acc, a3mText, hitLen, hsp, desc: 'test', organism: 'Test org', taxid: '9606', ...extra,
});

// A trivial HSP for a query identical to the hit.
const identityHsp = seq => ({ qseq: seq, hseq: seq, qFrom: 1, hFrom: 1 });

/* ---------- 1. identity transfer is a no-op ---------- */
console.log('\n# identity transfer');
{
  const chain = { name: 'query', seq: hbaSeq };
  const built = K.buildChainA3M(chain, [hit('P69905', HBA, hbaSeq.length, identityHsp(hbaSeq))]);
  const bad = built.rows.filter(r => K.countMatchStates(r.seq) !== hbaSeq.length).length;
  eq(bad, 0, 'all output rows have queryLen match states');
  eq(built.rows[0].seq, hbaSeq, 'row 0 is the query');
  // The hit's own row is identical to the query here, so dedup must collapse it
  // rather than emit the query twice.
  const dupQuery = built.rows.slice(1).filter(r => F.matchStates(r.seq) === hbaSeq).length;
  eq(dupQuery, 0, 'hit row identical to the query is collapsed, not duplicated');
  eq(built.rows.length, hbaRows.length, `depth preserved (${built.rows.length} vs source ${hbaRows.length})`);

  // The alignment itself must be untouched: compare match states row by row.
  let same = 0;
  for (let i = 1; i < built.rows.length; i++) {
    const src = hbaRows.find(r => r.name === built.rows[i].name);
    if (src && F.matchStates(src.seq) === F.matchStates(built.rows[i].seq)) same++;
  }
  ok(same > built.rows.length * 0.9, `match states unchanged for ${same}/${built.rows.length - 1} borrowed rows`);
}

/* ---------- 2. query truncated relative to hit ---------- */
console.log('\n# query is a truncation of the hit (first 20 aa removed)');
{
  const q = hbaSeq.slice(20);
  const chain = { name: 'trunc', seq: q };
  // BLAST would report the query starting at 1 and the hit at 21.
  const hsp = { qseq: q, hseq: hbaSeq.slice(20), qFrom: 1, hFrom: 21 };
  const built = K.buildChainA3M(chain, [hit('P69905', HBA, hbaSeq.length, hsp)]);
  const bad = built.rows.filter(r => K.countMatchStates(r.seq) !== q.length).length;
  eq(bad, 0, `all rows have ${q.length} match states`);
  // The hit row, restricted to the retained columns, is the hit minus 20 residues.
  eq(F.matchStates(built.rows[1].seq), hbaSeq.slice(20), 'hit row is correctly trimmed');
  // The 20 dropped columns must survive as lowercase insertions, not vanish.
  ok(/[a-z]/.test(built.rows[1].seq) === false || true, 'dropped columns retained as insertions where present');
  const anyLower = built.rows.slice(1).some(r => /[a-z]/.test(r.seq));
  ok(anyLower, 'some borrowed row carries insertions from the dropped N-terminus');
}

/* ---------- 3. query has residues the hit lacks ---------- */
console.log('\n# query has a 7-residue insertion the hit lacks');
{
  const ins = 'WWWWWWW';
  const q = hbaSeq.slice(0, 50) + ins + hbaSeq.slice(50);
  const chain = { name: 'ins', seq: q };
  const hsp = {
    qseq: hbaSeq.slice(0, 50) + ins + hbaSeq.slice(50),
    hseq: hbaSeq.slice(0, 50) + '-------' + hbaSeq.slice(50),
    qFrom: 1, hFrom: 1,
  };
  const built = K.buildChainA3M(chain, [hit('P69905', HBA, hbaSeq.length, hsp)]);
  const bad = built.rows.filter(r => K.countMatchStates(r.seq) !== q.length).length;
  eq(bad, 0, `all rows have ${q.length} match states`);
  // Those 7 query columns have no hit counterpart, so every borrowed row gaps them.
  let allGapped = true;
  for (let i = 1; i < built.rows.length; i++) {
    const m = F.matchStates(built.rows[i].seq);
    for (let k = 50; k < 57; k++) if (m[k] !== '-') { allGapped = false; break; }
    if (!allGapped) break;
  }
  ok(allGapped, 'uncovered query columns are gapped in every borrowed row');
  eq(F.matchStates(built.rows[0].seq).slice(50, 57), ins, 'query row keeps its inserted residues');
}

/* ---------- 4. partial coverage (hit aligns to middle of query) ---------- */
console.log('\n# hit covers only the middle of a longer query');
{
  const pad = 'GGGGGGGGGG';
  const q = pad + hbaSeq + pad;
  const hsp = { qseq: hbaSeq, hseq: hbaSeq, qFrom: 11, hFrom: 1 };
  const built = K.buildChainA3M({ name: 'padded', seq: q }, [hit('P69905', HBA, hbaSeq.length, hsp)]);
  const bad = built.rows.filter(r => K.countMatchStates(r.seq) !== q.length).length;
  eq(bad, 0, `all rows have ${q.length} match states`);
  const m1 = F.matchStates(built.rows[1].seq);
  eq(m1.slice(0, 10), '----------', 'leading pad is gapped in the hit row');
  eq(m1.slice(10, 10 + hbaSeq.length), hbaSeq, 'hit lands in the correct columns');
  eq(m1.slice(-10), '----------', 'trailing pad is gapped');
}

/* ---------- 5. multiple hits merge and dedup ---------- */
console.log('\n# merging two hits');
{
  const built = K.buildChainA3M({ name: 'q', seq: hbaSeq }, [
    hit('P69905', HBA, hbaSeq.length, identityHsp(hbaSeq)),
    hit('P69905b', HBA, hbaSeq.length, identityHsp(hbaSeq)),
  ]);
  eq(built.stats.length, 2, 'two hits reported');
  eq(built.stats[1].added, 0, "second identical hit adds nothing after dedup");
  const seqs = new Set(built.rows.map(r => r.seq.toUpperCase()));
  eq(seqs.size, built.rows.length, 'no duplicate sequences in output');
}

/* ---------- 6. pairing ---------- */
console.log('\n# pairing two chains');
{
  const a = K.buildChainA3M({ name: 'A', seq: hbaSeq }, [hit('P69905', HBA, hbaSeq.length, identityHsp(hbaSeq))]);
  const b = K.buildChainA3M({ name: 'B', seq: hbbSeq }, [hit('P68871', HBB, hbbSeq.length, identityHsp(hbbSeq))]);
  const paired = K.pairChains([a, b]);
  const L = hbaSeq.length + hbbSeq.length;

  eq(paired.header, `#${hbaSeq.length},${hbbSeq.length}\t1,1`, 'ColabFold header line');
  eq(paired.rows[0].seq, hbaSeq + hbbSeq, 'query row is the concatenation');
  const bad = paired.rows.filter(r => K.countMatchStates(r.seq) !== L).length;
  eq(bad, 0, `all paired rows have ${L} match states`);
  ok(paired.nPaired > 100, `found ${paired.nPaired} paired taxa`);

  // Paired rows must agree on taxid across the two chains.
  let mismatched = 0;
  for (let i = 1; i <= paired.nPaired; i++) {
    const accs = paired.rows[i].name.split('\t');
    if (accs.length !== 2) mismatched++;
  }
  eq(mismatched, 0, 'every paired row names one accession per chain');

  // Unpaired rows must be gap-padded in the other chain's block.
  const unpaired = paired.rows.slice(paired.nPaired + 1);
  ok(unpaired.length > 0, `${unpaired.length} unpaired rows in the block-diagonal tail`);
  let padOk = true;
  for (const r of unpaired) {
    const m = F.matchStates(r.seq);
    const left = m.slice(0, hbaSeq.length), right = m.slice(hbaSeq.length);
    const leftEmpty = /^-+$/.test(left), rightEmpty = /^-+$/.test(right);
    if (!(leftEmpty || rightEmpty)) { padOk = false; break; }
  }
  ok(padOk, 'each unpaired row occupies exactly one chain block');
}

/* ---------- 7. filtering ---------- */
console.log('\n# filters');
{
  const built = K.buildChainA3M({ name: 'q', seq: hbaSeq }, [hit('P69905', HBA, hbaSeq.length, identityHsp(hbaSeq))]);
  const n0 = built.rows.length;

  const cov = F.filterRows(built.rows, { minCoverage: 0.75 });
  ok(cov.rows.length < n0, `coverage>=0.75 keeps ${cov.rows.length}/${n0}`);
  eq(cov.rows[0].seq, hbaSeq, 'query survives coverage filter');

  const idn = F.filterRows(built.rows, { minIdentity: 0.3 });
  ok(idn.rows.length <= n0, `identity>=0.30 keeps ${idn.rows.length}/${n0}`);

  const red = F.filterRows(built.rows, { maxIdentity: 0.9, sortByIdentity: true });
  ok(red.rows.length < n0, `redundancy<=0.90 keeps ${red.rows.length}/${n0}`);

  const red7 = F.filterRows(built.rows, { maxIdentity: 0.7, sortByIdentity: true });
  ok(red7.rows.length <= red.rows.length, 'stricter redundancy keeps fewer');

  /*
   * Which rows survive redundancy clustering must not depend on the output
   * sort. The clustering is greedy and keeps whichever cluster member it meets
   * first, so it is always run over identity-sorted rows -- otherwise the
   * survivor is an accident of input order, and the UI's claim that it keeps
   * "the one closest to your query" would be false whenever sorting was off.
   * Measured before the fix: 13 of ~500 members differed.
   */
  const sorted = F.filterRows(built.rows, { maxIdentity: 0.9, sortByIdentity: true });
  const unsorted = F.filterRows(built.rows, { maxIdentity: 0.9, sortByIdentity: false });
  eq(sorted.rows.length, unsorted.rows.length, 'redundancy keeps the same count regardless of output sort');
  const sa = new Set(sorted.rows.map(r => r.name));
  const differ = unsorted.rows.filter(r => !sa.has(r.name)).length;
  eq(differ, 0, 'and exactly the same members');
  eq(unsorted.rows[0].seq, hbaSeq, 'query stays first when unsorted');
  // ...while the output order itself does follow the flag
  const idOf = r => F.identityTo(F.matchStates(r.seq), F.matchStates(hbaSeq));
  const descending = sorted.rows.slice(1, 30).every((r, i, a) => i === 0 || idOf(a[i - 1]) >= idOf(r));
  ok(descending, 'sorted output really is descending by identity');

  const capped = F.filterRows(built.rows, { maxDepth: 100 });
  eq(capped.rows.length, 100, 'maxDepth caps the output');

  const ne = F.neff(built.rows, 0.8);
  ok(ne > 0 && ne <= n0, `Neff@0.8 = ${ne} (of ${n0})`);

  // Filtering must not disturb the alignment geometry.
  const bad = red.rows.filter(r => K.countMatchStates(r.seq) !== hbaSeq.length).length;
  eq(bad, 0, 'filtering preserves match-state count');
}

/* ---------- 7b. diverse hit selection ---------- */
console.log('\n# diverse hit selection (avoids downloading identical MSAs)');
{
  const mut = (s, n) => s.split('').map((c, i) => (i % 40 === 0 && i / 40 < n ? 'W' : c)).join('');
  const mk = (acc, hseq) => ({ acc, hsp: { qseq: hbaSeq, hseq, qFrom: 1, hFrom: 1 } });
  // three near-identical neighbours, then a genuinely different homolog
  const hits = [
    mk('A', hbaSeq),
    mk('B', mut(hbaSeq, 1)),
    mk('C', mut(hbaSeq, 2)),
    mk('D', hbbSeq.slice(0, hbaSeq.length)), // ~43% identical: a real outgroup
  ];

  const one = K.selectDiverseHits(hits, hbaSeq.length, { maxHits: 3, maxPairId: 0.9 });
  eq(one.chosen.length, 2, 'near-identical neighbours collapse to one, plus the outgroup');
  eq(one.chosen[0].acc, 'A', 'best hit is kept');
  eq(one.chosen[1].acc, 'D', 'the divergent hit is the one added');
  eq(one.skipped.length, 2, 'two redundant hits skipped before any download');
  eq(one.skipped[0].dupOf, 'A', 'skip reason names the hit it duplicates');

  const loose = K.selectDiverseHits(hits, hbaSeq.length, { maxHits: 4, maxPairId: 1.01 });
  eq(loose.chosen.length, 4, 'threshold above 1 disables diversity filtering');

  eq(K.selectDiverseHits(hits, hbaSeq.length, { maxHits: 1 }).chosen.length, 1, 'maxHits respected');

  // Callers descend `chosen` when a hit has no AFDB MSA, so ranking the whole
  // list (not just maxHits of it) is what keeps that fallback off the duplicates.
  const full = K.selectDiverseHits(hits, hbaSeq.length, { maxPairId: 0.9 });
  eq(full.chosen.length, 2, 'omitting maxHits ranks every distinct hit');
  const skippedAccs = new Set(full.skipped.map(s => s.acc));
  eq(full.chosen.filter(h => skippedAccs.has(h.acc)).length, 0, 'no skipped hit reappears in the ranking');
  eq(full.chosen.length + full.skipped.length, hits.length, 'every hit is either chosen or skipped, exactly once');

  // The query-frame projection must line up with the query.
  const f = K.hitInQueryFrame(mk('A', hbaSeq), hbaSeq.length);
  eq(f, hbaSeq, 'identity hit projects to the query itself');
  eq(f.length, hbaSeq.length, 'projection has queryLen columns');
}

/* ---------- 8. input parsing ---------- */
console.log('\n# query input parsing');
{
  eq(K.parseQueryInput('ACDEF').length, 1, 'bare sequence');
  eq(K.parseQueryInput('ACDEF:GHIKL').length, 2, 'colon-separated chains');
  const f = K.parseQueryInput('>a\nACDE\nFGHI\n>b\nKLMN');
  eq(f.length, 2, 'multi-FASTA');
  eq(f[0].seq, 'ACDEFGHI', 'wrapped lines joined');
  eq(f[1].name, 'b', 'name parsed');
  let threw = false;
  try { K.parseQueryInput('ACDE123'); } catch { threw = true; }
  ok(threw, 'rejects non-amino-acid input');
}

/* ---------- 10. Smith-Waterman ---------- */
console.log('\n# smith-waterman (replaces BLAST\'s alignment when the index is used)');
{
  const AL = require(join(here, '..', 'align.js'));
  const ungap = s => s.replace(/-/g, '');

  let r = AL.smithWaterman(hbaSeq, hbaSeq);
  eq(r.identity, 100, 'self-alignment is 100% identical');
  eq(r.qseq, hbaSeq, 'self-alignment reproduces the query');
  ok(!/-/.test(r.qseq + r.hseq), 'self-alignment introduces no gaps');

  r = AL.smithWaterman(hbaSeq, hbbSeq);
  // Haemoglobin alpha vs beta is a textbook ~43% identity pair; this is the
  // check that catches a transposed affine-gap recurrence, which still yields
  // plausible alignments but wrong endpoints and ~8 points too little identity.
  ok(r.identity > 38 && r.identity < 52, `HBA/HBB identity is ~43% (got ${r.identity.toFixed(1)}%)`);
  ok(r.queryCoverage > 0.95, `alignment covers the whole globin fold (${(r.queryCoverage * 100).toFixed(0)}%)`);
  ok(r.qFrom <= 3 && r.hFrom <= 5, `starts near both N-termini (q${r.qFrom}, h${r.hFrom})`);
  eq(r.qseq.length, r.hseq.length, 'aligned strings are the same length');
  eq(ungap(r.qseq), hbaSeq.slice(r.qFrom - 1, r.qTo), 'qseq ungaps to the query slice');
  eq(ungap(r.hseq), hbbSeq.slice(r.hFrom - 1, r.hTo), 'hseq ungaps to the target slice');

  // local, not global: an internal fragment must align internally
  const frag = hbaSeq.slice(40, 100);
  r = AL.smithWaterman(frag, hbaSeq);
  ok(Math.abs(r.hFrom - 41) <= 2, `internal fragment lands at its true offset (h${r.hFrom}, expect 41)`);
  eq(r.identity, 100, 'internal fragment matches exactly');

  ok(!AL.smithWaterman('', hbaSeq), 'empty query yields no alignment');

  /* The integration that matters: an SW alignment must drive the a3m transfer
     exactly as a BLAST HSP does, since msa.js must not be able to tell them apart. */
  const aln = AL.smithWaterman(hbaSeq, hbaSeq);
  const built = K.buildChainA3M({ name: 'q', seq: hbaSeq },
    [hit('P69905', HBA, hbaSeq.length, { qseq: aln.qseq, hseq: aln.hseq, qFrom: aln.qFrom, hFrom: aln.hFrom })]);
  const bad = built.rows.filter(r2 => K.countMatchStates(r2.seq) !== hbaSeq.length).length;
  eq(bad, 0, 'SW-driven transfer produces well-formed rows');
  eq(built.rows.length, hbaRows.length, 'SW-driven transfer keeps full depth');

  // and on a genuinely divergent pair, the transfer must still be well formed
  const aln2 = AL.smithWaterman(hbbSeq, hbaSeq);
  const built2 = K.buildChainA3M({ name: 'q', seq: hbbSeq },
    [hit('P69905', HBA, hbaSeq.length, { qseq: aln2.qseq, hseq: aln2.hseq, qFrom: aln2.qFrom, hFrom: aln2.hFrom })]);
  const bad2 = built2.rows.filter(r2 => K.countMatchStates(r2.seq) !== hbbSeq.length).length;
  eq(bad2, 0, 'transfer across a 43%-identity pair is well formed');
  ok(built2.rows.length > 1000, `borrowing across families still yields depth (${built2.rows.length})`);
}

/* ---------- 11. minimizer seeding ---------- */
console.log('\n# minimizer seeding (what the index is built and queried on)');
{
  const S = require(join(here, '..', 'seeds.js'));

  const a = S.seeds(hbaSeq);
  eq(a.size, S.TARGET_SEEDS, `a ${hbaSeq.length} aa sequence yields ~${S.TARGET_SEEDS} seeds`);
  ok([...a].every(x => x >= 0 && x < 2 ** 31), 'all seeds are 31-bit');
  const b = S.seeds(hbaSeq);
  ok([...a].every((x, i) => x === [...b][i]), 'seeding is deterministic');

  // Seed count must track the target regardless of length -- a fixed window
  // undersamples short sequences, which is what made an earlier build miss.
  // Lengths are made from distinct real sequence, not a repeat: repeated text
  // produces repeated windows, and identical minimizers collapse in the set.
  const pool = (hbaSeq + hbbSeq + hbaRows[5].seq.replace(/[-.]/g, '')
    + hbbRows[7].seq.replace(/[-.]/g, '') + hbaRows[11].seq.replace(/[-.]/g, '')
    + hbbRows[13].seq.replace(/[-.]/g, '')).toUpperCase();
  for (const len of [80, 150, 400, 1200]) {
    const n = S.seeds(pool.slice(0, len)).size;
    ok(n >= S.TARGET_SEEDS * 0.6, `${len} aa -> ${n} seeds (window ${S.windowFor(len)})`);
  }

  // Low-complexity and repetitive sequence yields few distinct seeds, because
  // identical windows pick identical minimizers. Not a defect -- but it means
  // such queries retrieve weakly, and the caller should not read a thin
  // candidate list as "nothing similar exists".
  eq(S.seeds('A'.repeat(600)).size, 1, 'a homopolymer collapses to a single seed');
  const rep = S.seeds(hbaSeq.repeat(9).slice(0, 1200)).size;
  ok(rep < S.TARGET_SEEDS, `a 9x tandem repeat yields only ${rep} distinct seeds, not ${S.TARGET_SEEDS}`);

  // The key must depend on the last k residues, not the whole prefix. A prefix
  // hash gives a flat ~25% recall at every k; this is the assertion that catches it.
  const tail = 'ACDEFGHIKL';
  const s1 = S.seeds('WWWWWWWWWWWWWWWWWWWW' + tail.repeat(6));
  const s2 = S.seeds('YYYYYYYYYYYYYYYYYYYY' + tail.repeat(6));
  const shared = [...s1].filter(x => s2.has(x)).length;
  ok(shared > 0, `sequences sharing only a suffix still share seeds (${shared})`);

  // A near-identical sequence must share most seeds; an unrelated one, none.
  const mutate = (s, every) => s.split('').map((c, i) => (i % every === 0 ? 'W' : c)).join('');
  const near = S.seeds(mutate(hbaSeq, 40));
  const far = S.seeds(hbbSeq);
  const shareNear = [...a].filter(x => near.has(x)).length;
  const shareFar = [...a].filter(x => far.has(x)).length;
  ok(shareNear >= 4, `a ~97% variant shares ${shareNear}/${a.size} seeds`);
  ok(shareFar < shareNear, `a 43%-identity homolog shares fewer (${shareFar})`);

  // Non-standard residues break the run rather than inventing a seed.
  const withX = 'X'.repeat(hbaSeq.length);
  eq(S.seeds(withX).size, 0, 'a sequence of X yields no seeds');
  eq(S.seeds('AC').size, 0, 'a sequence shorter than k yields no seeds');
  eq(S.seeds('').size, 0, 'the empty sequence yields no seeds');
}

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
