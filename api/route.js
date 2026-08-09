const INDEX_URL = 'https://huggingface.co/datasets/sokrypton/afdb-msa-index/resolve/main';
const TOP_CANDIDATES = 20;

let sharedIndex = null;

function writeText(text) {
  document.body.textContent = text;
}

function jsonError(message) {
  writeText(JSON.stringify({ error: message }));
}

async function processChain(chain) {
  if (!sharedIndex) sharedIndex = new window.Search.MinimizerIndex(INDEX_URL);
  await sharedIndex.load();

  const cands = await sharedIndex.search(chain.seq, { topK: TOP_CANDIDATES });
  if (!cands.length) {
    throw new Error(
      `Nothing in AlphaFold DB resembles ${chain.name}. ` +
      `A sequence with no natural homologs has no alignment to borrow.`
    );
  }

  const info = await window.Api.candidateInfo(cands.map(c => c.acc));

  const scored = [];
  for (const c of cands) {
    const t = info.get(c.acc);
    if (!t) continue;
    const aln = window.Align.smithWaterman(chain.seq, t.seq);
    if (!aln) continue;
    scored.push({
      acc: c.acc,
      score: aln.score,
      identity: Number(aln.identity.toFixed(1)),
      queryCoverage: aln.queryCoverage,
      desc: t.desc,
      organism: t.organism,
      taxid: t.taxid,
      hsp: { qseq: aln.qseq, hseq: aln.hseq, qFrom: aln.qFrom, hFrom: aln.hFrom },
    });
  }

  scored.sort((a, b) => b.score - a.score);
  if (!scored.length) {
    throw new Error(`Could not retrieve sequences for any candidate of ${chain.name}.`);
  }

  let donor = null;
  for (const h of scored) {
    try {
      h.a3mText = await window.Api.fetchMsaCached(h.acc);
      donor = h;
      break;
    } catch (e) {
      if (e.name !== 'NoMsaError') throw e;
    }
  }

  if (!donor) {
    throw new Error(`None of the candidates for ${chain.name} has an AlphaFold DB MSA.`);
  }

  const built = window.MSAKit.buildChainA3M(chain, [donor]);
  return { chain, rows: built.rows, queryLen: chain.seq.length };
}

async function main() {
  const url = new URL(window.location.href);
  const seqText = (url.searchParams.get('seq') || '').trim();
  const format = (url.searchParams.get('format') || '').toLowerCase();

  if (!seqText) {
    jsonError('Missing sequence. Provide ?seq=<sequence>.');
    return;
  }

  let chains;
  try {
    chains = window.MSAKit.parseQueryInput(seqText);
  } catch (e) {
    jsonError(e.message);
    return;
  }

  if (!chains.length) {
    jsonError('No valid sequences found in input.');
    return;
  }

  let results;
  try {
    results = await Promise.all(chains.map(processChain));
  } catch (e) {
    jsonError(e.message || 'Pipeline failed.');
    return;
  }

  const wantsAll = format === 'all';
  const wantsPaired = format === 'paired';
  const wantsJson = format === 'json' || wantsAll || (!format && chains.length > 1);

  if (wantsJson) {
    const out = {};
    for (const r of results) out[r.chain.name] = window.MSAKit.formatA3M(r.rows);

    if (results.length > 1) {
      const chainData = results.map(r => ({ queryLen: r.queryLen, rows: r.rows }));
      const paired = window.MSAKit.pairChains(chainData);
      out['Paired Alignment'] = window.MSAKit.formatPairedA3M(paired);
    }

    writeText(JSON.stringify(out));
    return;
  }

  if (wantsPaired) {
    if (results.length < 2) {
      jsonError('format=paired requires at least two sequences/chains.');
      return;
    }
    const chainData = results.map(r => ({ queryLen: r.queryLen, rows: r.rows }));
    const paired = window.MSAKit.pairChains(chainData);
    writeText(window.MSAKit.formatPairedA3M(paired));
    return;
  }

  writeText(window.MSAKit.formatA3M(results[0].rows));
}

main();
