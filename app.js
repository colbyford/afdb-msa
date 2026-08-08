/* app.js -- UI wiring. All the real work lives in seeds/search/align/msa/filter. */

const $ = id => document.getElementById(id);

// The published index. There is nothing else to configure and nothing to fall
// back to: this is the search.
const INDEX_URL = 'https://huggingface.co/datasets/sokrypton/afdb-msa-index/resolve/main';
const TOP_CANDIDATES = 20;   // ranked by shared seeds; alignment picks the winner

const EXAMPLES = {
  hba: '>HBA_HUMAN\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR',
  hbab: '>HBA_HUMAN\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR\n>HBB_HUMAN\nMVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH',
  ubq: '>UBIQUITIN\nMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',
};

/* ---------- status ---------- */

let logLines = [];
function log(msg) {
  logLines.push(msg);
  $('log-panel').hidden = false;
  const el = $('log');
  el.textContent = logLines.join('\n');
  el.scrollTop = el.scrollHeight;
}
function status(msg, kind) {
  const el = $('status');
  el.textContent = msg || '';
  el.className = 'status' + (kind ? ' ' + kind : '');
}

function refreshSummary() {
  const el = $('seq-summary');
  try {
    const chains = MSAKit.parseQueryInput($('seq').value);
    if (!chains.length) { el.textContent = ''; return; }
    const parts = chains.map(c => `${c.name} (${c.seq.length} aa)`);
    el.textContent = `${chains.length} chain${chains.length > 1 ? 's' : ''}: ${parts.join(', ')}`
      + (chains.length > 1 ? ' — will be paired by species' : '');
    el.className = 'muted small';
  } catch (e) {
    el.textContent = e.message;
    el.className = 'small err';
  }
}

/* ---------- one chain ---------- */

let sharedIndex = null;

/*
 * The whole pipeline for one chain:
 *
 *   seeds -> index -> candidate accessions -> Smith-Waterman -> best donor
 *         -> that donor's AFDB MSA -> re-indexed onto this query
 *
 * There is no fallback. The index covers 239.6M AFDB entries and answers in
 * ~50 KB; when it finds nothing close, that is the answer, and it is reported
 * rather than papered over by a slower search that would not do better. An
 * earlier version tried a UniParc checksum shortcut first and EBI BLAST after,
 * and neither earned its complexity: the index already returns exact matches at
 * 100%, and below its range a borrowed alignment is not worth having.
 */
async function processChain(chain, cfg, signal, tag) {
  log(`\n[${tag}] ${chain.name}: ${chain.seq.length} aa`);

  if (!sharedIndex) sharedIndex = new Search.MinimizerIndex(INDEX_URL);
  status(`${tag}: searching AlphaFold DB…`);
  await sharedIndex.load();

  const t0 = Date.now();
  const cands = await sharedIndex.search(chain.seq, { topK: TOP_CANDIDATES, signal });
  log(`[${tag}] ${cands.length} candidates from ${cands.candidates} scored, `
    + `${cands.seedsFound}/${cands.seeds} seeds, ${(cands.bytesRead / 1024).toFixed(0)} KB, ${Date.now() - t0} ms`);

  if (!cands.length) {
    throw new Error(`Nothing in AlphaFold DB resembles ${chain.name}. A sequence with no natural `
      + `homologs has no alignment to borrow — designed proteins usually land here.`);
  }

  status(`${tag}: aligning ${cands.length} candidates…`);
  const info = await Api.candidateInfo(cands.map(c => c.acc), signal);

  const scored = [];
  for (const c of cands) {
    const t = info.get(c.acc);
    if (!t) continue;
    const aln = Align.smithWaterman(chain.seq, t.seq);
    if (!aln) continue;
    scored.push({
      acc: c.acc, hitLen: t.seq.length,
      identity: Number(aln.identity.toFixed(1)),
      queryCoverage: aln.queryCoverage,
      score: aln.score,
      desc: t.desc, organism: t.organism, taxid: t.taxid,
      hsp: { qseq: aln.qseq, hseq: aln.hseq, qFrom: aln.qFrom, hFrom: aln.hFrom },
    });
  }
  scored.sort((a, b) => b.score - a.score);
  if (!scored.length) throw new Error(`Could not retrieve sequences for any candidate of ${chain.name}.`);

  // Best donor that actually has an MSA. Nearly all do; a few (viral entries,
  // for one) 404, so walk down rather than fail on the first miss.
  let donor = null;
  for (const h of scored) {
    try {
      status(`${tag}: downloading MSA for ${h.acc}…`);
      h.a3mText = await Api.fetchMsaCached(h.acc, (got, total) => {
        status(`${tag}: downloading ${h.acc} — ${humanSize(got)}${total ? ` / ${humanSize(total)}` : ''}`);
      }, signal);
      donor = h;
      break;
    } catch (e) {
      if (e.name !== 'NoMsaError') throw e;
      log(`[${tag}] ${h.acc}: no MSA, trying the next candidate`);
    }
  }
  if (!donor) throw new Error(`None of the candidates for ${chain.name} has an AlphaFold DB MSA.`);
  log(`[${tag}] ${donor.acc} at ${donor.identity}% identity, `
    + `${(donor.queryCoverage * 100).toFixed(0)}% coverage — ${humanSize(donor.a3mText.length)}`);

  status(`${tag}: transferring the alignment onto your query…`);
  const built = MSAKit.buildChainA3M(chain, [donor]);
  const filtered = Filter.filterRows(built.rows, {
    minCoverage: cfg.minCov, minIdentity: cfg.minId, maxIdentity: cfg.maxId,
    maxDepth: cfg.maxDepth, sortByIdentity: cfg.sortById,
  });
  for (const w of filtered.warnings) log(`[${tag}]   ${w}`);
  log(`[${tag}] ${built.rows.length} sequences transferred, ${filtered.rows.length} after filters`);

  return { chain, donor, rows: filtered.rows, queryLen: chain.seq.length };
}

/* ---------- results ---------- */

function download(name, text) {
  const blob = new Blob([text], { type: 'text/plain' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = name;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

// "0.0 MB" for a small alignment reads like a bug, so scale the unit.
function humanSize(n) {
  if (n < 1024) return `${n} B`;
  if (n < 1024 * 1024) return `${(n / 1024).toFixed(0)} KB`;
  return `${(n / 1048576).toFixed(1)} MB`;
}

function el(tag, cls, text) {
  const n = document.createElement(tag);
  if (cls) n.className = cls;
  if (text !== undefined) n.textContent = text;
  return n;
}

function renderResults(results) {
  const root = $('results');
  root.hidden = false;
  root.textContent = '';

  for (const res of results) {
    const p = el('section', 'panel');
    p.appendChild(el('h2', null, res.chain.name));

    const d = res.donor;
    if (d.identity < 90) {
      p.appendChild(el('p', 'warn',
        `Closest match is ${d.identity}% identical. Seeding is reliable above ~90% and fades below ~70%, `
        + `so this alignment describes a relative of your sequence rather than your sequence itself.`));
    }
    if (d.queryCoverage < 0.8) {
      p.appendChild(el('p', 'warn',
        `The match covers only ${(d.queryCoverage * 100).toFixed(0)}% of your query; the rest is all-gap.`));
    }

    const t = el('table', 'hits');
    const head = el('tr');
    for (const h of ['Donor', 'Identity', 'Query cov.', 'Organism', 'Description']) head.appendChild(el('th', null, h));
    t.appendChild(head);
    const tr = el('tr');
    const a = el('a', null, d.acc);
    a.href = `https://alphafold.ebi.ac.uk/entry/${d.acc}`;
    a.target = '_blank'; a.rel = 'noopener';
    const td = el('td'); td.appendChild(a); tr.appendChild(td);
    tr.appendChild(el('td', 'num', `${d.identity}%`));
    tr.appendChild(el('td', 'num', `${(d.queryCoverage * 100).toFixed(0)}%`));
    tr.appendChild(el('td', null, d.organism || '—'));
    tr.appendChild(el('td', null, d.desc || '—'));
    t.appendChild(tr);
    p.appendChild(t);

    const neff = Filter.neff(res.rows, 0.8);
    const stats = el('p', 'stats');
    stats.innerHTML = `<b>${res.rows.length}</b> sequences &middot; <b>${res.queryLen}</b> columns &middot; `
      + `Neff@0.8 <b>${neff}</b> (${(neff / res.rows.length * 100).toFixed(0)}% effective)`;
    p.appendChild(stats);

    const text = MSAKit.formatA3M(res.rows);
    const b = el('button', 'primary', `Download ${res.chain.name}.a3m (${humanSize(text.length)})`);
    b.onclick = () => download(`${res.chain.name}.a3m`, text);
    p.appendChild(b);
    root.appendChild(p);
  }

  if (results.length > 1) {
    const p = el('section', 'panel');
    p.appendChild(el('h2', null, 'Paired alignment'));
    const chains = results.map(r => ({ queryLen: r.queryLen, rows: r.rows }));
    const paired = MSAKit.pairChains(chains);
    const total = paired.rows.length - 1;
    const stats = el('p', 'stats');
    stats.innerHTML = `<b>${paired.nPaired}</b> paired by species &middot; `
      + `<b>${total - paired.nPaired}</b> unpaired (block-diagonal) &middot; `
      + `<b>${paired.lens.reduce((a, b) => a + b, 0)}</b> columns`;
    p.appendChild(stats);
    p.appendChild(el('p', 'muted small',
      'ColabFold layout: a "#L1,L2  1,1" header, the concatenated query, species-paired rows, '
      + 'then each chain’s remaining sequences padded with gaps.'));
    const text = MSAKit.formatPairedA3M(paired);
    const b = el('button', 'primary', `Download paired.a3m (${humanSize(text.length)})`);
    b.onclick = () => download('paired.a3m', text);
    p.appendChild(b);
    root.appendChild(p);
  }
}

/* ---------- run ---------- */

let controller = null;

function readConfig() {
  const num = (id, dflt) => {
    const v = parseFloat($(id).value);
    return Number.isFinite(v) ? v : dflt;
  };
  return {
    minCov: num('mincov', 0),
    minId: num('minid', 0),
    maxId: num('maxid', 1),
    maxDepth: Math.max(0, num('maxdepth', 0)),
    sortById: $('sortid').checked,
  };
}

async function run() {
  logLines = [];
  $('log').textContent = '';
  $('results').hidden = true;
  status('');

  let chains;
  try {
    chains = MSAKit.parseQueryInput($('seq').value);
  } catch (e) {
    status(e.message, 'err');
    return;
  }
  if (!chains.length) { status('Paste a sequence first.', 'err'); return; }

  const cfg = readConfig();
  controller = new AbortController();
  $('run').disabled = true;
  $('cancel').hidden = false;

  try {
    const t0 = Date.now();
    const results = await Promise.all(chains.map((c, i) =>
      processChain(c, cfg, controller.signal, `chain ${i + 1}`)));
    renderResults(results);
    status(`Done in ${((Date.now() - t0) / 1000).toFixed(1)}s — ${results.map(r => r.rows.length).join(' + ')} sequences.`, 'ok');
    log('\nDone.');
  } catch (e) {
    if (e.name === 'AbortError' || /cancel/i.test(e.message)) {
      status('Cancelled.', '');
      log('Cancelled.');
    } else {
      status(e.message, 'err');
      log(`ERROR: ${e.message}`);
    }
  } finally {
    controller = null;
    $('run').disabled = false;
    $('cancel').hidden = true;
  }
}

document.addEventListener('DOMContentLoaded', () => {
  $('run').onclick = run;
  $('cancel').onclick = () => { if (controller) controller.abort(); };
  $('seq').addEventListener('input', refreshSummary);
  document.querySelectorAll('[data-example]').forEach(b => {
    b.onclick = () => {
      const k = b.dataset.example;
      $('seq').value = k === 'clear' ? '' : EXAMPLES[k];
      refreshSummary();
    };
  });
  refreshSummary();
});
