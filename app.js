/* app.js -- UI wiring. All the real work lives in msa.js, filter.js and api.js. */

const $ = id => document.getElementById(id);

const EXAMPLES = {
  hba: '>HBA_HUMAN\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR',
  hbab: '>HBA_HUMAN\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR\n>HBB_HUMAN\nMVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH',
  ubq: '>UBIQUITIN\nMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',
};

/* ---------- log & status ---------- */

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

/* ---------- query summary as you type ---------- */

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

/* ---------- one chain, end to end ---------- */

/*
 * Fast path. If the query is *exactly* a UniProt sequence -- the common case
 * when someone pastes an entry -- UniParc resolves it by CRC64 in under a
 * second, and no alignment is needed: the query IS the hit, so the transfer is
 * the identity map. This skips BLAST entirely (0.4 s instead of 8-250 s).
 *
 * Returns null when there is no exact match, and the caller falls back to BLAST.
 */
async function tryExactMatch(chain, cfg, signal, tag) {
  status(`${tag}: checking for an exact UniProt match…`);
  let accs;
  try {
    accs = await Api.uniparcAccessions(chain.seq, signal);
  } catch (e) {
    log(`[${tag}] UniParc lookup failed (${e.message}); falling back to BLAST`);
    return null;
  }
  if (!accs.length) {
    log(`[${tag}] no exact UniProt match — searching with BLAST`);
    return null;
  }
  log(`[${tag}] exact UniProt match: ${accs.slice(0, 6).join(', ')}${accs.length > 6 ? `, +${accs.length - 6} more` : ''}`);

  let got;
  try {
    status(`${tag}: downloading MSA…`);
    got = await Api.afdbFetchFirstAvailable(accs, (n, t) => {
      status(`${tag}: downloading ${humanSize(n)}${t ? ` / ${humanSize(t)}` : ''}`);
    }, signal, acc => log(`[${tag}] ${acc}: no AFDB MSA, trying the next accession`));
  } catch (e) {
    if (e.name !== 'NoMsaError') throw e;
    log(`[${tag}] none of the exactly-matching entries has an AFDB MSA — falling back to BLAST`);
    return null;
  }

  const meta = await Api.uniprotSummary(got.acc, signal);
  const hit = {
    acc: got.acc, db: 'exact', a3mText: got.text, hitLen: chain.seq.length,
    identity: 100, queryCoverage: 1, evalue: 0,
    desc: meta.desc || '', organism: meta.organism || '', taxid: meta.taxid || '',
    // query === hit, so the pairwise alignment is the identity map
    hsp: { qseq: chain.seq, hseq: chain.seq, qFrom: 1, hFrom: 1 },
  };
  log(`[${tag}] ${got.acc}: ${humanSize(got.text.length)} — BLAST skipped`);
  return { hits: [hit], all: [hit], skipped: [], exact: true };
}

/*
 * Middle path: the static minimizer index over AlphaFold DB. Replaces BLAST's
 * search with a few dozen ranged GETs, then reconstructs the one other thing
 * BLAST supplied -- the query/hit alignment -- with Smith-Waterman.
 *
 * Measured over 239,602,633 AFDB entries (99.4% of v6), against BLAST on full
 * UniProtKB as the reference:
 *
 *   hits at >=90% identity     index 100%   BLAST 81%
 *   median identity            index  95%   BLAST 95%
 *   time                       index   1ms  BLAST 287s
 *   read per query             ~50 KB of a 17 GB index
 *
 * The index wins outright where BLAST cannot: 63% of AFDB entries have been
 * deleted from current UniProtKB, so BLAST literally cannot return them and
 * settles for a distant relative. AFDB is its own authority here.
 *
 * Seeding is exact-k-mer, so sensitivity falls off below ~70% identity. That is
 * deliberate -- the target is close relatives worth borrowing an MSA from -- and
 * anything weaker falls through to BLAST.
 */
let sharedIndex = null;

async function tryIndex(chain, cfg, signal, tag) {
  if (!cfg.indexUrl) return null;
  try {
    if (!sharedIndex) sharedIndex = new Search.MinimizerIndex(cfg.indexUrl);
    status(`${tag}: searching the AlphaFold DB index…`);
    await sharedIndex.load();
  } catch (e) {
    log(`[${tag}] index unavailable (${e.message}); falling back to BLAST`);
    sharedIndex = null;
    return null;
  }

  const t0 = Date.now();
  let belowFloor = false;
  let cands;
  try {
    cands = await sharedIndex.search(chain.seq, { topK: cfg.topK, signal });
  } catch (e) {
    log(`[${tag}] index search failed (${e.message}); falling back to BLAST`);
    return null;
  }
  log(`[${tag}] index: ${cands.length} candidates from ${cands.candidates} scored, `
    + `${cands.seedsFound}/${cands.seeds} seeds found, `
    + `${(cands.bytesRead / 1024).toFixed(0)} KB read, ${Date.now() - t0} ms`);
  if (!cands.length) return null;

  // The index returns accessions; alignment needs residues.
  status(`${tag}: fetching ${cands.length} candidate sequences…`);
  const seqs = await Api.candidateSequences(cands.map(c => c.acc), signal,
    note => log(`[${tag}] ${note}`));
  if (!seqs.size) { log(`[${tag}] no candidate sequences retrievable; falling back to BLAST`); return null; }

  status(`${tag}: aligning candidates…`);
  const scored = [];
  for (const c of cands) {
    const target = seqs.get(c.acc);
    if (!target) continue;
    const aln = Align.smithWaterman(chain.seq, target);
    if (!aln) continue;
    scored.push({
      acc: c.acc, db: 'AFDB', hitLen: target.length,
      identity: Number(aln.identity.toFixed(1)),
      queryCoverage: aln.queryCoverage,
      evalue: null, score: aln.score,
      desc: '', organism: '', taxid: '',
      hsp: { qseq: aln.qseq, hseq: aln.hseq, qFrom: aln.qFrom, hFrom: aln.hFrom },
    });
  }
  scored.sort((a, b) => b.score - a.score);
  if (!scored.length) return null;

  const best = scored[0];
  log(`[${tag}] index best: ${best.acc} at ${best.identity}% identity, ${(best.queryCoverage * 100).toFixed(0)}% coverage`);
  if (best.identity < cfg.minIndexIdentity) {
    // Below the floor the index defers to BLAST, which searches more of UniProt
    // and is more sensitive. But only if BLAST can actually run: without a
    // contact email it cannot, and erroring out would discard a real hit --
    // an 86% match at full coverage still transfers a usable MSA. Take it, and
    // say plainly that it is under the bar.
    const blastPossible = /^[^@\s]+@[^@\s]+\.[^@\s]+$/.test(cfg.email);
    if (blastPossible) {
      log(`[${tag}] below the ${cfg.minIndexIdentity}% floor — falling back to BLAST, which searches far more of UniProt`);
      return null;
    }
    log(`[${tag}] best is ${best.identity}%, under the ${cfg.minIndexIdentity}% floor, and BLAST needs a contact email — using it anyway`);
    belowFloor = true;
  }

  // Same diversity rule as the BLAST path: merging near-identical hits buys
  // nothing, so screen before downloading any MSA.
  const sel = MSAKit.selectDiverseHits(scored, chain.seq.length, { maxPairId: cfg.maxPairId });
  const usable = [];
  const rejected = [];
  for (const h of sel.chosen) {
    if (usable.length >= cfg.maxHits) break;
    try {
      status(`${tag}: downloading MSA for ${h.acc}…`);
      h.a3mText = await Api.afdbFetchMsaCached(h.acc, (got, total) => {
        status(`${tag}: downloading ${h.acc} — ${humanSize(got)}${total ? ` / ${humanSize(total)}` : ''}`);
      }, signal);
      log(`[${tag}] ${h.acc}: ${humanSize(h.a3mText.length)}`);
      usable.push(h);
    } catch (e) {
      if (e.name !== 'NoMsaError') throw e;
      rejected.push(h.acc);
      log(`[${tag}] ${h.acc}: no AFDB MSA, trying the next candidate`);
    }
  }
  if (!usable.length) { log(`[${tag}] no index candidate has an AFDB MSA; falling back to BLAST`); return null; }

  return { hits: usable, all: scored, skipped: sel.skipped, rejected, exact: false, viaIndex: true, belowFloor };
}

/* Slow path: BLAST at EBI, then pick diverse hits. */
async function searchWithBlast(chain, cfg, signal, tag) {
  if (!/^[^@\s]+@[^@\s]+\.[^@\s]+$/.test(cfg.email)) {
    throw new Error('No exact UniProt match, so a BLAST search is needed — please enter a contact email (EBI requires one).');
  }
  status(`${tag}: submitting BLAST…`);

  const job = await Api.blastSubmit(chain.seq, {
    email: cfg.email, database: cfg.database, evalue: cfg.evalue, nHits: 50,
  });
  log(`[${tag}] BLAST job ${job}`);

  await Api.blastWait(job, (s, secs) => {
    status(`${tag}: BLAST ${s.toLowerCase()} (${secs}s)`);
  }, signal);

  const all = Api.parseBlastHits(await Api.blastResult(job));
  log(`[${tag}] ${all.length} hits`);
  if (!all.length) throw new Error(`No BLAST hits for ${chain.name}. Try a higher E-value or a different database.`);

  // Rank the hits for diversity BEFORE downloading anything -- see
  // selectDiverseHits. This ranks the whole list, so the availability loop below
  // can keep descending it rather than falling back on rejected duplicates.
  const sel = MSAKit.selectDiverseHits(all, chain.seq.length, { maxPairId: cfg.maxPairId });
  const dropped = sel.skipped.length;
  if (dropped) {
    log(`[${tag}] ${dropped} hit(s) >${cfg.maxPairId} identical to a better one — not downloaded ` +
        `(e.g. ${sel.skipped.slice(0, 3).map(s => `${s.acc}→${s.dupOf}`).join(', ')}${dropped > 3 ? ', …' : ''})`);
  }
  log(`[${tag}] ${sel.chosen.length} distinct hit(s) available to merge`);

  // Descend the diverse ranking, downloading as we go. There is no separate
  // existence probe: AFDB ignores Range and rejects HEAD under CORS, so probing
  // would download every file twice. A missing MSA surfaces as NoMsaError here.
  const usable = [];
  const rejected = [];
  for (const h of sel.chosen) {
    if (usable.length >= cfg.maxHits) break;
    try {
      status(`${tag}: downloading MSA for ${h.acc}…`);
      h.a3mText = await Api.afdbFetchMsaCached(h.acc, (got, total) => {
        status(`${tag}: downloading ${h.acc} — ${humanSize(got)}${total ? ` / ${humanSize(total)}` : ''}`);
      }, signal);
      log(`[${tag}] ${h.acc}: ${humanSize(h.a3mText.length)}, identity ${h.identity}%, coverage ${(h.queryCoverage * 100).toFixed(0)}%`);
      usable.push(h);
    } catch (e) {
      if (e.name !== 'NoMsaError') throw e;
      rejected.push(h.acc);
      log(`[${tag}] ${h.acc}: no AFDB MSA, trying the next distinct hit`);
    }
  }
  if (!usable.length) throw new Error(`None of the hits for ${chain.name} have an AlphaFold DB MSA.`);
  if (usable.length < cfg.maxHits) {
    log(`[${tag}] merging ${usable.length} hit(s), not ${cfg.maxHits} — no further hit is distinct enough to add anything`);
  }

  return { hits: usable, all, skipped: sel.skipped, rejected, exact: false };
}

/* ---------- one chain, end to end ---------- */

async function processChain(chain, cfg, signal, tag) {
  log(`\n[${tag}] ${chain.name}: ${chain.seq.length} aa`);

  const found = (cfg.tryExact && await tryExactMatch(chain, cfg, signal, tag))
    || (cfg.useIndex && await tryIndex(chain, cfg, signal, tag))
    || await searchWithBlast(chain, cfg, signal, tag);
  const { hits: usable, all, skipped, rejected, exact, viaIndex, belowFloor } = found;

  status(`${tag}: transferring alignment onto the query…`);
  const built = MSAKit.buildChainA3M(chain, usable);
  log(`[${tag}] merged depth ${built.rows.length}`);
  for (const s of built.stats) {
    const pct = s.total ? ((s.added / s.total) * 100).toFixed(0) : '0';
    log(`[${tag}]   ${s.acc}: +${s.added} new of ${s.total} (${pct}% new)`);
  }

  const filtered = Filter.filterRows(built.rows, {
    minCoverage: cfg.minCov, minIdentity: cfg.minId, maxIdentity: cfg.maxId,
    maxDepth: cfg.maxDepth, sortByIdentity: cfg.sortById,
  });
  for (const w of filtered.warnings) log(`[${tag}]   ${w}`);

  return { chain, hits: usable, allHits: all, skipped, rejected, exact, viaIndex, belowFloor, built, rows: filtered.rows, queryLen: chain.seq.length };
}

/* ---------- results rendering ---------- */

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

function hitTable(res) {
  const t = el('table', 'hits');
  const head = el('tr');
  for (const h of ['Accession', 'DB', 'Identity', 'Query cov.', 'E-value', 'Organism', 'Contributed']) {
    head.appendChild(el('th', null, h));
  }
  t.appendChild(head);
  res.hits.forEach((h, i) => {
    const st = res.built.stats[i] || {};
    const tr = el('tr');
    const a = el('a', null, h.acc);
    a.href = `https://alphafold.ebi.ac.uk/entry/${h.acc}`;
    a.target = '_blank'; a.rel = 'noopener';
    const td = el('td'); td.appendChild(a); tr.appendChild(td);
    tr.appendChild(el('td', null,
      h.db === 'exact' ? 'exact match' : h.db === 'AFDB' ? 'AFDB index' : h.db === 'SP' ? 'SwissProt' : 'TrEMBL'));
    tr.appendChild(el('td', 'num', `${h.identity}%`));
    tr.appendChild(el('td', 'num', `${(h.queryCoverage * 100).toFixed(0)}%`));
    tr.appendChild(el('td', 'num', h.evalue === null || h.evalue === undefined ? '—' : String(h.evalue)));
    tr.appendChild(el('td', null, h.organism || '—'));
    tr.appendChild(el('td', 'num', st.total ? `${st.added} / ${st.total}` : '—'));
    t.appendChild(tr);
  });
  return t;
}

function renderResults(results, cfg) {
  const root = $('results');
  root.hidden = false;
  root.textContent = '';

  for (const res of results) {
    const p = el('section', 'panel');
    p.appendChild(el('h2', null, res.chain.name));

    const best = res.hits[0];
    const warn = [];
    if (res.belowFloor) {
      warn.push(`Best hit is ${best.identity}% identical, below the ${'' + (res.floor || 90)}% floor you asked for. `
        + `BLAST could search further but needs a contact email. This alignment is still usable — judge it by the identity and coverage below.`);
    }
    if (best.identity < 30) warn.push(`Best hit is only ${best.identity}% identical — the borrowed alignment may be a poor fit.`);
    if (best.queryCoverage < 0.8) warn.push(`Best hit covers only ${(best.queryCoverage * 100).toFixed(0)}% of the query; uncovered columns are all-gap.`);
    for (const w of warn) p.appendChild(el('p', 'warn', w));

    if (res.exact) {
      p.appendChild(el('p', 'muted small',
        'Exact UniProt match — resolved by checksum in under a second, so no BLAST search was needed.'));
    } else if (res.viaIndex) {
      p.appendChild(el('p', 'muted small',
        'Found in the AlphaFold DB index \u2014 a few dozen ranged reads over 239.6M entries, no BLAST. '
        + 'The alignment below is Smith-Waterman, not BLAST\u2019s; identity and coverage are computed the same way.'));
    }
    p.appendChild(hitTable(res));
    if (res.skipped.length) {
      // There can be dozens of these; a few examples make the point.
      const eg = res.skipped.slice(0, 5).map(s => `${s.acc}→${s.dupOf}`).join(', ');
      p.appendChild(el('p', 'muted small',
        `Skipped ${res.skipped.length} hit(s) too similar to one already taken, before downloading them: `
        + eg + (res.skipped.length > 5 ? ', …' : '')));
    }

    const neff = Filter.neff(res.rows, 0.8);
    const stats = el('p', 'stats');
    stats.innerHTML =
      `<b>${res.rows.length}</b> sequences &middot; <b>${res.queryLen}</b> columns &middot; ` +
      `Neff@0.8 <b>${neff}</b> (${(neff / res.rows.length * 100).toFixed(0)}% effective)`;
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
    stats.innerHTML =
      `<b>${paired.nPaired}</b> paired by species &middot; <b>${total - paired.nPaired}</b> unpaired (block-diagonal) &middot; ` +
      `<b>${paired.lens.reduce((a, b) => a + b, 0)}</b> columns`;
    p.appendChild(stats);
    p.appendChild(el('p', 'muted small',
      'ColabFold layout: a "#L1,L2  1,1" header, the concatenated query, species-paired rows, then each chain’s remaining sequences padded with gaps.'));

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
    email: $('email').value.trim(),
    database: $('database').value,
    tryExact: $('tryexact').checked,
    useIndex: $('useindex').checked,
    indexUrl: $('indexurl').value.trim(),
    topK: Math.max(5, num('topk', 20)),
    minIndexIdentity: num('minindexid', 90),
    evalue: $('evalue').value.trim() || '1e-3',
    maxHits: Math.max(1, num('maxhits', 1)),
    maxPairId: num('pairid', 0.9),
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
  // Only the BLAST path needs an address; an exact-match run never contacts EBI.
  const emailOk = /^[^@\s]+@[^@\s]+\.[^@\s]+$/.test(cfg.email);
  if (!emailOk && !cfg.tryExact && !cfg.useIndex) {
    status('EBI requires a valid contact email on every BLAST job.', 'err');
    $('email').focus();
    return;
  }
  if (emailOk) localStorage.setItem('afdbmsa.email', cfg.email);

  controller = new AbortController();
  $('run').disabled = true;
  $('cancel').hidden = false;

  try {
    log(`Query: ${chains.length} chain(s), database=${cfg.database}, E<=${cfg.evalue}`);
    if (cfg.database === 'uniprotkb') log('Full UniProtKB typically takes 2-5 minutes per chain.');

    // Chains are independent; EBI tolerates a handful of concurrent jobs.
    const results = await Promise.all(chains.map((c, i) =>
      processChain(c, cfg, controller.signal, `chain ${i + 1}`)));

    renderResults(results, cfg);
    status(`Done — ${results.map(r => r.rows.length).join(' + ')} sequences.`, 'ok');
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

/* ---------- init ---------- */

document.addEventListener('DOMContentLoaded', () => {
  const saved = localStorage.getItem('afdbmsa.email');
  if (saved) $('email').value = saved;

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
