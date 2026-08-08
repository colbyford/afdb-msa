/*
 * node test/dom.mjs
 *
 * app.js is the one file the other tests never touch. There is no browser here,
 * so this stands up a minimal DOM shim -- just enough of document/element/Blob
 * for app.js to load, wire itself up, and render a result set -- and checks that
 * nothing throws. It is a smoke test for load-time and render-time errors, not a
 * substitute for clicking the real page.
 */
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import vm from 'node:vm';

const here = dirname(fileURLToPath(import.meta.url));
const root = join(here, '..');

let pass = 0, fail = 0;
const ok = (c, m) => { if (c) { pass++; console.log(`  ok   ${m}`); } else { fail++; console.log(`  FAIL ${m}`); } };

/* ---------- the shim ---------- */

function makeEl(tag) {
  const el = {
    tagName: (tag || 'div').toUpperCase(),
    children: [], attributes: {}, dataset: {}, style: {},
    className: '', textContent: '', value: '', checked: false, hidden: false,
    disabled: false, innerHTML: '',
    appendChild(c) { this.children.push(c); return c; },
    remove() {},
    click() { this._clicked = true; },
    focus() {},
    setAttribute(k, v) { this.attributes[k] = v; },
    addEventListener() {},
    querySelectorAll() { return []; },
    scrollTop: 0, scrollHeight: 0,
  };
  return el;
}

// Build the element table from index.html's id="..." attributes, so the shim
// stays in step with the real page instead of drifting from it.
const html = readFileSync(join(root, 'index.html'), 'utf8');
const ids = [...html.matchAll(/id="([^"]+)"/g)].map(m => m[1]);
const byId = new Map(ids.map(id => [id, makeEl('div')]));

const exampleButtons = [...html.matchAll(/data-example="([^"]+)"/g)].map(m => {
  const b = makeEl('button');
  b.dataset.example = m[1];
  return b;
});

const document = {
  getElementById: id => byId.get(id) || null,
  createElement: makeEl,
  body: makeEl('body'),
  addEventListener(type, fn) { if (type === 'DOMContentLoaded') this._ready = fn; },
  querySelectorAll: sel => (sel === '[data-example]' ? exampleButtons : []),
};

const store = new Map();
const sandbox = {
  document,
  window: {},
  localStorage: {
    getItem: k => (store.has(k) ? store.get(k) : null),
    setItem: (k, v) => store.set(k, String(v)),
  },
  Blob: class { constructor(parts) { this.parts = parts; } },
  URL: { createObjectURL: () => 'blob:x', revokeObjectURL() {} },
  AbortController: class { constructor() { this.signal = { aborted: false }; } abort() { this.signal.aborted = true; } },
  fetch: async () => { throw new Error('network disabled in this smoke test'); },
  setTimeout, clearTimeout, console,
};
sandbox.window = sandbox;
sandbox.globalThis = sandbox;
vm.createContext(sandbox);

/* ---------- load the page's scripts, in page order ---------- */

// index.html carries content-hash cache busting (app.js?v=d0262328), so strip
// the query before touching the filesystem.
const scripts = [...html.matchAll(/<script src="([^"?]+)(?:\?[^"]*)?"><\/script>/g)].map(m => m[1]);
console.log(`# loading ${scripts.join(', ')}`);
let loadError = null;
try {
  for (const s of scripts) {
    vm.runInContext(readFileSync(join(root, s), 'utf8'), sandbox, { filename: s });
  }
} catch (e) { loadError = e; }
ok(!loadError, `all scripts evaluate without throwing${loadError ? ` — ${loadError.message}` : ''}`);
// Assert the set, not a count -- a count breaks every time a module is added,
// and says nothing about whether the right ones are there.
const EXPECTED = ['seeds.js', 'align.js', 'search.js', 'msa.js', 'api.js', 'app.js'];
ok(EXPECTED.every(f => scripts.includes(f)),
  `index.html loads every module (missing: ${EXPECTED.filter(f => !scripts.includes(f)).join(', ') || 'none'})`);
// app.js last: it reads the globals the others publish.
ok(scripts[scripts.length - 1] === 'app.js', `app.js is loaded last (got ${scripts[scripts.length - 1]})`);

/*
 * Every asset must carry a content hash matching its bytes. Pages serves
 * cache-control: max-age=600, so without this a browser can hold a stale
 * app.js and pair it with a fresh index.html -- which crashed for real when
 * index.html dropped the filter controls and the cached app.js kept reading
 * them. Run tools/stamp.mjs after changing any asset.
 */
{
  const { createHash } = await import('node:crypto');
  const tags = [...html.matchAll(/(?:src|href)="([\w.-]+\.(?:js|css))\?v=([0-9a-f]+)"/g)];
  const untagged = [...html.matchAll(/(?:src|href)="([\w.-]+\.(?:js|css))"/g)].map(m => m[1]);
  ok(untagged.length === 0, `every asset is stamped${untagged.length ? ' -- missing: ' + untagged.join(', ') : ''}`);
  const stale = tags.filter(([, f, v]) =>
    createHash('sha256').update(readFileSync(join(root, f))).digest('hex').slice(0, 8) !== v).map(m => m[1]);
  ok(stale.length === 0, `all ${tags.length} hashes match the files${stale.length ? ' -- stale: ' + stale.join(', ') : ''}`);
}

console.log('\n# globals published');
for (const g of ['MSAKit', 'Api', 'Seeds', 'Align', 'Search']) ok(!!sandbox[g], `window.${g} is defined`);

/* ---------- init ---------- */

console.log('\n# DOMContentLoaded');

let initError = null;
try { document._ready(); } catch (e) { initError = e; }
ok(!initError, `init runs without throwing${initError ? ` — ${initError.message}` : ''}`);
ok(byId.get('seq') !== null, 'sequence box present');
ok(typeof byId.get('run').onclick === 'function', 'Run button is wired');
ok(typeof byId.get('cancel').onclick === 'function', 'Cancel button is wired');
ok(exampleButtons.every(b => typeof b.onclick === 'function'), `all ${exampleButtons.length} example buttons are wired`);

console.log('\n# example buttons');
const hbab = exampleButtons.find(b => b.dataset.example === 'hbab');
hbab.onclick();
ok(byId.get('seq').value.startsWith('>HBA_HUMAN'), 'paired example populates the textarea');
ok(/2 chains/.test(byId.get('seq-summary').textContent), `summary reports 2 chains — "${byId.get('seq-summary').textContent}"`);
const clear = exampleButtons.find(b => b.dataset.example === 'clear');
clear.onclick();
ok(byId.get('seq').value === '', 'clear empties the textarea');

console.log('\n# validation refuses to submit bad input');
byId.get('seq').value = '';
const run = byId.get('run').onclick;
await run();
ok(/paste a sequence/i.test(byId.get('status').textContent), `empty input is rejected — "${byId.get('status').textContent}"`);

byId.get('seq').value = 'MVLSPADK123';
await run();
ok(/amino-acid/i.test(byId.get('status').textContent), `bad residues are rejected — "${byId.get('status').textContent}"`);
ok(byId.get('run').disabled === false, 'Run is re-enabled after a rejected submit');

console.log('\n# the single search path');
ok(typeof sandbox.processChain === 'function', 'processChain exists');
// byId is a Map built from index.html's id attributes, so a removed control is
// `undefined`, not null.
ok(!byId.has('mincov') && !byId.has('maxid') && !byId.has('sortid'),
  'filter controls are not exposed in the page');
ok(sandbox.Search.MinimizerIndex.length >= 1, 'MinimizerIndex takes a base URL');
const appSrc = readFileSync(join(root, 'app.js'), 'utf8');
const urlMatch = /INDEX_URL\s*=\s*'([^']+)'/.exec(appSrc);
ok(urlMatch && urlMatch[1].startsWith('https://huggingface.co'),
  `index URL is baked in: ${urlMatch ? urlMatch[1].slice(0, 56) + '…' : 'NOT FOUND'}`);
ok(sandbox.Align.smithWaterman('MVLSPADKTNVK', 'MVLSPADKTNVK').identity === 100, 'aligner is live in the page context');
// No fallback: an unreachable index must surface as an error, not a silent
// downgrade to something slower.
byId.get('seq').value = 'MVLSPADKTNVKAAWGKVGAHAGEY';
await run();
ok(/error|fail|network/i.test(byId.get('status').textContent),
  `unreachable index reports an error rather than falling back — "${byId.get('status').textContent.slice(0, 60)}"`);

console.log('\n# rendering results');
// Feed renderResults a synthetic two-chain result and check it builds a DOM and
// a downloadable a3m without touching the network.
const mkRes = (name, seq, n) => {
  const rows = [{ name, seq }];
  for (let i = 0; i < n; i++) {
    rows.push({ name: `UniRef90_X${i} Tax=Org${i % 5} TaxID=${1000 + (i % 5)}`, seq });
  }
  return {
    chain: { name, seq }, queryLen: seq.length, rows,
    donor: {
      acc: 'A0A2J8INE6', identity: 98.6, queryCoverage: 1,
      organism: 'Pan troglodytes', desc: 'Hemoglobin subunit alpha',
    },
  };
};
let renderError = null;
try {
  sandbox.renderResults([mkRes('A', 'MVLSPADKTNVK', 20), mkRes('B', 'MVHLTPEEKSAV', 20)], {});
} catch (e) { renderError = e; }
ok(!renderError, `renderResults builds two chains + pairing${renderError ? ` — ${renderError.message}` : ''}`);
const results = byId.get('results');
ok(results.hidden === false, 'results section is revealed');
ok(results.children.length === 3, `two chain panels plus a paired panel (${results.children.length})`);

// The download buttons must produce real a3m text.
const panel = results.children[0];
const btn = panel.children.find(c => c.tagName === 'BUTTON');
ok(!!btn && /Download A\.a3m/.test(btn.textContent), `chain panel has a download button — "${btn && btn.textContent}"`);
let dlError = null;
try { btn.onclick(); } catch (e) { dlError = e; }
ok(!dlError, `download handler runs${dlError ? ` — ${dlError.message}` : ''}`);

const pairedPanel = results.children[2];
const pbtn = pairedPanel.children.find(c => c.tagName === 'BUTTON');
ok(!!pbtn && /paired\.a3m/.test(pbtn.textContent), `paired panel has a download button — "${pbtn && pbtn.textContent}"`);

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
