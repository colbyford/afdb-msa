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

const scripts = [...html.matchAll(/<script src="([^"]+)"><\/script>/g)].map(m => m[1]);
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
const EXPECTED = ['seeds.js', 'align.js', 'search.js', 'msa.js', 'filter.js', 'api.js', 'app.js'];
ok(EXPECTED.every(f => scripts.includes(f)),
  `index.html loads every module (missing: ${EXPECTED.filter(f => !scripts.includes(f)).join(', ') || 'none'})`);
// app.js last: it reads the globals the others publish.
ok(scripts[scripts.length - 1] === 'app.js', `app.js is loaded last (got ${scripts[scripts.length - 1]})`);

console.log('\n# globals published');
for (const g of ['MSAKit', 'Filter', 'Api', 'Seeds', 'Align', 'Search']) ok(!!sandbox[g], `window.${g} is defined`);

/* ---------- init ---------- */

console.log('\n# DOMContentLoaded');
store.set('afdbmsa.email', 'saved@example.org');
let initError = null;
try { document._ready(); } catch (e) { initError = e; }
ok(!initError, `init runs without throwing${initError ? ` — ${initError.message}` : ''}`);
ok(byId.get('email').value === 'saved@example.org', 'saved email is restored from localStorage');
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

byId.get('seq').value = 'MVLSPADKTNVK';
byId.get('email').value = 'not-an-email';
await run();
ok(/email/i.test(byId.get('status').textContent), `bad email is rejected — "${byId.get('status').textContent}"`);
ok(byId.get('run').disabled === false, 'Run is re-enabled after a rejected submit');

byId.get('seq').value = 'MVLSPADK123';
byId.get('email').value = 'a@b.co';
await run();
ok(/amino-acid/i.test(byId.get('status').textContent), `bad residues are rejected — "${byId.get('status').textContent}"`);

console.log('\n# index path wiring');
// The index must be optional: with no index URL configured, tryIndex has to bow
// out immediately rather than throw, so BLAST still runs.
ok(typeof sandbox.tryIndex === 'function', 'tryIndex exists');
byId.get('indexurl').value = '';
ok(await sandbox.tryIndex({ name: 'q', seq: 'MVLSPADK' }, { indexUrl: '' }, null, 't') === null,
  'no index URL configured -> falls through to BLAST');
byId.get('indexurl').value = 'data/index';
// fetch throws in this shim, standing in for a site published without an index.
ok(await sandbox.tryIndex({ name: 'q', seq: 'MVLSPADK' }, { indexUrl: 'data/index' }, null, 't') === null,
  'unreachable index -> falls through to BLAST rather than failing the run');
ok(sandbox.Search.MinimizerIndex.length >= 1, 'MinimizerIndex takes a base URL');
ok(sandbox.Align.smithWaterman('MVLSPADKTNVK', 'MVLSPADKTNVK').identity === 100, 'aligner is live in the page context');

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
    hits: [{ acc: 'P69905', db: 'SP', identity: 99.3, queryCoverage: 1, evalue: 1e-90, organism: 'Homo sapiens' }],
    built: { stats: [{ acc: 'P69905', added: n, total: n }] },
    skipped: [{ acc: 'P69907', dupOf: 'P69905' }],
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
