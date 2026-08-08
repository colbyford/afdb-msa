/*
 * node tools/stamp.mjs
 *
 * Rewrites the <script> and <link> tags in index.html to carry a content hash:
 *
 *   <script src="app.js?v=8f3c1a2b">
 *
 * GitHub Pages serves everything with cache-control: max-age=600, so after a
 * deploy a browser can hold a stale app.js for ten minutes -- and pair it with a
 * fresh index.html, which is worse than either alone. That mismatch produced a
 * real crash: index.html had dropped the filter controls while the cached app.js
 * was still reading #mincov, so run() threw on null.value.
 *
 * A hash in the query string makes the URL change whenever the bytes change, so
 * a stale asset can never be paired with a new page. Run before committing.
 */
import { readFileSync, writeFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const page = join(root, 'index.html');
let html = readFileSync(page, 'utf8');

const hash = f => {
  try {
    return createHash('sha256').update(readFileSync(join(root, f))).digest('hex').slice(0, 8);
  } catch (e) {
    // A tag pointing at a file that no longer exists -- almost always a removed
    // module whose <script> tag was missed, which the browser would report only
    // as a 404 and a missing global.
    console.error(`  index.html references ${f}, which does not exist`);
    process.exit(1);
  }
};

let n = 0;
html = html.replace(/(src|href)="([\w.-]+\.(?:js|css))(?:\?v=[0-9a-f]+)?"/g, (_, attr, file) => {
  n++;
  return `${attr}="${file}?v=${hash(file)}"`;
});

writeFileSync(page, html);
console.log(`  stamped ${n} asset(s):`);
for (const m of html.matchAll(/(?:src|href)="([\w.-]+\.(?:js|css)\?v=[0-9a-f]+)"/g)) console.log(`    ${m[1]}`);
