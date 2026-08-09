/*
 * node test/server-smoke.mjs
 *
 * Verifies the HTTP API starts and that input validation returns the expected
 * 400 JSON error for a missing sequence.
 */
import { spawn } from 'node:child_process';
import { setTimeout as delay } from 'node:timers/promises';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const here = dirname(fileURLToPath(import.meta.url));
const root = join(here, '..');
const PORT = 8080;
const BASE = `http://127.0.0.1:${PORT}`;

let pass = 0, fail = 0;
const ok = (c, m) => { if (c) { pass++; console.log(`  ok   ${m}`); } else { fail++; console.log(`  FAIL ${m}`); } };
const eq = (a, b, m) => ok(a === b, `${m} (got ${JSON.stringify(a)}, want ${JSON.stringify(b)})`);

async function waitForServer(url) {
  for (let i = 0; i < 20; i++) {
    if (server.exitCode !== null) throw new Error(`server exited early:\n${logs.trim()}`);
    try {
      const res = await fetch(url);
      await res.arrayBuffer();
      return;
    } catch {}
    await delay(200);
  }
  throw new Error('server did not start in time');
}

const server = spawn(process.execPath, ['server.mjs', '--port', String(PORT)], {
  cwd: root,
  stdio: ['ignore', 'pipe', 'pipe'],
});
const exited = new Promise(resolve => server.once('exit', resolve));

let logs = '';
server.stdout.on('data', chunk => { logs += chunk.toString('utf8'); });
server.stderr.on('data', chunk => { logs += chunk.toString('utf8'); });

try {
  await waitForServer(`${BASE}/`);

  console.log('# root endpoint');
  {
    const res = await fetch(`${BASE}/`);
    const text = await res.text();
    eq(res.status, 200, 'GET / returns 200');
    ok(text.includes('AFDB MSA API'), 'GET / returns the API help text');
  }

  console.log('\n# missing sequence validation');
  {
    const res = await fetch(`${BASE}/api?seq=`);
    const body = await res.json();
    eq(res.status, 400, 'empty seq returns 400');
    eq(body.error, 'Missing sequence. Provide ?seq= or POST the sequence in the request body.',
      'empty seq returns the documented error message');
  }
} finally {
  if (server.exitCode === null) server.kill('SIGTERM');
  await exited;
}

if (logs.trim()) console.log(`\n# server log\n${logs.trim()}`);
console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
