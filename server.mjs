/*
 * server.mjs -- HTTP API server for AFDB MSA.
 *
 * Exposes the same pipeline as the browser page via a simple HTTP endpoint.
 * Accepts a protein sequence query and returns a3m alignment file(s).
 *
 * Usage:
 *   node server.mjs [--port 8080]
 *
 * API:
 *   GET /api?seq=<sequence>
 *
 *   <sequence> may be:
 *     - A bare amino-acid sequence:       MVLSPADKTNVK...
 *     - Colon-separated multi-chain:      MVLSPA...:MVHLT...
 *     - FASTA (URL-encoded):              %3EHBA_HUMAN%0AMVLSPA...
 *     - POST with body in any of the above formats (Content-Type: text/plain)
 *
 *   Single chain  -> returns the chain's a3m (Content-Type: text/plain; name <chain>.a3m)
 *   Multi-chain   -> returns the paired a3m  (Content-Type: text/plain; name paired.a3m)
 *
 *   Query parameters:
 *     seq      required  Sequence(s) as described above
 *     format   optional  "single" | "paired" | "all"
 *                        single: return only the first/only chain's a3m (default for 1 chain)
 *                        paired: return the paired a3m (default for >1 chain)
 *                        all:    return a JSON object with one a3m per chain + the paired a3m
 *
 * Errors are returned as JSON: { error: "message" }
 */

import { createRequire } from 'node:module';
import { createServer } from 'node:http';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { request as httpRequest } from 'node:http';
import { request as httpsRequest } from 'node:https';

function abortError() {
  const e = new Error('The operation was aborted.');
  e.name = 'AbortError';
  return e;
}

function normalizeHeaders(headers = {}) {
  if (!headers) return {};
  if (typeof headers.entries === 'function') return Object.fromEntries(headers.entries());
  return { ...headers };
}

function makeCompatResponse(status, headersObj, data) {
  return {
    ok: status >= 200 && status < 300,
    status,
    headers: {
      get(name) {
        const v = headersObj[String(name).toLowerCase()];
        if (Array.isArray(v)) return v.join(', ');
        return v == null ? null : String(v);
      },
    },
    body: null,
    async text() { return data.toString('utf8'); },
    async arrayBuffer() {
      const view = data.buffer.slice(data.byteOffset, data.byteOffset + data.byteLength);
      return view;
    },
    async json() { return JSON.parse(data.toString('utf8')); },
  };
}

function fetchCompatNode(url, options = {}, redirects = 0) {
  return new Promise((resolve, reject) => {
    const u = new URL(url);
    const reqFn = u.protocol === 'https:' ? httpsRequest : httpRequest;
    const method = options.method || 'GET';
    const headers = normalizeHeaders(options.headers);
    const signal = options.signal;

    if (signal && signal.aborted) {
      reject(abortError());
      return;
    }

    const req = reqFn(
      u,
      { method, headers },
      res => {
        const status = res.statusCode || 0;
        const location = res.headers.location;
        if (location && [301, 302, 303, 307, 308].includes(status) && redirects < 5) {
          res.resume();
          const nextUrl = new URL(location, u).toString();
          const nextOptions = { ...options };
          if (status === 303) {
            nextOptions.method = 'GET';
            delete nextOptions.body;
          }
          resolve(fetchCompatNode(nextUrl, nextOptions, redirects + 1));
          return;
        }

        const chunks = [];
        res.on('data', c => chunks.push(c));
        res.on('end', () => {
          resolve(makeCompatResponse(status, res.headers, Buffer.concat(chunks)));
        });
      }
    );

    req.on('error', reject);

    if (signal) {
      const onAbort = () => req.destroy(abortError());
      signal.addEventListener('abort', onAbort, { once: true });
      req.on('close', () => signal.removeEventListener('abort', onAbort));
    }

    if (options.body != null) req.write(options.body);
    req.end();
  });
}

async function ensureFetch() {
  if (typeof globalThis.fetch === 'function') return;

  try {
    const undici = await import('undici');
    for (const k of ['fetch', 'Headers', 'Request', 'Response', 'FormData', 'File', 'Blob']) {
      if (typeof globalThis[k] === 'undefined' && undici[k]) {
        globalThis[k] = undici[k];
      }
    }
  } catch {
    globalThis.fetch = fetchCompatNode;
  }
}

await ensureFetch();

const here = dirname(fileURLToPath(import.meta.url));
const require = createRequire(import.meta.url);

// Load the shared modules (all support CommonJS exports)
const Seeds  = require(join(here, 'seeds.js'));
const Align  = require(join(here, 'align.js'));
const MSAKit = require(join(here, 'msa.js'));
const Api    = require(join(here, 'api.js'));
const Search = require(join(here, 'search.js'));

// Make globals available so the modules that check `typeof X` find them.
global.Seeds  = Seeds;
global.Align  = Align;
global.MSAKit = MSAKit;
global.Api    = Api;
global.Search = Search;

const INDEX_URL      = 'https://huggingface.co/datasets/sokrypton/afdb-msa-index/resolve/main';
const TOP_CANDIDATES = 20;

/* ---------- pipeline (mirrors app.js processChain) ---------- */

let sharedIndex = null;

async function processChain(chain, signal, tag) {
  if (!sharedIndex) sharedIndex = new Search.MinimizerIndex(INDEX_URL);
  await sharedIndex.load();

  const cands = await sharedIndex.search(chain.seq, { topK: TOP_CANDIDATES, signal });
  if (!cands.length) {
    throw new Error(
      `Nothing in AlphaFold DB resembles ${chain.name}. ` +
      `A sequence with no natural homologs has no alignment to borrow.`
    );
  }

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
  if (!scored.length) {
    throw new Error(`Could not retrieve sequences for any candidate of ${chain.name}.`);
  }

  let donor = null;
  for (const h of scored) {
    try {
      h.a3mText = await Api.fetchMsaCached(h.acc, null, signal);
      donor = h;
      break;
    } catch (e) {
      if (e.name !== 'NoMsaError') throw e;
    }
  }
  if (!donor) {
    throw new Error(`None of the candidates for ${chain.name} has an AlphaFold DB MSA.`);
  }

  const built = MSAKit.buildChainA3M(chain, [donor]);
  return { chain, donor, rows: built.rows, queryLen: chain.seq.length };
}

/* ---------- request body ---------- */

function readBody(req) {
  return new Promise((resolve, reject) => {
    const chunks = [];
    req.on('data', c => chunks.push(c));
    req.on('end', () => resolve(Buffer.concat(chunks).toString('utf8')));
    req.on('error', reject);
  });
}

/* ---------- HTTP handler ---------- */

function jsonError(res, status, message) {
  res.writeHead(status, { 'Content-Type': 'application/json' });
  res.end(JSON.stringify({ error: message }));
}

async function handleMsa(req, res) {
  const url = new URL(req.url, `http://${req.headers.host || 'localhost'}`);
  const ac = new AbortController();
  req.on('close', () => ac.abort());

  let seqText = url.searchParams.get('seq') || '';
  const format = url.searchParams.get('format') || '';  // "single" | "paired" | "all"

  // Allow POST with the sequence in the body
  if (req.method === 'POST' && !seqText) {
    seqText = await readBody(req);
  }

  seqText = seqText.trim();
  if (!seqText) {
    return jsonError(res, 400, 'Missing sequence. Provide ?seq= or POST the sequence in the request body.');
  }

  let chains;
  try {
    chains = MSAKit.parseQueryInput(seqText);
  } catch (e) {
    return jsonError(res, 400, e.message);
  }
  if (!chains.length) {
    return jsonError(res, 400, 'No valid sequences found in input.');
  }

  let results;
  try {
    results = await Promise.all(
      chains.map((c, i) => processChain(c, ac.signal, `chain ${i + 1}`))
    );
  } catch (e) {
    if (e.name === 'AbortError') {
      return jsonError(res, 499, 'Request cancelled.');
    }
    return jsonError(res, 502, e.message);
  }

  const wantsAll = format === 'all';
  const wantsPaired = format === 'paired';
  const wantsJson = wantsAll || (!format && chains.length > 1);

  if (wantsJson) {
    // Return JSON containing each chain's a3m and, when multi-chain, the paired a3m.
    const out = {};
    for (const r of results) {
      out[r.chain.name] = MSAKit.formatA3M(r.rows);
    }
    if (results.length > 1) {
      const chainData = results.map(r => ({ queryLen: r.queryLen, rows: r.rows }));
      const paired = MSAKit.pairChains(chainData);
      out['Paired Alignment'] = MSAKit.formatPairedA3M(paired);
    }
    res.writeHead(200, { 'Content-Type': 'application/json' });
    res.end(JSON.stringify(out));
  } else if (wantsPaired && results.length > 1) {
    const chainData = results.map(r => ({ queryLen: r.queryLen, rows: r.rows }));
    const paired = MSAKit.pairChains(chainData);
    const text = MSAKit.formatPairedA3M(paired);
    res.writeHead(200, {
      'Content-Type': 'text/plain; charset=utf-8',
      'Content-Disposition': 'attachment; filename="paired.a3m"',
    });
    res.end(text);
  } else {
    // Single chain (or forced single with format=single)
    const r = results[0];
    const text = MSAKit.formatA3M(r.rows);
    const filename = `${r.chain.name.replace(/[^a-zA-Z0-9_.-]/g, '_')}.a3m`;
    res.writeHead(200, {
      'Content-Type': 'text/plain; charset=utf-8',
      'Content-Disposition': `attachment; filename="${filename}"`,
    });
    res.end(text);
  }
}

/* ---------- router ---------- */

async function handler(req, res) {
  const path = new URL(req.url, 'http://localhost').pathname;

  if ((req.method === 'GET' || req.method === 'POST') && (path === '/api' || path === '/api/msa')) {
    return handleMsa(req, res);
  }

  if (req.method === 'GET' && path === '/') {
    res.writeHead(200, { 'Content-Type': 'text/plain' });
    res.end(
      'AFDB MSA API\n\n' +
      'GET  /api?seq=<sequence>\n' +
      'POST /api  (body = sequence, Content-Type: text/plain)\n\n' +
      '<sequence> formats:\n' +
      '  - Bare amino acids:         MVLSPADKTNVK...\n' +
      '  - Colon-separated chains:   MVLSPA...:MVHLT...\n' +
      '  - FASTA:                    >Name\\nMVLSPA...\\n>Name2\\nMVHLT...\n\n' +
      'Optional query params:\n' +
      '  format=single   return the first chain\'s a3m (default for 1 chain)\n' +
      '  format=paired   return only the species-paired a3m\n' +
      '  format=all      return JSON with every chain\'s a3m + paired\n' +
      '  default (>1 chain) returns JSON with one a3m per chain + paired\n'
    );
    return;
  }

  res.writeHead(404, { 'Content-Type': 'application/json' });
  res.end(JSON.stringify({ error: 'Not found. Try GET /api?seq=<sequence>' }));
}

/* ---------- main ---------- */

const args = process.argv.slice(2);
const portIdx = args.indexOf('--port');
const PORT = portIdx !== -1 ? Number(args[portIdx + 1]) : Number(process.env.PORT || 8080);

const server = createServer((req, res) => {
  handler(req, res).catch(err => {
    console.error('Unhandled error:', err);
    if (!res.headersSent) {
      res.writeHead(500, { 'Content-Type': 'application/json' });
      res.end(JSON.stringify({ error: 'Internal server error.' }));
    }
  });
});

server.listen(PORT, () => {
  console.log(`AFDB MSA API listening on http://localhost:${PORT}`);
  console.log(`  GET  http://localhost:${PORT}/api?seq=MVLSPADKTNVKAA...`);
  console.log(`  POST http://localhost:${PORT}/api  (body = sequence)`);
});
