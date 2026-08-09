/*
 * sw.js -- browser-only API shim for static hosting.
 *
 * Intercepts /api/msa requests and forwards them to an open client page,
 * which runs the same pipeline in app.js and returns a response payload.
 */

self.addEventListener('install', () => {
  self.skipWaiting();
});

self.addEventListener('activate', event => {
  event.waitUntil(self.clients.claim());
});

function isApiPath(url) {
  return url.pathname.endsWith('/api/msa');
}

function json(status, payload) {
  return new Response(JSON.stringify(payload), {
    status,
    headers: {
      'content-type': 'application/json; charset=utf-8',
      'access-control-allow-origin': '*',
      'access-control-allow-methods': 'GET,POST,OPTIONS',
      'access-control-allow-headers': 'content-type',
    },
  });
}

async function chooseClient(event) {
  if (event.clientId) {
    const client = await self.clients.get(event.clientId);
    if (client) return client;
  }
  const list = await self.clients.matchAll({ type: 'window', includeUncontrolled: true });
  return list[0] || null;
}

function callClient(client, request) {
  return new Promise((resolve, reject) => {
    const channel = new MessageChannel();
    const timeout = setTimeout(() => reject(new Error('Client timed out.')), 300000);

    channel.port1.onmessage = event => {
      clearTimeout(timeout);
      const msg = event.data || {};
      if (!msg.ok || !msg.response) {
        reject(new Error(msg.error || 'Client failed to produce a response.'));
        return;
      }
      resolve(msg.response);
    };

    client.postMessage({ type: 'afdb-api-request', request }, [channel.port2]);
  });
}

self.addEventListener('fetch', event => {
  const req = event.request;
  const url = new URL(req.url);
  if (url.origin !== self.location.origin || !isApiPath(url)) return;

  event.respondWith((async () => {
    const client = await chooseClient(event);
    if (!client) {
      return json(503, {
        error: 'Open the AFDB MSA page in a tab first. The browser API runs in-page on static hosting.',
      });
    }

    let body = '';
    if (req.method === 'POST') body = await req.text();

    try {
      const out = await callClient(client, {
        method: req.method,
        url: req.url,
        body,
      });
      return new Response(out.body || '', {
        status: out.status || 200,
        headers: out.headers || { 'content-type': 'text/plain; charset=utf-8' },
      });
    } catch (e) {
      return json(502, { error: e && e.message ? e.message : 'API bridge failed.' });
    }
  })());
});
