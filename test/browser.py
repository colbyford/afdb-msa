"""
python3 test/browser.py [indexDir]

The page in a real browser, driving the real index. Every other test runs the
same JavaScript under node against a DOM shim -- which cannot catch a CSS
mistake, a CORS preflight, a fetch that behaves differently under a real origin,
or the page simply not rendering.

Serves the site and the index from one origin with Range support, then drives
Chromium: fill in a sequence, press the button, wait for the a3m, and read back
what the user would actually see. Console errors and failed requests are fatal.
"""
import http.server, socketserver, threading, os, sys, re, json, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INDEX = sys.argv[1] if len(sys.argv) > 1 else None

class H(http.server.SimpleHTTPRequestHandler):
    def translate_path(self, path):
        p = path.split('?')[0].lstrip('/')
        if p.startswith('idx/') and INDEX:
            return os.path.join(INDEX, p[4:])
        return os.path.join(ROOT, p)

    def end_headers(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Accept-Ranges', 'bytes')
        super().end_headers()

    def do_GET(self):
        path = self.translate_path(self.path)
        rng = self.headers.get('Range')
        if not rng or not os.path.isfile(path):
            return super().do_GET()
        size = os.path.getsize(path)
        m = re.match(r'bytes=(\d+)-(\d*)', rng)
        a = int(m.group(1)); b = int(m.group(2)) if m.group(2) else size - 1
        b = min(b, size - 1); n = b - a + 1
        with open(path, 'rb') as f:
            f.seek(a); data = f.read(n)
        self.send_response(206)
        self.send_header('Content-Range', f'bytes {a}-{b}/{size}')
        self.send_header('Content-Length', str(len(data)))
        self.send_header('Content-Type', 'application/octet-stream')
        self.end_headers()
        self.wfile.write(data)

    def log_message(self, *a):
        pass

socketserver.TCPServer.allow_reuse_address = True
srv = socketserver.ThreadingTCPServer(('127.0.0.1', 0), H)
port = srv.server_address[1]
threading.Thread(target=srv.serve_forever, daemon=True).start()
base = f'http://127.0.0.1:{port}/'
print(f'serving {ROOT} at {base}' + (f' (index at /idx/)' if INDEX else ' (no index)'))

from playwright.sync_api import sync_playwright

HBA = ('MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNA'
       'LSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR')

# A variant, so the exact-match shortcut cannot fire and the index path is the
# one under test. Using HBA itself resolves by UniParc checksum in 1.5s and never
# touches the index at all.
HBA_VARIANT = ('MVLSAADKSNVKATWDKIGSHAGDYGGEALDRTFQSFPTTKTYFPHFDLSHGSAQVKAHGKKVAAALVEAVNHIDDIAGA'
               'LSKLSDLHAQKLRVDPVNFKLLGQCFLVVVAIHHPSALTPEVHASLDKFLCAVGNVLTAKYR')

# HBA with ~5% of positions changed: not an exact UniProt entry, so the checksum
# shortcut cannot fire, but close enough that a >=90% relative certainly exists.
HBA_MUT = ('MVLSPADKTNVKAAWGKVGAHAGEWGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNA'
           'LSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKFR')

npass = nfail = 0
def ok(cond, msg):
    global npass, nfail
    if cond: npass += 1; print(f'  ok   {msg}')
    else:    nfail += 1; print(f'  FAIL {msg}')

with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()
    errors, failed = [], []
    page.on('console', lambda m: errors.append(m.text) if m.type == 'error' else None)
    page.on('pageerror', lambda e: errors.append(str(e)))
    page.on('requestfailed', lambda r: failed.append(f'{r.url} {r.failure}'))

    print('\n# page loads')
    page.goto(base + 'index.html', wait_until='networkidle')
    ok(page.title() != '', f'title: "{page.title()}"')
    ok(not errors, f'no console errors{"" if not errors else ": " + errors[0][:120]}')
    ok(not failed, f'no failed requests{"" if not failed else ": " + failed[0][:120]}')

    print('\n# modules reached the page')
    for g in ['MSAKit', 'Filter', 'Api', 'Seeds', 'Align', 'Search']:
        ok(page.evaluate(f'typeof window.{g} !== "undefined"'), f'window.{g} defined')

    print('\n# it renders')
    ok(page.locator('#seq').is_visible(), 'sequence box visible')
    ok(page.locator('#run').is_visible(), 'Run button visible')
    h = page.evaluate('document.body.scrollHeight')
    ok(h > 400, f'page has real height ({h}px)')
    ok(page.evaluate('getComputedStyle(document.body).backgroundColor') != 'rgba(0, 0, 0, 0)',
       'stylesheet applied')
    w = page.evaluate('document.documentElement.scrollWidth <= window.innerWidth + 2')
    ok(w, 'no horizontal overflow')

    print('\n# examples wire up')
    page.click('button[data-example="hbab"]')
    ok('HBA_HUMAN' in page.input_value('#seq'), 'example fills the textarea')
    ok('2 chains' in page.text_content('#seq-summary'), f'summary: {page.text_content("#seq-summary")[:70]}')

    print('\n# validation')
    page.fill('#seq', '')
    page.click('#run'); page.wait_for_timeout(400)
    ok('aste a sequence' in page.text_content('#status'), f'empty rejected: "{page.text_content("#status")}"')
    page.fill('#seq', 'MVLSPADK123')
    page.click('#run'); page.wait_for_timeout(400)
    ok('amino-acid' in page.text_content('#status'), f'bad residues rejected')

    if INDEX:
        print('\n# exact-match path (no index needed)')
        page.fill('#indexurl', '')
        page.fill('#seq', HBA)
        page.fill('#email', '')
        t0 = time.time()
        page.click('#run')
        page.wait_for_selector('#results button.primary', timeout=180000)
        ok(True, f'exact match resolved in {time.time()-t0:.1f}s')
        ok('exact' in (page.text_content('#log') or '').lower(), 'log shows the exact-match path')

        print('\n# a real query through the INDEX (no exact match possible)')
        page.fill('#indexurl', base + 'idx')
        page.uncheck('#tryexact')
        page.fill('#seq', HBA_MUT)
        page.fill('#email', '')
        t0 = time.time()
        page.click('#run')
        try:
            page.wait_for_selector('#results button.primary', timeout=180000)
        except Exception:
            print('  status:', page.text_content('#status')[:160])
            print('  log:', (page.text_content('#log') or '')[:600])
            ok(False, 'a3m produced')
        else:
            dt = time.time() - t0
            label = page.text_content('#results button.primary')
            stats = page.text_content('#results .stats') or ''
            ok(True, f'a3m produced in {dt:.1f}s -- button reads "{label}"')
            ok('MB' in label or 'KB' in label, 'download button shows a real size')
            n = re.search(r'(\d+)\s*sequences', stats)
            ok(n and int(n.group(1)) > 100, f'stats: {stats.strip()[:90]}')
            log = page.text_content('#log') or ''
            ok('index:' in log.lower(), 'log shows the index path was used')
            ok('blast' not in log.lower() or 'falling back' not in log.lower(),
               'did not fall back to BLAST')
            m = re.search(r'KB read, (\d+) ms', log)
            if m: print(f'       index lookup: {m.group(0)}')

    ok(not errors, f'still no console errors{"" if not errors else ": " + errors[-1][:140]}')
    page.screenshot(path='/tmp/afdbmsa.png', full_page=True)
    print(f'\n  screenshot -> /tmp/afdbmsa.png')
    browser.close()

srv.shutdown()
print(f'\n{npass} passed, {nfail} failed')
sys.exit(1 if nfail else 0)
