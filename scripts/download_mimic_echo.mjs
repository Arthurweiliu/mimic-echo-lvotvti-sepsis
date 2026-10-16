#!/usr/bin/env node
/**
 * download_mimic_echo.mjs  (v3 - 稳健分块下载：断点续传 + 分块重试 + 探测大小)
 * 登录 PhysioNet（受控访问）下载 MIMIC-IV-ECHO 结构化 CSV。
 * 凭据：env PHYSIONET_USERNAME/PHYSIONET_PASSWORD 或同目录 _physionet_cred.json（绝不打印）。
 * 用法: node download_mimic_echo.mjs [输出目录] [并发=8]
 */
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const OUT_DIR = process.argv[2] || __dirname;
const CONNS = Math.max(1, parseInt(process.argv[3] || '8', 10));
const FILTER = process.argv[4] || '';
const BASE = 'https://physionet.org';
const FILES_URL = `${BASE}/files/mimic-iv-echo/1.0.1/`;
const LOGIN_URL = `${BASE}/login/`;
const UA = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36';
const REQ_TIMEOUT = 90000;

function loadCredentials() {
  let username = process.env.PHYSIONET_USERNAME, password = process.env.PHYSIONET_PASSWORD;
  if (!username || !password) {
    const p = path.join(__dirname, '_physionet_cred.json');
    if (fs.existsSync(p)) { try { const j = JSON.parse(fs.readFileSync(p, 'utf8')); username ||= j.username; password ||= j.password; } catch { } }
  }
  if (!username || !password) { console.error('未找到凭据（env 或 _physionet_cred.json）。'); process.exit(1); }
  return { username, password };
}
const jar = new Map();
const COOK = path.join(__dirname, '_physionet_cookies.json');
if (fs.existsSync(COOK)) { try { for (const [k, v] of Object.entries(JSON.parse(fs.readFileSync(COOK, 'utf8')))) jar.set(k, v); } catch { } }
function saveCookies() { try { fs.writeFileSync(COOK, JSON.stringify(Object.fromEntries(jar))); } catch { } }
function storeCookies(res) {
  const sc = res.headers.getSetCookie ? res.headers.getSetCookie() : [];
  for (const line of sc) { const [pair] = line.split(';'); const i = pair.indexOf('='); if (i > 0) jar.set(pair.slice(0, i).trim(), pair.slice(i + 1).trim()); }
}
const cookieHeader = () => [...jar.entries()].map(([k, v]) => `${k}=${v}`).join('; ');
async function req(url, opts = {}) {
  const headers = Object.assign({ 'User-Agent': UA, 'Referer': BASE + '/', ...(jar.size ? { Cookie: cookieHeader() } : {}) }, opts.headers || {});
  const ctrl = new AbortController();
  const t = setTimeout(() => ctrl.abort(), opts.timeout || REQ_TIMEOUT);
  try { return await fetch(url, { ...opts, headers, redirect: opts.redirect || 'manual', signal: ctrl.signal }); }
  finally { clearTimeout(t); }
}
async function login(username, password) {
  let last = 0;
  for (let a = 1; a <= 6; a++) {
    const r1 = await req(LOGIN_URL); const html = await r1.text();
    const m = html.match(/name="csrfmiddlewaretoken"\s+value="([^"]+)"/); const csrf = m ? m[1] : (jar.get('csrftoken') || '');
    if (!csrf) throw new Error('未取到 CSRF token。');
    const body = new URLSearchParams({ csrfmiddlewaretoken: csrf, username, password, next: '/content/mimic-iv-echo/1.0.1/' });
    const r = await req(LOGIN_URL, { method: 'POST', headers: { 'Content-Type': 'application/x-www-form-urlencoded', 'Referer': LOGIN_URL, 'Origin': BASE }, body });
    if (r.status === 302 || r.status === 301) return r;
    last = r.status;
    if (r.status === 403) { console.log(`  登录尝试 ${a} → 403（限流），等待 ${12 * a}s 重试...`); await new Promise(x => setTimeout(x, 12000 * a)); }
    else { console.log(`  登录尝试 ${a} → HTTP ${r.status}，等待重试...`); await new Promise(x => setTimeout(x, 5000 * a)); }
  }
  throw new Error('登录失败 HTTP ' + last);
}
async function assertAccess() { const r = await req(FILES_URL); if (r.status === 403) throw new Error('HTTP 403：账号未获授权/DUA 未生效。'); if (r.status >= 400) throw new Error('列目录 HTTP ' + r.status); return r.text(); }
async function probe(url) { const r = await req(url, { method: 'HEAD', redirect: 'follow' }); return { status: r.status, size: Number(r.headers.get('content-length') || 0), ranges: r.headers.get('accept-ranges') || '' }; }

async function fetchRange(url, part, s, e, need, total) {
  if (fs.existsSync(part) && fs.statSync(part).size === need) return;
  const res = await req(url, { headers: { Range: `bytes=${s}-${e}` }, redirect: 'follow' });
  if (res.status === 200 && need < total) throw new Error('server ignored Range');
  if (res.status !== 206 && res.status !== 200) throw new Error('HTTP ' + res.status);
  const tmp = part + '.tmp';
  const ws = fs.createWriteStream(tmp, { flags: 'w' });
  for await (const c of res.body) { if (!ws.write(c)) await new Promise(r => ws.once('drain', r)); }
  await new Promise(r => ws.end(r));
  if (fs.statSync(tmp).size !== need) { fs.rmSync(tmp, { force: true }); throw new Error('short read ' + fs.statSync(tmp).size + '/' + need); }
  fs.renameSync(tmp, part);
}
async function merge(parts, dest) {
  const out = fs.createWriteStream(dest, { flags: 'w' });
  for (const p of parts) { if (!fs.existsSync(p)) continue; await new Promise((res, rej) => { const rs = fs.createReadStream(p); rs.on('error', rej); rs.on('end', res); rs.pipe(out, { end: false }); }); }
  await new Promise(r => out.end(r));
}
async function rangedDownload(url, dest, size) {
  const partDir = dest + '.parts'; fs.mkdirSync(partDir, { recursive: true });
  const chunk = Math.ceil(size / CONNS); const tasks = [];
  for (let i = 0; i < CONNS; i++) { const s = i * chunk; if (s >= size) break; const e = Math.min(s + chunk - 1, size - 1); tasks.push({ s, e, need: e - s + 1, p: path.join(partDir, 'p' + i) }); }
  const ok = t => fs.existsSync(t.p) && fs.statSync(t.p).size === t.need;
  const bytes = () => tasks.reduce((a, t) => a + (fs.existsSync(t.p) ? fs.statSync(t.p).size : 0), 0);
  const base = path.basename(dest);
  const timer = setInterval(() => process.stdout.write(`\r  ${base}  ${(bytes() / size * 100).toFixed(1)}%  (${(bytes() / 1048576).toFixed(1)}/${(size / 1048576).toFixed(1)} MB, ${tasks.filter(ok).length}/${tasks.length} 块)   `), 3000);
  let round = 0;
  while (tasks.some(t => !ok(t)) && round < 60) {
    round++;
    const pending = tasks.filter(t => !ok(t)); let idx = 0;
    await Promise.all(Array.from({ length: Math.min(CONNS, pending.length) }, async () => {
      while (idx < pending.length) { const t = pending[idx++]; try { await fetchRange(url, t.p, t.s, t.e, t.need, size); } catch { /* 下一轮重试 */ } }
    }));
    if (tasks.some(t => !ok(t))) await new Promise(r => setTimeout(r, 2500));
  }
  clearInterval(timer); process.stdout.write('\n');
  if (tasks.some(t => !ok(t))) { console.log(`  未完成 ${tasks.filter(t => !ok(t)).length} 块（网络反复中断），已保留分块可续传，稍后重跑本命令即可。`); return false; }
  await merge(tasks.map(t => t.p), dest); fs.rmSync(partDir, { recursive: true, force: true });
  console.log(`完成: ${dest} (${(fs.statSync(dest).size / 1048576).toFixed(1)} MB)`);
  return true;
}
async function singleDownload(url, dest) {
  let start = fs.existsSync(dest) ? fs.statSync(dest).size : 0;
  const res = await req(url, { headers: start > 0 ? { Range: `bytes=${start}-` } : {}, redirect: 'follow' });
  if (res.status !== 200 && res.status !== 206) throw new Error('HTTP ' + res.status);
  const ws = fs.createWriteStream(dest, { flags: res.status === 206 ? 'a' : 'w' });
  for await (const c of res.body) { if (!ws.write(c)) await new Promise(r => ws.once('drain', r)); }
  await new Promise(r => ws.end(r));
  console.log(`完成: ${dest} (${(fs.statSync(dest).size / 1048576).toFixed(1)} MB)`);
  return true;
}
async function download(url, dest, size, ranges) {
  if (size > 4 * 1048576 && /bytes/i.test(ranges)) return rangedDownload(url, dest, size);
  return singleDownload(url, dest);
}
(async () => {
  const creds = loadCredentials();
  let indexHtml = null;
  if (jar.size) { const t = await req(FILES_URL); console.log('复用已存会话访问目录:', t.status); if (t.status === 200) indexHtml = await t.text(); }
  if (!indexHtml) {
    console.log(`登录中（用户 ${creds.username.replace(/(.).*(@.*)/, '$1***$2')}）...`);
    const r = await login(creds.username, creds.password);
    console.log(`登录响应: HTTP ${r.status}${r.status < 400 ? ' → ' + (r.headers.get('location') || '') : ''}`);
    saveCookies();
    indexHtml = await assertAccess();
    saveCookies();
  }
  const links = [...indexHtml.matchAll(/href="([^"]+\.(?:csv|csv\.gz))"/gi)].map(m => m[1]);
  const uniq = [...new Set(links.filter(h => /structured[-_]measurement|echo[-_]?(record|study)[-_]?list|record[-_]?list|study[-_]?list/i.test(h)))].filter(h => !FILTER || h.toLowerCase().includes(FILTER.toLowerCase()));
  if (!uniq.length) { console.log('未发现 CSV 链接。'); process.exit(3); }
  fs.mkdirSync(OUT_DIR, { recursive: true });
  const targets = [];
  for (const h of uniq) {
    const url = new URL(h, FILES_URL).href; const dest = path.join(OUT_DIR, decodeURIComponent(h.split('/').pop()));
    const p = await probe(url).catch(() => ({ size: 0, ranges: '' }));
    targets.push({ url, dest, size: p.size, ranges: p.ranges });
    console.log(`-> ${path.basename(dest)}  ≈ ${(p.size / 1048576).toFixed(1)} MB  分段=${/bytes/i.test(p.ranges)}`);
  }
  let all = true;
  for (const t of targets) { try { if (!await download(t.url, t.dest, t.size, t.ranges)) all = false; } catch (e) { all = false; console.log(`  ${path.basename(t.dest)} 出错: ${e.message}`); } }
  console.log(all ? '全部完成。' : '部分文件未完成（可重跑续传）。');
  process.exit(all ? 0 : 4);
})().catch(e => { console.error('错误: ' + e.message); process.exit(2); });
