#!/usr/bin/env node
/**
 * build_cohort_link.mjs — 把 TTE 表型表链接到 Sepsis-3 ICU stays（入科 ±72h）
 * 每个 stay 选一个 TTE：优先时间上最接近 ICU 入科者。
 * 用法: node build_cohort_link.mjs sepsis3_icustays.csv tte_pheno.csv linked_cohort.csv [windowHours=72]
 */
import fs from 'node:fs';
import readline from 'node:readline';

const [stayFile, tteFile, outFile, winArg] = process.argv.slice(2);
const WIN = (parseFloat(winArg) || 72) * 3600 * 1000;
if (!stayFile || !tteFile || !outFile) { console.error('usage: node build_cohort_link.mjs <sepsis3_icustays.csv> <tte_pheno.csv> <out.csv> [windowHours]'); process.exit(1); }

function parseCsvLine(line) {
  const out = []; let cur = '', q = false;
  for (let i = 0; i < line.length; i++) { const c = line[i];
    if (q) { if (c === '"') { if (line[i + 1] === '"') { cur += '"'; i++; } else q = false; } else cur += c; }
    else if (c === '"') q = true; else if (c === ',') { out.push(cur); cur = ''; } else cur += c; }
  out.push(cur); return out;
}
const t = (s) => { const d = new Date(String(s).replace(' ', 'T')); return isNaN(d) ? null : d.getTime(); };

async function readCsv(path) {
  const rl = readline.createInterface({ input: fs.createReadStream(path), crlfDelay: Infinity });
  let head = null; const rows = [];
  for await (const line of rl) { if (!line) continue; const f = parseCsvLine(line); if (!head) { head = f.map(x => x.replace(/^\uFEFF/, '')); continue; } const o = {}; head.forEach((h, i) => o[h] = f[i]); rows.push(o); }
  return { head, rows };
}

const stays = await readCsv(stayFile);
const ttes = await readCsv(tteFile);
const tteCols = ttes.head.filter(h => !['subject_id', 'measurement_id', 'measurement_datetime'].includes(h));

// 按 subject 分组 TTE
const bySubj = new Map();
for (const r of ttes.rows) { const k = r.subject_id; if (!bySubj.has(k)) bySubj.set(k, []); bySubj.get(k).push(r); }
for (const arr of bySubj.values()) arr.sort((a, b) => t(a.measurement_datetime) - t(b.measurement_datetime));

const out = ['stay_id,subject_id,hadm_id,intime,outtime,los,measurement_id,measurement_datetime,echo_dt_minus_intime_h,' + tteCols.join(',')];
let matched = 0, noTte = 0, noWindow = 0;
for (const s of stays.rows) {
  const list = bySubj.get(s.subject_id);
  const esc = (v) => { v = String(v ?? ''); return /[",\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v; };
  if (!list || !list.length) { noTte++; continue; }
  const it = t(s.intime), ot = t(s.outtime);
  const lo = it - WIN, hi = ot + WIN;
  const cands = list.filter(r => { const d = t(r.measurement_datetime); return d >= lo && d <= hi; });
  if (!cands.length) { noWindow++; continue; }
  // 选最接近入科的
  let best = cands[0], bestD = Math.abs(t(cands[0].measurement_datetime) - it);
  for (const c of cands) { const d = Math.abs(t(c.measurement_datetime) - it); if (d < bestD) { bestD = d; best = c; } }
  const dh = ((t(best.measurement_datetime) - it) / 3600000).toFixed(1);
  out.push([s.stay_id, s.subject_id, s.hadm_id, s.intime, s.outtime, s.los, best.measurement_id, best.measurement_datetime, dh, ...tteCols.map(c => esc(best[c]))].join(','));
  matched++;
}
fs.writeFileSync(outFile, '\ufeff' + out.join('\n') + '\n');
console.log(`stays 总数: ${stays.rows.length}`);
console.log(`链接成功: ${matched}`);
console.log(`该患者无任何 TTE: ${noTte}`);
console.log(`有 TTE 但不在 ±${WIN / 3600000}h 窗内: ${noWindow}`);
console.log(`输出: ${outFile}`);
