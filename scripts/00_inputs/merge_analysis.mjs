#!/usr/bin/env node
/** merge_analysis.mjs — 合并 linked_cohort.csv(超声) 与 cohort_covariates.csv(协变量) → analysis_data.csv */
import fs from 'node:fs';
import readline from 'node:readline';
import { fileURLToPath } from 'node:url';
const DIR = (process.argv[2] || fileURLToPath(new URL('.', import.meta.url))).replace(/[\\/]?$/, '/');
function parseCsvLine(line) { const out = []; let cur = '', q = false;
  for (let i = 0; i < line.length; i++) { const c = line[i];
    if (q) { if (c === '"') { if (line[i + 1] === '"') { cur += '"'; i++; } else q = false; } else cur += c; }
    else if (c === '"') q = true; else if (c === ',') { out.push(cur); cur = ''; } else cur += c; }
  out.push(cur); return out; }
async function readCsv(p) { const rl = readline.createInterface({ input: fs.createReadStream(p), crlfDelay: Infinity });
  let head = null; const rows = [];
  for await (const l of rl) { if (!l) continue; const f = parseCsvLine(l); if (!head) { head = f.map(x => x.replace(/^\uFEFF/, '')); continue; } const o = {}; head.forEach((h, i) => o[h] = f[i]); rows.push(o); } return { head, rows }; }
const t = s => { if (!s) return null; const d = new Date(String(s).replace(' ', 'T')); return isNaN(d) ? null : d; };

const echo = await readCsv(DIR + 'linked_cohort.csv');
const cov = await readCsv(DIR + 'cohort_covariates.csv');
const covById = new Map(cov.rows.map(r => [r.stay_id, r]));
const echoCols = echo.head.filter(h => !['stay_id', 'subject_id', 'hadm_id', 'intime', 'outtime', 'los'].includes(h));

const covCols = cov.head.filter(h => !['stay_id', 'subject_id', 'hadm_id'].includes(h));
const out = ['stay_id,subject_id,hadm_id,death28,death90,days_to_event,event_observed,' + covCols.join(',') + ',' + echoCols.join(',')];
const num = v => { if (v == null || v === '') return null; const x = Number(String(v).replace(/[^0-9.\-]/g, '')); return isNaN(x) ? null : x; };

let n = 0, d28 = 0, d90 = 0, lvefN = 0, tapseN = 0, vtiN = 0;
for (const e of echo.rows) {
  const c = covById.get(e.stay_id); if (!c) continue;
  const intime = t(c.icu_intime), disch = t(c.dischtime), dod = t(c.dod);
  // 存活分析：事件=死亡；如死亡时间晚于最后出院则用 dod，否则用出院时间截尾
  let death = 0, ev = t(''); // placeholder
  const hosp = c.hospital_expire_flag === '1' || c.hospital_expire_flag === 't';
  let deathDate = dod;
  if (!deathDate && hosp) deathDate = disch;
  let dt = null, observed = 0;
  if (deathDate && intime) { dt = (deathDate - intime) / 86400000; if (dt <= 28) d28++; if (dt <= 90) d90++; }
  if (dt != null && dt <= 90) observed = 1; // 90 天内观察到事件或截尾
  const death28 = (dt != null && dt <= 28) ? 1 : 0;
  const death90 = (dt != null && dt <= 90) ? 1 : 0;
  if (death28) { /* counted */ }
  if (num(e.lvef) != null) lvefN++;
  if (num(e.tapse) != null) tapseN++;
  if (num(e.lvot_vti) != null) vtiN++;
  const esc = v => { v = String(v ?? ''); return /[",\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v; };
  out.push([e.stay_id, e.subject_id, e.hadm_id, death28, death90, dt != null ? dt.toFixed(2) : '', observed].concat(covCols.map(k => esc(c[k])), echoCols.map(k => esc(e[k]))).join(','));
  n++;
}
fs.writeFileSync(DIR + 'analysis_data.csv', '\ufeff' + out.join('\n') + '\n');
console.log(`合并 N = ${n}`);
console.log(`28天死亡 = ${d28} (${(d28 / n * 100).toFixed(1)}%)`);
console.log(`90天死亡 = ${d90} (${(d90 / n * 100).toFixed(1)}%)`);
console.log(`有 LVEF 值: ${lvefN} (${(lvefN / n * 100).toFixed(1)}%) | TAPSE: ${tapseN} | LVOT VTI: ${vtiN}`);
