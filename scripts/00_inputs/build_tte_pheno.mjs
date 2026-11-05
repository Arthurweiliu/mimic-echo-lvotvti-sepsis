#!/usr/bin/env node
/**
 * build_tte_pheno.mjs — 从 structured-measurement.csv.gz 流式抽取 TTE 表型宽表
 * 只保留 test_type=tte 且 measurement 在白名单内的行，按 measurement_id 透视成一行。
 * 低内存：仅累积白名单变量的值。
 * 用法: node build_tte_pheno.mjs <structured-measurement.csv.gz> [out.csv]
 */
import fs from 'node:fs';
import zlib from 'node:zlib';
import readline from 'node:readline';

const input = process.argv[2];
const outPath = process.argv[3] || 'tte_pheno.csv';
if (!input) { console.error('usage: node build_tte_pheno.mjs <structured-measurement.csv.gz> [out.csv]'); process.exit(1); }

// 白名单：变量名 -> 输出列名（可直接用原名）
const VARS = new Set([
  'lvef','lvef_upper','biplane_lvef','lvef_3d','rest_lvef','lv_function1','lv_function2',
  'lvedd','lvesd','septal_thickness','inf_lat_thickness','lv_wall_thickness','lv_chamber_size','fs',
  'la_vol','la_dimen','la_size',
  'lvot_diam','lvot_vti','av_vti',
  'tapse','rv_diam','rv_function','rv_chamber_size','rv_wall_thickness',
  'tr_velocity','tr_mmhg','phtn_severity','tv_est_pa_press_min','pasp_underestimated_tr',
  'mv_peak_e','mv_peak_a','mv_peak_e_a','mv_e_decel',
  'sept_e_prime','lat_e_prime','e_e_prime','mean_E_to_eprime','septal_E_to_eprime','lateral_E_to_eprime',
  'diastolic_fcn','diastolic_grade',
  'ivc','ivc_diam','pericardial_effusion','tamponade','technical_quality','test_modifier',
  // measurement-validity variables (LVOT-VTI as an index of forward stroke volume)
  'aortic_stenosis','aortic_regurg','av_mean_grad','av_pk_grad','av_pk_vel','av_area_continuity',
  'av_mass','av_leaflets','rest_lvot','induc_lvot','lvot_pk_vel','lvot_mod_mid_apical','mv_valve_area',
  'body_surface_area','weight_kg','height_cm','gender'
]);

function parseCsvLine(line) {
  const out = []; let cur = '', q = false;
  for (let i = 0; i < line.length; i++) {
    const c = line[i];
    if (q) { if (c === '"') { if (line[i + 1] === '"') { cur += '"'; i++; } else q = false; } else cur += c; }
    else if (c === '"') q = true;
    else if (c === ',') { out.push(cur); cur = ''; }
    else cur += c;
  }
  out.push(cur); return out;
}

const rs = fs.createReadStream(input);
const stream = input.toLowerCase().endsWith('.gz') ? rs.pipe(zlib.createGunzip()) : rs;
const rl = readline.createInterface({ input: stream, crlfDelay: Infinity });

const map = new Map(); // measurement_id -> {subject_id, dt, vals:{}}
let header = null, n = 0, kept = 0;

for await (const line of rl) {
  if (!line) continue;
  const f = parseCsvLine(line);
  if (header === null) { header = f.map(x => x.trim().toLowerCase()); continue; }
  n++;
  const g = (name) => { const i = header.indexOf(name); return i >= 0 ? f[i] : ''; };
  if (String(g('test_type')).toLowerCase() !== 'tte') continue;
  const m = g('measurement');
  if (!VARS.has(m)) continue;
  const mid = g('measurement_id');
  const val = (g('result') ?? '').trim();
  let row = map.get(mid);
  if (!row) { row = { subject_id: g('subject_id'), dt: g('measurement_datetime'), vals: {} }; map.set(mid, row); }
  if (val !== '') row.vals[m] = val;      // 非空优先
  kept++;
}

const cols = [...VARS];
const out = ['subject_id,measurement_id,measurement_datetime,' + cols.join(',')];
const esc = (s) => { s = String(s ?? ''); return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; };
for (const [mid, r] of map) {
  out.push([r.subject_id, mid, r.dt, ...cols.map(c => esc(r.vals[c] ?? ''))].join(','));
}
fs.writeFileSync(outPath, '\ufeff' + out.join('\n') + '\n');
console.log(`扫描行数: ${n}`);
console.log(`TTE 表型检查数: ${map.size}`);
console.log(`白名单命中行: ${kept}`);
console.log(`输出: ${outPath}`);
