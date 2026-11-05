#!/usr/bin/env node
/**
 * build_variable_dictionary.mjs  (v2 — 低内存两遍法)
 * 从 MIMIC-IV-ECHO 的 structured_measurement.csv(.gz) 生成"完整变量字典"。
 *
 * v1 用 V8 Set 逐行累计去重 -> 28M 行时 OOM。v2 改为：
 *   Pass 1：收集全局 id 域（measurement_id / subject_id -> 稠密下标）、每个变量的行数/空值数、变量元数据
 *   Pass 2：为每个变量分配 bitset（按 id 域大小），标记成员，最后 popcount 得到精确去重计数
 * 零依赖（Node 内置模块）。输出列与原版一致，数值为精确值（非估计）。
 *
 * 用法:
 *   node build_variable_dictionary.mjs <structured_measurement.csv|.gz> [out.csv]
 */
import fs from 'node:fs';
import zlib from 'node:zlib';
import readline from 'node:readline';

const input = process.argv[2];
const outPath = process.argv[3] || 'variable_dictionary.csv';
if (!input) {
  console.error('usage: node build_variable_dictionary.mjs <structured_measurement.csv[.gz]> [out.csv]');
  process.exit(1);
}

function openStream(p) {
  const rs = fs.createReadStream(p, { highWaterMark: 1 << 20 });
  return p.toLowerCase().endsWith('.gz') ? rs.pipe(zlib.createGunzip()) : rs;
}

// 最简 RFC4180 行解析（支持引号内逗号与转义引号）
function parseCsvLine(line) {
  const out = [];
  let cur = '', q = false;
  for (let i = 0; i < line.length; i++) {
    const c = line[i];
    if (q) {
      if (c === '"') { if (line[i + 1] === '"') { cur += '"'; i++; } else q = false; }
      else cur += c;
    } else if (c === '"') q = true;
    else if (c === ',') { out.push(cur); cur = ''; }
    else cur += c;
  }
  out.push(cur);
  return out;
}

// 类别规则（顺序敏感：先具体后一般）
const RULES = [
  ['IVC/静脉回流', /(\bivc\b|inferior_vena|caval|下腔静脉)/],
  ['右心/右房 (RV/RA)', /(^rv|_rv\b|right_vent|right_atrial|\bra_|tapse|tricuspid_annulus|rvsp|rv_fac|rv_strain|fractional_area|右心|右房)/],
  ['三尖瓣/肺动脉压 (TR/PASP)', /(tr_velocity|tr_max|\btv_|tricuspid|\bpasp\b|pa_press|pulmonary_press|pulmonary_artery|mpap|三尖瓣|肺动脉压)/],
  ['主动脉/流出道 (Ao/LVOT)', /(aortic|\bav_|aorta|sinus_of|sinotubular|annulus|ascending|\baov\b|lvot|outflow|主动脉|流出道)/],
  ['左室尺寸/质量 (LV dims/mass)', /(lvedd|lvesd|lvid|lv_mass|lvmi|septal_thickness|inf_lat_thickness|ivsd|ivss|pwt|thickness|dimension|直径|内径|厚度)/],
  ['左房 (LA)', /(la_vol|left_atrial|la_diam|la_area|la_ap|\blai\b|la_strain|lars|左房)/],
  ['左室收缩功能 (LV systolic)', /(lvef|ejection_fraction|\bef\b|lvot_vti|stroke_volume|cardiac_output|cardiac_index|\bgls\b|strain|mapse|hyperdynamic|收缩功能|射血分数)/],
  ['二尖瓣/组织多普勒/充盈 (MV/TDI)', /(mv_peak|mv_e_|mv_a_|mitral|sept_e_prime|lat_e_prime|e_prime|eprime|e_e_prime|e_to_eprime|\bmva\b|pht|decel|\be_a\b|diastolic|ivrt|二尖瓣|舒张)/],
  ['基础/通用 (Basic)', /(body_surface|\bbsa\b|heart_rate|\bhr\b|weight|height|rhythm|test_type|性别|年龄|基础)/]
];

function categorize(m, d) {
  const s = `${m} ${d}`.toLowerCase();
  for (const [name, re] of RULES) if (re.test(s)) return name;
  return '其他/未归类 (Other)';
}

function eachRow(cb) {
  return new Promise((resolve, reject) => {
    const rl = readline.createInterface({ input: openStream(input), crlfDelay: Infinity });
    let header = null;
    rl.on('line', (line) => {
      if (!line) return;
      const f = parseCsvLine(line);
      if (header === null) {
        header = f.map((x) => x.trim().toLowerCase());
        if (!header.includes('measurement')) { console.error('表头缺 measurement 列，实际:', header.join(',')); process.exit(2); }
        return;
      }
      cb(f, header);
    });
    rl.on('close', resolve);
    rl.on('error', reject);
  });
}

(async () => {
  // ---------- Pass 1 ----------
  const agg = new Map();          // key -> {cat,m,desc,unit,rows,empty,nonempty}
  const studyIdx = new Map();     // measurement_id -> dense index
  const subjIdx = new Map();      // subject_id -> dense index
  let total = 0;

  await eachRow((f, header) => {
    total++;
    const g = (name) => { const i = header.indexOf(name); return i >= 0 ? f[i] : ''; };
    const measurement = g('measurement');
    const desc = g('measurement_description');
    const unit = g('unit');
    const sid = g('subject_id');
    const mid = g('measurement_id');
    const result = (g('result') ?? '').trim();

    const key = `${measurement}\t${desc}\t${unit}`;
    let a = agg.get(key);
    if (!a) { a = { cat: categorize(measurement, desc), m: measurement, desc, unit, rows: 0, empty: 0, nonempty: 0 }; agg.set(key, a); }
    a.rows++;
    if (result === '') a.empty++; else a.nonempty++;

    if (mid !== '' && !studyIdx.has(mid)) studyIdx.set(mid, studyIdx.size);
    if (sid !== '' && !subjIdx.has(sid)) subjIdx.set(sid, subjIdx.size);
  });

  const U_m = studyIdx.size, U_s = subjIdx.size;
  const b_m = (U_m + 7) >> 3, b_s = (U_s + 7) >> 3;

  // ---------- Pass 2 ----------
  const bits = new Map(); // key -> {s:Uint8Array, b:Uint8Array}
  for (const key of agg.keys()) bits.set(key, { s: new Uint8Array(b_m), b: new Uint8Array(b_s) });

  await eachRow((f, header) => {
    const g = (name) => { const i = header.indexOf(name); return i >= 0 ? f[i] : ''; };
    const key = `${g('measurement')}\t${g('measurement_description')}\t${g('unit')}`;
    const rec = bits.get(key); if (!rec) return;
    const mi = studyIdx.get(g('measurement_id')); if (mi !== undefined) rec.s[mi >> 3] |= (1 << (mi & 7));
    const si = subjIdx.get(g('subject_id')); if (si !== undefined) rec.b[si >> 3] |= (1 << (si & 7));
  });

  function popcount(arr) { let c = 0; for (let i = 0; i < arr.length; i++) { let v = arr[i]; while (v) { v &= v - 1; c++; } } return c; }

  const rows = [...agg.entries()].map(([key, a]) => {
    const rec = bits.get(key);
    return { ...a, studies: popcount(rec.s), subjects: popcount(rec.b) };
  }).sort((x, y) => x.cat.localeCompare(y.cat) || x.m.localeCompare(y.m));

  function csvEscape(s) { s = String(s ?? ''); return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; }
  const out = ['类别,measurement,measurement_description,unit,行数,检查数,患者数,检查覆盖率%,非空结果数,空结果数'];
  for (const a of rows) {
    const cov = U_m ? (a.studies / U_m * 100).toFixed(1) : '';
    out.push([a.cat, a.m, a.desc, a.unit, a.rows, a.studies, a.subjects, cov, a.nonempty, a.empty].map(csvEscape).join(','));
  }
  fs.writeFileSync(outPath, '\ufeff' + out.join('\n') + '\n');

  const byCat = new Map();
  for (const a of rows) byCat.set(a.cat, (byCat.get(a.cat) || 0) + 1);

  console.log(`处理行数: ${total}`);
  console.log(`唯一变量(measurement×desc×unit): ${rows.length}`);
  console.log(`全库检查数(去重 measurement_id): ${U_m}，患者数(去重 subject_id): ${U_s}`);
  for (const [c, n] of [...byCat].sort()) console.log(`  ${c}: ${n}`);
  console.log(`输出: ${outPath}`);
})().catch((e) => { console.error(e); process.exit(1); });
