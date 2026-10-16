#!/usr/bin/env node
/** summarize_congestion.mjs — 检查 analysis_data.csv 中淤血相关超声变量的覆盖率 */
import fs from 'node:fs';
import readline from 'node:readline';
import { fileURLToPath } from 'node:url';
const DIR = (process.argv[2] || fileURLToPath(new URL('.', import.meta.url))).replace(/[\\/]?$/, '/');
const F = DIR + 'analysis_data.csv';
const TARGETS = ['ivc_diam','ivc','tr_velocity','tr_mmhg','la_vol','la_dimen','la_size','rv_diam','rv_chamber_size','rv_function','body_surface_area','e_e_prime','mean_E_to_eprime','lvef','tapse','lvot_vti','lvot_diam'];
function parse(line){const o=[];let c='',q=false;for(let i=0;i<line.length;i++){const ch=line[i];if(q){if(ch==='"'){if(line[i+1]==='"'){c+='"';i++;}else q=false;}else c+=ch;}else if(ch==='"')q=true;else if(ch===','){o.push(c);c='';}else c+=ch;}o.push(c);return o;}
const rl = readline.createInterface({input: fs.createReadStream(F), crlfDelay: Infinity});
let head=null; const cnt={}, sum={};
for await (const line of rl){ if(!line) continue; const f=parse(line); if(!head){ head=f.map(x=>x.replace(/^\uFEFF/,'')); continue; } for(const t of TARGETS){ const i=head.indexOf(t); if(i<0) continue; const v=(f[i]??'').trim(); if(v!==''){ cnt[t]=(cnt[t]||0)+1; const n=parseFloat(v); if(!isNaN(n)) sum[t]=(sum[t]||0)+n; } } nRows++; }
var nRows=0;
console.log('总行数(不含表头):', nRows);
console.log('列是否存在:');
for(const t of TARGETS) console.log('  '+t+': '+(head.includes(t)?'YES':'NO'));
console.log('非空计数 / 非空均值:');
for(const t of TARGETS){ if(cnt[t]||sum[t]) console.log('  '+t+': n='+(cnt[t]||0)+' ('+((cnt[t]||0)/nRows*100).toFixed(1)+'%)  mean='+(sum[t]&&cnt[t]?(sum[t]/cnt[t]).toFixed(3):'NA')); }
