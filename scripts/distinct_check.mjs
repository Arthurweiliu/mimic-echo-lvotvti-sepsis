#!/usr/bin/env node
import fs from 'node:fs';
import readline from 'node:readline';
import { fileURLToPath } from 'node:url';
const DIR = (process.argv[2] || fileURLToPath(new URL('.', import.meta.url))).replace(/[\\/]?$/, '/');
const F = DIR + 'analysis_data.csv';
const COLS=['rv_function','rv_chamber_size','la_size','ivc','diastolic_fcn'];
function parse(l){const o=[];let c='',q=false;for(let i=0;i<l.length;i++){const ch=l[i];if(q){if(ch==='"'){if(l[i+1]==='"'){c+='"';i++;}else q=false;}else c+=ch;}else if(ch==='"')q=true;else if(ch===','){o.push(c);c='';}else c+=ch;}o.push(c);return o;}
const rl=readline.createInterface({input:fs.createReadStream(F),crlfDelay:Infinity});
let head=null; const maps={}; const seen=new Set();
for await(const line of rl){ if(!line) continue; const f=parse(line);
  if(!head){ head=f.map(x=>x.replace(/^\uFEFF/,'')); for(const c of COLS){ let n=c,k=0; while(seen.has(n)){n=c+'_'+(++k);} seen.add(n); maps[c]=head.indexOf(c); } continue; }
  for(const c of COLS){ const i=maps[c]; if(i<0) continue; const v=(f[i]??'').trim(); if(v==='') continue; maps['__'+c]=maps['__'+c]||{}; if(maps['__'+c]) {} (globalThis['m_'+c]=globalThis['m_'+c]||{})[v]=((globalThis['m_'+c]||{})[v]||0)+1; }
}
for(const c of COLS){ const m=globalThis['m_'+c]||{}; const arr=Object.entries(m).sort((a,b)=>b[1]-a[1]).slice(0,14);
  console.log('=== '+c+' ==='); for(const [k,v] of arr) console.log('  '+v+'  '+k); }
