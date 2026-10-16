#!/usr/bin/env node
/** build_congestion.mjs — 合并 AKI 结局 + 构建"右心/静脉淤血负担"评分 → congestion_analysis.csv
 *  用法: node build_congestion.mjs [数据目录]   (默认: 脚本所在目录) */
import fs from 'node:fs';
import readline from 'node:readline';
import { fileURLToPath } from 'node:url';
const DIR = (process.argv[2] || fileURLToPath(new URL('.', import.meta.url))).replace(/[\\/]?$/, '/');
function parse(l){const o=[];let c='',q=false;for(let i=0;i<l.length;i++){const ch=l[i];if(q){if(ch==='"'){if(l[i+1]==='"'){c+='"';i++;}else q=false;}else c+=ch;}else if(ch==='"')q=true;else if(ch===','){o.push(c);c='';}else c+=ch;}o.push(c);return o;}
async function readCsv(p){const rl=readline.createInterface({input:fs.createReadStream(p),crlfDelay:Infinity});let head=null;const rows=[];for await(const l of rl){if(!l)continue;const f=parse(l);if(!head){head=f.map(x=>x.replace(/^\uFEFF/,'').trim());continue;}const o={};head.forEach((h,i)=>o[h]=f[i]);rows.push({o,f});}return{head,rows};}
const num=v=>{if(v==null)return null;const s=String(v).replace(/[^0-9.\-]/g,'');if(s===''||s==='-')return null;const x=parseFloat(s);return isNaN(x)?null:x;};

const A=await readCsv(DIR+'analysis_data.csv');
const K=await readCsv(DIR+'aki_outcomes.csv');
const aki=new Map(K.rows.map(r=>[r.o.stay_id,r.o]));
// analysis_data 有重复 gender 列：用位置索引取第一个
const idxOf=(name)=>{let k=A.head.indexOf(name);return k;};
const I={};
for(const c of ['stay_id','admission_age','sofa_fd','sapsii','charlson_comorbidity_index','lactate_max_fd','creatinine_max','death28','death90','days_to_event','hospital_expire_flag','tr_velocity','tr_mmhg','la_size','la_dimen','la_vol','body_surface_area','rv_function','rv_chamber_size','rv_diam','ivc','ivc_diam','rrt_7d','invasive_vent_24h','norepi_eq_max_24h','lvef','gender','icu_intime','dischtime','dod']) I[c]=idxOf(c);
const T=s=>{if(!s)return null;const d=new Date(String(s).replace(' ','T'));return isNaN(d)?null:d.getTime();};

const RVP=['Mild global RV hypo','Moderate global RV hypo','Severe global hypo','RV function depressed','Apical free wall hypo',"Basal RV hypo (McConnell's sign)"];
const LApos=s=>/LAE|Elongated|Dilated|LAVI/i.test(s);
const RVsizePos=s=>/RVE/i.test(s);
const IVCpos=s=>/>2\.1|dilated/i.test(s);

const out=['stay_id,age,sex,sofa,sapsii,charlson,lactate,creat,death28,death90,days,fu_days,event90,hospital_death,aki_any,aki3,crrt,highTR,LAabn,RVabn,highIVC,score3,score4,tr_velocity,la_dimen,la_vol,rvi_bsa,rv_diam,lvef,norepi,vent24,rrt7d'];
let n=0, s3dist={0:0,1:0,2:0,3:0};
for(const {f} of A.rows){
  const g=c=>f[I[c]];
  const sid=g('stay_id'); const k=aki.get(sid)||{};
  const tr=num(g('tr_velocity')), trm=num(g('tr_mmhg')), lad=num(g('la_dimen')), lav=num(g('la_vol')), bsa=num(g('body_surface_area')), rvd=num(g('rv_diam'));
  const rvi = (lav!=null&&bsa)?lav/bsa:null;
  const highTR = (tr!=null&&tr>=2.8)||(trm!=null&&trm>=36) ? 1:0;
  const LAabn = (LApos(g('la_size'))||(lad!=null&&lad>=4.1)||(rvi!=null&&rvi>=34))?1:0;
  const RVabn = (RVP.includes(g('rv_function'))||RVsizePos(g('rv_chamber_size'))||(rvd!=null&&rvd>=4.2))?1:0;
  const highIVC = IVCpos(g('ivc')||'')?1:0;
  const score3 = highTR+LAabn+RVabn;
  const score4 = score3+highIVC;
  s3dist[score3]++;
  const intime=T(g('icu_intime')), disch=T(g('dischtime')), dodx=T(g('dod'));
  const hosp = (g('hospital_expire_flag')==='1'||g('hospital_expire_flag')==='t');
  let deathDate = dodx; if(!deathDate && hosp) deathDate = disch;
  // Survival time: time to death, else ADMINISTRATIVE censoring at day 90 (protocol §8).
  let fu=null, ev90=0;
  if(intime){
    if(deathDate){ const d0=(deathDate-intime)/86400000; fu = d0<0?0:d0; if(fu<=90) ev90=1; }
    else { fu = 90; }
  }
  const fuCap = fu==null?'':Math.min(fu,90).toFixed(2);
  const esc=v=>{v=String(v??'');return /[",\n]/.test(v)?'"'+v.replace(/"/g,'""')+'"':v;};
  out.push([sid, g('admission_age'), (f[I.gender]||''), g('sofa_fd'), g('sapsii'), g('charlson_comorbidity_index'), g('lactate_max_fd'), g('creatinine_max'),
    g('death28'), g('death90'), g('days_to_event'), fuCap, ev90, g('hospital_expire_flag'),
    (k.aki_any_7d==='t'?1:0), (k.aki3_7d==='t'?1:0), (k.crrt_7d==='t'?1:0),
    highTR, LAabn, RVabn, highIVC, score3, score4,
    tr??'', lad??'', lav??'', rvi!=null?rvi.toFixed(1):'', rvd??'', g('lvef'), g('norepi_eq_max_24h'), g('invasive_vent_24h'), g('rrt_7d')].map(esc).join(','));
  n++;
}
fs.writeFileSync(DIR+'congestion_analysis.csv','\uFEFF'+out.join('\n')+'\n');
console.log('N =',n);
console.log('score3 分布:',JSON.stringify(s3dist));
console.log('输出 congestion_analysis.csv');
