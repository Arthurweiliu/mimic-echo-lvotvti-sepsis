#!/usr/bin/env node
/** build_h1.mjs — 构建 H1 分析数据集 h1_data.csv */
import fs from 'node:fs';
import readline from 'node:readline';
import { fileURLToPath } from 'node:url';
// Data directory: 1st CLI argument, else this script's own directory.
const DIR = (process.argv[2] || fileURLToPath(new URL('.', import.meta.url))).replace(/[\\/]?$/, '/');
function parse(l){const o=[];let c='',q=false;for(let i=0;i<l.length;i++){const ch=l[i];if(q){if(ch==='"'){if(l[i+1]==='"'){c+='"';i++;}else q=false;}else c+=ch;}else if(ch==='"')q=true;else if(ch===','){o.push(c);c='';}else c+=ch;}o.push(c);return o;}
async function readCsv(p){const rl=readline.createInterface({input:fs.createReadStream(p),crlfDelay:Infinity});let head=null;const rows=[];for await(const l of rl){if(!l)continue;const f=parse(l);if(!head){head=f.map(x=>x.replace(/^\uFEFF/,'').trim());continue;}const o={};head.forEach((h,i)=>o[h]=f[i]);rows.push({o,f});}return{head,rows};}
const T=s=>{if(!s)return null;const d=new Date(String(s).replace(' ','T'));return isNaN(d)?null:d.getTime();};
const A=await readCsv(DIR+'analysis_data.csv');
const K=await readCsv(DIR+'aki_outcomes.csv'); const aki=new Map(K.rows.map(r=>[r.o.stay_id,r.o]));
const Y=await readCsv(DIR+'anchor_year_group.csv'); const yg=new Map(Y.rows.map(r=>[r.o.subject_id,r.o.anchor_year_group]));
const cols=['stay_id','subject_id','admission_age','sofa_fd','sapsii','charlson_comorbidity_index','lactate_max_fd','creatinine_max','mbp_min','heart_rate_mean','temperature_mean','spo2_min','bun_max','hemoglobin_min','platelets_min','wbc_max','bilirubin_total_max','inr_max','troponin_t_max_24h','norepi_eq_max_24h','invasive_vent_24h','rrt_7d','death28','death90','days_to_event','hospital_expire_flag','icu_intime','dischtime','dod','lvot_vti','lvot_diam','lvef','tr_velocity','la_dimen','rv_function','los_icu'];
const idx={}; A.head.forEach((h,i)=>{ if(cols.includes(h)&&idx[h]===undefined) idx[h]=i; });
const out=['stay_id,subject_id,age,sex,yg,sofa,sapsii,charlson,lactate,creat,mbp_min,hr_mean,temp_mean,spo2_min,bun,hgb,plt,wbc,bili,inr,troponin,norepi,vent,rrt,death28,death90,fu_days,event90,hospital_death,vti,vti_diam,lvef,trv,lad,rv_function,los_icu,aki_any,aki3,crrt'];
let n=0;
for(const {o,f} of A.rows){
  const g=c=>idx[c]!==undefined?f[idx[c]]:'';
  const sid=g('stay_id');
  const intime=T(g('icu_intime')), disch=T(g('dischtime')), dodx=T(g('dod'));
  const hosp='' + (g('hospital_expire_flag')==='1'||g('hospital_expire_flag')==='t'?1:0);
  let deathDate=dodx; if(!deathDate && hosp==='1') deathDate=disch;
  // Survival time = time to death, or ADMINISTRATIVE censoring at day 90 when no
  // death is recorded within 90 days (per protocol section 8: "censored at 90 days").
  // Do NOT censor at hospital discharge: that discards post-discharge follow-up.
  let fu='', ev90=0;
  if(intime){
    if(deathDate){
      let d=(deathDate-intime)/86400000; if(d<0)d=0;
      if(d<=90){ ev90=1; fu=Math.min(d,90).toFixed(2); } else { fu='90.00'; }
    } else { fu='90.00'; }
  }
  const k=aki.get(sid)||{};
  const esc=v=>{v=String(v??'');return /[",\n]/.test(v)?'"'+v.replace(/"/g,'""')+'"':v;};
  out.push([sid,g('subject_id'),g('admission_age'),f[A.head.indexOf('gender')]||'',yg.get(g('subject_id'))||'',g('sofa_fd'),g('sapsii'),g('charlson_comorbidity_index'),g('lactate_max_fd'),g('creatinine_max'),g('mbp_min'),g('heart_rate_mean'),g('temperature_mean'),g('spo2_min'),g('bun_max'),g('hemoglobin_min'),g('platelets_min'),g('wbc_max'),g('bilirubin_total_max'),g('inr_max'),g('troponin_t_max_24h'),g('norepi_eq_max_24h'),g('invasive_vent_24h'),g('rrt_7d'),g('death28'),g('death90'),fu,ev90,hosp,g('lvot_vti'),g('lvot_diam'),g('lvef'),g('tr_velocity'),g('la_dimen'),g('rv_function'),g('los_icu'),(k.aki_any_7d==='t'?1:0),(k.aki3_7d==='t'?1:0),(k.crrt_7d==='t'?1:0)].map(esc).join(','));
  n++;
}
fs.writeFileSync(DIR+'h1_data_allstays.csv','\uFEFF'+out.join('\n')+'\n');
console.log('N =',n,'→ h1_data_allstays.csv (all eligible stays); run 15_make_first_icu.R next');
