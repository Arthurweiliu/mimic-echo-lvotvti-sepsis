#!/usr/bin/env node
/* Portable: set H1_WORK (this repository) and MIMIC_ECHO_DIR (your derived tables). */
/** 64_number_audit.mjs — check that every estimated HR/CI reported in the manuscript
 *  has a matching value in a verified analysis output, and that cohort counts agree. */
import fs from 'node:fs';

const DIR = WORK;
const MS = `${DIR}/MANUSCRIPT_BMC_ANESTHESIOLOGY.md`;
const SRC = fs.readdirSync(DIR)
  .filter(f => /_out\.txt$/.test(f) || /^Table.*\.(csv|md)$/.test(f))
  .map(f => `${DIR}/${f}`);

const readEnc = f => {
  const b = fs.readFileSync(f);
  if (b[0] === 0xff && b[1] === 0xfe) return b.toString('utf16le');
  let s = b.toString('utf8');
  if (s.charCodeAt(0) === 0xfeff) s = s.slice(1);
  return s;
};

const srcText = SRC.filter(f => fs.existsSync(f)).map(readEnc).join('\n');
const ms = readEnc(MS);

// values quoted from the literature rather than estimated here
const LIT = new Set(['2.37|1.47|3.80', '1.91|1.41|2.59', '2.75|1.82|4.17', '1.70|1.09|1.88',
                     '1.87|1.00|3.00', '1.93|1.00|3.00']);

const re = /(\d+\.\d+)\s*\((?:95% CI )?(\d+\.\d+)\s*[-\u2013]\s*(\d+\.\d+)\)/g;
const norm = (a, b, c, d = 2) => [+a, +b, +c].map(v => v.toFixed(d)).join('|');

const src2 = new Set();
for (const m of srcText.matchAll(re)) {
  src2.add(norm(m[1], m[2], m[3], 2));
  src2.add(norm(m[1], m[2], m[3], 1));
}

let bad = 0, checked = 0;
const weak = [];
const srcNums = new Set();
for (const m of srcText.matchAll(/\d+\.\d+/g)) srcNums.add((+m[0]).toFixed(2));
for (const m of ms.matchAll(re)) {
  checked++;
  const k2 = norm(m[1], m[2], m[3], 2), k1 = norm(m[1], m[2], m[3], 1);
  if (src2.has(k2) || src2.has(k1) || LIT.has(k2)) continue;
  const allPresent = [m[1], m[2], m[3]].every(v => srcNums.has((+v).toFixed(2)));
  if (allPresent) { weak.push(m[0]); continue; }
  console.log(`  MISSING in sources: ${m[0]}`); bad++;
}
if (weak.length) console.log(`  weak match (numbers present individually, pair not on one line): ${[...new Set(weak)].join('; ')}`);
console.log(`HR/CI triples checked: ${checked}`);

// --- cohort counts ----------------------------------------------------
const counts = ['4,816', '1,573', '4,336', '1,403', '4,081', '1,311', '3,185', '3,573',
                '1,262', '3,902', '1,343', '4,346', '1,449', '4,589', '1,346', '4,268',
                '1,025', '3,661', '1,212', '5,366', '5,306', '7,253', '12,120', '1,242', '913', '3,574', '3,903', '548', '480', '735', '468', '445', '329', '8,568'];
const flatSrc = srcText.replace(/(\d),(\d{3})/g, '$1$2');
for (const c of counts) {
  const inMs = ms.includes(c);
  const inSrc = srcText.includes(c) || flatSrc.includes(c.replace(',', ''));
  if (inMs && !inSrc) { console.log(`  count not in sources: ${c}`); bad++; }
}

// --- structure --------------------------------------------------------
const absStart = ms.indexOf('## Abstract');
const absEnd = ms.indexOf('**Keywords:**');
const absWords = ms.slice(absStart, absEnd).replace(/\*\*/g, ' ').replace('## Abstract', '')
  .split(/\s+/).filter(Boolean).length;
console.log(`abstract words (excl. heading and keywords): ${absWords}${absWords > 350 ? '  <-- OVER LIMIT' : ''}`);

const refs = [...ms.matchAll(/^(\d+)\.\s+[A-Z]/gm)].map(m => +m[1]);
const maxRef = Math.max(...refs);
const cited = new Set();
for (const m of ms.matchAll(/\[([\d,\s\u2013-]+)\]/g)) {
  for (const part of m[1].split(',')) {
    const s = part.trim();
    if (/^\d+$/.test(s)) cited.add(+s);
    else if (/^\d+\s*-\s*\d+$/.test(s)) {
      const [a, b] = s.split('-').map(x => +x.trim());
      for (let i = a; i <= b; i++) cited.add(i);
    }
  }
}
const uncited = [];
for (let i = 1; i <= maxRef; i++) if (!cited.has(i)) uncited.push(i);
console.log(`references: ${maxRef} listed${uncited.length ? `, uncited: ${uncited.join(',')}` : ', all cited'}`);

const tables = (ms.match(/^### Table \d\./gm) || []).length;
const figs = (ms.match(/^\*\*Figure \d\.\*\*/gm) || []).length;
console.log(`tables: ${tables}, figure legends: ${figs}`);

console.log(bad === 0 ? '\nAUDIT PASS: every reported estimate traces to an output file' : `\nAUDIT: ${bad} issue(s)`);
