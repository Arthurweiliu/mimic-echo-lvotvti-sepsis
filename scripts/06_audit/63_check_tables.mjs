/* Portable: set H1_WORK (this repository) and MIMIC_ECHO_DIR (your derived tables). */
import fs from 'node:fs';
const M = process.argv[2] || file.path(WORK, "MANUSCRIPT_BMC_ANESTHESIOLOGY.md");
const md = fs.readFileSync(M, 'utf8').replace(/\r\n/g, '\n').split('\n');
let bad = 0, tables = 0, cur = null, expect = 0;
md.forEach((l, idx) => {
  if (/^\|/.test(l)) {
    const n = (l.match(/\|/g) || []).length;
    if (cur === null) { cur = idx; expect = n; tables++; }
    else if (n !== expect) { console.log(`line ${idx + 1}: ${n} pipes, expected ${expect} :: ${l.slice(0, 70)}`); bad++; }
  } else { cur = null; }
});
console.log(`tables: ${tables}, ragged rows: ${bad}`);
