// extract_echo_extra.mjs — stream the raw MIMIC-IV-ECHO structured measurements and pull
// only the measurements not in the extracted phenotype tables:
//   resting_hr (echo-time heart rate), mitral_regurg, tricuspid_regurg, pulm_regurg
// Output: echo_extra.csv  (measurement_id,measurement,result)
import fs from 'fs';
import zlib from 'zlib';
import readline from 'readline';
import { fileURLToPath } from 'node:url';

// Data directory: 1st CLI argument, else this script's own directory.
const DIR = (process.argv[2] || fileURLToPath(new URL('.', import.meta.url))).replace(/[\\/]?$/, '/');
const SRC = DIR + 'structured-measurement.csv.gz';
const OUT = DIR + 'echo_extra.csv';
const TARGET = new Set(['resting_hr', 'mitral_regurg', 'tricuspid_regurg', 'pulm_regurg', 'mitral_stenosis']);

function splitCSV(line) {
  const out = []; let f = '', q = false;
  for (let i = 0; i < line.length; i++) {
    const c = line[i];
    if (q) { if (c === '"') { if (line[i + 1] === '"') { f += '"'; i++; } else q = false; } else f += c; }
    else { if (c === '"') q = true; else if (c === ',') { out.push(f); f = ''; } else f += c; }
  }
  out.push(f); return out;
}

const rs = fs.createReadStream(SRC).pipe(zlib.createGunzip());
const rl = readline.createInterface({ input: rs, crlfDelay: Infinity });
const ws = fs.createWriteStream(OUT);
ws.write('measurement_id,measurement,result\n');

let n = 0, kept = 0, header = true;
rl.on('line', (line) => {
  n++;
  if (header) { header = false; return; }
  // quick filter before full parse
  const cheap = line.indexOf('resting_hr') >= 0 || line.indexOf('_regurg') >= 0 || line.indexOf('mitral_stenosis') >= 0;
  if (!cheap) return;
  const c = splitCSV(line);
  if (c.length < 7) return;
  const mid = c[1], meas = c[4], res = c[6];
  if (TARGET.has(meas)) { ws.write(`${mid},${meas},${res}\n`); kept++; }
});
rl.on('close', () => { ws.end(); console.log('lines read:', n, ' rows kept:', kept); });
