import fs from 'node:fs';
const lines = fs.readFileSync('variable_dictionary.csv', 'utf8').split('\n');
const H = lines[0];
const want = new Set(['lvef','lvef_upper','lvef_3d','biplane_lvef','rest_lvef','lv_function1','lv_function2','fs',
 'lvot_vti','lvot_diam','av_vti','tapse','rv_function','rv_chamber_size','rv_diam','rv_wall_thickness',
 'tr_velocity','tr_mmhg','tv_est_pa_press_min','pasp_underestimated_tr','e_e_prime','sept_e_prime','lat_e_prime',
 'mean_E_to_eprime','septal_E_to_eprime','lateral_E_to_eprime','mv_peak_e','mv_peak_a','mv_peak_e_a','mv_e_decel',
 'la_vol','la_dimen','la_size','ivc','ivc_diam','lv_dyssynchrony','lvedd','lvesd','septal_thickness',
 'inf_lat_thickness','lv_chamber_size','lv_wall_thickness','diastolic_fcn','diastolic_grade','phtn_severity',
 'pericardial_effusion','tamponade','technical_quality','test_modifier','body_surface_area','weight_kg','height_cm','gender']);
// dictionary rows have no embedded quotes (verified); plain comma split is safe here.
const out = [H]; const seen = new Set();
for (let i = 1; i < lines.length; i++) {
  if (!lines[i]) continue;
  const m = lines[i].split(',')[1];
  if (want.has(m) && !seen.has(m)) { seen.add(m); out.push(lines[i]); }
}
fs.writeFileSync('选题变量表.csv', '\ufeff' + out.join('\n') + '\n');
console.log('rows written:', out.length - 1);
console.log(out.join('\n'));
