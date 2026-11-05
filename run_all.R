## run_all.R — run every stage in order
## Usage:  H1_WORK=/path/to/repo MIMIC_ECHO_DIR=/path/to/derived Rscript run_all.R
source('config.R')

stages <- c(
  'scripts/01_cohort/51_lactate_rebuild_v2.R',
  'scripts/01_cohort/53_flow_mi.R',
  'scripts/01_cohort/66_counts.R',
  'scripts/02_primary/52_mi_primary_v2.R',
  'scripts/02_primary/55_ph_mi.R',
  'scripts/02_primary/58_ladder_mi.R',
  'scripts/02_primary/57_h1h2_mi.R',
  'scripts/03_headtohead/67_mi_headtohead.R',
  'scripts/03_headtohead/71_extra_mi.R',
  'scripts/03_headtohead/74_lrt_mi.R',
  'scripts/03_headtohead/75_incremental_lrt.R',
  'scripts/04_sensitivity/78_s4_mi.R',
  'scripts/04_sensitivity/79_s3_mi.R',
  'scripts/04_sensitivity/80_eco_mi.R',
  'scripts/05_tables_figures/50_table1.R',
  'scripts/05_tables_figures/54_cascade_table.R',
  'scripts/05_tables_figures/70_table34_mi.R',
  'scripts/05_tables_figures/69_figures_mi.R',
  'scripts/05_tables_figures/77_fig1_v7.R',
  'scripts/05_tables_figures/61_figures_v6.R'
)

## The two bootstrap jobs are left out of the default sequence because they take
## ~15 min each; run them separately (see README section 5) and then:
##   68_mi_delta.R with SET=H1 and SET=H2
##   73_optimism.R  with SET=H1 and SET=H2

root <- WORK
for (s in stages) {
  message('\n=== ', s, ' ===')
  t0 <- Sys.time()
  setwd(root)
  source(file.path(root, s), chdir = FALSE)
  message('    done in ', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), ' min')
}
message('\nAll stages finished; outputs are in ', OUT)
