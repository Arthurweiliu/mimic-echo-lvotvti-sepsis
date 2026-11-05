## config.R — paths for this reproduction package
##
## Two locations must be supplied:
##   MIMIC_ECHO_DIR : directory holding the derived tables exported from your own
##                    MIMIC-IV / MIMIC-IV-ECHO database (see README, "Inputs").
##   H1_WORK        : this repository (defaults to the working directory).
##
## All outputs are written to <H1_WORK>/results.

WORK <- Sys.getenv('H1_WORK', getwd())
DATA <- Sys.getenv('MIMIC_ECHO_DIR', file.path(WORK, '..', 'sepsis_echo_lit'))
OUT  <- file.path(WORK, 'results')

if (!dir.exists(DATA)) {
  stop('MIMIC_ECHO_DIR does not exist: ', DATA,
       '\nSet it with Sys.setenv(MIMIC_ECHO_DIR = "/path/to/derived/tables") before running.')
}
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setwd(OUT)

message('H1_WORK = ', WORK)
message('DATA    = ', DATA)
message('OUT     = ', OUT)
