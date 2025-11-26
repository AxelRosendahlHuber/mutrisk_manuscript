# Run manuscript analyses:
library(fs)
library(tictoc) # for speed testing

# 1. Run Mutrisk mutation rate estimates for all tissues (expected time ~1h)
#source("~/Nextcloud/Documents/mutrisk_manuscript/code/1_tissues/prepare_all_tissues.R")

# 2. Run resources scripts
resource_scripts = list.files("code/2_resources//",  full.names = TRUE)
for (script in resource_scripts) {
  print(script)
  tic()
  system(paste0("Rscript ", script), ignore.stdout = TRUE, ignore.stderr = TRUE)
  toc()
  gc()
}

# 3. Run Figure scripts
figure_scripts = list.files("code/3_figures/",  full.names = TRUE, recursive = TRUE)
for (script in figure_scripts) {
  print(script)
  tic()
  system(paste0("Rscript ", script), ignore.stdout = TRUE, ignore.stderr = TRUE)
  toc()
  gc()
}

