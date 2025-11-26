# script to systematically run all scripts
# script to run all scripts for the preparation
scripts = list.files("code/1_tissues/", recursive = TRUE, pattern = "prepare_", full.names = TRUE)

library(tictoc) # for speed testing
for (script in scripts) {
  print(script)
  tic()
  system(paste0("Rscript ", script), ignore.stdout = TRUE, ignore.stderr = TRUE)
  toc()
  gc()
}