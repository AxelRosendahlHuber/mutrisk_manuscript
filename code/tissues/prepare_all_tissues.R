# script to systematically run all scripts
# script to run all scripts for the preparation
preparation_scripts = list.files("code/tissues/", recursive = TRUE, pattern = "prepare_", full.names = TRUE)
scripts = preparation_scripts[!grepl("bladder|pancreas|prepare_all", preparation_scripts)]

library(tictoc) # for speed testing
for (script in scripts) {
  print(script)
  tic()
  system(paste0("Rscript ", script), ignore.stdout = TRUE, ignore.stderr = TRUE)
  toc()
  gc()
}