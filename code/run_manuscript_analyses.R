# Run manuscript analyses:
library(tictoc) # for speed testing

# # 0 (download data sources) automatically:
# if(!all(dir.exists("raw_data/blood/mutational_signatures_analysis/burden_all.txt/"),
#         dir.exists("raw_data/colon//"))) {
# source("code/01_download_data/Download_data.R")
# }

# additionally, download the data from BoostDM, BoostDM-CH, GENIE and COSMIC
# 1. Run Mutrisk mutation rate estimates for all tissues (expected time ~1h)
scripts = list.files("code/1_tissues/", recursive = TRUE, pattern = "pre", full.names = TRUE)
for (script in scripts) {
  print(script)
  tic()
  system(paste0("Rscript ", script))
  toc()
  gc()
}

# 2. Run resources scripts
resource_scripts = list.files("code/2_resources//",  full.names = TRUE)
for (script in resource_scripts) {
  print(script)
  tic()
  system(paste0("Rscript ", script))
  toc()
  gc()
}

# 3. Run Figure scripts (including Supplementary, and Supplementary Note figures)
figure_scripts = list.files("code/3_figures/",  full.names = TRUE, recursive = TRUE)
for (script in figure_scripts) {
  print(script)
  tic()
  system(paste0("Rscript ", script))
  toc()
  gc()
}

