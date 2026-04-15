library(httr2)
library(R.utils)


dirs = c("raw_data", "raw_data/blood", "raw_data/lung", "raw_data/colon","raw_data/lung",
         "processed_data", "processed_data/blood", "processed_data/lung", "processed_data/colon",
         "processed_data", "processed_data/blood", "processed_data/lung", "processed_data/colon",
         "plots", "plots/coverage_saturation", "plots/blood", "plots/lung", "plots/colon",
         "manuscript/", "manuscript/Figure_1/", "manuscript/Figure_2/", "manuscript/Figure_3/", "manuscript/Figure_4/",
         "manuscript/Figure_5/", "manuscript/figure_panels/",  "manuscript/Schematic_poster_presentations/",
         "manuscript/Supplementary_Figures/", "manuscript/Supplementary_Figures/Figure_S1/",
         "manuscript/Supplementary_Figures/Figure_S2/",
         "manuscript/Supplementary_Figures/Figure_S3/","manuscript/Supplementary_Figures/Figure_S4/",
         "manuscript/Supplementary_Figures/Figure_S5/", "manuscript/Supplementary_Figures/Figure_S6/",
         "manuscript/Supplementary_notes/", "manuscript/Supplementary_notes/Supplementary_Note_I/",
         "manuscript/Supplementary_notes/Supplementary_Note_III//",
         "manuscript/Supplementary_Tables/")

for (dir in dirs) {
  if(!dir.exists(dir)) {dir.create(dir)}
}


# tools to download the main mutation data from Mendeley data
urls <- c("https://data.mendeley.com/public-api/zip/np54zjkvxr/download/1",
          "https://data.mendeley.com/public-api/zip/x3vsxpspn4/download/2",
          "https://github.com/TimCoorens/Polymerase/archive/refs/heads/master.zip",
          "https://data.mendeley.com/public-api/zip/b53h2kwpyy/download/1")

dest <- c("raw_data/blood/blood.zip", "raw_data/colon/normal/colon_normal.zip", "raw_data/colon/hypermutated/colon_hm.zip", "raw_data/lung/lung.zip")

for (i in seq_along(urls)) {
  outdir = dirname(dest[i])
  if (!dir.exists(outdir)) { dir.create(outdir) }

  request(urls[i]) |>  req_perform(path = dest[i], verbosity = 0)
  unzip(zipfile = dest[i], exdir = outdir)
  # remove original zipfile
  file.remove(dest[[i]])

  # if other zipfiles in the outdir, unzip these too
  zips = list.files(outdir, pattern = "\\.zip$", full.names = TRUE, recursive = TRUE)
  for (z in zips) {
   unzip(z, exdir = outdir)
   file.remove(z)   # optional: remove the inner zip after extraction
  }

  # Detect and remove first-level directory
  top <- list.dirs(dirname(dest[i]), full.names = TRUE, recursive = FALSE)
  files <- list.files(top, full.names = TRUE)
  file.rename(files, file.path(dirname(dest[i]), basename(files)))
  unlink(top, recursive = TRUE)
}



# Add manual metadata files:
# download colon metadata:
# Loading metadata
zipfile = "raw_data/colon/normal/colon_normal_metadata.zip"
outdir = "raw_data/colon/normal/mmc.zip"
url = "https://www.cell.com/cms/10.1016/j.cell.2020.06.036/attachment/8e54fb1e-60b8-4c8a-ac3d-1ea905138889/mmc2.zip"
httr2::request(url) |>
  httr2::req_perform(path = outdir, verbosity = 0)
utils::unzip(zipfile = zipfile, exdir = "raw_data/colon/normal/")
file.remove(zipfile)


# TODO Download the lung organoid telomere data files:

