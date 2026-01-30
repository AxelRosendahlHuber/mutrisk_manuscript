# Figure S3
# summary of the signature refitting methods across all cohorts

# load metadata
md_files = list.files("processed_data/", recursive = TRUE, pattern = "_metadata",
                      full.names = TRUE)
names(md_files) = gsub("_.*", "", basename(md_files))
metadata = lapply(md_files, fread) |>
  rbindlist(idcol = "tissue", use.names = TRUE, fill = TRUE) |>
  dplyr::select(any_of(c("tissue", "sampleID", "category", "age", "donor"))) |>
  dplyr::distinct()



sig_contri_files = list.files("processed_data/", recursive = TRUE,
                                     pattern = "signature_contributions",
                                     full.names = TRUE)

names(sig_contri_files) = paste0(
  str_split_i(sig_contri_files,pattern = "/", i = c(3)), "_",
  str_split_i(sig_contri_files,pattern = "/", i = c(4)))



means_sig_contrib = sig_contribution_plot = list()
for (name in names(sig_contri_files)) {

  file = sig_contri_files[name]
  sig_per_sample = fread(file) |>
    left_join(metadata) |>
    group_by(donor) |>
    summarize(across(starts_with("SBS"), mean))

  mean_sig_contrib = sig_per_sample |>
    column_to_rownames("donor") |> as.matrix() |>
    prop.table(1) |>
    as.data.frame() |>
    rownames_to_column("donor") |>
    pivot_longer(starts_with("SBS"), names_to = "signature", values_to = "contribution")

  means_sig_contrib[[name]] = mean_sig_contrib |>
    group_by(signature) |>
    summarize(sd = sd(contribution),
              contribution = mean(contribution))

  sig_contribution_plot[[name]] = means_sig_contrib[[name]] |>
    ggplot(aes(x = signature, y = contribution)) +
    geom_col() +
    geom_point(data = mean_sig_contrib) +
    #geom_errorbar(aes(ymin = contribution, ymax = contribution + sd), width = 0.2) +
    labs(title = name, y = "relative signature contribution") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

sig_contribution_plot = wrap_plots(sig_contribution_plot)
ggsave("manuscript/Supplementary_Figures/Figure_S2/Figure_S2_all.png",
       sig_contribution_plot, width = 15, height = 10)


order(parse_number(means_sig_contrib[[2]]$signature))

mean_signatures = means_sig_contrib |>
  rbindlist(idcol = "name")

sigs = unique(mean_signatures$signature)
levels = sigs[order(parse_number(sigs), decreasing = FALSE)]
mean_signatures_ordered = mean_signatures |>
  mutate(name = gsub("_", " ", name),
         signature = factor(signature, levels = levels),
         name = factor(name, levels = c("blood normal", "colon normal", "colon IBD","colon POLD1",
                                        "colon POLE", "lung non-smoker", "lung ex-smoker", "lung smoker")))

sig_contribution_plot = ggplot(mean_signatures_ordered, aes(x = name, y = contribution, fill = signature)) +
  geom_col() +
  ggsci::scale_fill_igv() +
  labs(x = NULL) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(y = "relative signature contribution", title = "Signatures active across the 8 cohorts")
sig_contribution_plot
ggsave("manuscript/Supplementary_Figures/Figure_S2/Figure_S2.png" , sig_contribution_plot,
       width = 8, height = 5, bg = "white")

#  TODO use for the blood signature SBS1, SBS5 and SBS blood.