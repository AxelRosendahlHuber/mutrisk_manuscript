# Create overview figure given a 'cell muts' object.


plot_vaf = function(cell_muts, sample_name) {

  n_clones = cell_muts |> filter(donor %in% sample_name) |> n_distinct()
  number_samples = ifelse(n_clones >40, n_clones, 40) # if less than 40 clones in the sample, pick all of them

  set.seed(10)
  cell_muts |>
    filter(donor %in% sample_name) |>
    filter(sampleID %in% base::sample(unique(sampleID), number_samples)) |>
    mutate(y = as.numeric(as.factor(sampleID))) |>
    group_by(y) |>
    ggplot(aes(x = vaf, y = y, group = y)) +
    ggridges::geom_density_ridges() +
    geom_vline(xintercept = 0.5, linetype = "dashed") +
    theme_cowplot() +
    labs(title = paste0("40 clones from donor ",sample_name), y = "individual samples", x = "VAF of individiual mutations across sample")

}

create_vaf_overview = function(cell_muts, sample_names) {

  if (length(sample_names) != 2) {
    stop("Exactly 2 sample names must be provided")
  }

  VAF_per_sample = cell_muts |>
    group_by(donor, sampleID) |>
    summarize(mean_vaf = mean(vaf)) |>
    ggplot(aes(x = donor, y = mean_vaf)) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    ggbeeswarm::geom_quasirandom(size = 0.8) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_y_continuous(limits = c(0, 0.75),
                       breaks = scales::extended_breaks(4)) +
    labs(y = "average VAF per sample") |>
    prep_plot("A")


  vafplot1 = plot_vaf(cell_muts, sample_names[1]) |> prep_plot("B")
  vafplot2 = plot_vaf(cell_muts, sample_names[2]) |> prep_plot("C")


  VAF_per_sample / (vafplot1 + vafplot2)

}



SX001_plot = cell_muts |>
  filter(donor == "SX001") |>
  filter(sampleID %in% base::sample(unique(sampleID), 40)) |>
  mutate(y = as.numeric(as.factor(sampleID))) |>
  group_by(y) |>
  ggplot(aes(x = vaf, y = y, group = y)) +
  ggridges::geom_density_ridges() +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme_cowplot() +
  labs(title = "40 clones from donor AX001", y = "HSC/MPP clones", x = "VAF of individiual mutations across sample")

supplementary_note_plot_blood = VAF_per_sample / (AX001_plot | SX001_plot)
ggsave("manuscript/Supplementary_notes/Supplementary_Note_X/figure_blood.png", supplementary_note_plot_blood, width = 10, height = 10)
