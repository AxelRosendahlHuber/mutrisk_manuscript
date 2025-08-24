# use own 96-trinucleotide matrix functions
plot_96_profile2 = function (mut_matrix, relative = FALSE, free_y = TRUE, horizontal_labels = FALSE) {

  mut_matrix = as.data.frame(mut_matrix)
  number_sigs = dim(mut_matrix)[2]
  label_y = "Number of mutations"

  if (relative == TRUE) {

    label_y = "Relative contribution"

    if (number_sigs == 1 ) {
      norm_mut_matrix = mut_matrix/sum(mut_matrix)
      strand = sapply(rownames(mut_matrix), function(x) strsplit(x, "-")[[1]][2])
    } else if (number_sigs > 1) {
      norm_mut_matrix = prop.table(as.matrix(mut_matrix), 2)
      strand = sapply(rownames(mut_matrix), function(x) strsplit(x, "-")[[1]][2])
    }
    mut_matrix = norm_mut_matrix
  }

  tb <- mut_matrix %>% as.data.frame() %>% tibble::rownames_to_column("full_context") %>%
    dplyr::mutate(substitution = stringr::str_replace(full_context, "\\w\\[(.*)\\]\\w", "\\1"),
                  context = stringr::str_replace(full_context,"\\[.*\\]", "\\.")) %>%
    dplyr::select(-full_context) %>%
    tidyr::pivot_longer(c(-substitution, -context), names_to = "sample", values_to = "freq") %>%
    dplyr::mutate(sample = factor(sample, levels = unique(sample)))


  width <- 0.6

  scales = "free"
  if (free_y == TRUE) {
    scales = "free_y"
  }


  sample_df = tb |>
    dplyr::select(sample) |> distinct() |>
    mutate(substitution = "T>G")

  plt = ggplot(data = tb, aes(x = context, y = freq, fill = substitution, width = width)) +
    geom_bar(stat = "identity", linewidth = 0.2) +
    scale_y_continuous(expand =  expansion(mult = c(0, 0.1))) +
    ggh4x::facet_grid2(sample ~ substitution, scales = scales , axes = "x",remove_labels = "x",
                       strip = ggh4x::strip_themed(background_x = list(element_rect(fill = "#2EBAED"),
                          element_rect(fill = "#000000"),
                          element_rect(fill = "#DE1C14"),
                          element_rect(fill = "#D4D2D2"),
                          element_rect(fill = "#ADCC54"),
                          element_rect(fill = "#F0D0CE"))))  +
    scale_fill_manual(values = COLORS6) +
    guides(fill = "none") +
    theme_classic() +
     theme(
       panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_text(colour = "white", size = 12),
        strip.background = element_rect(colour = "white"),
        axis.line = element_line(linewidth = 0.3),
       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(x = NULL, y = label_y) +
    coord_cartesian(clip = 'off')

  if (horizontal_labels) {
    plt = plt +
      ggpp::geom_text_npc(data = sample_df, mapping = aes(npcx = 0.5, npcy = 0.95, label = sample), size = 4) +
      theme(strip.text.y = element_blank())
  }


  plt + scale_x_discrete(labels=gsub("\\[...\\]", ".", TRIPLETS_96))

  return(plt)
}
