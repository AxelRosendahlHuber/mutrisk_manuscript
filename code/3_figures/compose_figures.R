library(patchwork)
library(tidyverse)
library(figpatch)
source("code/0_functions/analysis_variables.R")

# Figure 1:
# load figure 1 graphic:
mg = 7 # set margins
figure_1A <- fig("manuscript/figure_panels/figure_1/Figure 1A.png") |> prep_plot(label = "A", all_margin = 0)
figure_1B = readRDS("manuscript/figure_panels/figure_1/figure_1B.rds") |> prep_plot(label = "B", all_margin = mg)
figure_1C = readRDS("manuscript/figure_panels/figure_1/figure_1C.rds") |> prep_plot(label = "C", all_margin = mg)
figure_1D = readRDS("manuscript/figure_panels/figure_1/figure_1D.rds") |> prep_plot(label = "D", all_margin = mg)

figure_1_middle = figure_1B + plot_spacer() + plot_layout(widths = c(2.5, 1))
figure_1_bottom = figure_1C + figure_1D + plot_layout(widths = c(2.5, 1))
figure_1  = figure_1A / figure_1_middle / figure_1_bottom + plot_layout(heights =  c(1,1, 1))

ggsave("manuscript/Figure_1/figure_1.png", figure_1, width = 15, height = 13)
ggsave("manuscript/Figure_1/figure_1.pdf", figure_1, width = 15, height = 13)

#### Figure 2
mg = 10 # set margins
figure_2A = readRDS("manuscript/figure_panels/figure_2/figure_2A.rds") |> prep_plot(label = "A", all_margin = mg)
figure_2B = fig("manuscript/figure_panels/figure_2/Figure_2B_small.png") |> prep_plot(label = "B", all_margin = 0)
figure_2C = readRDS("manuscript/figure_panels/figure_2/figure_2C.rds") |> prep_plot(label = "C", all_margin = mg)
figure_2D = readRDS("manuscript/figure_panels/figure_2/figure_2D.rds") |> prep_plot(label = "D", all_margin = mg)

figure_2_middle =  figure_2B + figure_2C
figure_2 = figure_2A / figure_2_middle / figure_2D
ggsave("manuscript/Figure_2/figure_2.png", figure_2, width = 12, height = 12)
ggsave("manuscript/Figure_2/figure_2.pdf", figure_2, width = 12, height = 12)

#### Figure 3
mg = 8

figure_3A = readRDS("manuscript/figure_panels/figure_3/figure_3A.rds") |> prep_plot(label = "A", all_margin = mg)
figure_3B = fig("manuscript/figure_panels/figure_3/figure_3B.png") |> prep_plot(label = "B", all_margin = mg)
figure_3_top = figure_3A + figure_3B + plot_layout(widths = c(2,1))

tissue_plots_raw = readRDS("manuscript/figure_panels/figure_3/figure_3CDE.rds")
tissue_plots = mapply(prep_plot, tissue_plots_raw, label = c("C", "D", "E"), all_margin = mg)

figure_3_bottom = wrap_plots(tissue_plots, nrow = 1, widths = c(4, 3, 1.2))
figure_3 = figure_3_top / figure_3_bottom + plot_layout(heights = c(1.2, 1))

ggsave("manuscript/Figure_3/figure_3.png", figure_3, width = 15, height = 10)
ggsave("manuscript/Figure_3/figure_3.pdf", figure_3, width = 15, height = 10)


##### Figure 4
mg = 5

# perform these operations later in the compose-figures part:
list_figure_4AB = readRDS("manuscript/figure_panels/figure_4/figures_AB.rds")
figure_4A = list_figure_4AB[[1]] |> prep_plot(label = 'A', all_margin = mg)
figure_4B = list_figure_4AB[[2]]  |> prep_plot(label = 'B', all_margin = mg)
figure_4_top = figure_4A + figure_4B

figs = readRDS("manuscript/figure_panels/figure_4/figures_C-G.rds")
annotated_figs = lapply(names(figs), \(x) prep_plot(figs[[x]], substr(x, 3,3), all_margin = mg))
figure_4_middle = wrap_plots(annotated_figs, nrow = 1)

figure_4H = readRDS("manuscript/figure_panels/figure_4/figures_adenoma.rds") |>
  prep_plot(label = 'H', all_margin = mg)

figures_4_CRC = readRDS("manuscript/figure_panels/figure_4/figure_4_CRC.rds")

annotated_figs = mapply(prep_plot, figures_4_CRC, c("I", "J", "K", "L", "M"), all_margin = 2)
figure_4_bottom = wrap_plots(c(figure_4H, annotated_figs), nrow = 1)

figure_4 = figure_4_top / figure_4_middle / figure_4_bottom# + plot_layout(heights = c(1.5,1, 1))

ggsave("manuscript/Figure_4/figure_4.png", figure_4, width = 22, height = 11)
ggsave("manuscript/Figure_4/figure_4.pdf", figure_4, width = 22, height = 11)

##### Figure 5
figure_5A = readRDS("manuscript/figure_panels/figure_5/figure_5A.rds") |> prep_plot("A", all_margin = 3)
figure_5B = readRDS("manuscript/figure_panels/figure_5/figure_5B.rds")
figure_5C = readRDS("manuscript/figure_panels/figure_5/figure_5C.rds")
figure_5D = readRDS("manuscript/figure_panels/figure_5/figure_5D.rds")

# save final completed plot
figure_5A / figure_5B / figure_5C / figure_5D
ggsave("manuscript/Figure_5/Figure_5.png", width = 13, height = 15)

#TODO check if it is neccesary to check extended figures;











