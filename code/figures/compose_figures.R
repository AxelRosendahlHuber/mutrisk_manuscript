library(patchwork)
library(tidyverse)
library(figpatch)
source("code/functions/analysis_variables.R")

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
figure_3 = figure_3_top / figure_3_bottom

ggsave("manuscript/Figure_3/figure_3.png", figure_3, width = 15, height = 10)
ggsave("manuscript/Figure_3/figure_3.pdf", figure_3, width = 15, height = 10)


##### Figure 4




##### Figure 5













