library(ggpubr)
prep_plot = function(plot, label, t = 5, r = 5 , l = 5, b = 5) {

  plot = plot + theme(plot.margin = margin(t,r,l, b, unit = "mm"))
  annotate_figure(plot, fig.lab = label, fig.lab.size = 20, fig.lab.face = "bold")
}