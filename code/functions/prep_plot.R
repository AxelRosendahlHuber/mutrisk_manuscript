library(ggpubr)
prep_plot = function(plot, label, all_margin, t = 5, r = 5 , l = 5, b = 5) {

  # if the 'all margin' comment is mentioned, this means that individual margin settings can be replaced with one
  if (!missing(all_margin)) {
    t = r = l = b = all_margin
  }

  plot = plot + theme(plot.margin = margin(t,r,l, b, unit = "mm"))
  annotate_figure(plot, fig.lab = label, fig.lab.size = 20, fig.lab.face = "bold")
}
