
library(scales)
plot_mrate = function(mrate_df, colors, title = NULL ) {
  mrate_df |> 
    ggplot(aes(x = age, y = mle, fill = category)) + 
    geom_pointrange(aes(ymin = cilow, ymax = cihigh, color = category)) + 
    scale_fill_manual(values = colors) + 
    scale_color_manual(values = colors) +
    scale_y_continuous(labels = comma) + 
    cowplot::theme_cowplot() + 
    labs(y = "number of mutated cells", subtitle = title,
         x = "Age (years)", color = NULL, fill = NULL) + 
    theme(legend.position = "inside", legend.position.inside = c(0.05,0.8 ), 
          plot.title = element_text(size = 11))
}

plot_mrate_sig = function(mrate_df, colors, title = NULL ) {
  
  mrate_df |> 
    ggplot(aes(x = age, y = mle, fill = signature)) + 
    geom_col() + 
    ggsci::scale_fill_igv() + 
    cowplot::theme_cowplot() + 
    scale_y_continuous(expand=expansion(mult=c(0,0.1)), labels = comma) + 
    labs(y = "number of mutated cells", subtitle = title, x = "Age (years)", fill = NULL) + 
    theme(legend.position = "inside", legend.position.inside = c(0.05,0.7), 
          legend.key.size = unit(4, "mm"))
}
