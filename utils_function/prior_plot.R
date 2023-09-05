prior_plot = function(data_to_plot){
  pp = ggplot(data_to_plot, aes(x=grid, y=AUC, group = model)) +
    geom_line(linewidth=0.8) + geom_point() +
    scale_y_continuous(n.breaks = 10,
                       limits = c(max(min(data_to_plot$AUC)-0.025, 0),
                                  min(max(data_to_plot$AUC)+0.025, 1))) +
    labs(title = "Prior robustness", x = "Prior values") +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10)) +
    theme_bw()
  
  pp
}