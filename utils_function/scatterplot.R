scatterplot = function(abundances){
  colnames(abundances) = c("x", "y")
  corr = cor(abundances)
  pp = ggplot(abundances, aes(x=x, y=y)) + 
    geom_point() + 
    labs(title = DATA_name, subtitle = glue("Correlation: {round(corr[1, 2], 2)}")) +
    theme_bw() +
    theme(plot.title = element_text(size = 15, hjust = 0.5),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.title = element_text(size = 11),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11))
  
  pp
}
