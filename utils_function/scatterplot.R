scatterplot = function(abundances){
  colnames(abundances) = c("x", "y")
  corr = cor(abundances)
  pp = ggplot(abundances, aes(x=x, y=y)) + 
    geom_point() + 
    labs(title = DATA, subtitle = glue("Correlation: {round(corr[1, 2], 2)}")) +
    theme_bw() +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          legend.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"))
  
  pp
}
