scatterplot = function(abundances){
  corr_log_abundance = cor(abundances)
  pp = ggplot(abundances, aes(x=log10_abundance, y=log10_abundance_validated)) + 
    geom_point() + 
    labs(title = "Abundance scatterplot", subtitle = glue("Correlation: {round(corr_log_abundance[1, 2], 2)}")) +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold")) +
    theme_bw()
  
  pp
}

