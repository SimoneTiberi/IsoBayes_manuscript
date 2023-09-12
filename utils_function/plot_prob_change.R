plot_prob_change = function(benchmark_df){
  pp = ggplot(benchmark_df, aes(class_Prob_prot_inc, Log2_FC_validation)) + 
    geom_violin(trim=FALSE, fill="plum") +
    #geom_boxplot(fill="plum", width = 0.1) +
    labs(title = "Change in protein and mRNA isoform relative abundances", 
         x = "P(Protein > Transcript)",
         y = "Log2FC of isoform relative abundances") +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10)) +
    theme_bw() +
    coord_cartesian(ylim = c(- 10, 10))
  
  pp
}
