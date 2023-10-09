plot_prob_change = function(benchmark_df, violin = FALSE){
  #vec_up = benchmark_df$Log2_FC_validation[benchmark_df$Prob_prot_inc > 0.8]
  #vec_down = benchmark_df$Log2_FC_validation[benchmark_df$Prob_prot_inc < 0.2]
  #iqr_up = IQR(vec_up)
  #iqr_down = IQR(vec_down)
  
  pp = ggplot(benchmark_df, aes(class_Prob_prot_inc, Log2_FC_validation))
  if(violin){
    pp = pp + geom_violin(fill="plum")
  }else{
    pp = pp + geom_boxplot(fill="plum", varwidth=T, outlier.shape = NA)
  }
  pp = pp + labs(title = DATA,
                 subtitle = "Change in protein and mRNA isoform relative abundances",
                 x = "P(Protein > Transcript)",
                 y = "Log2FC of isoform relative abundances") +
    theme_bw() +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10)) #+
    #coord_cartesian(ylim = c(quantile(vec_down, 0.25) - 1.5 * iqr_down,
     #                        quantile(vec_up, 0.75) + 1 * iqr_up)
    #)
  pp
}
