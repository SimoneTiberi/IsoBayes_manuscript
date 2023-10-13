plot_prob_change_group = function(benchmark_df, violin = FALSE){
  vec_up = benchmark_df$Log2_FC_validation[benchmark_df$class_Prob_prot_inc == "(0.8 ; 1]"]
  vec_down = benchmark_df$Log2_FC_validation[benchmark_df$Prob_prot_inc == "(0 ; 0.2]"]
  iqr_up = IQR(vec_up)
  iqr_down = IQR(vec_down)

  pp = ggplot(benchmark_df, aes(class_Prob_prot_inc, Log2_FC_validation, fill = Model))
  if(violin){
    pp = pp + geom_violin()
  }else{
    pp = pp + geom_boxplot(varwidth=T, outlier.shape = NA)
  }
  pp = pp + labs(title = DATA, 
         x = latex2exp::TeX("$Pr\\left(\\pi_p > \\pi_p^T\\right)$"),
         y = "Validated Log2FC") +
    theme_bw() +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10)) +
    scale_y_continuous(n.breaks = 8,
                       limits = c(quantile(vec_down, 0.25) - 2 * iqr_down,
                                  quantile(vec_up, 0.75) + 2 * iqr_up)
    )
  
  pp
}
