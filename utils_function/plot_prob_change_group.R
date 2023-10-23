plot_prob_change_group = function(benchmark_df, violin = FALSE){
  
  q25 = unlist(lapply(unique(benchmark_df$class_Prob_prot_inc), FUN = function(x){
    quantile(benchmark_df$Log2_FC_validation[benchmark_df$class_Prob_prot_inc == x], 0.25)
  }))
  Min = unlist(lapply(unique(benchmark_df$class_Prob_prot_inc), FUN = function(x){
    quantile(benchmark_df$Log2_FC_validation[benchmark_df$class_Prob_prot_inc == x], 0)
  }))
  q75 = unlist(lapply(unique(benchmark_df$class_Prob_prot_inc), FUN = function(x){
    quantile(benchmark_df$Log2_FC_validation[benchmark_df$class_Prob_prot_inc == x], 0.75)
  }))
  Max = unlist(lapply(unique(benchmark_df$class_Prob_prot_inc), FUN = function(x){
    quantile(benchmark_df$Log2_FC_validation[benchmark_df$class_Prob_prot_inc == x], 1)
  }))
  
  bottom_whisker = q25 - (q75 - q25) * 1.5
  top_whisker = q75 + (q75 - q25) * 1.5
  
  y_bottom = max(min(bottom_whisker, na.rm = TRUE), Min[which.min(bottom_whisker)], na.rm = TRUE)
  y_top = min(max(top_whisker, na.rm = TRUE), Max[which.max(top_whisker)], na.rm = TRUE)

  pp = ggplot(benchmark_df, aes(class_Prob_prot_inc, Log2_FC_validation, fill = Model)) +
    geom_boxplot(outlier.shape = NA)
  
  pp = pp + labs(title = DATA_name, 
         x = latex2exp::TeX("$Pr \\left(\\pi_p > \\pi_p^T \\right)$"),
         y = "Validated Log2-FC") +
    theme_bw() +
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          axis.title = element_text(size = 15),
          legend.title = element_blank(),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          legend.text = element_text(size = 11)) +
    scale_y_continuous(n.breaks = 8) +
    coord_cartesian(ylim = c(y_bottom, y_top))
  
  pp
}
