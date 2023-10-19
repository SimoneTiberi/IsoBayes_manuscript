plot_prob_change = function(benchmark_df, violin = FALSE){
  
  q25 = unlist(lapply(levels(benchmark_df$class_Prob_prot_inc), FUN = function(x){
    quantile(benchmark_df$Log2_FC_validation[benchmark_df$class_Prob_prot_inc == x], 0.25)
  }))
  Min = unlist(lapply(levels(benchmark_df$class_Prob_prot_inc), FUN = function(x){
    min(benchmark_df$Log2_FC_validation[benchmark_df$class_Prob_prot_inc == x])
  }))
  q75 = unlist(lapply(levels(benchmark_df$class_Prob_prot_inc), FUN = function(x){
    quantile(benchmark_df$Log2_FC_validation[benchmark_df$class_Prob_prot_inc == x], 0.75)
  }))
  Max = unlist(lapply(levels(benchmark_df$class_Prob_prot_inc), FUN = function(x){
    max(benchmark_df$Log2_FC_validation[benchmark_df$class_Prob_prot_inc == x])
  }))
  
  i_min = which.min(q25)
  i_max = which.max(q75)
  
  bottom_whisker = max(q25[i_min] - (q75[i_min] - q25[i_min]) * 1.5, Min[i_min])
  top_whisker = min(q75[i_max] + (q75[i_max] - q25[i_max]) * 1.5, Max[i_max])
  
  pp = ggplot(benchmark_df, aes(class_Prob_prot_inc, Log2_FC_validation))
  if(violin){
    pp = pp + geom_violin(fill="plum")
  }else{
    pp = pp + geom_boxplot(fill="plum", varwidth=T, outlier.shape = NA)
  }
  pp = pp + labs(title = DATA_name,
                 x = latex2exp::TeX("$Pr \\left(\\pi_p > \\pi_p^T \\right)$"),
                 y = "Validated Log2(FC)") +
    theme_bw() +
    theme(plot.title = element_text(size = 15, hjust = 0.5),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          legend.text = element_text(size = 11)) +
    scale_y_continuous(n.breaks = 8, limits = c(bottom_whisker, top_whisker))
  
  pp
}
