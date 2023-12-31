build_data_violin_plot = function(df, ths){
  #remove 0/0
  sel = df$tpm_validation == 0 & df$P_Y_validation == 0
  df = df[!sel, ]
  
  df$P_TPM = (df$tpm_validation+ths) / sum(df$tpm_validation+ths)
  df$P_Y_validation = (df$P_Y_validation+ths) / sum(df$P_Y_validation+ths)
  df$Log2_FC_validation = log2((df$P_Y_validation / df$P_TPM))
  
  df
}

convert_numeric_to_class = function(df, var, quantiles){
  quantiles[1] = -Inf
  class_Prob_prot_inc = rep("0", nrow(df))
  lab_levels = c()
  
  for (i in seq_len(length(quantiles)-1)) {
    sel = which(df[, var] > quantiles[i] & df[, var] <= quantiles[i+1])
    quantiles[1] = 0
    first_brack = "("
    if(i == 1){
      first_brack = "["
    }
    lab_level = glue("{first_brack}{round(quantiles[i], 2)}, {round(quantiles[i+1], 2)}]")
    class_Prob_prot_inc[sel] = lab_level
    lab_levels = c(lab_levels, lab_level)
  }
  df$class_Prob_prot_inc = class_Prob_prot_inc
  df$class_Prob_prot_inc = factor(df$class_Prob_prot_inc, levels = lab_levels)
  
  df
}

adjust_inf_log2FC = function(df, ths){
  p_tpm_adj = (df$TPM+ths)/sum(df$TPM+ths)
  pi = (df$Pi+ths)/sum(df$Pi+ths)
  df$Log2_FC_adj = log2(pi/p_tpm_adj)
  
  df
}
