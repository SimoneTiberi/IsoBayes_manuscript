build_data_violin_plot = function(df){
  ths_tpm=0#1e-16 #min(df$tpm_validation[df$tpm_validation>0])/10#1e-12
  ths_p = 0#1e-16 #min(df$P_Y_validation[df$P_Y_validation>0])/10#1e-12
  
  df$P_TPM = (df$tpm_validation+ths_tpm) / sum(df$tpm_validation+ths_tpm)
  df$P_Y_validation = (df$P_Y_validation+ths_p) / sum(df$P_Y_validation+ths_p)
  df$Log2_FC_validation = log2((df$P_Y_validation / df$P_TPM))
  
  df = na.omit(df) #remove 0/0
  df$Log2_FC_validation[df$Log2_FC_validation == Inf] = max(df$Log2_FC_validation[df$Log2_FC_validation<Inf])
  df$Log2_FC_validation[df$Log2_FC_validation == -Inf] = min(df$Log2_FC_validation[df$Log2_FC_validation>-Inf])
  
  df
}

convert_numeric_to_class = function(df, quantiles){
  quantiles[1] = -Inf
  class_Prob_prot_inc = rep("0", nrow(df))
  
  for (i in seq_len(length(quantiles)-1)) {
    sel = which(df$Prob_prot_inc > quantiles[i] & df$Prob_prot_inc <= quantiles[i+1])
    quantiles[1] = 0
    class_Prob_prot_inc[sel] = glue("({round(quantiles[i], 2)} ; {round(quantiles[i+1], 2)}]")
  }
  df$class_Prob_prot_inc = class_Prob_prot_inc
  
  df
}
