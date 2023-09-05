PATH_WD = commandArgs(trailingOnly = TRUE)[1]
DATA = commandArgs(trailingOnly = TRUE)[2]

# Set path, global variable and libraries
###########################################################################################
library(glue)
library(ggplot2)
source(glue("{PATH_WD}/utils_function/plot_prob_change.R"))
source(glue("{PATH_WD}/utils_function/plot_prob_change_group.R"))
source(glue("{PATH_WD}/utils_function/utils_benchmarking.R"))
source(glue("{PATH_WD}/utils_function/log_output.R"))

PATH_DATA = glue("{PATH_WD}/Data/{DATA}") 
PATH_RES_MODEL = glue("{PATH_WD}/Model_results/{DATA}")
PATH_RES_CHANGE = glue("{PATH_WD}/Change_protein_mRNA_isoform/{DATA}")
log_output(glue("change_prot_mRNA_iso_{DATA}"))


main = function(proteases){
  ######################################################
  # Changes in protein and mRNA isoform relative abundances (MM)
  ######################################################
  benchmark_df_all = list()
  
  for (protease in proteases) {
    # load validation dataset from metamorpheus
    load(glue("{PATH_DATA}/No{protease}/Validation_prot_psm"))
    # load model results
    load(glue("{PATH_RES_MODEL}/MM_psm_mRNA_PEP/{protease}/MM_psm_mRNA_PEP_MCMC.RData"))
    
    benchmark_df = merge(VALIDATION_DF_prot, res$isoform_results,
                         by.x = "proteins", by.y = "Isoform", all = T)
    
    # Eliminio isoforme non presenti nel validation set e in input
    benchmark_df = na.omit(benchmark_df)
    
    ths_tpm = min(benchmark_df$TPM[benchmark_df$TPM>0])
    benchmark_df$P_TPM = (benchmark_df$TPM+ths_tpm) / sum(benchmark_df$TPM+ths_tpm)
    benchmark_df$Log2_FC =log2((benchmark_df$Pi / benchmark_df$P_TPM)+1)
    benchmark_df_all = rbind(benchmark_df_all, benchmark_df)
    
    #quantiles = quantile(benchmark_df$Prob_prot_inc, seq(0, 1, 0.2))
    quantiles = seq(0, 1, 0.2)
    quantiles[1] = -Inf
    class_Prob_prot_inc = rep("0", nrow(benchmark_df))
    
    for (i in seq_len(length(quantiles)-1)) {
      sel = which(benchmark_df$Prob_prot_inc > quantiles[i] & benchmark_df$Prob_prot_inc <= quantiles[i+1])
      quantiles[1] = 0
      class_Prob_prot_inc[sel] = glue("({round(quantiles[i], 2)} ; {round(quantiles[i+1], 2)}]")
    }
    benchmark_df$class_Prob_prot_inc = class_Prob_prot_inc
    
    
    plot_change = plot_prob_change(benchmark_df)
    ggsave(glue("{PATH_RES_CHANGE}/{protease}/main_result.png"), plot = plot_change)
  }
  
  #quantiles = quantile(benchmark_df_all$Prob_prot_inc, seq(0, 1, 0.2))
  quantiles = seq(0, 1, 0.2)
  quantiles[1] = -Inf
  class_Prob_prot_inc = rep("0", nrow(benchmark_df_all))
  
  for (i in seq_len(length(quantiles)-1)) {
    sel = which(benchmark_df_all$Prob_prot_inc > quantiles[i] & benchmark_df_all$Prob_prot_inc <= quantiles[i+1])
    quantiles[1] = 0
    class_Prob_prot_inc[sel] = glue("({round(quantiles[i], 2)} ; {round(quantiles[i+1], 2)}]")
  }
  benchmark_df_all$class_Prob_prot_inc = class_Prob_prot_inc
  
  plot_change = plot_prob_change(benchmark_df_all)
  ggsave(glue("{PATH_RES_CHANGE}/main_result.png"), plot = plot_change)
  
  ######################################################
  # Changes in protein and mRNA isoform relative abundances - PEP-noPEP
  ######################################################
  benchmark_df_all = list()
  
  for (protease in proteases) {
    benchmark_df_models = list()
    for (model in c("_PEP", "")) {
      # load validation dataset from metamorpheus
      load(glue("{PATH_DATA}/No{protease}/Validation_prot_psm"))
      # load model results
      load(glue("{PATH_RES_MODEL}/MM_psm_mRNA{model}/{protease}/MM_psm_mRNA{model}_MCMC.RData"))
      
      benchmark_df = merge(VALIDATION_DF_prot, res$isoform_results[, c("Isoform", "Pi", "TPM", "Log2_FC", "Prob_prot_inc")],
                           by.x = "proteins", by.y = "Isoform", all = T)
      
      # Eliminio isoforme non presenti nel validation set e in input
      benchmark_df = na.omit(benchmark_df)
      
      ths_tpm = min(benchmark_df$TPM[benchmark_df$TPM>0])
      benchmark_df$P_TPM = (benchmark_df$TPM+ths_tpm) / sum(benchmark_df$TPM+ths_tpm)
      benchmark_df$Log2_FC = log2((benchmark_df$Pi / benchmark_df$P_TPM)+1)
      if(model == "_PEP"){
        label_model = "PEP"
      }else{
        label_model = "No PEP"
      }
      benchmark_df$Model = label_model
      
      #quantiles = quantile(benchmark_df$Prob_prot_inc, seq(0, 1, 0.2))
      quantiles = seq(0, 1, 0.2)
      quantiles[1] = -Inf
      class_Prob_prot_inc = rep("0", nrow(benchmark_df))
      
      for (i in seq_len(length(quantiles)-1)) {
        sel = which(benchmark_df$Prob_prot_inc > quantiles[i] & benchmark_df$Prob_prot_inc <= quantiles[i+1])
        quantiles[1] = 0
        class_Prob_prot_inc[sel] = glue("({round(quantiles[i], 2)} ; {round(quantiles[i+1], 2)}]")
      }
      benchmark_df$class_Prob_prot_inc = class_Prob_prot_inc
      
      benchmark_df_models = rbind(benchmark_df_models, benchmark_df)
    }
    plot_change = plot_prob_change_group(benchmark_df_models)
    ggsave(glue("{PATH_RES_CHANGE}/{protease}/PEP_no_PEP.png"), plot = plot_change)
    
    benchmark_df_all = rbind(benchmark_df_all, benchmark_df_models)
  }
  plot_change = plot_prob_change_group(benchmark_df_all)
  ggsave(glue("{PATH_RES_CHANGE}/PEP_no_PEP.png"), plot = plot_change)
  
  ######################################################
  # Changes in protein and mRNA isoform relative abundances - MM vs OpenMS
  ######################################################
  benchmark_df_all = list()
  
  for (protease in proteases) {
    benchmark_df_models = list()
    for (model in c("MM_psm", "OpenMS")) {
      # load validation dataset from metamorpheus
      load(glue("{PATH_DATA}/No{protease}/Validation_prot_psm"))
      # load model results
      load(glue("{PATH_RES_MODEL}/{model}_mRNA_PEP/{protease}/{model}_mRNA_PEP_MCMC.RData"))
      
      benchmark_df = merge(VALIDATION_DF_prot, res$isoform_results[, c("Isoform", "Pi", "TPM", "Log2_FC", "Prob_prot_inc")],
                           by.x = "proteins", by.y = "Isoform", all = T)
      
      if(model == "OpenMS"){
        # proteine presenti in input
        iso_input = get_score_from_idXML(glue("{PATH_DATA}/Only{protease}/merge_index_percolator_pep_switched.idXML"))
        benchmark_df = merge(benchmark_df, iso_input, by.x = "proteins", by.y = "Isoform", all = T)
        
        # Eliminio isoforme non presenti nel validation set
        benchmark_df = benchmark_df[!is.na(benchmark_df$Present), ]
        
        # Tengo quelle in input
        benchmark_df = benchmark_df[!is.na(benchmark_df$score), ]
        
        # 0 se modello non trova l'isoforma
        benchmark_df[is.na(benchmark_df)] = 0
        benchmark_df$score = NULL
      } else {
        # Eliminio isoforme non presenti nel validation set e in input
        benchmark_df = na.omit(benchmark_df)
      }
      
      ths_tpm = min(benchmark_df$TPM[benchmark_df$TPM>0])
      benchmark_df$P_TPM = (benchmark_df$TPM+ths_tpm) / sum(benchmark_df$TPM+ths_tpm)
      benchmark_df$Log2_FC = log2((benchmark_df$Pi / benchmark_df$P_TPM)+1)
      benchmark_df$Model = model
      
      #quantiles = quantile(benchmark_df$Prob_prot_inc, seq(0, 1, 0.2))
      quantiles = seq(0, 1, 0.2)
      quantiles[1] = -Inf
      class_Prob_prot_inc = rep("0", nrow(benchmark_df))
      
      for (i in seq_len(length(quantiles)-1)) {
        sel = which(benchmark_df$Prob_prot_inc > quantiles[i] & benchmark_df$Prob_prot_inc <= quantiles[i+1])
        quantiles[1] = 0
        class_Prob_prot_inc[sel] = glue("({round(quantiles[i], 2)} ; {round(quantiles[i+1], 2)}]")
      }
      benchmark_df$class_Prob_prot_inc = class_Prob_prot_inc
      
      benchmark_df_models = rbind(benchmark_df_models, benchmark_df)
    }
    plot_change = plot_prob_change_group(benchmark_df_models)
    ggsave(glue("{PATH_RES_CHANGE}/{protease}/OpenMS_vs_MM.png"), plot = plot_change)
    
    benchmark_df_all = rbind(benchmark_df_all, benchmark_df_models)
  }
  plot_change = plot_prob_change_group(benchmark_df_all)
  ggsave(glue("{PATH_RES_CHANGE}/OpenMS_vs_MM.png"), plot = plot_change)
}

main(proteases = list.dirs(PATH_RES_CHANGE, recursive = FALSE, full.names = FALSE))
