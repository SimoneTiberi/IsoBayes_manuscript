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
source(glue("{PATH_WD}/utils_function/scatterplot.R"))
source(glue("{PATH_WD}/utils_function/utils_change_iso_mrna.R"))

PATH_DATA = glue("{PATH_WD}/Data/{DATA}") 
PATH_RES_MODEL = glue("{PATH_WD}/Model_results/{DATA}")
PATH_RES_CHANGE = glue("{PATH_WD}/Change_protein_mRNA_isoform/{DATA}")
log_output(glue("change_prot_mRNA_iso_{DATA}"))

main = function(proteases, no_unique = FALSE){
  inputs = c("OpenMS", "MM_psm")
  quantiles = c(0, 0.2, 0.4, 0.6, 0.8, 1)
  
  for (input in inputs) {
    ######################################################
    # Changes in protein and mRNA isoform relative abundances (MM)
    ######################################################
    benchmark_df_all = list()
    
    for (protease in proteases) {
      # load validation data with res
      load(glue("{PATH_RES_MODEL}/{input}_mRNA_PEP/{protease}/Merged_validation_res_{input}_mRNA_PEP"))
      colnames(validation_dat)[grep("tpm", colnames(validation_dat))] = "tpm_validation"
      if(no_unique){
        validation_dat = validation_dat[validation_dat$Y_unique == 0, ]
        no_unique_nm = "_no_unique"
      }else{
        no_unique_nm = ""
      }
      validation_dat = build_data_violin_plot(validation_dat)
      benchmark_df_all = rbind(benchmark_df_all, validation_dat)
      validation_dat = convert_numeric_to_class(validation_dat, quantiles)
      validation_dat_extreme = convert_numeric_to_class(validation_dat, quantiles = c(0, 0.01, 0.99, 1))
      validation_dat_extreme = validation_dat_extreme[validation_dat_extreme$class_Prob_prot_inc != "c(0.01 ; 0.99]", ]
      
      plot_change = plot_prob_change(validation_dat_extreme)
      ggsave(glue("{PATH_RES_CHANGE}/{protease}/main_result_{input}{no_unique_nm}.png"), plot = plot_change)
      
      plot_change = plot_prob_change(validation_dat)
      ggsave(glue("{PATH_RES_CHANGE}/{protease}/main_result_extreme_{input}{no_unique_nm}.png"), plot = plot_change)
      
      ######################################################
      # Log2FC correlation
      ######################################################
      plot_scatter = scatterplot(validation_dat[, c("Log2_FC", "Log2_FC_validation")]) + 
        labs(x = "Log2_FC", y = "Log2_FC_validation")
      ggsave(glue("{PATH_RES_CHANGE}/{protease}/scatterplot_log2FC_{input}{no_unique_nm}.png"), plot = plot_scatter)
      
      plot_scatter = scatterplot(validation_dat[, c("Prob_prot_inc", "Log2_FC_validation")]) + 
        labs(x = "Prob_prot_inc", y = "Log2_FC_validation")
      ggsave(glue("{PATH_RES_CHANGE}/{protease}/scatterplot_probInc_log2FC_{input}{no_unique_nm}.png"), plot = plot_scatter)
    }
    benchmark_df_all = convert_numeric_to_class(benchmark_df_all, quantiles)
    benchmark_df_all_extreme = convert_numeric_to_class(benchmark_df_all, quantiles = c(0, 0.01, 0.99, 1))
    benchmark_df_all_extreme = benchmark_df_all_extreme[benchmark_df_all_extreme$class_Prob_prot_inc != "c(0.01 ; 0.99]", ]
    
    plot_change = plot_prob_change(benchmark_df_all)
    ggsave(glue("{PATH_RES_CHANGE}/main_result_{input}{no_unique_nm}.png"), plot = plot_change)
    
    plot_change = plot_prob_change(benchmark_df_all_extreme)
    ggsave(glue("{PATH_RES_CHANGE}/main_result_extreme_{input}{no_unique_nm}.png"), plot = plot_change)
    
    ######################################################
    # Log2FC correlation
    ######################################################
    plot_scatter = scatterplot(benchmark_df_all[, c("Log2_FC", "Log2_FC_validation")]) + 
      labs(x = "Log2_FC", y = "Log2_FC_validation")
    ggsave(glue("{PATH_RES_CHANGE}/scatterplot_log2FC_{input}{no_unique_nm}.png"), plot = plot_scatter)
    
    plot_scatter = scatterplot(benchmark_df_all[, c("Prob_prot_inc", "Log2_FC_validation")]) + 
      labs(x = "Prob_prot_inc", y = "Log2_FC_validation")
    ggsave(glue("{PATH_RES_CHANGE}/scatterplot_probInc_log2FC_{input}{no_unique_nm}.png"), plot = plot_scatter)
  }
  
  ######################################################
  # Changes in protein and mRNA isoform relative abundances - PEP-noPEP
  ######################################################
  for (input in inputs) {
    
    benchmark_df_all = list()
    
    for (protease in proteases) {
      benchmark_df_models = list()
      for (model in c("_PEP", "")) {
        # load validation data with res
        load(glue("{PATH_RES_MODEL}/{input}_mRNA{model}/{protease}/Merged_validation_res_{input}_mRNA{model}"))
        colnames(validation_dat)[grep("tpm", colnames(validation_dat))] = "tpm_validation"
        
        if(no_unique){
          validation_dat = validation_dat[validation_dat$Y_unique == 0, ]
          no_unique_nm = "_no_unique"
        }else{
          no_unique_nm = ""
        }
        
        validation_dat = build_data_violin_plot(validation_dat)
        validation_dat$Gene = NULL
        
        if(model == "_PEP"){
          label_model = "PEP"
        }else{
          label_model = "No PEP"
        }
        validation_dat$Model = label_model
        
        benchmark_df_all = rbind(benchmark_df_all, validation_dat)
        
        validation_dat = convert_numeric_to_class(validation_dat, quantiles)
        benchmark_df_models = rbind(benchmark_df_models, validation_dat)
      }
      
      sub_bench = benchmark_df_models[benchmark_df_models$Log2_FC_validation < Inf & benchmark_df_models$Log2_FC_validation > -Inf & !is.na(benchmark_df_models$Log2_FC_validation), ]
      
      plot_change = plot_prob_change_group(sub_bench)
      ggsave(glue("{PATH_RES_CHANGE}/{protease}/PEP_no_PEP_{input}{no_unique_nm}.png"), plot = plot_change)
    }
    benchmark_df_all = convert_numeric_to_class(benchmark_df_all, quantiles)
    sub_bench = benchmark_df_all[benchmark_df_all$Log2_FC_validation < Inf & benchmark_df_all$Log2_FC_validation > -Inf & !is.na(benchmark_df_all$Log2_FC_validation), ]
    plot_change = plot_prob_change_group(sub_bench)
    ggsave(glue("{PATH_RES_CHANGE}/PEP_no_PEP_{input}{no_unique_nm}.png"), plot = plot_change)
  }
  
  ######################################################
  # Changes in protein and mRNA isoform relative abundances - mRNA-no_mRNA
  ######################################################
  
  ######################################################
  # Changes in protein and mRNA isoform relative abundances - MM vs OpenMS
  ######################################################
  benchmark_df_all = list()
  
  for (protease in proteases) {
    benchmark_df_models = list()
    for (input in c("MM_psm", "OpenMS")) {
      # load validation data with res
      load(glue("{PATH_RES_MODEL}/{input}_mRNA_PEP/{protease}/Merged_validation_res_{input}_mRNA_PEP"))
      colnames(validation_dat)[grep("tpm", colnames(validation_dat))] = "tpm_validation"
      
      if(no_unique){
        validation_dat = validation_dat[validation_dat$Y_unique == 0, ]
        no_unique_nm = "_no_unique"
      }else{
        no_unique_nm = ""
      }
      validation_dat = build_data_violin_plot(validation_dat)
      validation_dat$Model = input
      benchmark_df_all = rbind(benchmark_df_all, validation_dat)
      validation_dat = convert_numeric_to_class(validation_dat, quantiles)
      
      benchmark_df_models = rbind(benchmark_df_models, validation_dat)
    }
    plot_change = plot_prob_change_group(benchmark_df_models)
    ggsave(glue("{PATH_RES_CHANGE}/{protease}/OpenMS_vs_MM{no_unique_nm}.png"), plot = plot_change)
  }
  benchmark_df_all = convert_numeric_to_class(benchmark_df_all, quantiles)
  plot_change = plot_prob_change_group(benchmark_df_all)
  ggsave(glue("{PATH_RES_CHANGE}/OpenMS_vs_MM{no_unique_nm}.png"), plot = plot_change)
}

main(proteases = list.dirs(PATH_RES_CHANGE, recursive = FALSE, full.names = FALSE))
main(proteases = list.dirs(PATH_RES_CHANGE, recursive = FALSE, full.names = FALSE),
     no_unique = TRUE)
