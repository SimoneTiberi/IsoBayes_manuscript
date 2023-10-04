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
source(glue("{PATH_WD}/utils_function/add_mrna.R"))

PATH_DATA = glue("{PATH_WD}/Data/{DATA}") 
PATH_RES_MODEL = glue("{PATH_WD}/Model_results/{DATA}")
PATH_RES_CHANGE = glue("{PATH_WD}/Change_protein_mRNA_isoform/{DATA}")
log_output(glue("change_prot_mRNA_iso_{DATA}"))

main = function(proteases, no_unique = FALSE){
  inputs = c("OpenMS", "MM_psm", "MM_intensities")
  quantiles = c(0, 0.2, 0.4, 0.6, 0.8, 1)
  
  for (input in inputs) {
    for (pep in c("", "_PEP")){
      for (mrna in c("", "_mRNA")) {
        ######################################################
        # Changes in protein and mRNA isoform relative abundances (MM)
        ######################################################
        benchmark_df_all = list()
        for (protease in proteases) {
          # load validation data with res
          load(glue("{PATH_RES_MODEL}/{input}{mrna}{pep}/{protease}/Merged_validation_res_{input}{mrna}{pep}"))
          if(mrna == ""){
            validation_dat = add_mrna(input, protease, pep, validation_dat)
          }
          colnames(validation_dat)[grep("tpm", colnames(validation_dat))] = "tpm_validation"
          if(no_unique){
            validation_dat = validation_dat[validation_dat$Y_unique == 0, ]
            no_unique_nm = "_no_unique"
          }else{
            no_unique_nm = ""
          }
          validation_dat = build_data_violin_plot(validation_dat, 1.5e-06)
          benchmark_df_all = rbind(benchmark_df_all, validation_dat)
          validation_dat = convert_numeric_to_class(validation_dat, quantiles)
          validation_dat_extreme = convert_numeric_to_class(validation_dat, quantiles = c(0, 0.05, 0.95, 1))
          validation_dat_extreme = validation_dat_extreme[validation_dat_extreme$class_Prob_prot_inc != "(0.05 ; 0.95]", ]
          
          plot_change = plot_prob_change(validation_dat)
          ggsave(glue("{PATH_RES_CHANGE}/{protease}/main_result_{input}{mrna}{pep}{no_unique_nm}.png"), plot = plot_change)
          
          plot_change = plot_prob_change(validation_dat_extreme)
          ggsave(glue("{PATH_RES_CHANGE}/{protease}/main_result_extreme_{input}{mrna}{pep}{no_unique_nm}.png"), plot = plot_change)
          
          ######################################################
          # Log2FC correlation
          ######################################################
          sel = validation_dat$TPM != 0 & validation_dat$tpm_validation != 0 & validation_dat$Y_validation != 0 & validation_dat$Pi != 0
          plot_scatter = scatterplot(validation_dat[sel, c("Log2_FC", "Log2_FC_validation")]) + 
            labs(x = "Log2_FC", y = "Log2_FC_validation")
          ggsave(glue("{PATH_RES_CHANGE}/{protease}/scatterplot_log2FC_{input}{mrna}{pep}{no_unique_nm}.png"), plot = plot_scatter)
          
          sel = validation_dat$tpm_validation != 0 & validation_dat$Y_validation != 0
          plot_scatter = scatterplot(validation_dat[sel, c("Prob_prot_inc", "Log2_FC_validation")]) + 
            labs(x = "Prob_prot_inc", y = "Log2_FC_validation")
          ggsave(glue("{PATH_RES_CHANGE}/{protease}/scatterplot_probInc_log2FC_{input}{mrna}{pep}{no_unique_nm}.png"), plot = plot_scatter)
        }
        benchmark_df_all = convert_numeric_to_class(benchmark_df_all, quantiles)
        benchmark_df_all_extreme = convert_numeric_to_class(benchmark_df_all, quantiles = c(0, 0.05, 0.95, 1))
        benchmark_df_all_extreme = benchmark_df_all_extreme[benchmark_df_all_extreme$class_Prob_prot_inc != "(0.05 ; 0.95]", ]
        
        plot_change = plot_prob_change(benchmark_df_all)
        ggsave(glue("{PATH_RES_CHANGE}/main_result_{input}{mrna}{pep}{no_unique_nm}.png"), plot = plot_change)
        
        plot_change = plot_prob_change(benchmark_df_all_extreme)
        ggsave(glue("{PATH_RES_CHANGE}/main_result_extreme_{input}{mrna}{pep}{no_unique_nm}.png"), plot = plot_change)
        
        ######################################################
        # Log2FC correlation
        ######################################################
        sel = benchmark_df_all$TPM != 0 & benchmark_df_all$tpm_validation != 0 & benchmark_df_all$Y_validation != 0 & benchmark_df_all$Pi != 0
        plot_scatter = scatterplot(benchmark_df_all[sel, c("Log2_FC", "Log2_FC_validation")]) + 
          labs(x = "Log2_FC", y = "Log2_FC_validation")
        ggsave(glue("{PATH_RES_CHANGE}/scatterplot_log2FC_{input}{mrna}{pep}{no_unique_nm}.png"), plot = plot_scatter)
        
        sel = benchmark_df_all$tpm_validation != 0 & benchmark_df_all$Y_validation != 0
        plot_scatter = scatterplot(benchmark_df_all[sel, c("Prob_prot_inc", "Log2_FC_validation")]) + 
          labs(x = "Prob_prot_inc", y = "Log2_FC_validation")
        ggsave(glue("{PATH_RES_CHANGE}/scatterplot_probInc_log2FC_{input}{mrna}{pep}{no_unique_nm}.png"), plot = plot_scatter)
      }
    }
  }
  
  ######################################################
  # Changes in protein and mRNA isoform relative abundances - PEP-noPEP
  ######################################################
  for (input in inputs) {
    for (mrna in c("", "_mRNA")) {
      
      benchmark_df_all = list()
      
      for (protease in proteases) {
        benchmark_df_models = list()
        for (model in c("_PEP", "")) {
          # load validation data with res
          load(glue("{PATH_RES_MODEL}/{input}{mrna}{model}/{protease}/Merged_validation_res_{input}{mrna}{model}"))
          if(mrna == ""){
            validation_dat = add_mrna(input, protease, validation_dat)
          }
          colnames(validation_dat)[grep("tpm", colnames(validation_dat))] = "tpm_validation"
          
          if(no_unique){
            validation_dat = validation_dat[validation_dat$Y_unique == 0, ]
            no_unique_nm = "_no_unique"
          }else{
            no_unique_nm = ""
          }
          
          validation_dat = build_data_violin_plot(validation_dat, 1.5e-06)
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
        plot_change = plot_prob_change_group(benchmark_df_models)
        ggsave(glue("{PATH_RES_CHANGE}/{protease}/PEP_no_PEP_{input}{mrna}{no_unique_nm}.png"), plot = plot_change)
      }
      benchmark_df_all = convert_numeric_to_class(benchmark_df_all, quantiles)
      plot_change = plot_prob_change_group(benchmark_df_all)
      ggsave(glue("{PATH_RES_CHANGE}/PEP_no_PEP_{input}{mrna}{no_unique_nm}.png"), plot = plot_change)
    }
  }
  
  ######################################################
  # Changes in protein and mRNA isoform relative abundances - MM vs OpenMS
  ######################################################
  for (mrna in c("", "_mRNA")) {
    benchmark_df_all = list()
    for (protease in proteases) {
      benchmark_df_models = list()
      for (input in c("MM_psm", "OpenMS")) {
        # load validation data with res
        load(glue("{PATH_RES_MODEL}/{input}{mrna}/{protease}/Merged_validation_res_{input}{mrna}"))
        if(mrna == ""){
          validation_dat = add_mrna(input, protease, validation_dat)
        }
        colnames(validation_dat)[grep("tpm", colnames(validation_dat))] = "tpm_validation"
        
        if(no_unique){
          validation_dat = validation_dat[validation_dat$Y_unique == 0, ]
          no_unique_nm = "_no_unique"
        }else{
          no_unique_nm = ""
        }
        validation_dat = build_data_violin_plot(validation_dat, 1.5e-06)
        validation_dat$Model = input
        benchmark_df_all = rbind(benchmark_df_all, validation_dat)
        validation_dat = convert_numeric_to_class(validation_dat, quantiles)
        
        benchmark_df_models = rbind(benchmark_df_models, validation_dat)
      }
      plot_change = plot_prob_change_group(benchmark_df_models)
      ggsave(glue("{PATH_RES_CHANGE}/{protease}/OpenMS_vs_MM{mrna}{no_unique_nm}.png"), plot = plot_change)
    }
    benchmark_df_all = convert_numeric_to_class(benchmark_df_all, quantiles)
    plot_change = plot_prob_change_group(benchmark_df_all)
    ggsave(glue("{PATH_RES_CHANGE}/OpenMS_vs_MM{mrna}{no_unique_nm}.png"), plot = plot_change)
  }
  
  ######################################################
  # Changes in protein and mRNA isoform relative abundances - MM_psm vs MM_intensities
  ######################################################
  for (mrna in c("", "_mRNA")) {
    benchmark_df_all = list()
    for (protease in proteases) {
      benchmark_df_models = list()
      for (input in c("MM_psm", "MM_intensities")) {
        # load validation data with res
        load(glue("{PATH_RES_MODEL}/{input}{mrna}/{protease}/Merged_validation_res_{input}{mrna}"))
        if(mrna == ""){
          validation_dat = add_mrna(input, protease, validation_dat)
        }
        colnames(validation_dat)[grep("tpm", colnames(validation_dat))] = "tpm_validation"
        
        if(no_unique){
          validation_dat = validation_dat[validation_dat$Y_unique == 0, ]
          no_unique_nm = "_no_unique"
        }else{
          no_unique_nm = ""
        }
        validation_dat = build_data_violin_plot(validation_dat, 1.5e-06)
        validation_dat$Model = input
        benchmark_df_all = rbind(benchmark_df_all, validation_dat)
        validation_dat = convert_numeric_to_class(validation_dat, quantiles)
        
        benchmark_df_models = rbind(benchmark_df_models, validation_dat)
      }
      plot_change = plot_prob_change_group(benchmark_df_models)
      ggsave(glue("{PATH_RES_CHANGE}/{protease}/MM_psm_vs_MM_intensities{mrna}{no_unique_nm}.png"), plot = plot_change)
    }
    benchmark_df_all = convert_numeric_to_class(benchmark_df_all, quantiles)
    plot_change = plot_prob_change_group(benchmark_df_all)
    ggsave(glue("{PATH_RES_CHANGE}/MM_psm_vs_MM_intensities{mrna}{no_unique_nm}.png"), plot = plot_change)
  }
  
  ######################################################
  # Changes in protein and mRNA isoform relative abundances - mRNA-no_mRNA
  ######################################################
  for (input in inputs) {
    
    benchmark_df_all = list()
    
    for (protease in proteases) {
      benchmark_df_models = list()
      for (model in c("_mRNA", "")) {
        # load validation data with res
        load(glue("{PATH_RES_MODEL}/{input}{model}/{protease}/Merged_validation_res_{input}{model}"))
        colnames(validation_dat)[grep("tpm", colnames(validation_dat))] = "tpm_validation"
        
        if(model == ""){
          load(glue("{PATH_RES_MODEL}/{input}{model}/{protease}/{input}{model}_MCMC.RData"))
          mrna_data = data.table::fread(glue("{PATH_DATA}/mrna_isoform.tsv"))
          res$isoform_results = merge(res$isoform_results, mrna_data[, c("isoname", "tpm")], by.x = "Isoform",
                                      by.y = "isoname", all.x = TRUE)
          res$isoform_results = res$isoform_results[!duplicated(res$isoform_results$Isoform), ]
          P_TPM = res$isoform_results$tpm/sum(res$isoform_results$tpm)
          res$isoform_results$Log2_FC = log2(res$isoform_results$Pi/P_TPM)
          
          res$isoform_results$Prob_prot_inc = vapply(seq_len(nrow(res$isoform_results)), function(i){
            mean(res$chain_Y[, i] > P_TPM[i])}, FUN.VALUE = numeric(1) )
          
          validation_dat = merge(validation_dat, res$isoform_results[, c("Isoform", "Prob_prot_inc")], by = "Isoform")
        }
        
        if(no_unique){
          validation_dat = validation_dat[validation_dat$Y_unique == 0, ]
          no_unique_nm = "_no_unique"
        }else{
          no_unique_nm = ""
        }
        
        validation_dat = build_data_violin_plot(validation_dat, 1.5e-06)
        validation_dat$Gene = NULL
        
        if(model == "_mRNA"){
          label_model = "mRNA"
        }else{
          label_model = "No mRNA"
        }
        validation_dat$Model = label_model
        
        benchmark_df_all = rbind(benchmark_df_all, validation_dat[, c("Log2_FC_validation", "Model", "Prob_prot_inc")])
        
        validation_dat = convert_numeric_to_class(validation_dat, quantiles)
        benchmark_df_models = rbind(benchmark_df_models, validation_dat[, c("Log2_FC_validation", "Model", "Prob_prot_inc")])
      }
      benchmark_df_models = convert_numeric_to_class(benchmark_df_models, quantiles)
      plot_change = plot_prob_change_group(benchmark_df_models)
      ggsave(glue("{PATH_RES_CHANGE}/{protease}/mRNA_no_mRNA_{input}{no_unique_nm}.png"), plot = plot_change)
    }
    benchmark_df_all = convert_numeric_to_class(benchmark_df_all, quantiles)
    plot_change = plot_prob_change_group(benchmark_df_all)
    ggsave(glue("{PATH_RES_CHANGE}/mRNA_no_mRNA_{input}{no_unique_nm}.png"), plot = plot_change)
  }
}

main(proteases = list.dirs(PATH_RES_CHANGE, recursive = FALSE, full.names = FALSE))
main(proteases = list.dirs(PATH_RES_CHANGE, recursive = FALSE, full.names = FALSE),
     no_unique = TRUE)
