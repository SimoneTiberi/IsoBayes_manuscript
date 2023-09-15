PATH_WD = commandArgs(trailingOnly = TRUE)[1]
DATA = commandArgs(trailingOnly = TRUE)[2]

# Set path, global variable and libraries
###########################################################################################
library(ggplot2)
library(glue)
source(glue("{PATH_WD}/utils_function/get_roc.R"))
source(glue("{PATH_WD}/utils_function/prior_plot.R"))
source(glue("{PATH_WD}/utils_function/utils_benchmarking.R"))
source(glue("{PATH_WD}/utils_function/concat_models.R"))
source(glue("{PATH_WD}/utils_function/log_output.R"))

PATH_TO_DATA = glue("{PATH_WD}/Data/{DATA}") 
PATH_RES = glue("{PATH_WD}/Model_results/{DATA}")
PATH_RES_roc = glue("{PATH_WD}/Robustness/{DATA}/")
load(glue("{PATH_WD}/utils_function/PALETTE_MODELS"))

log_output(glue("robustness_results_{DATA}"))

main = function(models, proteases){
  ###################################################################
  # PEP/no PEP and mRNA vs prot (OpenMS and MM)  
  ###################################################################
  selected_models = c("IsoBayes", "IsoBayes_fast", "IsoBayes_mRNA", "IsoBayes_fast_mRNA")
  
  for (input in c("OpenMS", "MM_psm", "MM_intensities")){
    benchmark_df_all = list()
    for(protease in proteases){
      benchmark_df = list()
      for (model in selected_models) {
        attribute_model = glue("{models[[model]][2]}{models[[model]][1]}")
        # load res and validation merged together
        load(glue("{PATH_RES}/{input}{attribute_model}/{protease}/Merged_validation_res_{input}{attribute_model}"))
        validation_dat = validation_dat[, c("Isoform", "Prob_present", "Present")]
        colnames(validation_dat)[2] = model
        colnames(validation_dat)[3] = glue("{colnames(validation_dat)[3]}_{model}")
        validation_dat = validation_dat[!duplicated(validation_dat$Isoform), ]
        benchmark_df = append(benchmark_df, list(validation_dat))
      }
      benchmark_df = concat_models(benchmark_df)
      
      plot_tab = get_roc(benchmark_df, selected_models)
      ggsave(glue("{PATH_RES_roc}/{protease}/ROC_PEP_vs_no_PEP_mRNA_vs_no_mRNA_{input}.png"), plot = plot_tab$gplot)
      write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_roc}/{protease}/SumTab_PEP_vs_no_PEP_mRNA_vs_no_mRNA_{input}.csv"), row.names = FALSE)
      benchmark_df_all = rbind(benchmark_df_all, benchmark_df)
    }
    plot_tab = get_roc(benchmark_df_all, selected_models)
    ggsave(glue("{PATH_RES_roc}/ROC_PEP_vs_no_PEP_mRNA_vs_no_mRNA_{input}.png"), plot = plot_tab$gplot)
    write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_roc}/SumTab_PEP_vs_no_PEP_mRNA_vs_no_mRNA_{input}.csv"), row.names = FALSE)
  }
  ###################################################################
  # OpenMS vs MM
  ###################################################################
  model = selected_models[1]
  benchmark_df_all = list()
  inputs = c("OpenMS", "MM_psm")
  
  for(protease in proteases){
    benchmark_df = list()
    for (input in inputs) {
      attribute_model = glue("{models[[model]][2]}{models[[model]][1]}") 
      # load res and validation merged together
      load(glue("{PATH_RES}/{input}{attribute_model}/{protease}/Merged_validation_res_{input}{attribute_model}"))
      validation_dat = validation_dat[, c("Isoform", "Prob_present", "Present")]
      colnames(validation_dat)[2] = glue("{model}_{input}")
      colnames(validation_dat)[3] = glue("{colnames(validation_dat)[3]}_{model}_{input}")
      validation_dat = validation_dat[!duplicated(validation_dat$Isoform), ]
      benchmark_df = append(benchmark_df, list(validation_dat))
    }
    benchmark_df = concat_models(benchmark_df)
    
    plot_tab = get_roc(benchmark_df, paste0(model, "_", inputs))
    ggsave(glue("{PATH_RES_roc}/{protease}/ROC_MM_vs_OpenMS.png"), plot = plot_tab$gplot)
    write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_roc}/{protease}/SumTab_MM_vs_OpenMS.csv"), row.names = FALSE)
    benchmark_df_all = rbind(benchmark_df_all, benchmark_df)
  }
  plot_tab = get_roc(benchmark_df_all, paste0(selected_models, "_", inputs))
  ggsave(glue("{PATH_RES_roc}/ROC_MM_vs_OpenMS.png"), plot = plot_tab$gplot)
  write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_roc}/SumTab_MM_vs_OpenMS.csv"), row.names = FALSE)
  
  ###################################################################
  # MM psm vs intensities
  ###################################################################
  model = selected_models[1]
  benchmark_df_all = list()
  inputs = c("MM_psm", "MM_intensities")
  
  for(protease in proteases){
    benchmark_df = list()
    for (input in inputs) {
      attribute_model = glue("{models[[model]][2]}{models[[model]][1]}") 
      # load res and validation merged together
      load(glue("{PATH_RES}/{input}{attribute_model}/{protease}/Merged_validation_res_{input}{attribute_model}"))
      validation_dat = validation_dat[, c("Isoform", "Prob_present", "Present")]
      colnames(validation_dat)[2] = glue("{model}_{input}")
      colnames(validation_dat)[3] = glue("{colnames(validation_dat)[3]}_{model}_{input}")
      validation_dat = validation_dat[!duplicated(validation_dat$Isoform), ]
      benchmark_df = append(benchmark_df, list(validation_dat))
    }
    benchmark_df = concat_models(benchmark_df)
    
    plot_tab = get_roc(benchmark_df, paste0(model, "_", inputs))
    ggsave(glue("{PATH_RES_roc}/{protease}/ROC_psm_vs_intensities.png"), plot = plot_tab$gplot)
    write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_roc}/{protease}/SumTab_psm_vs_intensities.csv"), row.names = FALSE)
    benchmark_df_all = rbind(benchmark_df_all, benchmark_df)
  }
  plot_tab = get_roc(benchmark_df_all, paste0(model, "_", inputs))
  ggsave(glue("{PATH_RES_roc}/ROC_psm_vs_intensities.png"), plot = plot_tab$gplot)
  write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_roc}/SumTab_psm_vs_intensities.csv"), row.names = FALSE)
  
  ###################################################################
  # OpenMS - prior robustness
  ###################################################################
  load(glue("{PATH_WD}/Model_results/prior_grid"))
  selected_models = rep("IsoBayes_mRNA", length(prior_grid))
  model = selected_models[1]
  benchmark_df_all = list()
  
  for(protease in proteases){
    benchmark_df = list()
    for (prior in prior_grid) {
      attribute_model = glue("{models[[model]][2]}{models[[model]][1]}_prior_{prior}")
      # load res and validation merged together
      load(glue("{PATH_RES}/OpenMS{attribute_model}/{protease}/Merged_validation_res_OpenMS{attribute_model}"))
      validation_dat = validation_dat[, c("Isoform", "Prob_present", "Present")]
      colnames(validation_dat)[2] = as.character(prior)
      colnames(validation_dat)[3] = glue("{colnames(validation_dat)[3]}_{prior}")
      validation_dat = validation_dat[!duplicated(validation_dat$Isoform), ]
      benchmark_df = append(benchmark_df, list(validation_dat))
    }
    benchmark_df = concat_models(benchmark_df)
    
    plot_tab = get_roc(benchmark_df, as.character(prior_grid))
    data_to_plot = plot_tab$sum_stat
    colnames(data_to_plot)[1] = "grid"
    
    pp = prior_plot(data_to_plot)
    
    ggsave(glue("{PATH_RES_roc}/{protease}/prior_plot.png"), plot = pp)
    write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_roc}/{protease}/SumTab_prior.csv"), row.names = FALSE)
    benchmark_df_all = rbind(benchmark_df_all, benchmark_df)
  }
  plot_tab = get_roc(benchmark_df_all, as.character(prior_grid))
  data_to_plot = plot_tab$sum_stat
  data_to_plot = plot_tab$sum_stat
  colnames(data_to_plot)[1] = "grid"
  
  pp = prior_plot(data_to_plot)
  
  ggsave(glue("{PATH_RES_roc}/prior_plot.png"), plot = pp)
  write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_roc}/SumTab_prior.csv"))
}

main(models = list(IsoBayes = c("_PEP", ""),
                   IsoBayes_fast = c("", ""),
                   IsoBayes_mRNA = c("_PEP", "_mRNA"),
                   IsoBayes_fast_mRNA = c("", "_mRNA")),
     proteases = list.dirs(glue("{PATH_WD}/Benchmark_results/{DATA}"), recursive = FALSE, full.names = FALSE)
     )
