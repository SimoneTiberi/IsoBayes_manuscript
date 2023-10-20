PATH_WD = commandArgs(trailingOnly = TRUE)[1]
DATA = commandArgs(trailingOnly = TRUE)[2]

if(DATA == "wtc11"){
  DATA_name = "WTC-11"
}else{
  DATA_name = "jurkat"
}

# Set path, global variable and libraries
###########################################################################################
library(ggplot2)
library(glue)
source(glue("{PATH_WD}/utils_function/get_roc.R"))
source(glue("{PATH_WD}/utils_function/prior_plot.R"))
source(glue("{PATH_WD}/utils_function/utils_benchmarking.R"))
source(glue("{PATH_WD}/utils_function/concat_models.R"))
source(glue("{PATH_WD}/utils_function/log_output.R"))
source(glue("{PATH_WD}/utils_function/scatterplot.R"))
source(glue("{PATH_WD}/utils_function/utils_change_iso_mrna.R"))
source(glue("{PATH_WD}/utils_function/plot_prob_change.R"))
source(glue("{PATH_WD}/utils_function/plot_prob_change_group.R"))
source(glue("{PATH_WD}/utils_function/compare_inputs.R"))
source(glue("{PATH_WD}/utils_function/compare_pep_fdr.R"))

PATH_TO_DATA = glue("{PATH_WD}/Data/{DATA}") 
PATH_RES = glue("{PATH_WD}/Model_results/{DATA}")
PATH_RES_roc = glue("{PATH_WD}/Robustness/{DATA}")
load(glue("{PATH_WD}/utils_function/PALETTE_MODELS"))

log_output(glue("robustness_results_{DATA}"))

main = function(models, proteases){
  ###################################################################
  # PEP/no PEP and mRNA vs prot (OpenMS and MM)  
  ###################################################################
  selected_models = c("IsoBayes", "IsoBayes_PEP", "IsoBayes_mRNA", "IsoBayes_mRNA_PEP")
  
  # for (input in c("OpenMS", "MM_psm", "MM_intensities")){
  #   benchmark_df_all = list()
  #   for(protease in proteases){
  #     benchmark_df = list()
  #     for (model in selected_models) {
  #       attribute_model = glue("{models[[model]][2]}{models[[model]][1]}")
  #       # load res and validation merged together
  #       load(glue("{PATH_RES}/{input}{attribute_model}/{protease}/Merged_validation_res_{input}{attribute_model}"))
  #       
  #       validation_dat = validation_dat[, c("Isoform", "Prob_present", "Present")]
  #       colnames(validation_dat)[2] = model
  #       colnames(validation_dat)[3] = glue("{colnames(validation_dat)[3]}_{model}")
  #       
  #       validation_dat = validation_dat[!duplicated(validation_dat$Isoform), ]
  #       benchmark_df = append(benchmark_df, list(validation_dat))
  #     }
  #     benchmark_df = concat_models(benchmark_df, union = TRUE)
  #   
  #     for (nm in selected_models) {
  #       benchmark_df[, paste0("Present_", nm)]
  #     }
  #     
  #     plot_tab = get_roc(benchmark_df, selected_models)
  #     ggsave(glue("{PATH_RES_roc}/{protease}/ROC_PEP_vs_no_PEP_mRNA_vs_no_mRNA_{input}.png"), plot = plot_tab$gplot)
  #     write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_roc}/{protease}/SumTab_PEP_vs_no_PEP_mRNA_vs_no_mRNA_{input}.csv"), row.names = FALSE)
  #     benchmark_df_all = rbind(benchmark_df_all, benchmark_df)
  #   }
  #   plot_tab = get_roc(benchmark_df_all, selected_models)
  #   ggsave(glue("{PATH_RES_roc}/ROC_PEP_vs_no_PEP_mRNA_vs_no_mRNA_{input}.png"), plot = plot_tab$gplot)
  #   write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_roc}/SumTab_PEP_vs_no_PEP_mRNA_vs_no_mRNA_{input}.csv"), row.names = FALSE)
  # }
  # OpenMS vs MM_psm
  #compare_inputs(inputs = c("OpenMS", "MM_psm"), selected_models[c(1, 3)], proteases, models)
  
  # MM_psm vs MM_intensities vs OpenMS
  compare_inputs(inputs = c("MM_psm", "MM_intensities", "OpenMS"), selected_models[c(1, 3)], proteases, models)
  
  # PEP vs FDR - MM_psm
  compare_pep_fdr("MM_psm", selected_models[c(1, 3)], proteases, models)
  
  # PEP vs FDR - MM_intensities
  compare_pep_fdr("MM_intensities", selected_models[c(1, 3)], proteases, models)
}

main(models = list(IsoBayes = c("", ""),
                   IsoBayes_PEP = c("_PEP", ""),
                   IsoBayes_mRNA = c("", "_mRNA"),
                   IsoBayes_mRNA_PEP = c("_PEP", "_mRNA")),
     proteases = list.dirs(glue("{PATH_WD}/Benchmark_results/{DATA}"), recursive = FALSE, full.names = FALSE)
     )
