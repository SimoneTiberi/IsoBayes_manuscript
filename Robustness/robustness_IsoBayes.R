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
dir.create(PATH_RES_roc)
load(glue("{PATH_WD}/utils_function/PALETTE_MODELS"))

log_output(glue("robustness_results_{DATA}"))

main = function(models, proteases){
  for (prot in proteases) {
    dir.create(paste0(PATH_RES_roc, "/", prot))
  }
  ###################################################################
  # PEP/no PEP and mRNA vs prot (OpenMS and MM)  
  ###################################################################
  selected_models = c("IsoBayes", "IsoBayes_PEP", "IsoBayes_mRNA", "IsoBayes_mRNA_PEP")
  
  # MM_psm vs MM_intensities vs OpenMS
  compare_inputs(inputs = c("MM_psm", "MM_intensities", "OpenMS"),
                 selected_models[c(1, 3)], proteases, models)
  
  # PEP vs FDR - MM_psm
  compare_pep_fdr("MM_psm", selected_models[c(1, 3)], proteases, models)
  
  # PEP vs FDR - MM_intensities
  compare_pep_fdr("MM_intensities", selected_models[c(1, 3)], proteases, models)
}

main(models = list(IsoBayes = c("", ""),
                   IsoBayes_PEP = c("_PEP", ""),
                   IsoBayes_mRNA = c("", "_mRNA"),
                   IsoBayes_mRNA_PEP = c("_PEP", "_mRNA")),
     proteases = list.dirs(glue("{PATH_WD}/Benchmark_results/{DATA}"),
                           recursive = FALSE, full.names = FALSE))
