PATH_WD = commandArgs(trailingOnly = TRUE)[1]
DATA = commandArgs(trailingOnly = TRUE)[2]

# Set path, global variable and libraries
###########################################################################################
library(glue)
library(IsoBayes)
library(ggplot2)

source(glue("{PATH_WD}/utils_function/validate_all_protease.R"))
source(glue("{PATH_WD}/utils_function/save_run_inferences.R"))
source(glue("{PATH_WD}/utils_function/plot_roc_model.R"))
source(glue("{PATH_WD}/utils_function/log_output.R"))

PATH_TO_DATA = glue("{PATH_WD}/Data/{DATA}")
PATH_TO_RES = glue("{PATH_WD}/Model_results/{DATA}")
load(glue("{PATH_WD}/utils_function/PALETTE_MODELS"))
set.seed(0)

log_output(glue("run_models_{DATA}"))
###########################################################################################

main = function(proteases, run_model = TRUE){
  # OpenMS data
  log_output(glue("run_models_{DATA}_OpenMS"))
  ###############################
  name = glue("OpenMS_mRNA_PEP")
  message(glue("---------- {name} ----------"))
  name_models = c("IsoBayes_mRNA", "TPM", "Baseline")
  ###############################
  for (protease in proteases) {
    message(protease)
    path_to_res_mod = glue("{PATH_TO_RES}/{name}/{protease}")
    if (run_model) {
      if(!dir.exists(path_to_res_mod)){dir.create(path_to_res_mod, recursive = TRUE)}
      path_to_peptides = paste0(PATH_TO_DATA, "/Only", protease, "/merge_index_percolator_pep_switched_0.1.idXML")
      
      data_loaded = load_data(path_to_peptides_psm = path_to_peptides,
                              path_to_tpm = paste0(PATH_TO_DATA, "/mrna_isoform.tsv"),
                              input_type = "openMS",
                              abundance_type = "psm",
                              PEP = TRUE,
                              FDR_thd = 0.1
      )
      map_iso_gene = glue("{PATH_WD}/Data/{DATA}/map_iso_gene_{DATA}.csv")
      save_run_inferences(data_loaded, path_to_res_mod, name,
                          map_iso_gene = map_iso_gene, save_chain = TRUE)
    }
    pp = plot_roc_model(path_to_res_mod, name, protease, name_models)
    ggsave(glue("{path_to_res_mod}/{name}.png"))
  }
  pp = validate_all_protease(proteases, name, name_models)
  ggsave(glue("{PATH_TO_RES}/{name}/{name}.png"))
  
  ###############################
  name = glue("OpenMS_mRNA")
  message(glue("---------- {name} ----------"))
  name_models = c("IsoBayes_fast_mRNA", "TPM", "Baseline")
  ###############################
  
  for (protease in proteases) {
    message(protease)
    path_to_res_mod = glue("{PATH_TO_RES}/{name}/{protease}")
    if (run_model) {
      if(!dir.exists(path_to_res_mod)){dir.create(path_to_res_mod, recursive = TRUE)}
      data_loaded = load_data(path_to_peptides_psm = paste0(PATH_TO_DATA, "/Only", protease, "/merge_index_percolator_pep_switched_0.01.idXML"),
                              path_to_tpm = paste0(PATH_TO_DATA, "/mrna_isoform.tsv"),
                              input_type = "openMS",
                              abundance_type = "psm",
                              PEP = FALSE,
                              FDR_thd = 0.01
      )
      map_iso_gene = glue("{PATH_WD}/Data/{DATA}/map_iso_gene_{DATA}.csv")
      save_run_inferences(data_loaded, path_to_res_mod, name, map_iso_gene = map_iso_gene, save_chain = TRUE)
    }
    pp = plot_roc_model(path_to_res_mod, name, protease, name_models)
    ggsave(glue("{path_to_res_mod}/{name}.png"))
  }
  pp = validate_all_protease(proteases, name, name_models)
  ggsave(glue("{PATH_TO_RES}/{name}/{name}.png"))
  
  ###############################
  name = glue("OpenMS_PEP")
  message(glue("---------- {name} ----------"))
  name_models = c("IsoBayes", "Baseline")
  ###############################
  
  for (protease in proteases) {
    message(protease)
    path_to_res_mod = glue("{PATH_TO_RES}/{name}/{protease}")
    if (run_model) {
      if(!dir.exists(path_to_res_mod)){dir.create(path_to_res_mod, recursive = TRUE)}
      
      data_loaded = load_data(path_to_peptides_psm = paste0(PATH_TO_DATA, "/Only", protease, "/merge_index_percolator_pep_switched_0.1.idXML"),
                              input_type = "openMS",
                              abundance_type = "psm",
                              PEP = TRUE,
                              FDR_thd = 0.1
      )
      save_run_inferences(data_loaded, path_to_res_mod, name, save_chain = TRUE)
    }
    pp = plot_roc_model(path_to_res_mod, name, protease, name_models)
    ggsave(glue("{path_to_res_mod}/{name}.png"))
  }
  pp = validate_all_protease(proteases, name, name_models)
  ggsave(glue("{PATH_TO_RES}/{name}/{name}.png"))
  
  ###############################
  name = glue("OpenMS")
  message(glue("---------- {name} ----------"))
  name_models = c("IsoBayes_fast", "Baseline")
  ###############################
  
  for (protease in proteases) {
    message(protease)
    path_to_res_mod = glue("{PATH_TO_RES}/{name}/{protease}")
    if (run_model) {
      if(!dir.exists(path_to_res_mod)){dir.create(path_to_res_mod, recursive = TRUE)}
      
      data_loaded = load_data(path_to_peptides_psm = paste0(PATH_TO_DATA, "/Only", protease, "/merge_index_percolator_pep_switched_0.01.idXML"),
                              input_type = "openMS",
                              abundance_type = "psm",
                              PEP = FALSE,
                              FDR_thd = 0.01
      )
      map_iso_gene = glue("{PATH_WD}/Data/{DATA}/map_iso_gene_{DATA}.csv")
      save_run_inferences(data_loaded, path_to_res_mod, name, map_iso_gene = map_iso_gene, save_chain = T)
    }
    pp = plot_roc_model(path_to_res_mod, name, protease, name_models)
    ggsave(glue("{path_to_res_mod}/{name}.png"))
  }
  pp = validate_all_protease(proteases, name, name_models)
  ggsave(glue("{PATH_TO_RES}/{name}/{name}.png"))

  ############################################################################################################
  # MetaMorpheus data (MM)
  ############################################################################################################
  for (abundance_type in c("psm", "intensities")) {
  
    log_output(glue("run_models_{DATA}_MM_{abundance_type}"))
    
    ###############################
    name = glue("MM_{abundance_type}_mRNA_PEP")
    message(glue("---------- {name} ----------"))
    name_models = c("IsoBayes_mRNA", "TPM", "Baseline")
    ###############################
    
    for (protease in proteases) {
      message(protease)
      path_to_res_mod = glue("{PATH_TO_RES}/{name}/{protease}")
      if (run_model) {
        if(!dir.exists(path_to_res_mod)){dir.create(path_to_res_mod, recursive = TRUE)}
        
        data_loaded = load_data(path_to_peptides_psm = paste0(PATH_TO_DATA, "/Only", protease, "/AllPeptides.psmtsv"),
                                path_to_peptides_intensities = paste0(PATH_TO_DATA, "/Only", protease, "/AllQuantifiedPeptides.tsv"),
                                path_to_tpm = paste0(PATH_TO_DATA, "/mrna_isoform.tsv"),
                                input_type = "metamorpheus",
                                abundance_type = abundance_type,
                                PEP = TRUE,
                                FDR_thd = 0.1
        )
        map_iso_gene = glue("{PATH_WD}/Data/{DATA}/map_iso_gene_{DATA}.csv")
        save_run_inferences(data_loaded, path_to_res_mod, name, map_iso_gene = map_iso_gene, save_chain = TRUE)
      }
      pp = plot_roc_model(path_to_res_mod, name, protease, name_models)
      ggsave(glue("{path_to_res_mod}/{name}.png"))
    }
    pp = validate_all_protease(proteases, name, name_models)
    ggsave(glue("{PATH_TO_RES}/{name}/{name}.png"))
    
    ###############################
    name = glue("MM_{abundance_type}_mRNA")
    message(glue("---------- {name} ----------"))
    name_models = c("IsoBayes_fast_mRNA", "TPM", "Baseline")
    ###############################
    
    for (protease in proteases) {
      message(protease)
      path_to_res_mod = glue("{PATH_TO_RES}/{name}/{protease}")
      if (run_model) {
        if(!dir.exists(path_to_res_mod)){dir.create(path_to_res_mod, recursive = TRUE)}
        data_loaded = load_data(path_to_peptides_psm = paste0(PATH_TO_DATA, "/Only", protease, "/AllPeptides.psmtsv"),
                                path_to_peptides_intensities = paste0(PATH_TO_DATA, "/Only", protease, "/AllQuantifiedPeptides.tsv"),
                                path_to_tpm = paste0(PATH_TO_DATA, "/mrna_isoform.tsv"),
                                input_type = "metamorpheus",
                                abundance_type = abundance_type,
                                PEP = FALSE,
                                FDR_thd = 0.01
        )
        save_run_inferences(data_loaded, path_to_res_mod, name, long_mcmc = FALSE, save_chain = TRUE)
      }
      pp = plot_roc_model(path_to_res_mod, name, protease, name_models)
      ggsave(glue("{path_to_res_mod}/{name}.png"))
    }
    pp = validate_all_protease(proteases, name, name_models)
    ggsave(glue("{PATH_TO_RES}/{name}/{name}.png"))
    
    ###############################
    name = glue("MM_{abundance_type}_PEP")
    message(glue("---------- {name} ----------"))
    name_models = c("IsoBayes", "Baseline")
    ###############################
    
    for (protease in proteases) {
      message(protease)
      path_to_res_mod = glue("{PATH_TO_RES}/{name}/{protease}")
      if (run_model) {
        if(!dir.exists(path_to_res_mod)){dir.create(path_to_res_mod, recursive = TRUE)}
        
        data_loaded = load_data(path_to_peptides_psm = paste0(PATH_TO_DATA, "/Only", protease, "/AllPeptides.psmtsv"),
                                path_to_peptides_intensities = paste0(PATH_TO_DATA, "/Only", protease, "/AllQuantifiedPeptides.tsv"),
                                input_type = "metamorpheus",
                                abundance_type = abundance_type,
                                PEP = TRUE,
                                FDR_thd = 0.1
        )
        save_run_inferences(data_loaded, path_to_res_mod, name, save_chain = TRUE)
      }
      pp = plot_roc_model(path_to_res_mod, name, protease, name_models)
      ggsave(glue("{path_to_res_mod}/{name}.png"))
    }
    pp = validate_all_protease(proteases, name, name_models)
    ggsave(glue("{PATH_TO_RES}/{name}/{name}.png"))
    
    ###############################
    name = glue("MM_{abundance_type}")
    message(glue("---------- {name} ----------"))
    ###############################
    name_models = c("IsoBayes_fast", "Baseline")
    
    for (protease in proteases) {
      message(protease)
      path_to_res_mod = glue("{PATH_TO_RES}/{name}/{protease}")
      if (run_model) {
        if(!dir.exists(path_to_res_mod)){dir.create(path_to_res_mod, recursive = TRUE)}
        
        data_loaded = load_data(path_to_peptides_psm = paste0(PATH_TO_DATA, "/Only", protease, "/AllPeptides.psmtsv"),
                                path_to_peptides_intensities = paste0(PATH_TO_DATA, "/Only", protease, "/AllQuantifiedPeptides.tsv"),
                                input_type = "metamorpheus",
                                abundance_type = abundance_type,
                                PEP = FALSE,
                                FDR_thd = 0.01
        )
        save_run_inferences(data_loaded, path_to_res_mod, name, save_chain = TRUE)
      }
      pp = plot_roc_model(path_to_res_mod, name, protease, name_models)
      ggsave(glue("{path_to_res_mod}/{name}.png"))
    }
    pp = validate_all_protease(proteases, name, name_models)
    ggsave(glue("{PATH_TO_RES}/{name}/{name}.png"))
  }
}

PROTEASES = list.dirs(PATH_TO_DATA, recursive = FALSE, full.names = FALSE)
PROTEASES = PROTEASES[grepl(pattern = "Only", PROTEASES)]
PROTEASES = gsub("Only", "", PROTEASES)

main(PROTEASES, run_model = TRUE)