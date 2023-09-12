PATH_WD = commandArgs(trailingOnly = TRUE)[1]
DATA = commandArgs(trailingOnly = TRUE)[2]

# Set path, global variable and libraries
###########################################################################################
library(glue)
library(IsoBayes)
library(ggplot2)
#source(glue("{PATH_WD}/utils_function/plot_results.R"))
source(glue("{PATH_WD}/utils_function/merge_validation.R"))
source(glue("{PATH_WD}/utils_function/validate_all_protease.R"))
source(glue("{PATH_WD}/utils_function/save_run_inferences_plot_validation.R"))
source(glue("{PATH_WD}/utils_function/get_roc.R"))
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
      data_loaded = load_data(path_to_peptides_psm = paste0(PATH_TO_DATA, "/Only", protease, "/merge_index_percolator_pep_switched.idXML"),
                              path_to_tpm = paste0(PATH_TO_DATA, "/mrna_isoform.tsv"),
                              input_type = "openMS",
                              abundance_type = "psm",
                              PEP = TRUE,
                              FDR_thd = 0.01
      )
      map_iso_gene = glue("{PATH_WD}/Data/{DATA}/map_iso_gene_{DATA}.csv")
      save_run_inferences_plot_validation(data_loaded, path_to_res_mod, name, map_iso_gene = map_iso_gene)
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
      data_loaded = load_data(path_to_peptides_psm = paste0(PATH_TO_DATA, "/Only", protease, "/merge_index_percolator_pep_switched.idXML"),
                              path_to_tpm = paste0(PATH_TO_DATA, "/mrna_isoform.tsv"),
                              input_type = "openMS",
                              abundance_type = "psm",
                              PEP = FALSE,
                              FDR_thd = 0.01
      )
      save_run_inferences_plot_validation(data_loaded, path_to_res_mod, name)
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
      
      data_loaded = load_data(path_to_peptides_psm = paste0(PATH_TO_DATA, "/Only", protease, "/merge_index_percolator_pep_switched.idXML"),
                              input_type = "openMS",
                              abundance_type = "psm",
                              PEP = TRUE,
                              FDR_thd = 0.01
      )
      save_run_inferences_plot_validation(data_loaded, path_to_res_mod, name)
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
      
      data_loaded = load_data(path_to_peptides_psm = paste0(PATH_TO_DATA, "/Only", protease, "/merge_index_percolator_pep_switched.idXML"),
                              input_type = "openMS",
                              abundance_type = "psm",
                              PEP = FALSE,
                              FDR_thd = 0.01
      )
      save_run_inferences_plot_validation(data_loaded, path_to_res_mod, name)
    }
    pp = plot_roc_model(path_to_res_mod, name, protease, name_models)
    ggsave(glue("{path_to_res_mod}/{name}.png"))
  }
  pp = validate_all_protease(proteases, name, name_models)
  ggsave(glue("{PATH_TO_RES}/{name}/{name}.png"))
  
  ############################################################################################################
  # MetaMorpheus data (MM)
  ############################################################################################################
  for (abundance_type in c("intensities", "psm")) {
    
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
                                FDR_thd = 0.01
                                )
        map_iso_gene = glue("{PATH_WD}/Data/{DATA}/map_iso_gene_{DATA}.csv")
        save_run_inferences_plot_validation(data_loaded, path_to_res_mod, name, map_iso_gene = map_iso_gene)
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
        save_run_inferences_plot_validation(data_loaded, path_to_res_mod, name, long_mcmc = FALSE)
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
                                FDR_thd = 0.01
        )
        save_run_inferences_plot_validation(data_loaded, path_to_res_mod, name)
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
        save_run_inferences_plot_validation(data_loaded, path_to_res_mod, name)
      }
      pp = plot_roc_model(path_to_res_mod, name, protease, name_models)
      ggsave(glue("{path_to_res_mod}/{name}.png"))
    }
    pp = validate_all_protease(proteases, name, name_models)
    ggsave(glue("{PATH_TO_RES}/{name}/{name}.png"))
    
    
  }
  ############################################################################################################
  # Prior robustness
  ############################################################################################################
  # Overwrite inference() with an extra argument to set own prior
  inference = function(loaded_data,
                       prior = 0.1,
                       map_iso_gene = NULL,
                       n_cores = 1,
                       K = 2000,
                       burn_in = 1000,
                       thin = 1) {
    
    if(is.null(map_iso_gene)){
      map_iso_gene = ""
    }
    
    input_check_inference(loaded_data, map_iso_gene, n_cores, K, burn_in, thin)
    
    if (is.null(loaded_data$PROTEIN_DF$TPM)) {
      message("Transcriptomics data not loaded. Inference will be based only on proteomics data.")
      loaded_data$prior = 0
    } else {
      loaded_data$prior = prior
    }
    if (file.exists(map_iso_gene)) {
      map_iso_gene_file = fread(map_iso_gene, header = FALSE)
    } 
    
    names(loaded_data) = formalArgs(set_MCMC_args)
    args_MCMC = do.call("set_MCMC_args", loaded_data)
    args_MCMC$params = list(n_cores = n_cores, K = K, burn_in = burn_in, thin = thin, PEP = loaded_data$PEP)
    sel_unique = loaded_data$PROTEIN_DF$Y_unique > 0
    rm(loaded_data)
    
    if (args_MCMC$params$PEP) {
      results_MCMC = do.call("run_MCMC_pep", args_MCMC)
    } else {
      results_MCMC = do.call("run_MCMC", args_MCMC)
    }
    
    if (args_MCMC$params$n_cores > 1) {
      old_order = unlist(lapply(results_MCMC$groups, function(x){x$proteins})
      )
      old_order = c(old_order, results_MCMC$one_pept_one_prot)
      old_order = sort(old_order, index.return = TRUE)$ix
      results_MCMC$PI = results_MCMC$PI[, old_order]
      results_MCMC$Y = results_MCMC$Y[, old_order]
      #results_MCMC$isoform_results = results_MCMC$isoform_results[old_order, ]
    }
    
    results_MCMC$isoform_results = stat_from_MCMC_Y(results_MCMC$Y)
    results_MCMC = get_res_MCMC(results_MCMC, args_MCMC$prot_df$protein_name)
    
    if (!is.null(args_MCMC$prot_df$TPM)) {
      results_MCMC$isoform_results = stat_from_TPM(results_MCMC$isoform_results, args_MCMC$prot_df$TPM, results_MCMC$PI)
      reorder_col = c("Isoform", "Prob_present", "Abundance", "CI_LB", "CI_UB",
                      "Pi", "Pi_CI_LB", "Pi_CI_UB", "TPM", "Log2_FC", "Prob_prot_inc")
    } else {
      reorder_col = c("Isoform", "Prob_present", "Abundance", "CI_LB", "CI_UB",
                      "Pi", "Pi_CI_LB", "Pi_CI_UB")
    }
    
    if (file.exists(map_iso_gene)) {
      results_MCMC$isoform_results = merge(results_MCMC$isoform_results, map_iso_gene_file, by.x = "Isoform", by.y = "V1")
      colnames(results_MCMC$isoform_results)[ncol(results_MCMC$isoform_results)] = "Gene"
      
      res_norm = normalize_by_gene(results_MCMC, tpm = !is.null(args_MCMC$prot_df$TPM))
      reorder_col = c("Gene", reorder_col)
      
      results_MCMC$Y = aggregate.data.frame(t(results_MCMC$Y), by = list(results_MCMC$isoform_results$Gene), FUN = sum)
      Gene = results_MCMC$Y[, 1]
      results_MCMC$Y = t(results_MCMC$Y[, -1])
      # 0.95 CI for protein abundance:
      CI = hdi(results_MCMC$Y, credMass = 0.95)
      gene_abundance = data.frame(Gene = Gene,
                                  Abundance = colMeans(results_MCMC$Y),
                                  CI_LB = CI[1, ],
                                  CI_UB = CI[2, ]
      )
    } else {
      res_norm = NULL
      gene_abundance = NULL
    }
    
    # add a small threshold to isoform with unique peptides
    results_MCMC$isoform_results$Prob_present[sel_unique] = results_MCMC$isoform_results$Prob_present[sel_unique] + 0.0001
    
    list(isoform_results = results_MCMC$isoform_results[, reorder_col],
         normalized_isoform_results = res_norm,
         gene_abundance = gene_abundance
    )
  }
  environment(inference) = environment(load_data)
  
  prior_grid = c(0, 0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.4, 0.8, 1)
  save(prior_grid, file = glue("{PATH_WD}/Model_results/prior_grid"))
  log_output(glue("run_models_{DATA}_MM_prior_grid"))
  
  for (prior in prior_grid) {
    
    name = glue("OpenMS_mRNA_PEP_prior_{prior}")
    message(glue("---------- {name} ----------"))
    name_models = c("IsoBayes_mRNA", "TPM", "Baseline")
    
    for (protease in proteases) {
      message(protease)
      path_to_res_mod = glue("{PATH_TO_RES}/{name}/{protease}")
      if (run_model) {
        if(!dir.exists(path_to_res_mod)){dir.create(path_to_res_mod, recursive = TRUE)}
        data_loaded = load_data(path_to_peptides_psm = paste0(PATH_TO_DATA, "/Only", protease, "/merge_index_percolator_pep_switched.idXML"),
                                path_to_tpm = paste0(PATH_TO_DATA, "/mrna_isoform.tsv"),
                                input_type = "openMS",
                                abundance_type = "psm",
                                PEP = TRUE,
                                FDR_thd = 0.01
        )
        save(data_loaded, file = glue("{path_to_res_mod}/{name}_data_loaded.RData"))
        
        res = inference(data_loaded, prior, map_iso_gene = NULL, n_cores = 8, K = 2000, thin = 1)
        save(res, file = glue("{path_to_res_mod}/{name}_MCMC.RData"))
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
