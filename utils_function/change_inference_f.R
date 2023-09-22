change_inference_f = function(){
  env_iso = environment(load_data)
  inference_mod = function (loaded_data, map_iso_gene = NULL, n_cores = 1, K = 2000, 
                        burn_in = 1000, thin = 1) 
  {
    if (is.null(map_iso_gene)) {
      map_iso_gene = ""
    }
    input_check_inference(loaded_data, map_iso_gene, n_cores, 
                          K, burn_in, thin)
    if (is.null(loaded_data$PROTEIN_DF$TPM)) {
      message("Transcriptomics data not loaded. Inference will be based only on proteomics data.")
      loaded_data$prior = 0
    }
    else {
      loaded_data$prior = 0.1
    }
    if (file.exists(map_iso_gene)) {
      map_iso_gene_file = fread(map_iso_gene, header = FALSE)
    }
    names(loaded_data) = formalArgs(set_MCMC_args)
    args_MCMC = do.call("set_MCMC_args", loaded_data)
    args_MCMC$params = list(n_cores = n_cores, K = K, burn_in = burn_in, 
                            thin = thin, PEP = loaded_data$PEP)
    sel_unique = loaded_data$prot_df$Y_unique > 0
    rm(loaded_data)
    if (args_MCMC$params$PEP) {
      results_MCMC = do.call("run_MCMC_pep", args_MCMC)
    }
    else {
      results_MCMC = do.call("run_MCMC", args_MCMC)
    }
    if (args_MCMC$params$n_cores > 1) {
      old_order = unlist(lapply(results_MCMC$groups, function(x) {
        x$proteins
      }))
      old_order = c(old_order, results_MCMC$one_pept_one_prot)
      old_order = sort(old_order, index.return = TRUE)$ix
      results_MCMC$PI = results_MCMC$PI[, old_order]
      results_MCMC$Y = results_MCMC$Y[, old_order]
    }
    chain_Y = results_MCMC$Y
    results_MCMC$isoform_results = stat_from_MCMC_Y(results_MCMC$Y)
    results_MCMC = get_res_MCMC(results_MCMC, args_MCMC$prot_df$protein_name)
    if (!is.null(args_MCMC$prot_df$TPM)) {
      results_MCMC$isoform_results = stat_from_TPM(results_MCMC$isoform_results, 
                                                   args_MCMC$prot_df$TPM, results_MCMC$PI)
      reorder_col = c("Isoform", "Prob_present", "Abundance", 
                      "CI_LB", "CI_UB", "Pi", "Pi_CI_LB", "Pi_CI_UB", "TPM", 
                      "Log2_FC", "Prob_prot_inc")
    }
    else {
      reorder_col = c("Isoform", "Prob_present", "Abundance", 
                      "CI_LB", "CI_UB", "Pi", "Pi_CI_LB", "Pi_CI_UB")
    }
    if (file.exists(map_iso_gene)) {
      results_MCMC = map_isoform_to_gene(results_MCMC, map_iso_gene_file)
      res_norm = normalize_by_gene(results_MCMC, tpm = !is.null(args_MCMC$prot_df$TPM))
      reorder_col = c("Gene", reorder_col)
      results_MCMC$Y = aggregate.data.frame(t(results_MCMC$Y), 
                                            by = list(results_MCMC$isoform_results$Gene), FUN = sum)
      Gene = results_MCMC$Y[, 1]
      results_MCMC$Y = t(results_MCMC$Y[, -1])
      CI = hdi(results_MCMC$Y, credMass = 0.95)
      gene_abundance = data.frame(Gene = Gene, Abundance = colMeans(results_MCMC$Y), 
                                  CI_LB = CI[1, ], CI_UB = CI[2, ])
    }
    else {
      res_norm = NULL
      gene_abundance = NULL
    }
    results_MCMC$isoform_results$Prob_present[sel_unique] = results_MCMC$isoform_results$Prob_present[sel_unique] + 
      1e-04
    list(isoform_results = results_MCMC$isoform_results[, reorder_col], 
         normalized_isoform_results = res_norm, gene_abundance = gene_abundance,
         chain_Y = chain_Y)
  }
  environment(inference_mod) = env_iso
  
  inference_mod
}
