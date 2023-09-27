save_run_inferences = function(data_loaded, path_to_res_mod, name, map_iso_gene = NULL, long_mcmc = FALSE, inference_function = NULL){
  save(data_loaded, file = glue("{path_to_res_mod}/{name}_data_loaded.RData"))
  
  if(!is.null(inference_function)){
    inference = inference_function
  }
  
  res = inference(data_loaded, map_iso_gene = map_iso_gene, n_cores = 8, K = 2000, thin = 1)
  save(res, file = glue("{path_to_res_mod}/{name}_MCMC.RData"))
  
  if(long_mcmc){
    res = inference(data_loaded, map_iso_gene = map_iso_gene, n_cores = 8, K = 10000, thin = 5)
    save(res, file = glue("{path_to_res_mod}/{name}_MCMC_10000.RData"))
  }
}
