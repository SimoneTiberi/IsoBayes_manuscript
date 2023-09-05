plot_roc_model = function(path_to_res_mod, name, protease, name_models){
  
  validation_dat_path = glue("{PATH_WD}/Model_results/{DATA}/{name}/{protease}/Merged_validation_res_{name}")
  if(!file.exists(validation_dat_path)){
    merge_validation(protease, name,
                     data_loaded_file = glue("{path_to_res_mod}/{name}_data_loaded.RData"),
                     MCMC_file = glue("{path_to_res_mod}/{name}_MCMC.RData"),
                     validation_file = glue("{PATH_TO_DATA}/No{protease}/Validation_prot_psm")
                     )
    }
  load(validation_dat_path)
  colnames(validation_dat)[colnames(validation_dat) == "Prob_present"] = name_models[1]
  
  get_roc(validation_dat, name_models)
}
