validate_all_protease = function(proteases, name, name_models){
  all_validation_res = cbind()
  
  for (protease in proteases) {
    path_to_res_mod = paste0(PATH_TO_RES, "/", name, "/", protease)
    validation_dat_path = glue("{PATH_WD}/Model_results/{DATA}/{name}/{protease}/Merged_validation_res_{name}.RData")
    
    if(!file.exists(validation_dat_path)){
      merge_validation(protease, name,
                       data_loaded_file = glue("{path_to_res_mod}/{name}_data_loaded.RData"),
                       MCMC_file = glue("{path_to_res_mod}/{name}_MCMC.RData"),
                       validation_file = glue("{PATH_TO_DATA}/No{protease}/Validation_prot_psm")
      )
    }
  }
}