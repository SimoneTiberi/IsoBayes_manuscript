merge_validation = function(protease, name, data_loaded_file, MCMC_file, validation_file){
  load(data_loaded_file)
  load(MCMC_file)
  load(validation_file)
  validation_dat = merge(res$isoform_results, VALIDATION_DF_prot, by.y = "proteins", by.x = "Isoform")
  validation_dat = merge(validation_dat, data_loaded$PROTEIN_DF[, c("protein_name", "Y_unique")],
                         by.y = "protein_name", by.x = "Isoform", all.x = TRUE)
  validation_dat$Baseline = as.numeric(validation_dat$Y_unique > 0)
  
  save(validation_dat, file = glue("{PATH_WD}/Model_results/{DATA}/{name}/{protease}/Merged_validation_res_{name}.RData"))
}
