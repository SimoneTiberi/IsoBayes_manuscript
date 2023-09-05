plot_results = function(res, VALIDATION_DF_prot, data_loaded, path_to_data, protease, name){
  
  validation_dat = merge(res$isoform_results, VALIDATION_DF_prot, by.y = "proteins", by.x = "Isoform")
  validation_dat = merge(validation_dat, data_loaded$PROTEIN_DF[, c("protein_name", "Y_unique")],
                         by.y = "protein_name", by.x = "Isoform", all.x = TRUE)
  validation_dat_no_unique = validation_dat[validation_dat$Y_unique == 0, ]
  
  rocobj = pROC::roc(validation_dat$Present, validation_dat$Prob_present, quiet = TRUE)
  if (!is.null(validation_dat$TPM)) {
    rocobj_tpm = pROC::roc(validation_dat$Present, validation_dat$TPM, quiet = TRUE)
  } else {
    rocobj_tpm = NULL
  }
  rocobj_baseline = pROC::roc(validation_dat$Present, ifelse(validation_dat$Y_unique > 0, 1, 0), quiet = TRUE) # detect an isoform in the original data IF it has a unique isoform
  rocobj2 = pROC::roc(validation_dat_no_unique$Present, validation_dat_no_unique$Prob_present, quiet = TRUE)
  if (!is.null(validation_dat$TPM)) {
    rocobj2_tpm = pROC::roc(validation_dat_no_unique$Present, validation_dat_no_unique$TPM, quiet = TRUE)
  } else {
    rocobj2_tpm = NULL
  }
  
  if(is.null(rocobj_tpm)){
    pROC::ggroc(list(our_model = rocobj, our_model2 = rocobj2, unique_only = rocobj_baseline)) +
      geom_abline(slope = 1, intercept = 1) +
      scale_colour_discrete(name="Models",
                            breaks = c("our_model", "our_model2", "unique_only"),
                            labels=c(paste0("Model_with_UP\nAUC: ", round(rocobj$auc, 3)),
                                     paste0("Model_without_UP\nAUC: ", round(as.numeric(rocobj2$auc), 3)),
                                     paste0("Baseline\nAUC: ", round(rocobj_baseline$auc, 3)))
      )
    # `comparison` is a global variable
    ggsave(glue("{path_to_data}/ROC_{name}.png"))
    
  }else{
    pROC::ggroc(list(our_model = rocobj, our_model2 = rocobj2, tpm_only = rocobj_tpm,
                     tpm_only_2 = rocobj2_tpm, unique_only = rocobj_baseline)) +
      geom_abline(slope = 1, intercept = 1) +
      scale_colour_discrete(name="Models",
                            breaks = c("our_model", "our_model2", "tpm_only", "tpm_only_2", "unique_only"),
                            labels=c(paste0("Model_with_UP\nAUC: ", round(rocobj$auc, 3)),
                                     paste0("Model_without_UP\nAUC: ", round(as.numeric(rocobj2$auc), 3)),
                                     paste0("Tpm_only\nAUC: ", round(as.numeric(rocobj_tpm$auc), 3)),
                                     paste0("Tpm_only_without_UP\nAUC: ", round(as.numeric(rocobj2_tpm$auc), 3)),
                                     paste0("Baseline\nAUC: ", round(rocobj_baseline$auc, 3)))
      )
    ggsave(glue("{path_to_data}/ROC_{name}.png"))
  }
}
