add_mrna = function(input, protease, pep, validation_dat){
  load(glue("{PATH_RES_MODEL}/{input}{pep}/{protease}/{input}{pep}_MCMC.RData"))
  mrna_data = data.table::fread(glue("{PATH_DATA}/mrna_isoform.tsv"))
  colnames(mrna_data)[grep("tpm", colnames(mrna_data))] = "TPM"
  res$isoform_results = merge(res$isoform_results, mrna_data[, c("isoname", "TPM")], by.x = "Isoform",
                              by.y = "isoname", all.x = TRUE)
  res$isoform_results = res$isoform_results[!duplicated(res$isoform_results$Isoform), ]
  P_TPM = res$isoform_results$TPM/sum(res$isoform_results$TPM)
  res$isoform_results$Log2_FC = log2(res$isoform_results$Pi/P_TPM)
  
  res$isoform_results$Prob_prot_inc = vapply(seq_len(nrow(res$isoform_results)), function(i){
    mean(res$chain_Y[, i] > P_TPM[i])}, FUN.VALUE = numeric(1) )
  
  validation_dat = merge(validation_dat, res$isoform_results[, c("Isoform", "Prob_prot_inc", "Log2_FC", "TPM")], by = "Isoform")
  
  validation_dat
}
