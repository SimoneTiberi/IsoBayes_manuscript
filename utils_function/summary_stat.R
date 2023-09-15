df = rbind()
for (DATA in c("jurkat", "wtc11")) {
  PATH_TO_DATA = glue("{PATH_WD}/Data/{DATA}")
  if(DATA == "jurkat"){
    proteases = c("ArgC", "AspN", "Chym", "GluC", "LysC", "Trypsin")
  }else{
    proteases = c("AspN", "Chym", "LysC", "Trypsin")
  }
  for (protease in proteases) {
    
    data_loaded = load_data(path_to_peptides_psm = paste0(PATH_TO_DATA, "/Only", protease, "/AllPeptides.psmtsv"),
                            path_to_peptides_intensities = paste0(PATH_TO_DATA, "/Only", protease, "/AllQuantifiedPeptides.tsv"),
                            path_to_tpm = paste0(PATH_TO_DATA, "/mrna_isoform.tsv"),
                            input_type = "metamorpheus",
                            abundance_type = "psm",
                            PEP = FALSE,
                            FDR_thd = 0.01
    )
    vec = c(DATA, protease, round(sum(data_loaded$PROTEIN_DF$Y_unique > 0) / length(data_loaded$PROTEIN_DF$Y_unique), 2))
    df = rbind(df, vec)
  }
}

