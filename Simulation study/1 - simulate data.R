rm(list =ls())

setwd("~/Desktop/Proteomics/SIMULATION STUDY")

# loop over datasets:
datasets = c("jurkat", "WTC-11")
for(data in datasets){
  
  path = paste0("data/", data, "/")
  names = list.files(path)  
  
  # LOOP OVER proteases:
  for(id in 1:length(names) ){
    name = names[id]
    
    # set parameters:
    path_to_peptides_psm = paste0(path, name, "/AllPeptides.psmtsv")
    
    # set path to TPMs:
    if(data == "jurkat"){
      tpm_path = paste0("TPM/jurkat_isoform_kallisto.tsv")
    }else{
      tpm_path = paste0("TPM/wtc-11_mrna_isoform.tsv")
    }
    
    library(IsoBayes)
    SE = generate_SE(path_to_peptides_psm = path_to_peptides_psm,
                     PEP = FALSE, 
                     input_type = "metamorpheus",
                     abundance_type = "psm",
                     FDR_thd = 0.01)
    
    PSM_data_loaded = input_data(SE, path_to_tpm = tpm_path)
    
    set.seed(169612)
    res = inference(PSM_data_loaded)[[1]]
    
    res$Abundance = round(res$Abundance)
    
    ground_truth = data.frame(Isoform = res$Isoform, 
                              presence = res$Abundance > 0.5,
                              abundance = res$Abundance,
                              TPM = res$TPM)
    N_proteins = nrow(ground_truth); N_proteins
    
    source("0 - compute_ECs.R")
    
    # get ECs:
    EC_numeric = obtain_EC_numeric(path_to_peptides_psm = path_to_peptides_psm,
                                   path_to_tpm = tpm_path,
                                   PEP = FALSE,
                                   FDR_thd = 0.01)
    min(unlist(EC_numeric))
    max(unlist(EC_numeric))
    
    # NOW, we compute the opposite: peptides associated to each protein.
    N_peptides = length(EC_numeric)
    EC_numeric = lapply(EC_numeric, unique)
    n_proteins_per_peptide = sapply(EC_numeric, length)
    peptides = rep(1:N_peptides, n_proteins_per_peptide)
    
    DF_protein_peptide = data.frame(Protein = unlist(EC_numeric), Peptide = peptides )
    head(DF_protein_peptide)
    Peptides_per_protein = split(DF_protein_peptide$Peptide, f =  DF_protein_peptide$Protein)
    head(Peptides_per_protein);
    
    N_proteins; length(Peptides_per_protein)
    
    # match protein ID
    X_protein = round(ground_truth$abundance)
    
    # sample peptide abundance:
    X_peptide = rep(0, N_peptides)
    
    set.seed(169612)
    for(i in 1:length(Peptides_per_protein)){
      x = X_protein[i]
      if(x > 0){
        index = Peptides_per_protein[[i]]
        len = length(index)
        
        #X_peptide[index] = X_peptide[index] + rpois(len, x/len)
        if(len > 1){
          # consider the PEP here?
          X_peptide[index] = X_peptide[index] + c(rmultinom(n = 1, size = x, prob = rep(1, len) ))
        }else{
          X_peptide[index] = X_peptide[index] + x
        }
      }
    }
    sum(X_protein); sum(X_peptide)
    
    N_peptides_per_protein = sapply(Peptides_per_protein, length)
    
    # simulated data in:
    head(X_peptide); length(X_peptide) # peptide abundance
    head(EC_numeric); length(EC_numeric) # proteins each peptide is associated to
    N_proteins; N_peptides
    
    TPM = ground_truth$TPM
    
    # ground truth:
    head(ground_truth); dim(ground_truth)
    
    # store data in the same format provided by load_data:
    PEPTIDE_DF = data.frame(Y = X_peptide)
    PEPTIDE_DF$EC_numeric = EC_numeric
    
    PROTEIN_DF = data.frame(protein_name = ground_truth$Isoform,
                            TPM = ground_truth$TPM,
                            protein_length = N_peptides_per_protein,
                            Y_unique = 0)
    
    DATA = list(PEPTIDE_DF = PEPTIDE_DF,
                PROTEIN_DF_unique = NULL,
                PROTEIN_DF = PROTEIN_DF,
                PEP = FALSE)
    
    filename = paste0("Simulation study - MN/simulated data/", data, "_", name, ".RData")
    
    save(DATA, ground_truth, 
         file = filename)
  }
}
