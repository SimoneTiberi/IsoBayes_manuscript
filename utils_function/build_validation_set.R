PATH_WD = commandArgs(trailingOnly = TRUE)[1]
DATA = commandArgs(trailingOnly = TRUE)[2]

# Set path, global variable and libraries
###########################################################################################
library(glue)
library(IsoBayes)

source(glue("{PATH_WD}/utils_function/load_AllPeptides.R"))
source(glue("{PATH_WD}/utils_function/build_validation.R"))
source(glue("{PATH_WD}/utils_function/log_output.R"))

PATH_TO_DATA = glue("{PATH_WD}/Data/{DATA}")
set.seed(0)
log_output(glue("building_validation_set_{DATA}"))

# build validation set
###########################################################################################
main = function(no_proteases){
  for (model in c("psm", "intensities")) {
    
    message(glue("Isoform Validation with MetaMorpheus data ({model})"))
    no_proteases = no_proteases[grepl(pattern = "No", no_proteases)]
    
    for (no_protease in no_proteases) {
      message(rev(unlist(strsplit(no_protease, "/")))[1])
      variables = c("Protein Accession", "QValue", "Decoy/Contaminant/Target", "Base Sequence")
      data_list = list(x = load_AllPeptides(no_protease, model), y = NULL)
      if(model == "intensities"){
        data_list = IsoBayes:::build_intensity(data_list$x, glue("{no_protease}/AllPeptides.psmtsv"), variables)
      }else{
        data_list$y = data_list$x$`PSM Count (unambiguous, <0.01 q-value)`
      }
      PEPTIDE_DF = IsoBayes:::build_peptide_df(data_list$y, data_list$x) ; rm(data_list)
      PEPTIDE_DF$n_proteins_per_peptide = sapply(strsplit(PEPTIDE_DF$EC, split = "\\|"), length)
      PEPTIDE_DF_prot = PEPTIDE_DF[PEPTIDE_DF$n_proteins_per_peptide == 1, ] ; rm(PEPTIDE_DF)
      PEPTIDE_DF_prot$Y[PEPTIDE_DF_prot$QValue > 0.01] = 0 # FDR > 1%
      
      VALIDATION_DF_prot = build_validation(Y = PEPTIDE_DF_prot$Y, var = PEPTIDE_DF_prot$EC, var_name = "proteins")
      
      if(model == "intensities"){
        VALIDATION_DF_prot$Y_validation = VALIDATION_DF_prot$Y_validation/sum(VALIDATION_DF_prot$Y_validation) * 10^5
        # round intensities to closest integer, BUT we add 0.5 so that very small intensities (between 0 and 0.5) are rounded to 1.
        VALIDATION_DF_prot$Y_validation = round(VALIDATION_DF_prot$Y_validation + 0.5)
      }
      
      tpm_df = data.table::fread(glue("{PATH_TO_DATA}/mrna_isoform.tsv"))
      VALIDATION_DF_prot = merge(VALIDATION_DF_prot, tpm_df, by.x = "proteins", by.y = "isoname", all.x = TRUE)
      VALIDATION_DF_prot$tpm[is.na(VALIDATION_DF_prot$tpm)] = 0
      #VALIDATION_DF_prot = na.omit(VALIDATION_DF_prot)
      
      save(VALIDATION_DF_prot, file = glue("{no_protease}/Validation_prot_{model}"))
    }
    rm(PEPTIDE_DF_prot) ; rm(VALIDATION_DF_prot)
    
    message(glue("Gene Validation with MetaMorpheus data ({model})"))
    for (no_protease in no_proteases) {
      message(rev(unlist(strsplit(no_protease, "/")))[1])
      variables = c("Protein Accession", "QValue", "Decoy/Contaminant/Target", "Base Sequence")
      data_list = list(x = load_AllPeptides(no_protease, model), y = NULL)
      if(model == "intensities"){
        data_list = IsoBayes:::build_intensity(data_list$x, glue("{no_protease}/AllPeptides.psmtsv"), variables)
      }else{
        data_list$y = data_list$x$`PSM Count (unambiguous, <0.01 q-value)`
      }
      PEPTIDE_DF = IsoBayes:::build_peptide_df(data_list$y, data_list$x) ; rm(data_list)
      EC_gene = lapply(strsplit(PEPTIDE_DF$EC, split = "\\|"), function(x){gsub("-.*", "", x)})
      EC_gene = lapply(EC_gene, function(x){unique(x)})
      PEPTIDE_DF$EC_gene = EC_gene
      
      PEPTIDE_DF$n_proteins_per_gene = sapply(PEPTIDE_DF$EC_gene, length)
      PEPTIDE_DF_prot = PEPTIDE_DF[PEPTIDE_DF$n_proteins_per_gene == 1, ] ; rm(PEPTIDE_DF)
      PEPTIDE_DF_prot$Y[PEPTIDE_DF_prot$QValue > 0.01] = 0 #FDR > 0.01
      PEPTIDE_DF_prot$EC_gene = unlist(PEPTIDE_DF_prot$EC_gene)
      
      VALIDATION_DF_prot = build_validation(Y = PEPTIDE_DF_prot$Y, var = PEPTIDE_DF_prot$EC_gene, var_name = "Gene")
      
      if(model == "intensities"){
        VALIDATION_DF_prot$Y_validation = VALIDATION_DF_prot$Y_validation/sum(VALIDATION_DF_prot$Y_validation) * 10^5
        # round intensities to closest integer, BUT we add 0.5 so that very small intensities (between 0 and 0.5) are rounded to 1.
        VALIDATION_DF_prot$Y_validation = round(VALIDATION_DF_prot$Y_validation + 0.5)
      }
      save(VALIDATION_DF_prot, file = glue("{no_protease}/Validation_gene_{model}"))
    }
    rm(PEPTIDE_DF_prot) ; rm(VALIDATION_DF_prot)
  }
}

main(no_proteases = list.dirs(PATH_TO_DATA, recursive = FALSE))