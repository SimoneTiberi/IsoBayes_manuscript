################## ################## ################## ##################
# Study multi-mapping peptides in MM PSM counts:
################## ################## ################## ##################
rm(list =ls())

setwd("~/Desktop/Proteomics")

# loop over datasets:
datasets = c("jurkat", "wtc11")
MULTI_jurkat = matrix(NA, nrow = 4, ncol = 6)
MULTI_WTC11 = matrix(NA, nrow = 4, ncol = 4)

for(data in datasets){
  
  path = paste0("Data/", data, "/")
  names = list.files(path)  
  names = names[grep("Only", names)]
  if(data == "jurkat"){
    colnames(MULTI_jurkat) = names
  }else{
    colnames(MULTI_WTC11) = names
  }
  
  # LOOP OVER proteases:
  for(id in 1:length(names) ){
    name = names[id]
    
    # set parameters:
    path_to_peptides_psm = paste0(path, name, "/AllPeptides.psmtsv")
    
    library(IsoBayes)
    SE = generate_SE(path_to_peptides_psm = path_to_peptides_psm,
                     PEP = FALSE, 
                     input_type = "metamorpheus",
                     abundance_type = "psm",
                     FDR_thd = 0.01)
    
    PSM_data_loaded = input_data(SE)
    
    EC_numeric = unique(IsoBayes:::convert_EC_to_num(strsplit(SE$EC, split = "\\|"),
                                                     unique(unlist(strsplit(SE$EC, split = "\\|")))))
    
    n_proteins = sapply(EC_numeric, length)
    avg_multi = mean(n_proteins > 1)
    n_multi = mean(n_proteins[n_proteins > 1])
    
    multi = sum(PSM_data_loaded$PEPTIDE_DF$Y)
    uni = sum(PSM_data_loaded$PROTEIN_DF$Y_unique)
    
    # compute abundance associated to each gene
    isoform = PSM_data_loaded$PROTEIN_DF$protein_name
    gene = sapply(strsplit(isoform, "-"), function(x) x[[1]])
    
    head(PSM_data_loaded$PEPTIDE_DF)
    PSM_data_loaded$PEPTIDE_DF$n_genes = sapply(PSM_data_loaded$PEPTIDE_DF$EC_numeric, function(x){
      length(unique(gene[x]))
    })
    
    if(data == "jurkat"){
      MULTI_jurkat[1,id] = multi/(multi+uni)
      MULTI_jurkat[2,id] = avg_multi
      MULTI_jurkat[3,id] = n_multi
      MULTI_jurkat[4,id] = sum(PSM_data_loaded$PEPTIDE_DF$Y[ PSM_data_loaded$PEPTIDE_DF$n_genes > 1 ])/(multi+uni)
    }else{
      MULTI_WTC11[1,id] = multi/(multi+uni)
      MULTI_WTC11[2,id] = avg_multi
      MULTI_WTC11[3,id] = n_multi
      MULTI_WTC11[4,id] = sum(PSM_data_loaded$PEPTIDE_DF$Y[ PSM_data_loaded$PEPTIDE_DF$n_genes > 1 ])/(multi+uni)
    }
  }
}

rownames(MULTI_jurkat) = rownames(MULTI_WTC11) = c("Multi counts",
                                                   "Multi peptides",
                                                   "N protein per multi peptide",
                                                   "Gene multi counts")
# add mean:
MULTI = rbind( t(MULTI_jurkat), t(MULTI_WTC11))
MULTI = rbind(MULTI,  Average = colMeans(MULTI) )

\begin{table}[ht]
\centering
\begin{tabular}{rrrrr}
\hline
& Multi counts & Multi peptides & N protein per multi peptide & Gene multi counts \\ 
\hline
OnlyArgC & 0.54 & 0.60 & 3.50 & 0.10 \\ 
OnlyAspN & 0.45 & 0.58 & 3.44 & 0.06 \\ 
OnlyChym & 0.51 & 0.64 & 3.64 & 0.07 \\ 
OnlyGluC & 0.45 & 0.59 & 3.40 & 0.06 \\ 
OnlyLysC & 0.42 & 0.55 & 3.37 & 0.07 \\ 
OnlyTrypsin & 0.43 & 0.56 & 3.44 & 0.08 \\ 
OnlyAspN.1 & 0.60 & 0.59 & 3.52 & 0.47 \\ 
OnlyChym.1 & 0.63 & 0.65 & 3.84 & 0.49 \\ 
OnlyLysC.1 & 0.51 & 0.52 & 3.31 & 0.41 \\ 
OnlyTrypsin.1 & 0.47 & 0.47 & 3.13 & 0.37 \\ 
Average & 0.50 & 0.58 & 3.46 & 0.22 \\ 
\hline
\end{tabular}
\end{table}
