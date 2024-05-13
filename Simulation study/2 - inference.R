rm(list =ls())

setwd("~/Desktop/Proteomics/SIMULATION STUDY/Simulation study - MN")

names = list.files("simulated data/")

PERFORMANCE = data.frame(AUC = rep(0,length(names)), 
                         log10corr = rep(0,length(names)), 
                         HPDCI_coverage = rep(0,length(names)),
                         abundance_absent = rep(0,length(names)),
                         abundance_present = rep(0,length(names)),
                         prob_absent = rep(0,length(names)),
                         prob_present = rep(0,length(names)))

rownames(PERFORMANCE) = substring(names, 1, nchar(names) - 6)

# LOOP OVER proteases:
for(id in 1:length(names) ){
  name = names[id]
  
  filename = paste0("simulated data/", name)
  load(filename)
  
  # DO NOT Account for protein length:
  DATA$PROTEIN_DF$protein_length = 1
  
  # remove TPMs
  DATA$PROTEIN_DF$TPM = NULL

  library(IsoBayes)
  set.seed(169612)
  res = inference(DATA)[[1]]
  
  RES = merge(res, ground_truth, all.y = TRUE, by = "Isoform")
  rm(res); rm(ground_truth)
  
  head(RES)
  
  # AUC:
  library(pROC)
  PERFORMANCE[id, 1] = auc(RES$presence, 1-RES$Prob_present)
  
  # log10-corr and scatterplot:
  if(FALSE){
    plot( log10(1 + RES$Abundance), log10(1 + RES$abundance))
    abline(0,1)
  }
  PERFORMANCE[id, 2] = cor( log10(1 + RES$Abundance), log10(1 + RES$abundance))
  
  # 0.95 level HPD CI coverage:
  PERFORMANCE[id, 3] = mean((RES$abundance >= RES$CI_LB) & (RES$abundance <= RES$CI_UB))
  
  if(FALSE){
    # ROC presence:
    library(iCOBRA)
    DF_COBRA <- COBRAData(
      pval = data.frame(IsoBayes = 1-RES$Prob_present),
      truth = data.frame(status = RES$presence))
    perf <- calculate_performance(DF_COBRA, binary_truth = "status")
    cobra_plot <- prepare_data_for_plot(perf, facetted = TRUE)
    # plot ROC curve
    plot_roc(cobra_plot)
  }
  
  # average abundance for present and absent isoforms:
  PERFORMANCE[id,4:5] =  tapply(RES$Abundance, RES$presence, mean)
  
  # average Prob for present and absent isoforms:
  PERFORMANCE[id,6:7] = tapply(RES$Prob_present, RES$presence > 0, mean)
  
  print(PERFORMANCE[id, ])
}
round(PERFORMANCE,2)
colMeans(PERFORMANCE)

save(PERFORMANCE, file = "results/PERFORMANCE.RData")


# AUC , logCORR high
# coverage > 0.95
# abundance present = 51 * abundance absent
# prob present = 3.2 * prob absent

PERFORMANCE = rbind(PERFORMANCE, colMeans(PERFORMANCE))

library(xtable)
xtable(PERFORMANCE)

\begin{table}[ht]
\centering
\begin{tabular}{rrrrrrrr}
\hline
& AUC & log2corr & HPDCI\_coverage & abundance\_absent & abundance\_present & prob\_absent & prob\_present \\ 
\hline
jurkat\_ArgC & 0.91 & 0.86 & 0.98 & 0.68 & 5.73 & 0.27 & 0.78 \\ 
jurkat\_AspN & 0.93 & 0.89 & 0.98 & 0.58 & 8.74 & 0.26 & 0.79 \\ 
jurkat\_Chym & 0.87 & 0.83 & 0.98 & 0.70 & 4.42 & 0.31 & 0.73 \\ 
jurkat\_GluC & 0.93 & 0.90 & 0.98 & 0.60 & 9.14 & 0.25 & 0.79 \\ 
jurkat\_LysC & 0.95 & 0.92 & 0.97 & 0.46 & 11.85 & 0.23 & 0.82 \\ 
jurkat\_Trypsin & 0.95 & 0.92 & 0.98 & 0.46 & 11.31 & 0.22 & 0.82 \\ 
WTC-11\_AspN & 0.92 & 0.82 & 0.96 & 0.79 & 7.30 & 0.29 & 0.81 \\ 
WTC-11\_Chymo & 0.89 & 0.78 & 0.94 & 1.15 & 9.63 & 0.35 & 0.80 \\ 
WTC-11\_LysC & 0.93 & 0.86 & 0.97 & 0.62 & 7.44 & 0.28 & 0.82 \\ 
WTC-11\_Trypsin & 0.94 & 0.88 & 0.97 & 0.65 & 7.28 & 0.28 & 0.85 \\ 
11 & 0.92 & 0.86 & 0.97 & 0.67 & 8.28 & 0.27 & 0.80 \\ 
\hline
\end{tabular}
\end{table}