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

save(PERFORMANCE, file = "results/PERFORMANCE_mRNA.RData")

# AUC , logCORR very high
# coverage > 0.95
# abundance present = 185 * abundance absent
# prob present = 6.9 * prob absent


PERFORMANCE = rbind(PERFORMANCE, colMeans(PERFORMANCE))

library(xtable)
xtable(PERFORMANCE)

\begin{table}[ht]
\centering
\begin{tabular}{rrrrrrrr}
\hline
& AUC & log2corr & HPDCI\_coverage & abundance\_absent & abundance\_present & prob\_absent & prob\_present \\ 
\hline
jurkat\_ArgC & 0.96 & 0.95 & 0.99 & 0.24 & 6.01 & 0.09 & 0.68 \\ 
jurkat\_AspN & 0.97 & 0.97 & 0.99 & 0.16 & 8.97 & 0.10 & 0.70 \\ 
jurkat\_Chym & 0.96 & 0.96 & 0.98 & 0.16 & 4.91 & 0.09 & 0.68 \\ 
jurkat\_GluC & 0.98 & 0.97 & 0.99 & 0.15 & 9.40 & 0.09 & 0.71 \\ 
jurkat\_LysC & 0.98 & 0.98 & 0.98 & 0.14 & 11.98 & 0.10 & 0.72 \\ 
jurkat\_Trypsin & 0.98 & 0.98 & 0.99 & 0.13 & 11.46 & 0.09 & 0.72 \\ 
WTC-11\_AspN & 0.97 & 0.96 & 0.99 & 0.17 & 7.62 & 0.10 & 0.68 \\ 
WTC-11\_Chymo & 0.97 & 0.96 & 0.99 & 0.18 & 10.14 & 0.11 & 0.64 \\ 
WTC-11\_LysC & 0.96 & 0.96 & 0.98 & 0.16 & 7.68 & 0.11 & 0.73 \\ 
WTC-11\_Trypsin & 0.96 & 0.96 & 0.97 & 0.16 & 7.47 & 0.12 & 0.76 \\ 
11 & 0.97 & 0.97 & 0.98 & 0.16 & 8.56 & 0.10 & 0.70 \\ 
\hline
\end{tabular}
\end{table}