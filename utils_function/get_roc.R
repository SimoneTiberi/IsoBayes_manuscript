get_roc = function(benchmark, name_models){
  name_models = sort(name_models)
  scores = list()
  labels = list()
  for (name_model in name_models) {
    scores = append(scores, list(benchmark[, name_model]))
    labels = append(labels, list(benchmark$Present))
  }
  roc_pr_curves = precrec::mmdata(scores, labels, name_models)
  roc_pr_curves = precrec::evalmod(roc_pr_curves)
  
  auc_roc = round(precrec::auc(roc_pr_curves)[precrec::auc(roc_pr_curves)[, "curvetypes"] == "ROC", "aucs"], 3)
  #auc_pr = round(precrec::auc(roc_pr_curves)[precrec::auc(roc_pr_curves)[, "curvetypes"] == "PRC", "aucs"], 3)
  
  sum_stat = data.frame(Model = name_models, AUC = auc_roc)
  
  labs_roc = paste0(name_models, "\nAUC: ", auc_roc)
  #labs_pr = paste0(name_models, "\nAUC: ", auc_pr)
  
  palette_models = PALETTE_MODELS[name_models]
  
  pp = autoplot(roc_pr_curves, "ROC") +
    scale_colour_manual("Model", values = palette_models, labels = labs_roc) +
    geom_line(linewidth=0.8) +
    labs(title = "ROC curves") +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10)
    )
  
  list(gplot = pp, sum_stat = sum_stat)
  #autoplot(roc_pr_curves, "PR") + scale_colour_discrete(name="Models", breaks = name_models, labels=labs_pr) +
  #ggtitle(paste0("PR curve ", protease))
}
