get_roc = function(benchmark, name_models, protease = ""){
  normalize = function(x) {
    return((x- min(x)) /(max(x)-min(x)))
  }
  
  name_models = sort(name_models)
  scores = list()
  labels = list()
  for (name_model in name_models) {
    id_na = !is.na(benchmark[, name_model])
    scores = append(scores, list(benchmark[id_na, name_model]))
    labels = append(labels, list(benchmark[id_na, glue("Present_{name_model}")]))
  }
  roc_pr_curves = precrec::mmdata(scores, labels, name_models, dsids = 1:length(name_models))
  roc_pr_curves = precrec::evalmod(roc_pr_curves)
  
  auc_roc = round(precrec::auc(roc_pr_curves)[precrec::auc(roc_pr_curves)[, "curvetypes"] == "ROC", "aucs"], 3)
  sum_stat = data.frame(Model = name_models, AUC = auc_roc)
  labs_roc = name_models#paste0(name_models, "\nAUC: ", auc_roc)
  palette_models = PALETTE_MODELS[name_models]
  
  pp = autoplot(roc_pr_curves, "ROC")
  
  check = grepl("no_unique", name_models)
  if(any(check)){
    pp$mapping = aes(x = x, y = y, colour = modname, linetype = modname)
    linetype_models = palette_models
    linetype_models = rep("solid", length(linetype_models))
    linetype_models[check] = "dotted"
  }
  
  pp = pp + scale_colour_manual("", values = palette_models, labels = labs_roc)
  
  if(any(check)){
    pp = pp + scale_linetype_manual("", values = linetype_models, labels = labs_roc)
  }else{
    pp = pp + geom_line(linewidth = 1) + labs(title = glue("{DATA_name}{protease}")) +
      theme(plot.title = element_text(size = 11, hjust = 0.5),
            axis.title = element_text(size = 10),
            #legend.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(size = 8, angle = 45, vjust = 0.8),
            axis.text.y = element_text(size = 8),
            legend.text = element_text(size = 11)
      )
  }
  
  pp = pp + xlab("FPR") + ylab("TPR")
  
  list(gplot = pp, sum_stat = sum_stat)
}
