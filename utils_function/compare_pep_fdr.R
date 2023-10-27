compare_pep_fdr = function(input, sub_selected_models, proteases, models){
  benchmark_df_all = list()
  
  for(protease in proteases){
    benchmark_df = list()
    
    for (model in sub_selected_models) {
      for (pep in c("_PEP", "")) {
        attribute_model = glue("{models[[model]][2]}{models[[model]][1]}") 
        # load res and validation merged together
        load(glue("{PATH_RES}/{input}{attribute_model}{pep}/{protease}/Merged_validation_res_{input}{attribute_model}{pep}.RData"))
        
        if(!grepl("_mRNA", model)){
          load(glue("{PATH_RES}/{input}{attribute_model}{pep}/{protease}/{input}{attribute_model}{pep}_MCMC.RData"))
          mrna_data = data.table::fread(glue("{PATH_TO_DATA}/mrna_isoform.tsv"))
          colnames(mrna_data)[grep("tpm", colnames(mrna_data))] = "TPM"
          
          validation_dat = merge(validation_dat,
                                 mrna_data[, c("isoname", "TPM")],
                                 by.x = "Isoform",
                                 by.y = "isoname",
                                 all.x = TRUE)
          
          validation_dat[is.na(validation_dat)] = 0
          validation_dat = validation_dat[!duplicated(validation_dat$Isoform), ]
          P_TPM = validation_dat$TPM/sum(validation_dat$TPM)
          validation_dat$Log2_FC = log2(validation_dat$Pi/P_TPM)
          
          validation_dat$Prob_prot_inc = vapply(seq_len(nrow(validation_dat)), function(i){
            mean(res$chains$PI[, i] > P_TPM[i])}, FUN.VALUE = numeric(1))
        }
        
        colnames(validation_dat)[grep("tpm", colnames(validation_dat))] = "tpm_validation"
        colnames(validation_dat) = glue("{colnames(validation_dat)}_IsoBayes{attribute_model}{pep}")
        colnames(validation_dat)[grep(glue("Prob_present_IsoBayes{attribute_model}{pep}"), colnames(validation_dat))] = glue("IsoBayes{attribute_model}{pep}")
        colnames(validation_dat)[grep(glue("Present_IsoBayes{attribute_model}{pep}"), colnames(validation_dat))] = glue("Present_IsoBayes{attribute_model}{pep}")
        colnames(validation_dat)[grep(glue("Isoform_IsoBayes{attribute_model}{pep}"), colnames(validation_dat))] = "Isoform"
        validation_dat = validation_dat[!duplicated(validation_dat$Isoform), ]
        benchmark_df = append(benchmark_df, list(validation_dat))
      }
    }
    benchmark_df = concat_models(benchmark_df, union = FALSE)
    
    plot_tab = get_roc(benchmark_df, paste0(sub_selected_models, rep(c("_PEP", ""), each=2)), protease = glue("- {protease}"))
    ggsave(glue("{PATH_RES_roc}/{protease}/ROC_{input}_pep_vs_no_pep.png"), plot = plot_tab$gplot)
    save(plot_tab, file = glue("{PATH_RES_roc}/{protease}/ROC_{input}_pep_vs_no_pep.rdata"))
    write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_roc}/{protease}/SumTab_{input}_pep_vs_no_pep.csv"), row.names = FALSE)
    benchmark_df_all = rbind(benchmark_df_all, benchmark_df)
  }
  plot_tab = get_roc(benchmark_df_all, paste0(sub_selected_models, rep(c("_PEP", ""), each=2)))
  ggsave(glue("{PATH_RES_roc}/ROC_{input}_pep_vs_no_pep.png"), plot = plot_tab$gplot)
  save(plot_tab, file = glue("{PATH_RES_roc}/ROC_{input}_pep_vs_no_pep.rdata"))
  write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_roc}/SumTab_{input}_pep_vs_no_pep.csv"), row.names = FALSE)
  
  ## abundance of OpensMS vs MM
  plot_prob_change_data = rbind()
  plot_prob_change_data_extreme = rbind()
  for (pep in c("_PEP", "")) {
    for (mrna in c("", "_mRNA")) {
      scat_bench = scatterplot(log10(benchmark_df_all[, c(glue("Abundance_IsoBayes{mrna}{pep}"),
                                                          glue("Y_validation_IsoBayes{mrna}{pep}"))] + 1))  + 
        labs(x = "Log10(Estimated Abundance + 1)", y = "Log10(Validated Abundance + 1)") 
      
      ggsave(glue("{PATH_RES_roc}/scatterplot_{input}{pep}{mrna}_pep_vs_no_pep.png"), plot = scat_bench)
      save(scat_bench, file = glue("{PATH_RES_roc}/scatterplot_{input}{pep}{mrna}_pep_vs_no_pep.rdata"))
      
      # change mrna prot
      sub_data = benchmark_df_all[, c(glue("Prob_prot_inc_IsoBayes{mrna}{pep}"),
                                      glue("TPM_IsoBayes{mrna}{pep}"),
                                      glue("tpm_validation_IsoBayes{mrna}{pep}"),
                                      glue("P_Y_validation_IsoBayes{mrna}{pep}"),
                                      glue("Log2_FC_IsoBayes{mrna}{pep}"),
                                      glue("Y_validation_IsoBayes{mrna}{pep}"),
                                      glue("Pi_IsoBayes{mrna}{pep}")
      )
      ]
      sub_data = build_data_violin_plot(sub_data, 1.5e-06)
      sub_data = convert_numeric_to_class(sub_data, glue("Prob_prot_inc_IsoBayes{mrna}{pep}"),
                                          seq(0, 1, 0.2))
      sub_data_extreme = convert_numeric_to_class(sub_data, glue("Prob_prot_inc_IsoBayes{mrna}{pep}"),
                                                  quantiles = c(0, 0.01, 0.99, 1))
      sub_data_extreme = sub_data_extreme[sub_data_extreme$class_Prob_prot_inc != "(0.01, 0.99]", ]
      
      if(mrna == ""){
        sub_data$mrna = "No mRNA"
        sub_data_extreme$mrna = "No mRNA"
      }else{
        sub_data$mrna = "mRNA"
        sub_data_extreme$mrna = "mRNA"
      }
      if(pep == ""){
        sub_data$Model = "No PEP"
        sub_data_extreme$Model =  "No PEP"
      }else{
        sub_data$Model = "PEP"
        sub_data_extreme$Model =  "PEP"
      }
      
      plot_prob_change_data = rbind(plot_prob_change_data,
                                    sub_data[, c("class_Prob_prot_inc",
                                                 "Log2_FC_validation",
                                                 "Model", "mrna")]
      )
      plot_prob_change_data_extreme = rbind(plot_prob_change_data_extreme,
                                            sub_data_extreme[, c("class_Prob_prot_inc",
                                                                 "Log2_FC_validation",
                                                                 "Model", "mrna")]
                                            )
      ths = 1.5e-06
      p_tpm_adj = (sub_data[, glue("TPM_IsoBayes{mrna}{pep}")]+ths)/sum(sub_data[, glue("TPM_IsoBayes{mrna}{pep}")]+ths)
      pi = (sub_data[, glue("Pi_IsoBayes{mrna}{pep}")]+ths)/sum(sub_data[, glue("Pi_IsoBayes{mrna}{pep}")]+ths)
      sub_data$Log2_FC_adj = log2(pi/p_tpm_adj)
      
      scat_bench = scatterplot(sub_data[, c("Log2_FC_adj", "Log2_FC_validation")])  + 
        labs(x = "Log2-FC", y = "Validated Log2-FC")
      
      ggsave(glue("{PATH_RES_roc}/scatterplot_log2fc_{input}{pep}{mrna}_pep_vs_no_pep.png"), plot = scat_bench)
      save(scat_bench, file = glue("{PATH_RES_roc}/scatterplot_log2fc_{input}{pep}{mrna}_pep_vs_no_pep.rdata"))
    }
  }
  # plot prob change
  plot_change = plot_prob_change_group(plot_prob_change_data[plot_prob_change_data$mrna == "mRNA", ])
  ggsave(glue("{PATH_RES_roc}/change_mrna_prot_{input}_mRNA_pep_vs_no_pep.pdf"), plot = plot_change)
  save(plot_change, file = glue("{PATH_RES_roc}/change_mrna_prot_{input}_mRNA_pep_vs_no_pep.rdata"))
  
  plot_change = plot_prob_change_group(plot_prob_change_data[plot_prob_change_data$mrna != "mRNA", ])
  ggsave(glue("{PATH_RES_roc}/change_mrna_prot_{input}_pep_vs_no_pep.pdf"), plot = plot_change)
  save(plot_change, file = glue("{PATH_RES_roc}/change_mrna_prot_{input}_pep_vs_no_pep.rdata"))
  
  # plot prob change extreme
  plot_change = plot_prob_change_group(plot_prob_change_data_extreme[plot_prob_change_data_extreme$mrna == "mRNA", ])
  ggsave(glue("{PATH_RES_roc}/change_mrna_prot_extreme_{input}_mRNA_pep_vs_no_pep.pdf"), plot = plot_change)
  save(plot_change, file = glue("{PATH_RES_roc}/change_mrna_prot_extreme_{input}_mRNA_pep_vs_no_pep.rdata"))
  
  plot_change = plot_prob_change_group(plot_prob_change_data_extreme[plot_prob_change_data_extreme$mrna != "mRNA", ])
  ggsave(glue("{PATH_RES_roc}/change_mrna_prot_extreme_{input}_pep_vs_no_pep.pdf"), plot = plot_change)
  save(plot_change, file = glue("{PATH_RES_roc}/change_mrna_prot_extreme_{input}_pep_vs_no_pep.rdata"))
}