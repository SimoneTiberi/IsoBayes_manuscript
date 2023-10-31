PATH_WD = commandArgs(trailingOnly = TRUE)[1]
DATA = commandArgs(trailingOnly = TRUE)[2]

if(DATA == "wtc11"){
  DATA_name = "WTC-11"
}else{
  DATA_name = "jurkat"
}

# Set path, global variable and libraries
###########################################################################################
library(ggplot2)
library(readr)
library(glue)
source(glue("{PATH_WD}/utils_function/utils_benchmarking.R"))
source(glue("{PATH_WD}/utils_function/get_roc.R"))
source(glue("{PATH_WD}/utils_function/memory_plot.R"))
source(glue("{PATH_WD}/utils_function/run_time_plot.R"))
source(glue("{PATH_WD}/utils_function/log_output.R"))
source(glue("{PATH_WD}/utils_function/scatterplot.R"))
source(glue("{PATH_WD}/utils_function/utils_change_iso_mrna.R"))
source(glue("{PATH_WD}/utils_function/plot_prob_change.R"))

PATH_RES_COMPETITORS = glue("{PATH_WD}/Benchmark_results/{DATA}")
PATH_TO_DATA = glue("{PATH_WD}/Data/{DATA}")
PATH_TO_RES = glue("{PATH_WD}/Model_results/{DATA}")
load(glue("{PATH_WD}/utils_function/PALETTE_MODELS"))

log_output(glue("benchmarking_results_{DATA}"))

main = function(models, proteases){
    # COMPETITORS vs IsoBayes_openMS vs IsoBayes_openMS_mRNA [PERFORMANCE]
    benchmark_df_all = list()
    selected_models = c("IsoBayes", "IsoBayes_mRNA")
    
    for(protease in proteases){
      # load validation dataset from metamorpheus
      load(glue("{PATH_TO_DATA}/No{protease}/Validation_prot_psm"))
      
      benchmark_df = list(VALIDATION_DF_prot)
      colnames(benchmark_df[[1]])[grep("tpm", colnames(benchmark_df[[1]]))] = "tpm_validation"
      for (model in selected_models) {
        attribute_model = glue("{models[[model]][2]}{models[[model]][1]}") 
        load(glue("{PATH_TO_RES}/OpenMS{attribute_model}/{protease}/OpenMS{attribute_model}_MCMC.RData"))
        load(glue("{PATH_TO_RES}/OpenMS{attribute_model}/{protease}/OpenMS{attribute_model}_data_loaded.RData"))
        res$isoform_results = merge(res$isoform_results, data_loaded$PROTEIN_DF[, c("protein_name", "Y_unique")],
                                    by.x = "Isoform", by.y = "protein_name")
        if(grepl("_mRNA", model)){
          benchmark_df = append(benchmark_df, list(res$isoform_results))
        }else{
          mrna_data = data.table::fread(glue("{PATH_TO_DATA}/mrna_isoform.tsv"))
          colnames(mrna_data)[grep("tpm", colnames(mrna_data))] = "TPM"
          
          res$isoform_results = merge(res$isoform_results,
                                      mrna_data[, c("isoname", "TPM")],
                                      by.x = "Isoform",
                                      by.y = "isoname",
                                      all.x = TRUE)
          
          res$isoform_results[is.na(res$isoform_results)] = 0
          res$isoform_results = res$isoform_results[!duplicated(res$isoform_results$Isoform), ]
          P_TPM = res$isoform_results$TPM/sum(res$isoform_results$TPM)
          res$isoform_results$Log2_FC = log2(res$isoform_results$Pi/P_TPM)
          
          res$isoform_results$Prob_prot_inc = vapply(seq_len(nrow(res$isoform_results)), function(i){
            mean(res$chains$PI[, i] > P_TPM[i])}, FUN.VALUE = numeric(1) )
          
          benchmark_df = append(benchmark_df, list(res$isoform_results))
        }
      }
      
      for (i in seq_len(length(benchmark_df)-1)) {
        colnames(benchmark_df[[2]]) = paste0(colnames(benchmark_df[[2]]), "_", selected_models[i])
        benchmark_df[[1]] = merge(benchmark_df[[1]], benchmark_df[[2]],
                                  by.x = "proteins", by.y = paste0("Isoform_", selected_models[i]), all = T)
        benchmark_df[[2]] = NULL
      }
      # EPIFANY Epifany
      res = get_score_from_idXML(paste0(PATH_RES_COMPETITORS, "/",  protease, "/epifany.idXML"))
      colnames(res) = paste0(colnames(res), "_EPIFANY")
      benchmark_df = merge(benchmark_df, res, by.x = "proteins", by.y = "Isoform_EPIFANY", all = T)
      
      # FIDO
      res = get_score_from_idXML(paste0(PATH_RES_COMPETITORS, "/", protease, "/fido.idXML"))
      colnames(res) = paste0(colnames(res), "_Fido")
      benchmark_df = merge(benchmark_df, res, by.x = "proteins", by.y = "Isoform_Fido", all = T)
      
      # PIA
      path_file = paste0(PATH_RES_COMPETITORS, "/", protease, "/pia_results.mzTab")
      system(paste0("fgrep -v MTD ", path_file, " > ",  path_file, ".red"))
      res = read_delim(paste0(path_file, ".red"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
      res = as.data.frame(res[res$PRH == "PRT", c("accession", "best_search_engine_score[1]")])
      colnames(res) = c("Isoform_PIA", "score_PIA")
      res$score_PIA = as.numeric(res$score_PIA)
      benchmark_df = merge(benchmark_df, res, by.x = "proteins", by.y = "Isoform_PIA", all = T)
      
      # proteine presenti in input
      iso_input = get_score_from_idXML(glue("{PATH_TO_DATA}/Only{protease}/merge_index_percolator_pep_switched_0.01.idXML"))
      benchmark_df = merge(benchmark_df, iso_input, by.x = "proteins", by.y = "Isoform", all = T)
      
      # Eliminio isoforme non presenti nel validation set
      benchmark_df = benchmark_df[!is.na(benchmark_df$Present), ]
      
      # Tengo quelle in input
      benchmark_df = benchmark_df[!is.na(benchmark_df$score), ]
      
      # 0 se modello non trova l'isoforma
      benchmark_df[is.na(benchmark_df)] = 0
      benchmark_df$score = NULL
      
      colnames(benchmark_df) = gsub(".*Prob_present_", "", colnames(benchmark_df))
      colnames(benchmark_df) = gsub(".*score_", "", colnames(benchmark_df))
      
      for (nm in c(selected_models, "EPIFANY", "Fido", "PIA")) {
        benchmark_df[, paste0("Present_", nm)] = benchmark_df$Present
      }
      
      plot_tab = get_roc(benchmark_df,
                         c(selected_models, "EPIFANY", "Fido", "PIA"),
                         protease = glue(" - {protease}"))
      ggsave(glue("{PATH_RES_COMPETITORS}/{protease}/ROC_main_result.png"), plot = plot_tab$gplot)
      save(plot_tab, file = glue("{PATH_RES_COMPETITORS}/{protease}/ROC_main_result.rdata"))
      shared_vs_all_auc = plot_tab$sum_stat
      
      write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_COMPETITORS}/{protease}/SumTab_main_result.csv"),
                row.names = FALSE)
      benchmark_df_all = rbind(benchmark_df_all, benchmark_df)
      
      # Focus on validation without isoform with Unique Peptide (UP)
      plot_tab = get_roc(benchmark_df[benchmark_df$Y_unique_IsoBayes == 0, ],
                         c(selected_models, "EPIFANY", "Fido", "PIA"),
                         protease = glue(" - {protease}")
      )
      ggsave(glue("{PATH_RES_COMPETITORS}/{protease}/ROC_main_result_no_UP.png"), plot = plot_tab$gplot)
      save(plot_tab, file = glue("{PATH_RES_COMPETITORS}/{protease}/no_UP_ROC_main_result.rdata"))
      shared_vs_all_auc = cbind(shared_vs_all_auc, plot_tab$sum_stat$AUC)
      colnames(shared_vs_all_auc) = c("Model", "AUC_all", "AUC_only_shared")
      
      write.csv(shared_vs_all_auc,
                file = glue("{PATH_RES_COMPETITORS}/{protease}/SumTab_main_result_no_UP.csv"),
                row.names = FALSE)
    }
    
    for (noUP in c("", "no_UP_")) {
      if(noUP == "no_UP_"){
        sel = benchmark_df_all$Y_unique_IsoBayes == 0 # equal to Y_unique_IsoBayes_mRNA
      }else{
        sel = benchmark_df_all$Y_unique_IsoBayes > -Inf
      }
      plot_tab = get_roc(benchmark_df_all[sel, ], c(selected_models, "EPIFANY", "Fido", "PIA"))
      ggsave(glue("{PATH_RES_COMPETITORS}/{noUP}ROC_main_result.png"), plot = plot_tab$gplot)
      save(plot_tab, file = glue("{PATH_RES_COMPETITORS}/{noUP}ROC_main_result.rdata"))
      write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_COMPETITORS}/{noUP}SumTab_main_result.csv"), row.names = FALSE)
      
      ## abundance of main model (we consider the validation set used to benchmark all models)
      for (mrna in c("", "_mRNA")) {
        scat_bench = scatterplot(log10(benchmark_df_all[sel, c(glue("Abundance_IsoBayes{mrna}"), "Y_validation")] + 1))  + 
          labs(x = "Log10(Estimated abundance + 1)", y = "Log10(Validated Abundance + 1)") 
        
        ggsave(glue("{PATH_RES_COMPETITORS}/{noUP}scatterplot_benchmark{mrna}.png"), plot = scat_bench)
        save(scat_bench, file = glue("{PATH_RES_COMPETITORS}/{noUP}scatterplot_benchmark{mrna}.rdata"))
      }
      
      # change mrna prot
      for (mrna in c("", "_mRNA")) {
        sub_data = benchmark_df_all[sel, c(glue("Prob_prot_inc_IsoBayes{mrna}"),
                                           glue("TPM_IsoBayes{mrna}"),
                                           "tpm_validation", "P_Y_validation",
                                           glue("Log2_FC_IsoBayes{mrna}"),
                                           "Y_validation", glue("Pi_IsoBayes{mrna}"))
        ]
        sub_data = build_data_violin_plot(sub_data, 1.5e-06)
        sub_data = convert_numeric_to_class(sub_data, glue("Prob_prot_inc_IsoBayes{mrna}"),
                                            seq(0, 1, 0.2))
        sub_data_extreme = convert_numeric_to_class(sub_data, glue("Prob_prot_inc_IsoBayes{mrna}"),
                                                    quantiles = c(0, 0.01, 0.99, 1))
        sub_data_extreme = sub_data_extreme[sub_data_extreme$class_Prob_prot_inc != "(0.01, 0.99]", ]
        
        plot_change = plot_prob_change(sub_data) 
        ggsave(glue("{PATH_RES_COMPETITORS}/{noUP}change_mrna_prot{mrna}.pdf"), plot = plot_change)
        save(plot_change, file = glue("{PATH_RES_COMPETITORS}/{noUP}change_mrna_prot{mrna}.rdata"))
        
        plot_change = plot_prob_change(sub_data_extreme) 
        ggsave(glue("{PATH_RES_COMPETITORS}/{noUP}change_mrna_prot_extreme{mrna}.pdf"), plot = plot_change)
        save(plot_change, file = glue("{PATH_RES_COMPETITORS}/{noUP}change_mrna_prot_extreme{mrna}.rdata"))
        
        ths = 1.5e-06
        p_tpm_adj = (sub_data[, glue("TPM_IsoBayes{mrna}")]+ths)/sum(sub_data[, glue("TPM_IsoBayes{mrna}")]+ths)
        pi = (sub_data[, glue("Pi_IsoBayes{mrna}")]+ths)/sum(sub_data[, glue("Pi_IsoBayes{mrna}")]+ths)
        sub_data$Log2_FC_adj = log2(pi/p_tpm_adj)

        scat_bench = scatterplot(sub_data[, c(glue("Log2_FC_adj"), "Log2_FC_validation")])  + 
          labs(x = "Log2-FC", y = "Validated Log2-FC")

        ggsave(glue("{PATH_RES_COMPETITORS}/{noUP}scatterplot_log2fc{mrna}.png"), plot = scat_bench)
        save(scat_bench, file = glue("{PATH_RES_COMPETITORS}/{noUP}scatterplot_log2fc{mrna}.rdata"))
      }
      
      ### GENE CORRELATION
      # load validation dataset from metamorpheus
      load(glue("{PATH_TO_DATA}/No{protease}/Validation_gene_psm"))
      
      for (mrna in c("", "_mRNA")) {
        bench_gene = aggregate(benchmark_df_all[sel, glue("Abundance_IsoBayes{mrna}")],
                               by = list(benchmark_df_all[sel, glue("Gene_IsoBayes{mrna}")]),
                               FUN = sum)
        
        bench_gene = merge(VALIDATION_DF_prot, bench_gene,
                           by.x= "Gene", by.y = "Group.1", all = T)
        bench_gene = na.omit(bench_gene)
        
        scat_bench = scatterplot(log10(bench_gene[, c("x", "Y_validation")] + 1))  + 
          labs(x = "Log10(Gene Abundance + 1)", y = "Log10(Validated Gene Abundance + 1)")
        
        ggsave(glue("{PATH_RES_COMPETITORS}/{noUP}scatterplot_gene{mrna}.png"), plot = scat_bench)
        save(scat_bench, file = glue("{PATH_RES_COMPETITORS}/{noUP}scatterplot_gene{mrna}.rdata"))
      }
    }

  
  #################################################################################################
  # RUN_TIME and RAM
  #################################################################################################
  Data = rbind()
  
  for (ff in list.files(glue("{PATH_WD}/Benchmark_results/results"), pattern = "res_used")) {
    res_used = readLines(glue("{PATH_WD}/Benchmark_results/results/{ff}"))
    
    i_name = grep('job name', res_used)
    job_name = gsub(".*job name = ", "", res_used[i_name])
    job_name = gsub(", queue.*", "", job_name)
    
    i_res_used = grep('walltime', res_used)
    ramMB = round(as.numeric(gsub("kb resources_used.ncp.*", "",
                                  gsub(".*used.mem=", "",
                                       res_used[i_res_used]))) * 0.0009765625)
    
    wallTime = gsub(".*walltime=", "", res_used[i_res_used])
    wallTime = as.POSIXlt(wallTime,format="%H:%M:%S")
    wallTime = unclass(wallTime)
    wallTime = round(wallTime$min + wallTime$sec / 60, 1)
    
    info_model = strsplit(job_name, "_")[[1]]
    
    if(info_model[3] == "IsoBayes" & info_model[[length(info_model)]] == "mRNA"){
      info_model[3] = glue("{info_model[3]}_mRNA")
    }
    if(substr(info_model[3], 1, 8) == "IsoBayes" & info_model[5] == "TRUE"){
      info_model[3] = glue("{info_model[3]}_PEP")
    }
    if(substr(info_model[3], 1, 8) == "IsoBayes"){
      df = data.frame(Model = info_model[3], data = info_model[2], protease = info_model[7], RunTime = wallTime, RAM = ramMB)
    }else{
      df = data.frame(Model = info_model[3], data = info_model[2], protease = info_model[4], RunTime = wallTime, RAM = ramMB)
    }
    
    Data = rbind(Data, df)
  }
  Data$Model[grep("PiaTot", Data$Model)] = "PIA"
  Data$Model[grep("Epifany", Data$Model)] = "EPIFANY"
  
  data_protease_all = NULL
  data_protease_all_RAM = NULL
  
  for (protease in proteases) {
    # select model to plot
    data_protease = Data[Data$data == DATA & Data$protease == protease, ]
    data_protease = data_protease[!grepl("_PEP", data_protease$Model), ]
    data_protease_mean = aggregate.data.frame(data_protease[, c("RunTime", "RAM")], by = list(data_protease$Model), FUN = mean)
    data_protease_sd = aggregate.data.frame(data_protease[, c("RunTime", "RAM")], by = list(data_protease$Model), FUN = sd)
    colnames(data_protease_sd) = paste0(colnames(data_protease_sd), "_sd")
    data_protease = merge(data_protease_mean, data_protease_sd, by.x = "Group.1", by.y = "Group.1_sd")
    colnames(data_protease)[1] = "Model"
    
    # RUN TIME
    data_protease$RunTime = round(data_protease$RunTime, 1)
    id_sort = sort(data_protease$RunTime, decreasing = TRUE, index.return = TRUE)$ix
    data_protease$Model = factor(data_protease$Model, levels = data_protease$Model[id_sort])
    
    pp = run_time_plot(data_protease, title = glue("Runtime - {protease} - {DATA_name}"))
    ggsave(glue("{PATH_RES_COMPETITORS}/{protease}/Run-Time.png"), plot = pp)
    
    if(is.null(data_protease_all)){
      data_protease_all = data_protease
    }else{
      data_protease_all$RunTime = data_protease_all$RunTime + data_protease$RunTime
    }
    
    # Memory usage (RAM MB)
    data_protease$RAM = round(data_protease$RAM)
    id_sort = sort(data_protease$RAM, decreasing = TRUE, index.return = TRUE)$ix
    data_protease$Model = factor(data_protease$Model, levels = data_protease$Model[id_sort])
    
    pp = memory_plot(data_protease, title = glue("RAM - {protease} - {DATA_name}"))
    ggsave(glue("{PATH_RES_COMPETITORS}/{protease}/Memory_usage.png"), plot = pp)
    
    if(is.null(data_protease_all_RAM)){
      data_protease_all_RAM = data_protease
    }else{
      data_protease_all_RAM$RAM = data_protease_all_RAM$RAM + data_protease$RAM
    }
  }
  #############################################
  # average time
  #############################################
  data_protease_all$RunTime = round(data_protease_all$RunTime / length(proteases), 1)
  id_sort = sort(data_protease_all$RunTime, decreasing = TRUE, index.return = TRUE)$ix
  data_protease_all$Model = factor(data_protease_all$Model, levels = data_protease_all$Model[id_sort])
  
  pp = run_time_plot(data_protease_all, title = glue("{DATA_name}"))
  write.csv(data_protease_all_RAM, glue("{PATH_RES_COMPETITORS}/Average_Run-Time.csv"))
  ggsave(glue("{PATH_RES_COMPETITORS}/Average_Run-Time.png"), plot = pp)
  save(pp, file = glue("{PATH_RES_COMPETITORS}/Average_Run-Time.rdata"))
  
  #############################################
  # average RAM
  #############################################
  data_protease_all_RAM$RAM = round(data_protease_all_RAM$RAM / length(proteases) / 1000, 1)
  id_sort = sort(data_protease_all_RAM$RAM, decreasing = TRUE, index.return = TRUE)$ix
  data_protease_all_RAM$Model = factor(data_protease_all_RAM$Model, levels = data_protease_all_RAM$Model[id_sort])
  
  pp = memory_plot(data_protease_all_RAM, title = glue("{DATA_name}")) + ylab("Memory (GB)")
  write.csv(data_protease_all_RAM, glue("{PATH_RES_COMPETITORS}/Average_Memory_usage.csv"))
  ggsave(glue("{PATH_RES_COMPETITORS}/Average_Memory_usage.png"), plot = pp)
  save(pp, file = glue("{PATH_RES_COMPETITORS}/Average_Memory_usage.rdata"))
  
}

main(proteases = list.dirs(glue("{PATH_WD}/Benchmark_results/{DATA}"), recursive = FALSE, full.names = FALSE),
     models = list(IsoBayes = c("", ""),
                   IsoBayes_PEP = c("_PEP", ""),
                   IsoBayes_mRNA = c("", "_mRNA"),
                   IsoBayes_mRNA_PEP = c("_PEP", "_mRNA")
     )
)