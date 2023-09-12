PATH_WD = commandArgs(trailingOnly = TRUE)[1]
DATA = commandArgs(trailingOnly = TRUE)[2]

# Set path, global variable and libraries
###########################################################################################
library(ggplot2)
library(readr)
library(glue)
source(glue("{PATH_WD}/utils_function/utils_benchmarking.R"))
source(glue("{PATH_WD}/utils_function/get_roc.R"))
source(glue("{PATH_WD}/utils_function/prior_plot.R"))
source(glue("{PATH_WD}/utils_function/log_output.R"))

PATH_RES_COMPETITORS = glue("{PATH_WD}/Benchmark_results/{DATA}")
PATH_TO_DATA = glue("{PATH_WD}/Data/{DATA}")
PATH_TO_RES = glue("{PATH_WD}/Model_results/{DATA}/")
load(glue("{PATH_WD}/utils_function/PALETTE_MODELS"))

log_output(glue("benchmarking_results_{DATA}"))

main = function(models, proteases){
  # COMPETITORS vs IsoBayes_openMS vs IsoBayes_openMS_mRNA [PERFORMANCE]
  benchmark_df_all = list()
  selected_models = c("IsoBayes", "IsoBayes_mRNA")
  
  for(protease in proteases){
    # load validation dataset from metamorpheus
    load(glue("{PATH_TO_DATA}/No{protease}/Validation_prot_psm"))
    
    benchmark_df = list(VALIDATION_DF_prot[, c("proteins", "Present")])
    for (model in selected_models) {
      attribute_model = glue("{models[[model]][2]}{models[[model]][1]}") 
      load(glue("{PATH_TO_RES}/OpenMS{attribute_model}/{protease}/OpenMS{attribute_model}_MCMC.RData"))
      load(glue("{PATH_TO_RES}/OpenMS{attribute_model}/{protease}/OpenMS{attribute_model}_data_loaded.RData"))
      res$isoform_results = merge(res$isoform_results, data_loaded$PROTEIN_DF,
                                  by.x = "Isoform", by.y = "protein_name")
      benchmark_df = append(benchmark_df, list(res$isoform_results[, c("Isoform", "Prob_present", "Y_unique")]))
    }
    
    for (i in seq_len(length(benchmark_df)-1)) {
      colnames(benchmark_df[[2]]) = paste0(colnames(benchmark_df[[2]]), "_", selected_models[i])
      benchmark_df[[1]] = merge(benchmark_df[[1]], benchmark_df[[2]],
                                by.x = "proteins", by.y = paste0("Isoform_", selected_models[i]), all = T)
      benchmark_df[[2]] = NULL
    }
    # EPIFANY
    res = get_score_from_idXML(paste0(PATH_RES_COMPETITORS, "/",  protease, "/epifany.idXML"))
    colnames(res) = paste0(colnames(res), "_Epifany")
    benchmark_df = merge(benchmark_df, res, by.x = "proteins", by.y = "Isoform_Epifany", all = T)
    
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
    iso_input = get_score_from_idXML(glue("{PATH_TO_DATA}/Only{protease}/merge_index_percolator_pep_switched.idXML"))
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
    
    benchmark_df$IsoBayes_mRNA_no_unique = benchmark_df$IsoBayes_mRNA
    benchmark_df$IsoBayes_mRNA_no_unique[benchmark_df$Y_unique_IsoBayes != 0] = NA
    
    benchmark_df$IsoBayes_no_unique = benchmark_df$IsoBayes
    benchmark_df$IsoBayes_no_unique[benchmark_df$Y_unique_IsoBayes != 0] = NA
    
    plot_tab = get_roc(benchmark_df, c(selected_models, "Epifany", "Fido", "PIA", "IsoBayes_mRNA_no_unique", "IsoBayes_no_unique"))
    ggsave(glue("{PATH_RES_COMPETITORS}/{protease}/ROC_main_result.png"), plot = plot_tab$gplot)
    write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_COMPETITORS}/{protease}/SumTab_main_result.csv"), row.names = FALSE)
    benchmark_df_all = rbind(benchmark_df_all, benchmark_df)
    
    # Focus on validation without isoform with Unique Peptide (UP)
    plot_tab = get_roc(benchmark_df[benchmark_df$Y_unique_IsoBayes == 0, ], c(selected_models, "Epifany", "Fido", "PIA"))
    ggsave(glue("{PATH_RES_COMPETITORS}/{protease}/ROC_main_result_no_UP.png"), plot = plot_tab$gplot)
    shared_vs_all_auc = plot_tab$sum_stat
    
    plot_tab = get_roc(benchmark_df, c(selected_models, "Epifany", "Fido", "PIA"))
    shared_vs_all_auc = cbind(shared_vs_all_auc, plot_tab$sum_stat$AUC)
    colnames(shared_vs_all_auc) = c("Model", "AUC_only_shared", "AUC_all")
    write.csv(shared_vs_all_auc, file = glue("{PATH_RES_COMPETITORS}/{protease}/SumTab_main_result_no_UP.csv"), row.names = FALSE)
  }
  plot_tab = get_roc(benchmark_df_all, c(selected_models, "Epifany", "Fido", "PIA"))
  ggsave(glue("{PATH_RES_COMPETITORS}/ROC_main_result.png"), plot = plot_tab$gplot)
  write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_COMPETITORS}/SumTab_main_result.csv"), row.names = FALSE)
  
  # Focus on validation without isoform with Unique Peptide (UP)
  plot_tab = get_roc(benchmark_df_all[benchmark_df_all$Y_unique_IsoBayes == 0, ], c(selected_models, "Epifany", "Fido", "PIA"))
  ggsave(glue("{PATH_RES_COMPETITORS}/ROC_main_result_no_UP.png"), plot = plot_tab$gplot)
  write.csv(plot_tab$sum_stat, file = glue("{PATH_RES_COMPETITORS}/SumTab_main_result_no_UP.csv"), row.names = FALSE)
  
  #################################################################################################
  # RUN_TIME and RAM
  #################################################################################################
  Data = rbind()
  
  for (ff in list.files(glue("{PATH_WD}/Benchmark_results"), pattern = "res_used")) {
    res_used = readLines(glue("{PATH_WD}/Benchmark_results/{ff}"))
    
    i_name = grep('job name', res_used)
    job_name = gsub(".*job name = ", "", res_used[i_name])
    job_name = gsub(", queue.*", "", job_name)
    
    i_res_used = grep('walltime', res_used)
    ramMB = round(as.numeric(gsub("kb resources_used.ncp.*", "", gsub(".*used.mem=", "", res_used[i_res_used]))) * 0.0009765625)
    
    wallTime = gsub(".*walltime=", "", res_used[i_res_used])
    wallTime = as.POSIXlt(wallTime,format="%H:%M:%S")
    wallTime = unclass(wallTime)
    wallTime = round(wallTime$min + wallTime$sec / 60, 2)
    
    info_model = strsplit(job_name, "_")[[1]]
    if(info_model[3] == "IsoBayes" & info_model[5] == "FALSE"){
      info_model[3] = glue("{info_model[3]}_fast")
    }
    df = data.frame(Model = info_model[3], data = info_model[2], protease = rev(info_model)[1], RunTime = wallTime, RAM = ramMB)
    Data = rbind(Data, df)
  }
  Data$Model[grep("PiaTot", Data$Model)] = "PIA"
  Data$Model[Data$Model == "IsoBayes"] = "IsoBayes_mRNA"
  Data$Model[Data$Model == "IsoBayes_fast"] = "IsoBayes_fast_mRNA"
  
  data_protease_all = NULL
  data_protease_all_RAM = NULL
  
  for (protease in proteases) {
    # select model to plot
    data_protease = Data[Data$data == DATA & Data$protease == protease, ]
    
    # RUN TIME
    id_sort = sort(data_protease$RunTime, decreasing = TRUE, index.return = TRUE)$ix
    data_protease$Model = factor(data_protease$Model, levels = data_protease$Model[id_sort])
    
    ggplot(data_protease, aes(x = Model, y = RunTime, label = RunTime, fill = Model)) +
      geom_bar(stat = "identity") +
      scale_fill_manual("Model", values = PALETTE_MODELS) +
      geom_text(size = 4.25, colour = "black",  vjust=-0.25) +
      labs(y = "Run-Time (Min)",
           x = "Tool",
           title = glue("Run-Time - {protease} - {DATA}")) +
      theme_classic() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 11, face = "bold"),
            legend.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(size = 10, face = "bold", angle = 10, margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
            axis.text.y = element_text(size = 10, face = "bold"),
            legend.position = "none") 
    ggsave(glue("{PATH_RES_COMPETITORS}/{protease}/Run-Time.png"))
    
    if(is.null(data_protease_all)){
      data_protease_all = data_protease
    }else{
      data_protease_all$RunTime = data_protease_all$RunTime + data_protease$RunTime
    }
    
    # Memory usage (RAM MB)
    id_sort = sort(data_protease$RAM, decreasing = TRUE, index.return = TRUE)$ix
    data_protease$Model = factor(data_protease$Model, levels = data_protease$Model[id_sort])
    
    ggplot(data_protease, aes(x = Model, y = RAM, label = RAM, fill = Model)) +
      geom_bar(stat = "identity") +
      scale_fill_manual("Model", values = PALETTE_MODELS) +
      geom_text(size = 4.25, colour = "black",  vjust=-0.25) +
      labs(y = "RAM (MB)",
           x = "Tool",
           title = glue("RAM - {protease} - {DATA}")) +
      theme_classic() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 11, face = "bold"),
            legend.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(size = 10, face = "bold", angle = 10, margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
            axis.text.y = element_text(size = 10, face = "bold"),
            legend.position = "none")

    ggsave(glue("{PATH_RES_COMPETITORS}/{protease}/Memory_usage.png"))
    
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
  
  ggplot(data_protease_all, aes(x = Model, y = RunTime, label = RunTime, fill = Model)) +
    geom_bar(stat = "identity") +
    scale_fill_manual("Model", values = PALETTE_MODELS) +
    geom_text(size = 4.25, colour = "black",  vjust=-0.25) +
    labs(y = "Run-Time (Min)",
         x = "Tool",
         title = glue("Run-Time - {DATA}")) +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold", angle = 10, margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
          axis.text.y = element_text(size = 10, face = "bold"),
          legend.position = "none") 
  ggsave(glue("{PATH_RES_COMPETITORS}/Average_Run-Time.png"))
  
  #############################################
  # average RAM
  #############################################
  
  data_protease_all_RAM$RAM = round(data_protease_all_RAM$RAM/length(proteases)/2)
  id_sort = sort(data_protease_all_RAM$RAM, decreasing = TRUE, index.return = TRUE)$ix
  data_protease_all_RAM$Model = factor(data_protease_all_RAM$Model, levels = data_protease_all_RAM$Model[id_sort])
  
  ggplot(data_protease_all_RAM, aes(x = Model, y = RAM, label = RAM, fill = Model)) +
    geom_bar(stat = "identity") +
    scale_fill_manual("Model", values = PALETTE_MODELS) +
    geom_text(size = 4.25, colour = "black",  vjust=-0.25) +
    labs(y = "RAM (MB)",
         x = "Tool",
         title = glue("RAM - {DATA}")) +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold", angle = 10, margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
          axis.text.y = element_text(size = 10, face = "bold"),
          legend.position = "none") 
  ggsave(glue("{PATH_RES_COMPETITORS}/Average_Memory_usage.png"))
}

main(proteases = list.dirs(glue("{PATH_WD}/Benchmark_results/{DATA}"), recursive = FALSE, full.names = FALSE),
     models = list(IsoBayes = c("_PEP", ""),
                   IsoBayes_fast = c("", ""),
                   IsoBayes_mRNA = c("_PEP", "_mRNA"),
                   IsoBayes_fast_mRNA = c("", "_mRNA")
                   )
)

