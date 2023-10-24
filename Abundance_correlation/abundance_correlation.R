PATH_WD = commandArgs(trailingOnly = TRUE)[1]
DATA = commandArgs(trailingOnly = TRUE)[2]

# Set path, global variable and libraries
###########################################################################################
library(glue)
library(ggplot2)
source(glue("{PATH_WD}/Benchmark_results/utils_benchmarking.R"))
source(glue("{PATH_WD}/utils_function/scatterplot.R"))
source(glue("{PATH_WD}/utils_function/build_OpenMS_validation.R"))
source(glue("{PATH_WD}/utils_function/prior_plot.R"))
source(glue("{PATH_WD}/utils_function/log_output.R"))

PATH_RES_COMPETITORS = glue("{PATH_WD}/Benchmark_results/{DATA}")
PATH_TO_DATA = glue("{PATH_WD}/Data/{DATA}")
PATH_RES_MODEL = glue("{PATH_WD}/Model_results/{DATA}")
PATH_RES_ABUNDANCE = glue("{PATH_WD}/Abundance_correlation/{DATA}" )
log_output(glue("abundance_correlation_{DATA}"))

main = function(proteases){
  ######################################################
  # Scatterplot: log10 abundance vs log10 abundance validated - OpenMS, MM_psm, MM_intensity
  ######################################################
  inputs = c("OpenMS", "MM_psm", "MM_intensities")
  
  for (input in inputs) {
    for (pep in c("", "_PEP")){
      for (mrna in c("", "_mRNA")) {
        benchmark_df_all = list()
        SumTab_corr = rbind() 
        SumTab_corr_no_unique = rbind()
        for (protease in proteases) {
          # load validation dataset merged with model results
          load(glue("{PATH_WD}/Model_results/{DATA}/{input}{mrna}{pep}/{protease}/Merged_validation_res_{input}{mrna}{pep}.RData"))
          benchmark_df = validation_dat ; rm(validation_dat)
          
          benchmark_df_all = rbind(benchmark_df_all, benchmark_df)
          abundances = data.frame(log10_abundance = log10(benchmark_df$Abundance + 1),
                                  log10_abundance_validated = log10(benchmark_df$Y_validation + 1)
          )
          plot_corr = scatterplot(abundances) + 
            labs(x = "Correlation Log10(Abundance)", y = "Correlation Log10(Validated Abundance)")
          ggsave(glue("{PATH_RES_ABUNDANCE}/{protease}/scatterplot_{input}{mrna}.png"), plot = plot_corr)
          SumTab_corr = rbind(SumTab_corr, c(protease, cor(abundances$log10_abundance, abundances$log10_abundance_validated)))
          
          abundances = abundances[benchmark_df$Y_unique == 0, ]
          plot_corr = scatterplot(abundances) + 
            labs(x = "Correlation Log10(Abundance)", y = "Correlation Log10(Validated Abundance)")
          ggsave(glue("{PATH_RES_ABUNDANCE}/{protease}/scatterplot_no_unique_{input}{mrna}{pep}.png"), plot = plot_corr)
          SumTab_corr_no_unique = rbind(SumTab_corr_no_unique,
                                        c(protease, cor(abundances$log10_abundance,
                                                        abundances$log10_abundance_validated)))
        }
        abundances = data.frame(log10_abundance = log10(benchmark_df_all$Abundance + 1),
                                log10_abundance_validated = log10(benchmark_df_all$Y_validation + 1)
        )
        plot_corr = scatterplot(abundances)+ 
          labs(x = "Correlation Log10(Abundance)", y = "Correlation Log10(Validated Abundance)")
        ggsave(glue("{PATH_RES_ABUNDANCE}/scatterplot_{input}{mrna}{pep}.png"), plot = plot_corr)
        
        SumTab_corr = rbind(SumTab_corr, c("All", cor(abundances$log10_abundance, abundances$log10_abundance_validated)))
        SumTab_corr = as.data.frame(SumTab_corr)
        colnames(SumTab_corr) = c("Protease", "Log10_correlation")
        SumTab_corr$Log10_correlation = round(as.numeric(SumTab_corr$Log10_correlation), 2)
        
        write.csv(SumTab_corr, file = glue("{PATH_RES_ABUNDANCE}/correlation_{input}{mrna}{pep}.csv"), row.names = FALSE)
        
        # focus on shared 
        abundances = abundances[benchmark_df_all$Y_unique == 0, ]
        plot_corr = scatterplot(abundances)+ 
          labs(x = "Correlation Log10(Abundance)", y = "Correlation Log10(Validated Abundance)")
        ggsave(glue("{PATH_RES_ABUNDANCE}/scatterplot_no_unique_{input}{mrna}{pep}.png"), plot = plot_corr)
        
        SumTab_corr_no_unique = rbind(SumTab_corr_no_unique, c("All", cor(abundances$log10_abundance, abundances$log10_abundance_validated)))
        SumTab_corr_no_unique = as.data.frame(SumTab_corr_no_unique)
        colnames(SumTab_corr_no_unique) = c("Protease", "Log10_correlation")
        SumTab_corr_no_unique$Log10_correlation = round(as.numeric(SumTab_corr_no_unique$Log10_correlation), 2)
        
        write.csv(SumTab_corr_no_unique, file = glue("{PATH_RES_ABUNDANCE}/correlation_no_unique_{input}{mrna}{pep}.csv"), row.names = FALSE)
      } 
    }
  }
  
  ######################################################
  # Log10 abundance vs log10 abundance validated - PEP vs no_PEP (MM with mRNA)
  ######################################################
  SumTab_corr_all = data.frame(matrix(nrow = length(proteases)+1, ncol = 0))
  for (model in c("_PEP", "")) {
    benchmark_df_all = list()
    SumTab_corr = rbind()
    for (protease in proteases) {
      # load validation dataset merged with model results
      load(glue("{PATH_RES_MODEL}/MM_psm_mRNA{model}/{protease}/Merged_validation_res_MM_psm_mRNA{model}.RData"))
      benchmark_df = validation_dat ; rm(validation_dat)
      benchmark_df_all = rbind(benchmark_df_all, benchmark_df)
      
      SumTab_corr = rbind(SumTab_corr, c(protease, cor(log10(benchmark_df$Abundance + 1),
                                                       log10(benchmark_df$Y_validation + 1)))
      )
    }
    SumTab_corr = rbind(SumTab_corr, c("All", cor(log10(benchmark_df_all$Abundance + 1),
                                                  log10(benchmark_df_all$Y_validation + 1)))
    )
    SumTab_corr = as.data.frame(SumTab_corr)
    colnames(SumTab_corr) = c("Protease", glue("Log10_correlation{model}"))
    SumTab_corr[, 2] = round(as.numeric(SumTab_corr$Log10_correlation), 3)
    SumTab_corr_all = cbind(SumTab_corr_all, SumTab_corr)
  }
  write.csv(SumTab_corr_all, file = glue("{PATH_RES_ABUNDANCE}/correlation_PEP_noPEP.csv"), row.names = FALSE)
  
  ######################################################
  # Log10 abundance vs log10 abundance validated - psm vs intensities (MM)
  ######################################################
  SumTab_corr_all = data.frame(matrix(nrow = length(proteases)+1, ncol = 0))
  for (input in c("MM_psm", "MM_intensities")){
    tab = read.csv(glue("{PATH_RES_ABUNDANCE}/correlation_{input}.csv"))
    colnames(tab) = c("Protease", glue("Log10_correlation_{input}"))
    SumTab_corr_all = cbind(SumTab_corr_all, tab)
  }
  write.csv(SumTab_corr_all, file = glue("{PATH_RES_ABUNDANCE}/correlation_psm_intensities.csv"), row.names = FALSE)
  
  ######################################################
  # Log10 abundance vs log10 abundance validated - MM vs OpenMS
  ######################################################
  SumTab_corr_all = data.frame(matrix(nrow = length(proteases)+1, ncol = 0))
  for (input in c("MM_psm", "OpenMS")){
    tab = read.csv(glue("{PATH_RES_ABUNDANCE}/correlation_{input}.csv"))
    colnames(tab) = c("Protease", glue("Log10_correlation_{input}"))
    SumTab_corr_all = cbind(SumTab_corr_all, tab)
  }
  write.csv(SumTab_corr_all, file = glue("{PATH_RES_ABUNDANCE}/correlation_MM_OpenMS.csv"), row.names = FALSE)
  
  ######################################################
  # Log10 abundance vs log10 abundance validated - mRNA vs no_mRNA (MM)
  ######################################################
  for (input in inputs) {
    SumTab_corr_all = data.frame(matrix(nrow = length(proteases)+1, ncol = 0))
    SumTab_corr_all_no_unique = data.frame(matrix(nrow = length(proteases)+1, ncol = 0))
    for (model in c("", "_mRNA")) {
      benchmark_df_all = list()
      SumTab_corr = rbind()
      SumTab_corr_no_unique = rbind()
      
      for (protease in proteases) {
        # load validation dataset merged with model results
        load(glue("{PATH_RES_MODEL}/{input}{model}_PEP/{protease}/Merged_validation_res_{input}{model}_PEP.RData"))
        benchmark_df = validation_dat ; rm(validation_dat)
        benchmark_df_no_unique = benchmark_df[benchmark_df$Y_unique == 0, ]
        benchmark_df_all = rbind(benchmark_df_all, benchmark_df)
        
        SumTab_corr = rbind(SumTab_corr, c(protease, cor(log10(benchmark_df$Abundance + 1),
                                                         log10(benchmark_df$Y_validation + 1)))
        )
        
        SumTab_corr_no_unique = rbind(SumTab_corr_no_unique, c(protease, cor(log10(benchmark_df_no_unique$Abundance + 1),
                                                                             log10(benchmark_df_no_unique$Y_validation + 1)))
        )
      }
      SumTab_corr = rbind(SumTab_corr, c("All", cor(log10(benchmark_df_all$Abundance + 1),
                                                    log10(benchmark_df_all$Y_validation + 1)))
      )
      SumTab_corr = as.data.frame(SumTab_corr)
      colnames(SumTab_corr) = c("Protease", glue("Log10_correlation{model}"))
      SumTab_corr[, 2] = round(as.numeric(SumTab_corr$Log10_correlation), 3)
      SumTab_corr_all = cbind(SumTab_corr_all, SumTab_corr)
      
      # focus on shared
      benchmark_df_all = benchmark_df_all[benchmark_df_all$Y_unique == 0, ]
      SumTab_corr_no_unique = rbind(SumTab_corr_no_unique, c("All", cor(log10(benchmark_df_all$Abundance + 1),
                                                                        log10(benchmark_df_all$Y_validation + 1)))
      )
      SumTab_corr_no_unique = as.data.frame(SumTab_corr_no_unique)
      colnames(SumTab_corr_no_unique) = c("Protease", glue("Log10_correlation{model}"))
      SumTab_corr_no_unique[, 2] = round(as.numeric(SumTab_corr_no_unique$Log10_correlation), 3)
      
      SumTab_corr_all_no_unique = cbind(SumTab_corr_all_no_unique, SumTab_corr_no_unique)
    }
    write.csv(SumTab_corr_all, file = glue("{PATH_RES_ABUNDANCE}/correlation_mRNA_no_mRNA_{input}.csv"), row.names = FALSE)
    write.csv(SumTab_corr_no_unique, file = glue("{PATH_RES_ABUNDANCE}/correlation_mRNA_no_mRNA_no_unique_{input}.csv"), row.names = FALSE)
  }
  
  ######################################################
  # Log10 abundance vs log10 abundance validated - prior grid (OpenMS with mRNA)
  ######################################################
  load(glue("{PATH_WD}/Model_results/prior_grid"))
  SumTab_corr_all = data.frame(matrix(nrow = length(proteases)+1, ncol = 0))
  
  for (prior in prior_grid) {
    benchmark_df_all = list()
    SumTab_corr = rbind()
    for(protease in proteases){
      # load validation dataset from metamorpheus
      load(glue("{PATH_TO_DATA}/No{protease}/Validation_prot_psm"))
      load(glue("{PATH_RES_MODEL}/OpenMS_mRNA_PEP_prior_{prior}/{protease}/OpenMS_mRNA_PEP_prior_{prior}_MCMC.RData"))
      
      benchmark_df = VALIDATION_DF_prot[, c("proteins", "Y_validation", "Present")]
      benchmark_df = merge(benchmark_df, res$isoform_results[, c("Isoform", "Abundance")],
                           by.x = "proteins", by.y = "Isoform", all = T)
      
      colnames(benchmark_df)[1] = "Isoform"
      benchmark_df = build_OpenMS_validation(benchmark_df, protease)
      benchmark_df_all = rbind(benchmark_df_all, benchmark_df)
      
      SumTab_corr = rbind(SumTab_corr, c(protease, cor(log10(benchmark_df$Abundance + 1),
                                                       log10(benchmark_df$Y_validation + 1)))
      )
    }
    SumTab_corr = rbind(SumTab_corr, c("All", cor(log10(benchmark_df_all$Abundance + 1),
                                                  log10(benchmark_df_all$Y_validation + 1)))
    )
    SumTab_corr = as.data.frame(SumTab_corr)
    colnames(SumTab_corr) = c("Protease", glue("Log10_correlation_prior_{prior}"))
    SumTab_corr[, 2] = round(as.numeric(SumTab_corr$Log10_correlation), 3)
    SumTab_corr_all = cbind(SumTab_corr_all, SumTab_corr)
  }
  write.csv(SumTab_corr_all, file = glue("{PATH_RES_ABUNDANCE}/correlation_prior.csv"), row.names = FALSE)  
  
  data_to_plot = SumTab_corr_all
  data_to_plot = data_to_plot[length(proteases)+1, grep("Log", colnames(data_to_plot))]
  
  data_to_plot = data.frame(grid = as.character(prior_grid),
                            AUC = as.numeric(data_to_plot[1, ]),
                            model = rep("a", length(prior_grid))
  )
  pp = prior_plot(data_to_plot) +
    ylab("Correlation - Log10(Abundance)")
  
  ggsave(glue("{PATH_RES_ABUNDANCE}/correlation_prior_plot.png"), plot = pp)
  
  ######################################################
  # Gene correlation
  ######################################################
  for (input in c("MM_psm", "OpenMS")){
    SumTab_corr_all = data.frame(matrix(nrow = length(proteases)+1, ncol = 0))
    benchmark_df_all = list()
    SumTab_corr = rbind()
    for (protease in proteases) {
      # load validation dataset from metamorpheus
      load(glue("{PATH_TO_DATA}/No{protease}/Validation_gene_psm"))
      # load model results
      load(glue("{PATH_RES_MODEL}/{input}_mRNA_PEP/{protease}/{input}_mRNA_PEP_MCMC.RData"))
      
      benchmark_df = merge(VALIDATION_DF_prot, res$gene_abundance[, c("Gene", "Abundance")],
                           by= "Gene", all = T)
      
      # Eliminio isoforme non presenti nel validation set e in input
      benchmark_df = na.omit(benchmark_df)
      benchmark_df_all = rbind(benchmark_df_all, benchmark_df)
      
      SumTab_corr = rbind(SumTab_corr, c(protease, cor(log10(benchmark_df$Abundance + 1),
                                                       log10(benchmark_df$Y_validation + 1)))
      )
    }
    SumTab_corr = rbind(SumTab_corr, c("All", cor(log10(benchmark_df_all$Abundance + 1),
                                                  log10(benchmark_df_all$Y_validation + 1)))
    )
    SumTab_corr = as.data.frame(SumTab_corr)
    colnames(SumTab_corr) = c("Protease", glue("Log10_correlation_gene"))
    SumTab_corr[, 2] = round(as.numeric(SumTab_corr$Log10_correlation), 3)
    SumTab_corr_all = cbind(SumTab_corr_all, SumTab_corr)
    
    write.csv(SumTab_corr_all, file = glue("{PATH_RES_ABUNDANCE}/{input}_correlation_gene.csv"), row.names = FALSE)
  }
}

main(proteases = list.dirs(PATH_RES_ABUNDANCE, recursive = FALSE, full.names = FALSE))