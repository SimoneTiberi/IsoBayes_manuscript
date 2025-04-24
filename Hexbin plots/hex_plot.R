library(glue) ; library(ggplot2)

# TODO:
# MM intensities data

# add gene-level plot!!!

setwd("~/Desktop/IsoBayes_manuscript-main")

source(glue("./utils_function/utils_benchmarking.R"))
source(glue("./utils_function/get_roc.R"))
source(glue("./utils_function/memory_plot.R"))
source(glue("./utils_function/run_time_plot.R"))
source(glue("./utils_function/log_output.R"))
source(glue("./utils_function/scatterplot.R"))
source(glue("./utils_function/utils_change_iso_mrna.R"))
source(glue("./utils_function/plot_prob_change.R"))

hex = function(abundances){
  colnames(abundances) = c("x", "y")
  corr = cor(abundances)
  abundances = abundances[ abundances$x > 0.31,] # remove 0 and 1
  abundances = abundances[ abundances$y > 0.31,] # remove 0 and 1
  pp =  ggplot(abundances, aes(x = x, y = y)) +
    geom_hex(bins = 15) +
    labs(title = DATA_name, subtitle = glue("Correlation: {round(corr[1, 2], 2)}")) +
    theme_bw() +
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.title = element_text(size = 11),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11))
  
  pp
}

hist_plot = function(df_sub){
  colnames(df_sub) = c("x")
  ggplot(df_sub, aes(x = x)) +
    geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
    labs(title = DATA_name,
         x = "Validated Log2-FC") +
    theme_bw() +
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          axis.title = element_text(size = 15),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          legend.text = element_text(size = 11))
}

CORR = matrix(NA, 2, 2)
colnames(CORR) = c("IsoBayes", "IsoBayes_mRNA")
rownames(CORR) = c("jurkat", "WTC-11")
for (DATA in c("jurkat", "wtc11")){
  benchmark_df_all = read.csv(glue("Hexbin plots/df_{DATA}.csv"))
  DATA_name = DATA
  PATH_TO_DATA = glue("./Data/{DATA}")
  PATH_RES_COMPETITORS = glue("./new_res/{DATA}")
  
  for (noUP in c("", "no_UP_", "gene_w_more_ISO_")) { 
    if(noUP == "no_UP_"){
      sel = benchmark_df_all$Y_unique_IsoBayes == 0 # equal to Y_unique_IsoBayes_mRNA
      
      i = ifelse(DATA == "jurkat", 1, 2)
      CORR[i , 1] = print(cor( log10(1+benchmark_df_all$tpm_validation), log10(1+benchmark_df_all$Abundance_IsoBayes)))
      CORR[i , 2] = print(cor( log10(1+benchmark_df_all$tpm_validation), log10(1+benchmark_df_all$Abundance_IsoBayes_mRNA)))
      
      for (mrna in c("", "_mRNA")) {
        scat_bench = hex(log10(benchmark_df_all[sel, c(glue("Abundance_IsoBayes{mrna}"), "tpm_validation")] + 1))  + 
          labs(x = "Log10(Estimated abundance + 1)", y = "Log10(TMP + 1)") 
        
        ggsave(glue("{PATH_RES_COMPETITORS}/{noUP}scatterplot_TPM{mrna}.png"), plot = scat_bench)
        save(scat_bench, file = glue("{PATH_RES_COMPETITORS}/{noUP}scatterplot_benchmark{mrna}.rdata"))
      }
    }else if (noUP == ""){
      sel = benchmark_df_all$Y_unique_IsoBayes > -Inf
    }else{
      sel = benchmark_df_all$count_iso_gene > 1
    }
    
    ## abundance of main model (we consider the validation set used to benchmark all models)
    for (mrna in c("", "_mRNA")) {
      scat_bench = hex(log10(benchmark_df_all[sel, c(glue("Abundance_IsoBayes{mrna}"), "Y_validation")] + 1))  + 
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
      
      sel_sub = sub_data[, glue("TPM_IsoBayes{mrna}")] == 0 | sub_data[, glue("Pi_IsoBayes{mrna}")] == 0
      sub_data = sub_data[sel_sub,]
      ths = 1.5e-06
      p_tpm_adj = (sub_data[, glue("TPM_IsoBayes{mrna}")]+ths)/sum(sub_data[, glue("TPM_IsoBayes{mrna}")]+ths)
      pi = (sub_data[, glue("Pi_IsoBayes{mrna}")]+ths)/sum(sub_data[, glue("Pi_IsoBayes{mrna}")]+ths)
      sub_data$Log2_FC_adj = log2(pi/p_tpm_adj)
      
      plot_change = hist_plot(data.frame(sub_data[,"Log2_FC_adj"])) 
      ggsave(glue("{PATH_RES_COMPETITORS}/{noUP}hist{mrna}.pdf"), plot = plot_change)
      save(plot_change, file = glue("{PATH_RES_COMPETITORS}/{noUP}hist{mrna}.rdata"))
    }
    
    for (mrna in c("", "_mRNA")) {
      ### GENE CORRELATION
      # load validation dataset from metamorpheus
      bench_gene = aggregate(benchmark_df_all[sel, glue("Abundance_IsoBayes{mrna}")],
                             by = list(benchmark_df_all[sel, glue("Gene_IsoBayes{mrna}")]),
                             FUN = sum)
      colnames(bench_gene) = c("Gene", "x")
      
      VALIDATION_gene = read.csv(glue("Hexbin plots/gene_{DATA}.csv"))[,c(2,3)]
      colnames(VALIDATION_gene) = c("Gene", "Y_validation")
      
      bench_gene = merge(VALIDATION_gene, bench_gene,
                         by.x= "Gene")
      bench_gene = na.omit(bench_gene)
      bench_gene = bench_gene[bench_gene$Gene != "0",]
      
      scat_bench = hex(log10(bench_gene[, c("x", "Y_validation")] + 1))  +
        labs(x = "Log10(Gene Abundance + 1)", y = "Log10(Validated Gene Abundance + 1)")
      
      ggsave(glue("{PATH_RES_COMPETITORS}/{noUP}scatterplot_gene{mrna}.png"), plot = scat_bench)
      save(scat_bench, file = glue("{PATH_RES_COMPETITORS}/{noUP}scatterplot_gene{mrna}.rdata"))
    }
  }
}

if(FALSE){
  round(CORR, 2)
  library(xtable)
  xtable(t(CORR))
  if(FALSE){
    \begin{table}[ht]
    \centering
    \begin{tabular}{rrr}
    \hline
    & jurkat & WTC-11 \\ 
    \hline
    IsoBayes & 0.52 & 0.45 \\ 
    IsoBayes\_mRNA & 0.65 & 0.66 \\ 
    \hline
    \end{tabular}
    \caption{Correlation between log10 mRNA and estimated protein isoform abundances (i.e., log10(abundance + 1)). In each
      cell line, we considered results from all protesease}
    \end{table}
  }
}
############ MERGE PLOT ########################################
library(gridExtra)
library(grid)

size = 10 ; scale = 0.5
source(glue("./utils_function/grid_arrange_shared_legend.R"))


for (noUP in c("", "no_UP_", "gene_w_more_ISO_")) {
  list_plot_bench_abundance = list()
  list_plot_bench_change = list()
  list_plot_bench_time = list()
  list_plot_bench_memory = list()
  
  for (dat in c("jurkat", "wtc11")) {
    
    for (mrna in c("", "_mRNA")) {
      load(glue("new_res/{dat}/{noUP}scatterplot_benchmark{mrna}.rdata"))
      list_plot_bench_abundance = append(list_plot_bench_abundance, list(scat_bench))
      
      load(glue("new_res/{dat}/{noUP}scatterplot_gene{mrna}.rdata"))
      list_plot_bench_abundance = append(list_plot_bench_abundance, list(scat_bench))
      
      load(glue("new_res/{dat}/{noUP}hist{mrna}.rdata"))
      list_plot_bench_change = append(list_plot_bench_change, list(plot_change))
      
    }
  }
  
  a = grid.arrange(list_plot_bench_abundance[[1]], list_plot_bench_abundance[[5]]+ ylab(NULL), nrow =1)
  ggsave(glue("new_res/{noUP}scatterplot_abundance_benchmark.png"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_abundance[[2]], list_plot_bench_abundance[[6]]+ ylab(NULL), nrow =1)
  ggsave(glue("new_res/{noUP}scatterplot_abundance_gene_benchmark.png"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_abundance[[3]], list_plot_bench_abundance[[7]]+ ylab(NULL), nrow =1)
  ggsave(glue("new_res/{noUP}scatterplot_abundance_benchmark_mRNA.png"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_abundance[[4]], list_plot_bench_abundance[[8]]+ ylab(NULL), nrow =1)
  ggsave(glue("new_res/{noUP}scatterplot_abundance_gene_benchmark_mRNA.png"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_change[[1]],
                   list_plot_bench_change[[3]] +
                     ylab(NULL),
                   nrow =1) 
  ggsave(glue("new_res/{noUP}change_mrna_prot.pdf"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_change[[2]],
                   list_plot_bench_change[[4]]+ ylab(NULL), nrow =1)
  ggsave(glue("new_res/{noUP}change_mrna_prot_mRNA.pdf"), plot = a,
         height = size * scale, width = size)
  
}
