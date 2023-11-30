PATH_WD = commandArgs(trailingOnly = TRUE)[1]

library(glue)
library(gridExtra)
library(grid)
library(ggplot2)

size = 10 ; scale = 0.5
source(glue("{PATH_WD}/utils_function/grid_arrange_shared_legend.R"))
list_plot_bench = list()

for (noUP in c("", "no_UP_", "gene_w_more_ISO_")) {

  list_plot_bench_abundance = list()
  list_plot_bench_change = list()
  list_plot_bench_time = list()
  list_plot_bench_memory = list()
  
  for (dat in c("jurkat", "wtc11")) {
    proteases = list.dirs(glue("{PATH_WD}/Benchmark_results/{dat}"),
                          recursive = FALSE,
                          full.names = FALSE
                          )
    for (prot in proteases) {
      load(glue("{PATH_WD}/Benchmark_results/{dat}/{prot}/{noUP}ROC_main_result.rdata"))
      list_plot_bench = append(list_plot_bench, list(plot_tab))
    }
    
    for (mrna in c("", "_mRNA")) {
      load(glue("{PATH_WD}/Benchmark_results/{dat}/{noUP}scatterplot_benchmark{mrna}.rdata"))
      list_plot_bench_abundance = append(list_plot_bench_abundance, list(scat_bench))
      
      load(glue("{PATH_WD}/Benchmark_results/{dat}/{noUP}scatterplot_gene{mrna}.rdata"))
      list_plot_bench_abundance = append(list_plot_bench_abundance, list(scat_bench))
      
      load(glue("{PATH_WD}/Benchmark_results/{dat}/{noUP}change_mrna_prot{mrna}.rdata"))
      list_plot_bench_change = append(list_plot_bench_change, list(plot_change))
      
      load(glue("{PATH_WD}/Benchmark_results/{dat}/{noUP}change_mrna_prot_extreme{mrna}.rdata"))
      list_plot_bench_change = append(list_plot_bench_change, list(plot_change))
      
      load(glue("{PATH_WD}/Benchmark_results/{dat}/{noUP}scatterplot_log2fc{mrna}.rdata"))
      list_plot_bench_change = append(list_plot_bench_change, list(scat_bench))
    }
    load(glue("{PATH_WD}/Benchmark_results/{dat}/Average_Run-Time.rdata"))
    list_plot_bench_time = append(list_plot_bench_time, list(pp))
    
    load( glue("{PATH_WD}/Benchmark_results/{dat}/Average_Memory_usage.rdata"))
    list_plot_bench_memory = append(list_plot_bench_memory, list(pp))
  }
  
  a = grid.arrange(list_plot_bench_time[[1]], list_plot_bench_time[[2]] + ylab(NULL), nrow =1)
  ggsave(glue("{PATH_WD}/Benchmark_results/Average_Run-Time.png"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_memory[[1]], list_plot_bench_memory[[2]] + ylab(NULL), nrow =1)
  ggsave(glue("{PATH_WD}/Benchmark_results/Average_Memory_usage.png"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_abundance[[1]], list_plot_bench_abundance[[5]]+ ylab(NULL), nrow =1)
  ggsave(glue("{PATH_WD}/Benchmark_results/{noUP}scatterplot_abundance_benchmark.png"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_abundance[[2]], list_plot_bench_abundance[[6]]+ ylab(NULL), nrow =1)
  ggsave(glue("{PATH_WD}/Benchmark_results/{noUP}scatterplot_abundance_gene_benchmark.png"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_abundance[[3]], list_plot_bench_abundance[[7]]+ ylab(NULL), nrow =1)
  ggsave(glue("{PATH_WD}/Benchmark_results/{noUP}scatterplot_abundance_benchmark_mRNA.png"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_abundance[[4]], list_plot_bench_abundance[[8]]+ ylab(NULL), nrow =1)
  ggsave(glue("{PATH_WD}/Benchmark_results/{noUP}scatterplot_abundance_gene_benchmark_mRNA.png"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_change[[1]],
                   list_plot_bench_change[[7]] +
                     ylab(NULL),
                   nrow =1) 
  ggsave(glue("{PATH_WD}/Benchmark_results/{noUP}change_mrna_prot.pdf"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_change[[2]],
                   list_plot_bench_change[[8]]+ ylab(NULL), nrow =1)
  ggsave(glue("{PATH_WD}/Benchmark_results/{noUP}change_mrna_prot_extreme.pdf"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_change[[3]], list_plot_bench_change[[9]]+ ylab(NULL), nrow =1)
  ggsave(glue("{PATH_WD}/Benchmark_results/{noUP}scatterplot_log2fc.png"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_change[[4]],
                   list_plot_bench_change[[10]]+ ylab(NULL), nrow =1)
  ggsave(glue("{PATH_WD}/Benchmark_results/{noUP}change_mrna_prot_mRNA.pdf"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_change[[5]],
                   list_plot_bench_change[[11]] + ylab(NULL), nrow =1)
  ggsave(glue("{PATH_WD}/Benchmark_results/{noUP}change_mrna_prot_extreme_mRNA.pdf"), plot = a,
         height = size * scale, width = size)
  
  a = grid.arrange(list_plot_bench_change[[6]], list_plot_bench_change[[12]]+ ylab(NULL), nrow =1)
  ggsave(glue("{PATH_WD}/Benchmark_results/{noUP}scatterplot_log2fc_mRNA.png"), plot = a,
         height = size * scale, width = size)
}

a = ggpubr::ggarrange(list_plot_bench[[1]]$gplot + xlab(NULL),
                  list_plot_bench[[2]]$gplot + xlab(NULL) + ylab(NULL),
                  list_plot_bench[[3]]$gplot + xlab(NULL) + ylab(NULL),
                  list_plot_bench[[4]]$gplot + xlab(NULL) + ylab(NULL),
                  list_plot_bench[[5]]$gplot + xlab(NULL) + ylab(NULL),
                  list_plot_bench[[6]]$gplot,
                  list_plot_bench[[7]]$gplot + ylab(NULL),
                  list_plot_bench[[8]]$gplot + ylab(NULL),
                  list_plot_bench[[9]]$gplot + ylab(NULL),
                  list_plot_bench[[10]]$gplot + ylab(NULL),
                  nrow = 2, ncol = 5, common.legend = TRUE, legend="bottom")

ggsave(glue("{PATH_WD}/Benchmark_results/ROC_main_result.png"), plot = a,
       height = size * scale, width = size)

a = ggpubr::ggarrange(list_plot_bench[[11]]$gplot + xlab(NULL),
                               list_plot_bench[[12]]$gplot + xlab(NULL) + ylab(NULL),
                               list_plot_bench[[13]]$gplot + xlab(NULL) + ylab(NULL),
                               list_plot_bench[[14]]$gplot + xlab(NULL) + ylab(NULL),
                               list_plot_bench[[15]]$gplot + xlab(NULL) + ylab(NULL),
                               list_plot_bench[[16]]$gplot,
                               list_plot_bench[[17]]$gplot + ylab(NULL),
                               list_plot_bench[[18]]$gplot + ylab(NULL),
                               list_plot_bench[[19]]$gplot + ylab(NULL),
                               list_plot_bench[[20]]$gplot + ylab(NULL),
                      nrow = 2, ncol = 5, common.legend = TRUE, legend="bottom")

ggsave(glue("{PATH_WD}/Benchmark_results/no_UP_ROC_main_result.png"), plot = a,
       height = size * scale, width = size)

a = ggpubr::ggarrange(list_plot_bench[[21]]$gplot + xlab(NULL),
                      list_plot_bench[[22]]$gplot + xlab(NULL) + ylab(NULL),
                      list_plot_bench[[23]]$gplot + xlab(NULL) + ylab(NULL),
                      list_plot_bench[[24]]$gplot + xlab(NULL) + ylab(NULL),
                      list_plot_bench[[25]]$gplot + xlab(NULL) + ylab(NULL),
                      list_plot_bench[[26]]$gplot,
                      list_plot_bench[[27]]$gplot + ylab(NULL),
                      list_plot_bench[[28]]$gplot + ylab(NULL),
                      list_plot_bench[[29]]$gplot + ylab(NULL),
                      list_plot_bench[[30]]$gplot + ylab(NULL),
                      nrow = 2, ncol = 5, common.legend = TRUE, legend="bottom")

ggsave(glue("{PATH_WD}/Benchmark_results/gene_w_more_ISO_ROC_main_result.png"), plot = a,
       height = size * scale, width = size)

############################# 3.1 ###################################################
list_plot_roc = list()
list_plot_abundance = list()
list_plot_change = list()

for (dat in c("jurkat", "wtc11")) {
  proteases = list.dirs(glue("{PATH_WD}/Benchmark_results/{dat}"),
                        recursive = FALSE,
                        full.names = FALSE
                        )
  for (prot in proteases) {
    load(glue("{PATH_WD}/Robustness/{dat}/{prot}/ROC_MM_psm_vs_MM_intensities_vs_OpenMS.rdata"))
    list_plot_roc = append(list_plot_roc, list(plot_tab))
  }
  for (mrna in c("", "_mRNA")) {
    
    load(glue("{PATH_WD}/Robustness/{dat}/change_mrna_prot{mrna}_MM_psm_vs_MM_intensities_vs_OpenMS.rdata"))
    plot_change$name = glue("change_mrna_prot{mrna}_MM_psm_vs_MM_intensities_vs_OpenMS")
    list_plot_change = append(list_plot_change, list(plot_change))
    
    load(glue("{PATH_WD}/Robustness/{dat}/change_mrna_prot_extreme{mrna}_MM_psm_vs_MM_intensities_vs_OpenMS.rdata"))
    plot_change$name = glue("change_mrna_prot_extreme{mrna}_MM_psm_vs_MM_intensities_vs_OpenMS")
    list_plot_change = append(list_plot_change, list(plot_change))
    
    for (input in c("OpenMS", "MM_psm")) {
      load(glue("{PATH_WD}/Robustness/{dat}/scatterplot_{input}{mrna}_MM_psm_vs_MM_intensities_vs_OpenMS.rdata"))
      scat_bench$name = glue("scatterplot_{input}{mrna}_MM_psm_vs_MM_intensities_vs_OpenMS")
      list_plot_abundance = append(list_plot_abundance, list(scat_bench))
      
      load(glue("{PATH_WD}/Robustness/{dat}/scatterplot_log2fc_{input}{mrna}_MM_psm_vs_MM_intensities_vs_OpenMS.rdata"))
      scat_bench$name = glue("scatterplot_log2fc_{input}{mrna}_MM_psm_vs_MM_intensities_vs_OpenMS")
      list_plot_change = append(list_plot_change, list(scat_bench))
    }
  }
}

a = ggpubr::ggarrange(list_plot_roc[[1]]$gplot + xlab(NULL),
                      list_plot_roc[[2]]$gplot + xlab(NULL) + ylab(NULL),
                      list_plot_roc[[3]]$gplot + xlab(NULL) + ylab(NULL),
                      list_plot_roc[[4]]$gplot + xlab(NULL) + ylab(NULL),
                      list_plot_roc[[5]]$gplot + xlab(NULL) + ylab(NULL),
                      list_plot_roc[[6]]$gplot,
                      list_plot_roc[[7]]$gplot + ylab(NULL),
                      list_plot_roc[[8]]$gplot + ylab(NULL),
                      list_plot_roc[[9]]$gplot + ylab(NULL),
                      list_plot_roc[[10]]$gplot + ylab(NULL),
                      nrow = 2, ncol = 5, common.legend = TRUE, legend="bottom")

ggsave(glue("{PATH_WD}/Robustness/ROC_MM_psm_vs_MM_intensities_vs_OpenMS.png"), plot = a,
       height = size * scale, width = size)

half_length = length(list_plot_abundance)/2
for (i in 1:half_length) {
  a = grid.arrange(list_plot_abundance[[i]], list_plot_abundance[[i+half_length]] + ylab(NULL), nrow =1)
  ggsave(glue("{PATH_WD}/Robustness/{list_plot_abundance[[i]]$name}.png"), plot = a,
         height = size * scale, width = size)
}

half_length = length(list_plot_change)/2
for (i in 1:half_length) {
  ext = "png"
  if(grepl("change", list_plot_change[[i]]$name)){
    a = grid_arrange_shared_legend(list_plot_change[[i]], list_plot_change[[i+half_length]] + ylab(NULL), nrow =1)
    ext = "pdf"
  }else{
    a = grid.arrange(list_plot_change[[i]], list_plot_change[[i+half_length]] + ylab(""), nrow =1)
  }
  ggsave(glue("{PATH_WD}/Robustness/{list_plot_change[[i]]$name}.{ext}"), plot = a,
         height = size * scale, width = size)
}

############################# 3.2 ###################################################
list_plot_roc = list()
list_plot_abundance = list()
list_plot_change = list()

for (dat in c("jurkat", "wtc11")) {
  for (mrna in c("", "_mRNA")) {
    for (input in c("MM_intensities", "MM_psm")) {
      load(glue("{PATH_WD}/Robustness/{dat}/scatterplot_{input}{mrna}_MM_psm_vs_MM_intensities_vs_OpenMS.rdata"))
      scat_bench$name = glue("scatterplot_{input}{mrna}_MM_psm_vs_MM_intensities_vs_OpenMS")
      list_plot_abundance = append(list_plot_abundance, list(scat_bench))
      
      load(glue("{PATH_WD}/Robustness/{dat}/scatterplot_log2fc_{input}{mrna}_MM_psm_vs_MM_intensities_vs_OpenMS.rdata"))
      scat_bench$name = glue("scatterplot_log2fc_{input}{mrna}_MM_psm_vs_MM_intensities_vs_OpenMS")
      list_plot_change = append(list_plot_change, list(scat_bench))
    }
  }
}

half_length = length(list_plot_abundance)/2
for (i in 1:half_length) {
  ext = "png"
  a = grid.arrange(list_plot_abundance[[i]],
                   list_plot_abundance[[i+half_length]] + ylab(NULL),
                   nrow =1)
  ggsave(glue("{PATH_WD}/Robustness/{list_plot_abundance[[i]]$name}.{ext}"), plot = a,
         height = size * scale, width = size)
}

half_length = length(list_plot_change)/2
for (i in 1:half_length) {
  ext = "png"
  if(grepl("change", list_plot_change[[i]]$name)){
    ext = "pdf"
    a = grid_arrange_shared_legend(list_plot_change[[i]],
                                   list_plot_change[[i+half_length]] + ylab(NULL),
                                   nrow =1)
  }else{
    a = grid.arrange(list_plot_change[[i]],
                     list_plot_change[[i+half_length]] + ylab(NULL),
                     nrow =1)
  }
  ggsave(glue("{PATH_WD}/Robustness/{list_plot_change[[i]]$name}.{ext}"), plot = a,
         height = size * scale, width = size)
}

############################# 4 ###################################################
for (input in c("MM_psm", "MM_intensities")) {
  list_plot_roc = list()
  list_plot_abundance = list()
  list_plot_change = list()
  for (dat in c("jurkat", "wtc11")) {
    proteases = list.dirs(glue("{PATH_WD}/Benchmark_results/{dat}"), recursive = FALSE, full.names = FALSE)
    for (prot in proteases) {
      load(glue("{PATH_WD}/Robustness/{dat}/{prot}/ROC_{input}_pep_vs_no_pep.rdata"))
      list_plot_roc = append(list_plot_roc, list(plot_tab))
    }
    for (mrna in c("", "_mRNA")) {
      for (pep in c("_PEP", "")) {
        load(glue("{PATH_WD}/Robustness/{dat}/scatterplot_{input}{pep}{mrna}_pep_vs_no_pep.rdata"))
        scat_bench$name = glue("scatterplot_{input}{pep}{mrna}_pep_vs_no_pep")
        list_plot_abundance = append(list_plot_abundance, list(scat_bench))
        
        load(glue("{PATH_WD}/Robustness/{dat}/change_mrna_prot_{input}{mrna}_pep_vs_no_pep.rdata"))
        plot_change$name = glue("change_mrna_prot_{input}{mrna}_pep_vs_no_pep")
        list_plot_change = append(list_plot_change, list(plot_change))
        
        load(glue("{PATH_WD}/Robustness/{dat}/change_mrna_prot_extreme_{input}{mrna}_pep_vs_no_pep.rdata"))
        plot_change$name = glue("change_mrna_prot_extreme_{input}{mrna}_pep_vs_no_pep")
        list_plot_change = append(list_plot_change, list(plot_change))
        
        load(glue("{PATH_WD}/Robustness/{dat}/scatterplot_log2fc_{input}{pep}{mrna}_pep_vs_no_pep.rdata"))
        scat_bench$name = glue("scatterplot_log2fc_{input}{pep}{mrna}_pep_vs_no_pep")
        list_plot_change = append(list_plot_change, list(scat_bench))
      }
    }
  }
  a = ggpubr::ggarrange(list_plot_roc[[1]]$gplot + xlab(NULL),
                        list_plot_roc[[2]]$gplot + xlab(NULL) + ylab(NULL),
                        list_plot_roc[[3]]$gplot + xlab(NULL) + ylab(NULL),
                        list_plot_roc[[4]]$gplot + xlab(NULL) + ylab(NULL),
                        list_plot_roc[[5]]$gplot + xlab(NULL) + ylab(NULL),
                        list_plot_roc[[6]]$gplot,
                        list_plot_roc[[7]]$gplot + ylab(NULL),
                        list_plot_roc[[8]]$gplot + ylab(NULL),
                        list_plot_roc[[9]]$gplot + ylab(NULL),
                        list_plot_roc[[10]]$gplot + ylab(NULL),
                        nrow = 2, ncol = 5,
                        common.legend = TRUE,
                        legend="bottom"
                        )
  ggsave(glue("{PATH_WD}/Robustness/ROC_{input}_pep_vs_no_pep.png"), plot = a,
         height = size * scale, width = size)
  
  half_length = length(list_plot_abundance)/2
  for (i in 1:half_length) {
    ext = "png"
    a = grid.arrange(list_plot_abundance[[i]],
                     list_plot_abundance[[i+half_length]] + ylab(NULL),
                     nrow =1)
    ggsave(glue("{PATH_WD}/Robustness/{list_plot_abundance[[i]]$name}.{ext}"), plot = a,
           height = size * scale, width = size)
  }
  
  half_length = length(list_plot_change)/2
  for (i in 1:half_length) {
    ext = "png"
    if(grepl("change", list_plot_change[[i]]$name)){
      ext = "pdf"
      a = grid_arrange_shared_legend(list_plot_change[[i]],
                                     list_plot_change[[i+half_length]] + ylab(NULL),
                                     nrow =1)
    }else{
      a = grid.arrange(list_plot_change[[i]],
                       list_plot_change[[i+half_length]] + ylab(NULL),
                       nrow =1)
    }
    
    ggsave(glue("{PATH_WD}/Robustness/{list_plot_change[[i]]$name}.{ext}"), plot = a,
           height = size * scale, width = size)
  }
}
