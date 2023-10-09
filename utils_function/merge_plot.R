library(glue)
library(gridExtra)
library(grid)

PATH_WD = "/home/jbollon/prot_iso_mrna_dev/IsoBayes_paper"
size = 10 ; scale = 0.5
source(glue("{PATH_WD}/utils_function/grid_arrange_shared_legend.R"))

list_plot_bench = list()
list_plot_bench_abundance = list()
list_plot_bench_change = list()
list_plot_bench_time = list()
list_plot_bench_memory = list()

for (dat in c("jurkat", "wtc11")) {
  load(glue("{PATH_WD}/Benchmark_results/{dat}/ROC_main_result"))
  list_plot_bench = append(list_plot_bench, list(plot_tab))
  
  for (mrna in c("", "_mRNA")) {
    load(glue("{PATH_WD}/Benchmark_results/{dat}/scatterplot_benchmark{mrna}"))
    list_plot_bench_abundance = append(list_plot_bench_abundance, list(scat_bench))
    
    load(glue("{PATH_WD}/Benchmark_results/{dat}/change_mrna_prot{mrna}"))
    list_plot_bench_change = append(list_plot_bench_change, list(plot_change))
    
    load(glue("{PATH_WD}/Benchmark_results/{dat}/change_mrna_prot_extreme{mrna}"))
    list_plot_bench_change = append(list_plot_bench_change, list(plot_change))
    
    load(glue("{PATH_WD}/Benchmark_results/{dat}/scatterplot_log2fc{mrna}"))
    list_plot_bench_change = append(list_plot_bench_change, list(scat_bench))
  }
  load(glue("{PATH_WD}/Benchmark_results/{dat}/Average_Run-Time"))
  list_plot_bench_time = append(list_plot_bench_time, list(pp))
  
  load( glue("{PATH_WD}/Benchmark_results/{dat}/Average_Memory_usage"))
  list_plot_bench_memory = append(list_plot_bench_memory, list(pp))
}

a = grid_arrange_shared_legend(list_plot_bench[[1]]$gplot, list_plot_bench[[2]]$gplot, nrow =1)
ggsave(glue("{PATH_WD}/Benchmark_results/ROC_main_result.png"), plot = a,
       height = size * scale, width = size)

a = grid.arrange(list_plot_bench_time[[1]], list_plot_bench_time[[2]], nrow =1)
ggsave(glue("{PATH_WD}/Benchmark_results/Average_Run-Time.png"), plot = a,
       height = size * scale, width = size)

a = grid.arrange(list_plot_bench_memory[[1]], list_plot_bench_memory[[2]], nrow =1)
ggsave(glue("{PATH_WD}/Benchmark_results/Average_Memory_usage.png"), plot = a,
       height = size * scale, width = size)

a = grid.arrange(list_plot_bench_abundance[[1]], list_plot_bench_abundance[[3]], nrow =1)
ggsave(glue("{PATH_WD}/Benchmark_results/scatterplot_abundance_benchmark.png"), plot = a,
       height = size * scale, width = size)

a = grid.arrange(list_plot_bench_abundance[[2]], list_plot_bench_abundance[[4]], nrow =1)
ggsave(glue("{PATH_WD}/Benchmark_results/scatterplot_abundance_benchmark_mRNA.png"), plot = a,
       height = size * scale, width = size)

a = grid.arrange(list_plot_bench_change[[1]], list_plot_bench_change[[7]], nrow =1) 
ggsave(glue("{PATH_WD}/Benchmark_results/change_mrna_prot.png"), plot = a,
       height = size * scale, width = size)

a = grid.arrange(list_plot_bench_change[[2]], list_plot_bench_change[[8]], nrow =1)
ggsave(glue("{PATH_WD}/Benchmark_results/change_mrna_prot_extreme.png"), plot = a,
       height = size * scale, width = size)

a = grid.arrange(list_plot_bench_change[[3]], list_plot_bench_change[[9]], nrow =1)
ggsave(glue("{PATH_WD}/Benchmark_results/scatterplot_log2fc.png"), plot = a,
       height = size * scale, width = size)

a = grid.arrange(list_plot_bench_change[[4]], list_plot_bench_change[[10]], nrow =1)
ggsave(glue("{PATH_WD}/Benchmark_results/change_mrna_prot_mRNA.png"), plot = a,
       height = size * scale, width = size)

a = grid.arrange(list_plot_bench_change[[5]], list_plot_bench_change[[11]], nrow =1)
ggsave(glue("{PATH_WD}/Benchmark_results/change_mrna_prot_extreme_mRNA.png"), plot = a,
       height = size * scale, width = size)

a = grid.arrange(list_plot_bench_change[[6]], list_plot_bench_change[[12]], nrow =1)
ggsave(glue("{PATH_WD}/Benchmark_results/scatterplot_log2fc_mRNA.png"), plot = a,
       height = size * scale, width = size)
