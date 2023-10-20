library(xtable)
library(glue)
PATH_WD = "/home/jbollon/prot_iso_mrna_dev/IsoBayes_paper"

final_tab = rbind()
vec_dat = c()
vec_prot = c()

for (dat in c("jurkat", "wtc11")) {
  if(dat == "wtc11"){
    DATA_name = "WTC-11"
  }else{
    DATA_name = "jurkat"
  }
  proteases = list.dirs(glue("{PATH_WD}/Benchmark_results/{dat}"), recursive = FALSE, full.names = FALSE)
  for (prot in proteases) {
    vec_dat = c(vec_dat, DATA_name)
    vec_prot = c(vec_prot, prot)
    tab = data.table::fread(glue("{PATH_WD}/Benchmark_results/{dat}/{prot}/SumTab_main_result.csv"))
    i_ord = match(c("IsoBayes_mRNA", "IsoBayes", "EPIFANY", "PIA", "Fido"), tab$Model)
    final_tab = rbind(final_tab, t(tab[i_ord, "AUC"]))
  }
}

final_tab = rbind(final_tab, colMeans(final_tab))
final_tab = round(final_tab, 2)
final_tab = cbind(c(vec_dat, NA), c(vec_prot, NA), final_tab)

print(xtable(final_tab, type = "latex"), include.rownames=FALSE)

##################################################################################################
# no up
##################################################################################################
final_tab = rbind()
vec_dat = c()
vec_prot = c()

for (dat in c("jurkat", "wtc11")) {
  if(dat == "wtc11"){
    DATA_name = "WTC-11"
  }else{
    DATA_name = "jurkat"
  }
  proteases = list.dirs(glue("{PATH_WD}/Benchmark_results/{dat}"), recursive = FALSE, full.names = FALSE)
  for (prot in proteases) {
    vec_dat = c(vec_dat, DATA_name)
    vec_prot = c(vec_prot, prot)
    tab = data.table::fread(glue("{PATH_WD}/Benchmark_results/{dat}/{prot}/SumTab_main_result_no_UP.csv"))
    i_ord = match(c("IsoBayes_mRNA", "IsoBayes", "EPIFANY", "PIA", "Fido"), tab$Model)
    final_tab = rbind(final_tab, t(tab[i_ord, "AUC_only_shared"]))
  }
}

final_tab = rbind(final_tab, colMeans(final_tab))
final_tab = round(final_tab, 2)
final_tab = cbind(c(vec_dat, NA), c(vec_prot, NA), final_tab)

print(xtable(final_tab, type = "latex"), include.rownames=FALSE)

##################################################################################################
# Correlation openms
##################################################################################################
final_tab = rbind()
vec_method = c("IsoBayes_mRNA", "IsoBayes")

for (mrna in c("_mRNA", "")) {
  row_val = c()
  for (noUP in c("", "no_UP_")) {
    for (dat in c("jurkat", "wtc11")) {
      load(glue("{PATH_WD}/Benchmark_results/{dat}/{noUP}scatterplot_benchmark{mrna}.rdata"))
      correlation = strsplit(scat_bench$labels$subtitle, " ")[[1]][2]
      
      row_val = c(row_val, correlation)
    }
  }
  final_tab = rbind(final_tab, row_val)
}
final_tab = cbind(vec_method, final_tab)

print(xtable(final_tab, type = "latex"), include.rownames=FALSE)

##################################################################################################
# Correlation Gene Abundance
##################################################################################################
final_tab = rbind()
vec_method = c("IsoBayes_mRNA", "IsoBayes")

for (mrna in c("_mRNA", "")) {
  row_val = c()
  for (dat in c("jurkat", "wtc11")) {
    load(glue("{PATH_WD}/Benchmark_results/{dat}/scatterplot_gene{mrna}.rdata"))
    correlation = strsplit(scat_bench$labels$subtitle, " ")[[1]][2]
    
    row_val = c(row_val, correlation)
  }
  
  final_tab = rbind(final_tab, row_val)
}
final_tab = cbind(vec_method, final_tab)

print(xtable(final_tab, type = "latex"), include.rownames=FALSE)

##################################################################################################
# Correlation log2fc openms
##################################################################################################
final_tab = rbind()
vec_method = c("IsoBayes_mRNA", "IsoBayes")

for (mrna in c("_mRNA", "")) {
  row_val = c()
  for (noUP in c("", "no_UP_")) {
    for (dat in c("jurkat", "wtc11")) {
      load(glue("{PATH_WD}/Benchmark_results/{dat}/{noUP}scatterplot_log2fc{mrna}.rdata"))
      correlation = strsplit(scat_bench$labels$subtitle, " ")[[1]][2]
      
      row_val = c(row_val, correlation)
    }
  }
  final_tab = rbind(final_tab, row_val)
}
final_tab = cbind(vec_method, final_tab)

print(xtable(final_tab, type = "latex"), include.rownames=FALSE)


