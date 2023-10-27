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
  proteases = list.dirs(glue("{PATH_WD}/Benchmark_results/{dat}"),
                        recursive = FALSE,
                        full.names = FALSE
                        )
  for (prot in proteases) {
    vec_dat = c(vec_dat, DATA_name)
    vec_prot = c(vec_prot, prot)
    tab = data.table::fread(glue("{PATH_WD}/Benchmark_results/{dat}/{prot}/SumTab_main_result.csv"))
    i_ord = match(c("IsoBayes_mRNA", "IsoBayes", "EPIFANY", "PIA", "Fido"),
                  tab$Model
                  )
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
  proteases = list.dirs(glue("{PATH_WD}/Benchmark_results/{dat}"),
                        recursive = FALSE,
                        full.names = FALSE
                        )
  for (prot in proteases) {
    vec_dat = c(vec_dat, DATA_name)
    vec_prot = c(vec_prot, prot)
    tab = data.table::fread(glue("{PATH_WD}/Benchmark_results/{dat}/{prot}/SumTab_main_result_no_UP.csv"))
    i_ord = match(c("IsoBayes_mRNA", "IsoBayes", "EPIFANY", "PIA", "Fido"),
                  tab$Model
                  )
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
# Correlation openms log2fc
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

####################################################################################
# OpenMS - MM_psm - MM_intensities
####################################################################################
for (dat in c("jurkat", "wtc11")) {
  final_tab = rbind()
  vec_dat = c()
  vec_prot = c()
  if(dat == "wtc11"){
    DATA_name = "WTC-11"
  }else{
    DATA_name = "jurkat"
  }
  proteases = list.dirs(glue("{PATH_WD}/Robustness/{dat}"),
                        recursive = FALSE,
                        full.names = FALSE
                        )
  for (prot in proteases) {
    vec_dat = c(vec_dat, DATA_name)
    vec_prot = c(vec_prot, prot)
    tab = data.table::fread(glue("{PATH_WD}/Robustness/{dat}/{prot}/SumTab_MM_psm_vs_MM_intensities_vs_OpenMS.csv"))
    i_ord = match(c("IsoBayes_OpenMS", "IsoBayes_MM_psm",
                    "IsoBayes_MM_intensities", "IsoBayes_mRNA_OpenMS",
                    "IsoBayes_mRNA_MM_psm", "IsoBayes_mRNA_MM_intensities"),
                  tab$Model
                  )
    final_tab = rbind(final_tab, t(tab[i_ord, "AUC"]))
  }
  final_tab = rbind(final_tab, colMeans(final_tab))
  final_tab = round(final_tab, 2)
  final_tab = cbind(c(vec_dat, NA), c(vec_prot, NA), final_tab)
  
  print(xtable(final_tab, type = "latex"), include.rownames=FALSE)
}

##################################################################################################
# Correlation OpenMS - MM_psm - MM_intensities
##################################################################################################
final_tab = rbind()

for (dat in c("jurkat", "wtc11")) {
  row_val = c()
  for (mrna in c("", "_mRNA")) {
    for (input in c("OpenMS", "MM_psm", "MM_intensities")) {
      load(glue("{PATH_WD}/Robustness/{dat}/scatterplot_{input}{mrna}_MM_psm_vs_MM_intensities_vs_OpenMS.rdata"))
      correlation = strsplit(scat_bench$labels$subtitle, " ")[[1]][2]
      row_val = c(row_val, correlation)
    }
  }
  final_tab = rbind(final_tab, row_val)
}
print(xtable(final_tab, type = "latex"), include.rownames=FALSE)

####################################################################################
# Pep vs no Pep -- MM_psm - MM_intensities
####################################################################################
for (input in c("MM_psm", "MM_intensities")) {
  for (dat in c("jurkat", "wtc11")) {
    final_tab = rbind()
    vec_dat = c()
    vec_prot = c()
    if(dat == "wtc11"){
      DATA_name = "WTC-11"
    }else{
      DATA_name = "jurkat"
    }
    proteases = list.dirs(glue("{PATH_WD}/Robustness/{dat}"),
                          recursive = FALSE,
                          full.names = FALSE
                          )
    for (prot in proteases) {
      vec_dat = c(vec_dat, DATA_name)
      vec_prot = c(vec_prot, prot)
      tab = data.table::fread(glue("{PATH_WD}/Robustness/{dat}/{prot}/SumTab_{input}_pep_vs_no_pep.csv"))
      i_ord = match(c("IsoBayes_PEP", "IsoBayes",
                      "IsoBayes_mRNA_PEP", "IsoBayes_mRNA"),
                    tab$Model
                    )
      final_tab = rbind(final_tab, t(tab[i_ord, "AUC"]))
    }
    final_tab = rbind(final_tab, colMeans(final_tab))
    final_tab = round(final_tab, 2)
    final_tab = cbind(c(vec_dat, NA), c(vec_prot, NA), final_tab)
    
    print(input)
    print(xtable(final_tab, type = "latex"), include.rownames=FALSE)
  }
}

##################################################################################################
# Correlation (pep vs no pep) MM_psm - MM_intensities 
##################################################################################################
for (input in c("MM_psm", "MM_intensities")) {
  for (dat in c("jurkat", "wtc11")) {
    row_val = c()
    for (mrna in c("", "_mRNA")) {
      for (pep in c("_PEP", "")) {
        load(glue("{PATH_WD}/Robustness/{dat}/scatterplot_{input}{pep}{mrna}_pep_vs_no_pep.rdata"))
        correlation = strsplit(scat_bench$labels$subtitle, " ")[[1]][2]
        row_val = c(row_val, correlation)
      }
    }
    print(c(input, dat))
    print(xtable(t(as.matrix(row_val)), type = "latex"), include.rownames=FALSE)
  }
}

####################################################################################################
# MM - Run time & memory
####################################################################################################
Data = rbind()

for (ff in list.files(glue("{PATH_WD}/Benchmark_results/internal_results"), pattern = "res_used")) {
  res_used = readLines(glue("{PATH_WD}/Benchmark_results/internal_results/{ff}"))
  
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
  df = data.frame(Model = info_model[3], data = info_model[2],
                  protease = info_model[7], mRNA = info_model[8],
                  pep = info_model[5],
                  RunTime = wallTime, RAM = ramMB
                  )
  
  Data = rbind(Data, df)
}

Data[is.na(Data)] = ""

for (model in c("psm", "intensities")) {
  for (dat in c("jurkat", "wtc11")){
    row_dataset = cbind()
    for (mrna in c("", "mRNA")) {
      for (pep in c("TRUE", "FALSE")) {
        data_protease_all = cbind()
        for (protease in proteases){
          data_protease = Data[Data$Model == model & Data$data == dat & Data$mRNA == mrna & Data$pep == pep & Data$protease == protease, ]
          data_protease_mean = colMeans(data_protease[, c("RunTime", "RAM")])
          data_protease_all = cbind(data_protease_all, as.matrix(data_protease_mean))
        }
        row_dataset = cbind(row_dataset, as.matrix(rowMeans(data_protease_all)))
      }
    }
    print(model) ; print(dat)
    row_dataset[2, ] = round(row_dataset[2, ]/1000, 1)
    print(xtable(row_dataset, type = "latex", digits = 1), include.rownames=FALSE)
  }
}