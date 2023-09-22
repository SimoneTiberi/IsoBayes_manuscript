library(ggplot2)
library(glue)
theme_set(theme_classic())
dataset = "jurkat"
path_std = glue("/home/jbollon/prot_iso_mrna_dev/IsoBayes/Benchmark_results/Benchmarking_{dataset}")
path_no_PEP = glue("/home/jbollon/prot_iso_mrna_dev/IsoBayes/Benchmark_results/Benchmarking_{dataset}_no_PEP")
data_plot_all_RAM = NULL
data_plot_all = NULL
list_data_enzyme_time = list()
list_data_enzyme_RAM = list()
#PROTEASES = c("AspN", "Chym", "LysC", "Trypsin", "ArgC", "GluC")

palette_models = c("Fido" = "#87bc45",
                   "Epifany" = "#b33dc6",
                   "IsoBayes_mRNA" = "#ea5545",
                   "IsoBayes_mRNA\nNo PEP" = "#f46a9b",
                   "PIA" = "#27aeef")

for (enzyme in list.dirs(path_std, recursive = F, full.names = F)) {
  Data = rbind()
  
  for (path in c(path_std, path_no_PEP)) {
    for (ff in list.files(glue("{path}/{enzyme}"), pattern = "res_used")) {
      res_used = readLines(glue("{path}/{enzyme}/{ff}"))
      i_res_used = grep('walltime', res_used)
      
      ramMB = round(as.numeric(gsub("kb resources_used.ncp.*", "", gsub(".*used.mem=", "", res_used[i_res_used]))) * 0.0009765625)
      
      wallTime = gsub(".*walltime=", "", res_used[i_res_used])
      wallTime = as.POSIXlt(wallTime,format="%H:%M:%S")
      wallTime = unclass(wallTime)
      wallTime = round(wallTime$min + wallTime$sec / 60, 2)
      
      if(path == path_no_PEP){
        ff = gsub("IsoBayes", "IsoBayesnoPEP", ff)
      }
      
      df = data.frame(tool = strsplit(ff, "_")[[1]][2], RunTime = wallTime, RAM = ramMB)
      Data = rbind(Data, df)
    }
  }
  
  Data$tool[grep("PiaTot", Data$tool)] = "PIA"
  Data$tool[grep("epifany", Data$tool)] = "Epifany"
  Data$tool[grep("fido", Data$tool)] = "Fido"
  Data$tool[grep("IsoBayesnoPEP", Data$tool)] = "IsoBayes\nNo PEP"
  
  # select model to plot
  data_plot = Data
  
  # RUN TIME
  data_plot = aggregate(data_plot$RunTime, by = list(data_plot$tool), FUN = sum)
  colnames(data_plot) = c("Model", "RunTime")
  id_sort = sort(data_plot$RunTime, decreasing = TRUE, index.return = TRUE)$ix
  data_plot$Model = factor(data_plot$Model, levels = data_plot$Model[id_sort])
  
  ggplot(data_plot, aes(x = Model, y = RunTime, label = RunTime, fill = Model)) +
    geom_bar(stat = "identity") +
    scale_fill_manual("Model", values = palette_models) +
    geom_text(size = 4.25, colour = "black",  vjust=-0.25) +
    labs(y = "Run-Time (Min)",
         x = "Tool",
         title = glue("Run-Time - {enzyme} - {dataset}")) +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10)
    )
  ggsave(glue("{path_std}/{enzyme}/Run-Time.png"))
  
  if(is.null(data_plot_all)){
    data_plot_all = data_plot
  }else{
    data_plot_all$RunTime = data_plot_all$RunTime + data_plot$RunTime
  }
  list_data_enzyme_time = append(list_data_enzyme_time, list(data_plot, enzyme))
  
  #Resident Set Size: how much memory a process is consuming in our physical RAM
  data_plot = Data
  Data_max = aggregate(data_plot$RAM, by = list(data_plot$tool), FUN = max)
  colnames(Data_max) = c("Model", "RAM")
  
  id_sort = sort(Data_max$RAM, decreasing = TRUE, index.return = TRUE)$ix
  Data_max$Model = factor(Data_max$Model, levels = Data_max$Model[id_sort])
  
  ggplot(Data_max, aes(x = Model, y = RAM, label = RAM, fill = Model)) +
    geom_bar(stat = "identity") +
    scale_fill_manual("Model", values = palette_models) +
    geom_text(size = 4.25, colour = "black",  vjust=-0.25) +
    labs(y = "RAM (MB)",
         x = "Tool",
         title = glue("RAM - {enzyme} - {dataset}")) +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10)
    )
  ggsave(glue("{path_std}/{enzyme}/Memory_usage.png"))
  
  if(is.null(data_plot_all_RAM)){
    data_plot_all_RAM = Data_max
  }else{
    data_plot_all_RAM$RAM = data_plot_all_RAM$RAM + Data_max$RAM
  }
  list_data_enzyme_RAM = append(list_data_enzyme_RAM, list(Data_max, enzyme))
}

#############################################
# average time
#############################################

data_plot_all$RunTime = round(data_plot_all$RunTime / length(list_data_enzyme_time)/2, 1)
id_sort = sort(data_plot_all$RunTime, decreasing = TRUE, index.return = TRUE)$ix
data_plot_all$Model = factor(data_plot_all$Model, levels = data_plot_all$Model[id_sort])

ggplot(data_plot_all, aes(x = Model, y = RunTime, label = RunTime, fill = Model)) +
  geom_bar(stat = "identity") +
  scale_fill_manual("Model", values = palette_models) +
  geom_text(size = 4.25, colour = "black",  vjust=-0.25) +
  labs(y = "Run-Time (Min)",
       x = "Tool",
       title = glue("Run-Time - {enzyme} - {dataset}")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10)
  )
ggsave(glue("{path_std}/Average_Run-Time.png"))

#############################################
# average RAM
#############################################

data_plot_all_RAM$RAM = round(data_plot_all_RAM$RAM/length(list_data_enzyme_RAM)/2)
id_sort = sort(data_plot_all_RAM$RAM, decreasing = TRUE, index.return = TRUE)$ix
data_plot_all_RAM$Model = factor(data_plot_all_RAM$Model, levels = data_plot_all_RAM$Model[id_sort])

ggplot(data_plot_all_RAM, aes(x = Model, y = RAM, label = RAM, fill = Model)) +
  geom_bar(stat = "identity") +
  scale_fill_manual("Model", values = palette_models) +
  geom_text(size = 4.25, colour = "black",  vjust=-0.25) +
  labs(y = "RAM (MB)",
       x = "Tool",
       title = glue("RAM - {enzyme} - {dataset}")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10)
  )
ggsave(glue("{path_std}/Average_Memory_usage.png"))

