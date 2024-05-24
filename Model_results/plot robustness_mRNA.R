rm(list = ls())
setwd("~/Desktop/codice_simone")
PATH_WD = "~/Desktop/codice_simone"
id = 0
gg = list()
for(DATA in c("jurkat", "wtc11")){
  for(PEP in c(FALSE, TRUE)){
    id = id + 1
    name = paste0("AUC_", DATA, "_PEP_", PEP,".RData")
    load(name)
    
    library(ggplot2)
    DF = data.frame(AUC = AUC[nrow(AUC),],
                    c = colnames(AUC))
    
    if(DATA == "wtc11"){
      data = "WTC-11"
    }else{
      data = DATA
    }
    
    if(PEP){
      title = paste0(data, " - PEP mode")
    }else{
      title = paste0(data, " - FDR mode")
    }
    
    gg[[id]] = ggplot( aes(x=c, y=AUC), data = DF) +
      #geom_line( color="grey") +
      geom_point(color="black", 
                 size=2) + 
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
  }
}


AA = ggpubr::ggarrange( gg[[1]] + xlab(""),
                        gg[[2]] + xlab("") + ylab(""),
                        gg[[3]], 
                        gg[[4]] + ylab(""),
                          ncol = 2, nrow = 2 )

ggsave(filename = "c_robustness_mRNA.pdf",
       plot = AA,
       device = "pdf",
       width = 8,
       height = 8,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

