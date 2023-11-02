PATH_WD = commandArgs(trailingOnly = TRUE)[1]

# Set path, global variable and libraries
################################################################################
library(glue)
source(glue("{PATH_WD}/utils_function/log_output.R"))

PATH_TO_FASTA = glue("{PATH_WD}/Data/jurkat/isoform.fasta")
log_output("get_gene_jurkat")

# get jurkat gene
################################################################################
main = function(){
  fasta = data.table::fread(PATH_TO_FASTA, header = FALSE)
  fasta = fasta[grep(">", fasta$V1), ]
  fasta$V1 = gsub(">", "", fasta$V1)
  fasta$isoform = gsub(" OS.*", "", fasta$V1)
  fasta$gene = gsub(".*GN=", "", fasta$V1)
  fasta$V1 = NULL
  
  print("Number of genes: {nrow(fasta)}")
  
  data.table::fwrite(fasta, glue("{PATH_WD}/Data/jurkat/map_iso_gene_jurkat.csv"))
}

main()