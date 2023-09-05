PATH_WD = commandArgs(trailingOnly = TRUE)[1]

# Set path, global variable and libraries
###########################################################################################
library(glue)
source(glue("{PATH_WD}/utils_function/log_output.R"))

PATH_TO_FASTA = glue("{PATH_WD}/Data/wtc11/isoform.fasta")
PATH_TO_FASTA_PB = glue("{PATH_WD}/Data/wtc11/wtc11_protein_refined.fasta")
log_output("get_gene_wtc11")

# get jurkat gene
###########################################################################################
main = function(){
  fasta = data.table::fread(PATH_TO_FASTA, header = FALSE)
  fasta = fasta[grep(">", fasta$V1), ]
  fasta$isoform = gsub(">", "", fasta$V1)
  fasta$gene = gsub("-.*", "", fasta$isoform)
  fasta$V1 = NULL
  
  fasta_PB = data.table::fread(PATH_TO_FASTA_PB, header = FALSE)
  fasta_PB = fasta_PB[grep(">", fasta_PB$V1), ]
  fasta_PB$V1 = gsub(">pb\\|", "", fasta_PB$V1)
  fasta_PB$isoform = gsub("\\|.*", "", fasta_PB$V1)
  fasta_PB$gene_pb = gsub(".*GN=", "", fasta_PB$V1)
  fasta_PB$V1 = NULL
  
  fasta = merge(fasta, fasta_PB, by = "isoform", all.x = TRUE)
  sel = grep("PB\\.", fasta$isoform)
  fasta$gene[sel] = fasta$gene_pb[sel]
  fasta$gene_pb = NULL
  
  print(glue("Number of genes: {nrow(fasta)}"))
  
  data.table::fwrite(fasta, glue("{PATH_WD}/Data/wtc11/map_iso_gene_wtc11.csv"))
}

main()