# CPM FILE
mrna = data.table::fread("wtc11/wtc11_orf_refined.tsv")

# conversion file
conv = data.table::fread("accession_map_gencode_uniprot_pacbio.tsv")

# merge and delete
mrna_pac_name = merge(mrna, conv, by.x = "base_acc", by.y = "pacbio_acc")
mrna_pac_name = mrna_pac_name[, c("base_acc", "CPM", "gencode_acc")]

# concatenate gene name under pacbio
sub = mrna_pac_name[mrna_pac_name$gencode_acc != "", 2:3]
sub = sub[, c(2, 1)]
colnames(sub) = colnames(mrna_pac_name)[1:2]
mrna_pac_name$gencode_acc = NULL
mrna_pac_name = rbind(mrna_pac_name, sub)

colnames(mrna_pac_name) = c("isoname", "tpm")
data.table::fwrite(mrna_pac_name, "wtc11/mrna_isoform.tsv")