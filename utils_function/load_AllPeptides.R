load_AllPeptides = function(path, model){
  files_in_protease = list.files(path, recursive = TRUE)
  if(model == "psm"){
    path_file = files_in_protease[grep("AllPeptides.psmtsv", files_in_protease)]
  }else if(model == "intensities"){
    path_file = files_in_protease[grep("AllQuantifiedPeptides.tsv", files_in_protease)]
  }
  
  as.data.frame(data.table::fread(paste(path, path_file, sep = "/")))
}
