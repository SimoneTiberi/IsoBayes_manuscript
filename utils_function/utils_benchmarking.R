get_attributes = function(record, begins, ends, isNumeric = FALSE){
  out = vapply(record, function(x){gsub(paste0(".*", begins), "", x)},
               FUN.VALUE = character(1))
  if(isNumeric){
    out = vapply(out, function(x){as.numeric(gsub(paste0(ends, ".*"), "", x))},
                 FUN.VALUE = numeric(1))
  }else{
    out = vapply(out, function(x){gsub(paste0(ends, ".*"), "", x)},
                 FUN.VALUE = character(1))
  }
  names(out) = NULL
  
  out
}

get_score_from_idXML = function(path_to_xml){
  idXML = data.table::fread(path_to_xml, sep = NULL, header = FALSE)
  
  keep = rep(FALSE, nrow(idXML))
  keep = keep + (substr(idXML$V1, 1, 14) == "\t\t\t<ProteinHit")
  idXML = idXML[keep > 0, ]
  
  Isoform = get_attributes(idXML$V1, "accession=\"", "\" score=")
  score = get_attributes(idXML$V1, "score=\"", "\" sequence", isNumeric = TRUE)
  
  data.frame(Isoform, score)
}



load_and_merge = function(path, name_model){
  load(path)
  results$isoform_results = results$isoform_results[, c("proteins", "Probability_present")]
  colnames(results$isoform_results) = paste0(colnames(results$isoform_results), name_model)
  
  merge(benchmark_intersect, results$isoform_results, by.x = "proteins",
        by.y = paste0("proteins", name_model), all = T)
}