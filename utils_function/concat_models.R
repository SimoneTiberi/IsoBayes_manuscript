concat_models = function(df){
  for (i in seq_len(length(df)-1)) {
    df[[1]] = merge(df[[1]], df[[2]], by = "Isoform", all = T)
    df[[2]] = NULL
  }

  na.omit(df[[1]])
}
