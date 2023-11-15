concat_models = function(df, union = FALSE){
  for (i in seq_len(length(df)-1)) {
    df[[1]] = merge(df[[1]], df[[2]], by = "Isoform", all = T)
    if(union){
      pres_col = grep("Present", colnames(df[[1]]))[i:(i+1)]
      
      df[[1]][is.na(df[[1]][, pres_col[1]]), pres_col[1]] = df[[1]][is.na(df[[1]][, pres_col[1]]), pres_col[2]]
      df[[1]][is.na(df[[1]][, pres_col[2]]), pres_col[2]] = df[[1]][is.na(df[[1]][, pres_col[2]]), pres_col[1]]
      
      df[[1]][is.na(df[[1]])] = 0
    }else{
      df[[1]] = na.omit(df[[1]])
    }
    
    df[[2]] = NULL
  }
  
  df[[1]]
}
