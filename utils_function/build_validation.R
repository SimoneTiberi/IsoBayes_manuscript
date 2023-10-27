build_validation = function(Y, var, var_name = "unit"){
  
  var_unique = sort(unique(var))
  message(glue("Number of {var_name}: {length(var_unique)}"))
  
  Y_unique = sapply(var_unique, function(id){sum(Y[var == id])})
  
  VALIDATION_DF = data.frame(var_name = var_unique, Y_validation = Y_unique,
                             Present = Y_unique > 0,
                             P_Y_validation = Y_unique/sum(Y_unique)
                             )
  colnames(VALIDATION_DF)[1] = var_name
  
  VALIDATION_DF
}