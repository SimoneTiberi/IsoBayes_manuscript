log_output = function(name){
  file.remove(glue("{PATH_WD}/Log_files/{name}.txt"))
  file.create(glue("{PATH_WD}/Log_files/{name}.txt"), overwrite = TRUE)
  zz = file(glue("{PATH_WD}/Log_files/{name}.txt"), open = "w")
  sink(zz, append = TRUE, type = "output")
  sink(zz, append = TRUE, type = "message")
}