build_OpenMS_validation = function(benchmark_df, protease){
  # proteine presenti in input
  iso_input = get_score_from_idXML(glue("{PATH_TO_DATA}/Only{protease}/merge_index_percolator_pep_switched.idXML"))
  benchmark_df = merge(benchmark_df, iso_input, by = "Isoform", all = T)
  
  # Eliminio isoforme non presenti nel validation set
  benchmark_df = benchmark_df[!is.na(benchmark_df$Present), ]
  
  # Tengo quelle in input
  benchmark_df = benchmark_df[!is.na(benchmark_df$score), ]
  
  # 0 se modello non trova l'isoforma
  benchmark_df[is.na(benchmark_df)] = 0
  benchmark_df$score = NULL
  
  benchmark_df
}
