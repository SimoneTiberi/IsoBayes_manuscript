library(IsoBayes)

PEP = as.logical(commandArgs(trailingOnly = TRUE)[6])
FDR_thd = as.numeric(commandArgs(trailingOnly = TRUE)[7])

if(PEP){
  FDR_thd = 0.1
}else{
  FDR_thd = 0.01
}

set.seed(0)
system.time(data_loaded <- load_data(path_to_peptides_psm = commandArgs(trailingOnly = TRUE)[1],
                        path_to_peptides_intensities = commandArgs(trailingOnly = TRUE)[2], 
                        path_to_tpm = commandArgs(trailingOnly = TRUE)[3],
                        input_type = commandArgs(trailingOnly = TRUE)[4], 
                        abundance_type = commandArgs(trailingOnly = TRUE)[5],
                        PEP = PEP,
                        FDR_thd = FDR_thd
			)
)

system.time(res_pack <- inference(data_loaded,
		     K = as.numeric(commandArgs(trailingOnly = TRUE)[8]),
                     n_cores = as.numeric(commandArgs(trailingOnly = TRUE)[9]),
		     thin = as.numeric(commandArgs(trailingOnly = TRUE)[10])
		     )
)
