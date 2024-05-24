rm(list = ls())
setwd("~/Desktop/codice_simone")
PATH_WD = "~/Desktop/codice_simone"

# Set path, global variable and libraries
###########################################################################################
library(glue)
library(IsoBayes)
library(ggplot2)

source(glue("{PATH_WD}/utils_function/merge_validation.R"))
source(glue("{PATH_WD}/utils_function/validate_all_protease.R"))
source(glue("{PATH_WD}/utils_function/save_run_inferences.R"))
source(glue("{PATH_WD}/utils_function/get_roc.R"))
source(glue("{PATH_WD}/utils_function/plot_roc_model.R"))
source(glue("{PATH_WD}/utils_function/log_output.R"))

load(glue("{PATH_WD}/utils_function/PALETTE_MODELS"))

###########################################################################################

DATA = "wtc11" # jurkat or wtc11
PEP = TRUE # TRUE or FALSE

PATH_TO_DATA = glue("{PATH_WD}/Data/{DATA}")
PATH_TO_RES = glue("{PATH_WD}/Model_results/{DATA}")
path_to_tpm = paste0(PATH_TO_DATA, "/mrna_isoform.tsv")
proteases = list.dirs(PATH_TO_DATA, recursive = FALSE, full.names = FALSE)
proteases = proteases[grepl(pattern = "Only", proteases)]
proteases = gsub("Only", "", proteases)
proteases

############################################################################################################
# Prior robustness
############################################################################################################
prior_grid = seq(0, 1, 0.1)
save(prior_grid, file = glue("{PATH_WD}/Model_results/prior_grid"))

AUC = matrix(NA, 
             nrow = length(proteases),
             ncol = length(prior_grid))

for (protease in proteases) {
  message(protease)
  
  for (prior in prior_grid){
    message(prior)
    
    if(PEP){
      name = glue("OpenMS_mRNA_PEP_prior_{prior}")
    }else{
      name = glue("OpenMS_mRNA_NoPEP_prior_{prior}")
    }
    message(glue("---------- {name} ----------"))
    name_models = c("IsoBayes_mRNA")
    
    path_to_res_mod = glue("{PATH_TO_RES}/{name}/{protease}")
    if(!dir.exists(path_to_res_mod)){dir.create(path_to_res_mod, recursive = TRUE)}
    SE = IsoBayes::generate_SE(
      path_to_peptides_psm = paste0(PATH_TO_DATA, "/Only",
                                    protease,
                                    "/merge_index_percolator_pep_switched_0.01.idXML"),
      input_type = "openMS",
      abundance_type = "psm",
      PEP = PEP,
      FDR_thd = 0.01)
    
    data_loaded = IsoBayes::input_data(SE, 
                                       path_to_tpm = path_to_tpm)
    
    save(data_loaded, file = glue("{path_to_res_mod}/{name}_data_loaded.RData"))
    
    
    set.seed(123)
    res = IsoBayes::inference(data_loaded,
                              map_iso_gene = NULL, n_cores = 8,
                              prior = prior)
    
    save(res, file = glue("{path_to_res_mod}/{name}_MCMC.RData"))
    
    AUC[which( proteases == protease),
        which( prior_grid == prior)] = plot_roc_model(path_to_res_mod, name, protease, name_models)$AUC
    print(AUC)
  }
}
AUC

AUC = round(rbind(AUC, mean = colMeans(AUC)), 3)

colnames(AUC) = prior_grid
rownames(AUC) = c(proteases, "average")

name = paste0("AUC_", DATA, "_PEP_", PEP,".RData")

save(AUC, file = name)

round(rbind(AUC, mean = colMeans(AUC)), 3)
# complete run: AUC_wtc11_PEP_TRUE.RData
name

# make table of AUCs:
# rows = protease 
# cols = prior
# add jurkat
# try with and without PEP (2 tables)


# [1] "AUC_jurkat_PEP_FALSE.RData"
          0   0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9     1
ArgC    0.800 0.845 0.846 0.846 0.846 0.845 0.845 0.844 0.844 0.844 0.843
AspN    0.821 0.860 0.861 0.860 0.860 0.859 0.858 0.857 0.857 0.857 0.856
Chym    0.785 0.850 0.849 0.850 0.850 0.850 0.849 0.849 0.849 0.849 0.848
GluC    0.816 0.855 0.855 0.855 0.853 0.853 0.852 0.852 0.851 0.850 0.849
LysC    0.842 0.859 0.860 0.860 0.860 0.859 0.859 0.859 0.858 0.857 0.857
Trypsin 0.836 0.861 0.861 0.861 0.860 0.860 0.859 0.859 0.858 0.858 0.858
average 0.817 0.855 0.855 0.855 0.855 0.854 0.854 0.853 0.853 0.852 0.852
mean    0.817 0.855 0.855 0.855 0.855 0.854 0.854 0.853 0.853 0.852 0.852

[1] "AUC_wtc11_PEP_FALSE.RData"
          0   0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9     1
AspN    0.808 0.888 0.889 0.888 0.888 0.887 0.886 0.885 0.884 0.884 0.883
Chym    0.754 0.854 0.854 0.853 0.852 0.851 0.849 0.848 0.847 0.847 0.846
LysC    0.834 0.875 0.876 0.877 0.877 0.876 0.876 0.876 0.875 0.875 0.874
Trypsin 0.841 0.863 0.866 0.867 0.867 0.867 0.867 0.867 0.867 0.867 0.867
average 0.809 0.870 0.871 0.871 0.871 0.870 0.869 0.869 0.868 0.868 0.868

[1] "AUC_jurkat_PEP_TRUE.RData"
          0   0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9     1
ArgC    0.810 0.857 0.859 0.858 0.859 0.859 0.858 0.858 0.857 0.856 0.856
AspN    0.832 0.872 0.873 0.873 0.872 0.871 0.870 0.869 0.868 0.867 0.866
Chym    0.789 0.853 0.855 0.856 0.856 0.856 0.855 0.855 0.854 0.853 0.853
GluC    0.825 0.864 0.865 0.864 0.863 0.862 0.861 0.860 0.859 0.858 0.858
LysC    0.859 0.878 0.880 0.879 0.877 0.877 0.877 0.875 0.875 0.874 0.874
Trypsin 0.850 0.876 0.877 0.877 0.875 0.874 0.874 0.872 0.872 0.871 0.871
average 0.828 0.867 0.868 0.868 0.867 0.866 0.866 0.865 0.864 0.863 0.863
mean    0.828 0.867 0.868 0.868 0.867 0.866 0.866 0.865 0.864 0.863 0.863
          
# [1] "AUC_wtc11_PEP_TRUE.RData"
          0   0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9     1
AspN    0.809 0.887 0.888 0.887 0.885 0.884 0.883 0.882 0.881 0.880 0.879
Chym    0.759 0.860 0.860 0.858 0.857 0.855 0.854 0.853 0.852 0.850 0.850
LysC    0.860 0.899 0.899 0.898 0.897 0.896 0.896 0.894 0.893 0.893 0.892
Trypsin 0.871 0.892 0.893 0.893 0.893 0.892 0.892 0.891 0.891 0.891 0.890
average 0.825 0.884 0.885 0.884 0.883 0.882 0.881 0.880 0.879 0.879 0.878
mean    0.825 0.884 0.885 0.884 0.883 0.882 0.881 0.880 0.879 0.879 0.878