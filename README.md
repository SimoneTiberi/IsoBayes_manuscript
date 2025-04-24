This repository gathers all the scripts required to reproduce the results (plots and tables) reported in the paper "IsoBayes: a Bayesian approach for single-isoform proteomics inference".
Input data are stored in figshare [here](https://figshare.com/account/home#/projects/183988).

The code is organized in 8 folder:

- Benchmark_results, which contains the R code to process the benchmark results. At the end of the analysis, the folder will store the benchmark results;
- Chapters, which at the end of the analysis will store the main results reported in the aforementioned paper;
- Containers, that stores all the Singularity Recipes to build the containers required to execute the entire analysis;
- Hexbin plots, where are the script and data to obtain the hexbin plots in the manuscript;
- Model_results, which contains the R code to run the IsoBayes model with different parameters and datasets. At the end of the analysis, the folder will store the model results;
- Robustness, that includes the R code to perform the robustness analysis. At the end of the analysis, the folder will store the robustness results;
- Simulation study, which contains the R code to run the simulation study;
- sh_files, which contains the bash files to execute the sequential tasks required by the benchmark and data analysis;
- utils_function, that includes several R files with functions called during the entire analysis.


### Reproducibility
R and IsoBayes library, used for benchmarking and data analysis, have been packaged up in a Singularity container in order to be portable and reproducible. We built the container with IsoBayes version 1.0.0 from [Bioconductor repository](doi:10.18129/B9.bioc.IsoBayes).
Please, install [singularity](https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps) before proceeding.
Once Singularity has been installed, all the results (plots and tables) in the IsoBayes paper can be reproduce by following the next four steps:

1 - download the singularity containers and the input data required to perform the entire analysis:
```shell
git clone https://github.com/SimoneTiberi/IsoBayes_manuscript.git
cd IsoBayes_manuscript
sh sh_files/download_data.sh
```
2 - preprocess the input data in order to get all the required files (.idXML files, mRNA dataset, validation set) to run benchmark and validation analysis:
```shell
mkdir Log_files
sh sh_files/pre_process_data.sh
```
3 - once data have been downloaded and preprocessed, to run the entire benchmark analysis, move the IsoBayes\_manuscript folder on a HPC infrastructure endowed with OpenPBS software and then execute the next two lines:
```shell
sh sh_files/run_benchmark.sh
sh sh_files/internal_main_benchmark.sh
```
4 - finally, with the benchmarked results saved in the Benchmark\_results folder, move the IsoBayes\_manuscript folder on your local operating system and execute:
```shell
sh sh_files/get_models_and_results.sh
```

