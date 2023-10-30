This repository gathers all the scripts required to reproduce the results (plots and tables) reported in the paper "IsoBayes: a Bayesian approach for single-isoform proteomics inference".
Input data are stored [here](https://figshare.com/account/home#/projects/183988).

### Reproducibility
R and IsoBayes library, used for benchmarking and data analysis, have been packaged up in a Singularity container in order to be portable and reproducible.
Please, install [singularity](https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps) before proceeding.
Once Singularity has been installed, all the results (plots and tables) in the IsoBayes paper can be reproduce by following the next four steps:

1 - download the input data required to perform the entire analysis:
```shell
wget https://github.com/SimoneTiberi/IsoBayes_manuscript.git
cd IsoBayes_manuscript
sh sh_files/download_data.sh
```
2 - preprocess the input data in order to get all the required files (.idXML files, mRNA dataset, validation set) to run benchmark and validation analysis:
```shell
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

