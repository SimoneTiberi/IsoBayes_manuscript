This repository gathers all the scripts required to reproduce the results (plots and tables) reported in the paper "IsoBayes: a Bayesian approach for single-isoform proteomics inference".
Input data are stored here: https://figshare.com/account/home#/projects/183988.

# Reproducibility
R and IsoBayes library, used for benchmarking and data analysis, have been packaged up in a Singularity container in order to be portable and reproducible.
Please, install singularity (https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps) before proceeding.
Once Singularity has been installed, all the results (plots and tables) in the IsoBayes paper can be reproduce by running the main.sh file:

```shell
sh main.sh
```

The main.sh file executes five sequential tasks:

- sh\_files/download\_data.sh to download the input data required to perform the entire analysis reported in IsoBayes paper;
- sh sh\_files/pre\_process\_data.sh to preprocess the input data in order to get all the required files (.idXML files, mRNA dataset, validation set) to run benchmark and validation analysis;
- sh sh\_files/run\_benchmark.sh to compare IsoBayes performance with respect to the major competitors (EPIFANY, Fido, PIA);
- sh sh\_files/internal\_main\_benchmark.sh to benchmark IsoBayes with different inputs and parameters;
- sh sh\_files/get\_models\_and\_results.sh to return the final plots and tables reported in the IsoBayes paper.
 

Please, note thate we carried out the benchmark analysis on a HPC endowed with OpenPBS software for job scheduling and workload management.
The run\_benchmark.sh and internal\_main\_benchmark.sh files must be run on a HPC with OpenPBS commands.
