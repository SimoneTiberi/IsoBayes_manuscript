# Reproducibility

R and IsoBayes library, used for benchmarking and data analysis, have been packaged up in a Singularity container in order to be portable and reproducible.
All the results (plots and tables) in the paper "x" can be reproduce by launching the following code

```shell
sh main.sh
sh main_benchmark.sh
```

The main.sh file automatically downloads the input data from figshare(link) and returns the results reported in the paper.
Then, the main_benchmark.sh file performs the benchmark analysis.

Please, note thate we carried out the benchmark analysis on a HPC endowed with OpenPBS for job scheduling and workload management.
The main_benchmark.sh file must be run on a HPC with OpenPBS commands.
