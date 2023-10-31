#!/bin/bash

# Chapter 1
mkdir 1
cp ../Benchmark_results/ROC_main_result.png 1/
cp ../Benchmark_results/Average_Memory_usage.png 1/
cp ../Benchmark_results/Average_Run-Time.png 1/
cp ../Benchmark_results/scatterplot_abundance_*.png 1/
cp ../Benchmark_results/change_mrna_prot*.pdf 1/
cp ../Benchmark_results/scatterplot_log2fc*.png 1/

# Chapter 2
mkdir 2
cp ../Benchmark_results/no_UP_ROC_main_result.png 2/
cp ../Benchmark_results/no_UP_scatterplot_abundance_*.png 2/
cp ../Benchmark_results/no_UP_change_mrna_prot*.pdf 2/
cp ../Benchmark_results/no_UP_scatterplot_log2fc*.png 2/


# Chapter 3
mkdir 3
cp ../Robustness/ROC_MM_psm_vs_MM_intensities_vs_OpenMS.png 3/
cp ../Robustness/scatterplot_*MM_psm_vs_MM_intensities_vs_OpenMS.png 3/
cp ../Robustness/change_mrna_prot_*MM_psm_vs_MM_intensities_vs_OpenMS.pdf 3/

# Chapter 4
mkdir 4
cp ../Robustness/ROC_*pep_vs_no_pep.png 4/
cp ../Robustness/scatterplot_*pep_vs_no_pep.png 4/
cp ../Robustness/change_mrna_prot_*pep_vs_no_pep.pdf 4/
