#!/bin/bash
#rm Data/wtc11/mrna_isoform.tsv
#rm Data/*/map_iso_gene_*
#rm Data/*/*/Validation_*
#rm Data/*/*/merge_index*
#rm Data/*/*/pia*
rm Log_files/*
rm -r Model_results/*/*
rm Model_results/prior_grid
rm Benchmark_results/*/*/*.png
rm Benchmark_results/*/*/*.txt
rm Benchmark_results/*/*/*.csv
rm Benchmark_results/*/*.png
rm Benchmark_results/*/*.csv

rm Change_protein_mRNA_isoform/*/*.png
rm Change_protein_mRNA_isoform/*/*/*.png
rm Change_protein_mRNA_isoform/*/*.csv
rm Change_protein_mRNA_isoform/*/*/*.csv


rm Robustness/*/*.csv
rm Robustness/*/*.png
rm Robustness/*/*/*.csv
rm Robustness/*/*/*.png

rm Abundance_correlation/*/*.csv
rm Abundance_correlation/*/*.png
rm Abundance_correlation/*/*/*.png
