#!/bin/bash
rm Data/wtc11/mrna_isoform.tsv
rm Data/*/map_iso_gene_*
rm Data/*/*/Validation_*
rm Data/*/*/merge_index*
rm Data/*/*/pia*

rm Log_files/*

rm -r Model_results/jurkat
rm -r Model_results/wtc11

rm Benchmark_results/*/*/*.png
rm Benchmark_results/*/*/*.txt
rm Benchmark_results/*/*/*.csv
rm Benchmark_results/*/*/*ROC*
rm Benchmark_results/*/*

rm -r Robustness/jurkat
rm -r Robustness/wtc11
rm Robustness/*.png
rm Robustness/*.pdf

rm -r Chapters/1
rm -r Chapters/2
rm -r Chapters/3
rm -r Chapters/4


