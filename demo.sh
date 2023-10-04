#!/bin/bash

check=0
n_models=50
echo $check
while [ $check != $n_models ]
do
	echo $check
	sleep 5
	check=$(wc -l Benchmark_results/list_jobid.txt | cut -f1 -d' ')
done

while read -r line
do
	tracejob "$line" -n 2 > Benchmark_results/"$line"'_res_used.txt'
done < Benchmark_results/list_jobid.txt


