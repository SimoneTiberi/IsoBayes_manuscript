#!/bin/bash

a=$(find ./ -name "competitors_*")
echo -n > Benchmark_results/internal_list_jobid.txt
for file in $a
do
	id=${file##*.o} # retain the part after ".o"
	echo $id >> Benchmark_results/internal_list_jobid.txt
done
