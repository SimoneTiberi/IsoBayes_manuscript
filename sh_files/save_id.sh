#!/bin/bash

a=$(find ./ -maxdepth 1 -name "competitors_*")
echo -n > Benchmark_results/list_jobid.txt
for file in $a
do
	id=${file##*.o} # retain the part after ".o"
	echo $id >> Benchmark_results/list_jobid.txt
done
