#!/bin/bash

export MAIN_PATH=$(pwd)
export NTHREADS=8
n_models=80
tot_models=0

echo -n > Benchmark_results/internal_list_jobid.txt

for n_run in {1..10}
do
	echo $n_run
	for dt in jurkat wtc11
	do
        	export data=$dt
		export path_to_data=$MAIN_PATH/Data/$data
        	export path_to_res=$MAIN_PATH/Benchmark_results/$data
        	mkdir $path_to_res

        	for file in $path_to_data/Only*
        	do
                	export f=$file
			export NAME_PROTEASE=${f##*Only} # retain the part after "Only"
                	export out=$path_to_res/$NAME_PROTEASE
                	mkdir $out
			sh sh_files/run_internal_IsoBayes.sh
		done
	done
	
	check=0
	tot_models=$(( $n_run * $n_models - 20 ))
	while [ $check -lt $tot_models ]
	do
		echo $tot_models
		sleep 30
		sh sh_files/save_internal_id.sh
		check=$(wc -l Benchmark_results/internal_list_jobid.txt | cut -f1 -d' ')
		echo $check
	done
done

check=0
tot_models=$(( $n_run * $n_models ))
while [ $check != $tot_models ]
do
   sleep 30
   sh sh_files/save_internal_id.sh
   check=$(wc -l Benchmark_results/internal_list_jobid.txt | cut -f1 -d' ')
   echo $check
   echo $tot_models
done

mkdir Benchmark_results/internal_results
while read -r line
do
    tracejob -n 5 "$line" > Benchmark_results/internal_results/"$line"'_res_used.txt'
done < Benchmark_results/internal_list_jobid.txt

mkdir Benchmark_results/internal_results/tmp
mv competitors* Benchmark_results/internal_results/tmp
