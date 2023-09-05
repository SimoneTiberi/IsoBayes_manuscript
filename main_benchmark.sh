#!/bin/bash

export MAIN_PATH=$(pwd)
export NTHREADS=8
n_models=50

echo -n > Benchmark_results/list_jobid.txt

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
		
		sh sh_files/run_fido.sh
		sh sh_files/run_pia.sh
		sh sh_files/run_IsoBayes.sh
		sh sh_files/run_epifany.sh
	done
done

check=0

while [ $check != $n_models ]
do
   sleep 5
   check=$(wc -l Benchmark_results/list_jobid.txt | cut -f1 -d' ')
done

while read -r line
do
    tracejob "$line" > Benchmark_results/"$line"'_res_used.txt'
done < Benchmark_results/list_jobid.txt

