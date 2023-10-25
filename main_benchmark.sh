#!/bin/bash

# wget figshare
echo "Downloading Data from.."
path_data=/homenfs/jbollon/RNA_PROT/IsoBayes_paper/Data

export MAIN_PATH=$(pwd)
export NTHREADS=8
n_models=70

echo -n > Benchmark_results/list_jobid.txt


for n_run in {1..2}
do
	for dt in jurkat wtc11
	do
		export data=$dt
		export path_to_data=$path_data/$data
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
	tot_models=$(( $n_run * $n_models - 20 ))
	echo "Number of models executed"
	echo $tot_models

	while [ $check -lt $tot_models ]
	do
		sleep 15
		sh save_id.sh
		check=$(wc -l Benchmark_results/list_jobid.txt | cut -f1 -d' ')
	done
done
			      
check=0
tot_models=$(( $n_run * $n_models ))
while [ $check != $tot_models ]
do
	sleep 30
	sh save_id.sh
	check=$(wc -l Benchmark_results/internal_list_jobid.txt | cut -f1 -d' ')	
done

mkdir Benchmark_results/results

while read -r line
do
    tracejob -n 3 "$line" > Benchmark_results/results/"$line"'_res_used.txt'
done < Benchmark_results/list_jobid.txt
