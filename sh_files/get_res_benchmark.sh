#!/bin/bash

for data in jurkat wtc11
do
	path_to_res=/homenfs/jbollon/RNA_PROT/Benchmarking_$data #jurkat #wtc11
	#path_to_res=/homenfs/jbollon/RNA_PROT/Benchmarking_$data'_no_PEP'
	
	for f in $path_to_res/*
	do
		for file in $(find $f -name *hpcfe.txt)
		do
			model=${file%_*}
			id_job=${file##*_} # retain the part after "Only"
			id_job=${id_job%.5*}
			tracejob $id_job -n 3 > $model'_res_used.txt'
		done 
	done
done
