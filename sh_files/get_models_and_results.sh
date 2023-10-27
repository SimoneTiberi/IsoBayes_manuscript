#!/bin/bash
MAIN_PATH=$(pwd)

for data in 'jurkat' 'wtc11'
do
	echo '--- Run IsoBayes models ('$data') ---'
	singularity exec Containers/IsoBayes.img Rscript Model_results/run_models.R $MAIN_PATH $data
	
	if [ -d "Benchmark_results/results" ]; then
		if [ -d "Benchmark_results/internal_results" ]; then
			echo '--- Get benchmark results ---'
        		singularity exec Containers/getResults.img Rscript Benchmark_results/benchmarking_plot.R $MAIN_PATH $data
		fi
	fi

	echo '--- Get robustness results ---'
	singularity exec Containers/getResults.img Rscript Robustness/robustness_IsoBayes.R $MAIN_PATH $data      
done

singularity exec Containers/getResults.img Rscript utils_function/merge_plot.R
