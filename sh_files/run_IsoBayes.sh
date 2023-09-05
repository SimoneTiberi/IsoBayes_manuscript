#!/bin/bash

nmcmc=2000
	
for pep in 'FALSE' 'TRUE'
do
	echo "#!/bin/bash" >> temp.sh
       	echo "cd "$MAIN_PATH >> temp.sh
       	echo "singularity exec --bind $f Containers/IsoBayes.img Rscript sh_files/run_IsoBayes.R  $f/merge_index_percolator_pep_switched.idXML '' $path_to_data/mrna_isoform.tsv 'openMS' 'psm' $pep 0.01 $nmcmc $NTHREADS 1" >> temp.sh
	echo "echo \$PBS_JOBID >> Benchmark_results/list_jobid.txt" >> temp.sh
       
       	qsub -N competitors_$data'_IsoBayes'_pep_$pep'_'$nmcmc'_'$NAME_PROTEASE -j oe -l walltime=01:00:00 -l select=1:ncpus=10:mpiprocs=10:mem=30gb -q fatnodes -m ae temp.sh
       	rm temp.sh
done	
