#!/bin/bash

nmcmc=2000
path_psm=$f/AllPeptides.psmtsv
path_quant=$f/AllQuantifiedPeptides.tsv

for abundance in 'psm' 'intensities'
do
	if [ "$abundance" = "intensities" ]
	then
		abundance='intensities'
	else
		abundance='psm'
	fi
	for mrna in '_mRNA' ''
	do
		if [ "$mrna" = "_mRNA" ]
		then
			path_mrna=$path_to_data/mrna_isoform.tsv
		else
			path_mrna='NULL'
		fi

		for pep in 'FALSE' 'TRUE'
		do
			if [ "$pep" = "TRUE" ]
			then
				FDR=0.1
			else
				FDR=0.01
			fi	
			echo "#!/bin/bash" >> temp.sh
       			echo "cd "$MAIN_PATH >> temp.sh
       			echo "singularity exec --bind $f Containers/IsoBayes.img Rscript sh_files/run_IsoBayes.R  $path_psm $path_quant $path_mrna 'metamorpheus' $abundance $pep $FDR $nmcmc $NTHREADS 1" >> temp.sh
			echo "echo \$PBS_JOBID >> Benchmark_results/internal_list_jobid.txt" >> temp.sh
       
       			qsub -N competitors_$data'_'$abundance'_pep_'$pep'_'$nmcmc'_'$NAME_PROTEASE''$mrna -j oe -l walltime=00:30:00 -l select=1:ncpus=8:mem=20gb -q cpunodes -m ae temp.sh
       			rm temp.sh
		done	
	done
done
