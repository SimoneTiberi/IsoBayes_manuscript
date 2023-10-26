#!/bin/bash

echo "#!/bin/bash" >> temp.sh
echo "cd $MAIN_PATH" >> temp.sh
echo "singularity exec --bind $f Containers/OpenMS.img Epifany -in $f/merge_index_percolator_pep_switched.idXML -threads $NTHREADS -out $out/epifany.idXML > $out/logEpifany.txt" >> temp.sh
echo "echo \$PBS_JOBID >> Benchmark_results/list_jobid.txt" >> temp.sh
qsub -N jobs_log/competitors_$data'_Epifany'_$NAME_PROTEASE -j oe -l walltime=01:00:00 -l select=1:ncpus=10:mpiprocs=10:mem=30gb -q cpunodes -m ae temp.sh
rm temp.sh

