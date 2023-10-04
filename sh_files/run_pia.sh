#!/bin/bash

echo "#!/bin/bash" >> temp.sh
echo "cd $MAIN_PATH" >> temp.sh
echo "singularity run --bind $f Containers/PIA_tool/pia.sif -c -o $f/pia-compilation.xml $f/merge_index_percolator_pep_switched.idXML > $out/logPiaComp.txt" >> temp.sh
echo "cp Containers/PIA_tool/pia.json $f/pia.json" >> temp.sh
echo "sed -i "s+./piaExport-proteins.mzTab+$out/pia_results.mzTab+g" $f/pia.json" >> temp.sh
echo "singularity run --bind $f Containers/PIA_tool/pia.sif $f/pia.json $f/pia-compilation.xml > $out/logPia.txt" >> temp.sh
echo "echo \$PBS_JOBID >> Benchmark_results/list_jobid.txt" >> temp.sh

qsub -N competitors_$data'_PiaTot'_$NAME_PROTEASE -j oe -l walltime=01:00:00 -l select=1:ncpus=10:mpiprocs=10:mem=30gb -q cpunodes -m ae temp.sh
rm temp.sh
