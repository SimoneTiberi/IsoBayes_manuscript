#!/bin/bash

NTHREADS=4
FASTA=isoform.fasta
DECOY_STR='DECOY_'
#FDR=0.01

for FDR in 0.1 # 0.01
do
	for f in ./Only*
	do
        if [ -d $f ]; then
                echo '------------------------------------'
                echo 'Processing '$f' ...'
                echo '------------------------------------'

                for mz in $f/Individual_File_Results/*.mzID
                do
                        IDFileConverter -in $mz -threads $NTHREADS -out $mz.idXML
                done

                IDMerger -in $f/Individual_File_Results/*.idXML -threads $NTHREADS -merge_proteins_add_PSMs -out $f/merge.idXML
                rm $f/Individual_File_Results/*.idXML

                ENZYME=${f##*/}  # retain the part after the last "/"
                ENZYME="$(grep -Po $ENZYME',\K[^/]*' ../conversion_enzyme_index.csv)"
                PeptideIndexer -in $f/merge.idXML -enzyme:name $ENZYME -threads $NTHREADS -decoy_string_position prefix -decoy_string $DECOY_STR -fasta $FASTA -out $f/merge_index.idXML
                rm $f/merge.idXML

                ENZYME=${f##*/}  # retain the part after the last "/"
                ENZYME="$(grep -Po $ENZYME',\K[^/]*' ../conversion_enzyme.csv)"
                PercolatorAdapter -in $f/merge_index.idXML -enzyme $ENZYME -threads $NTHREADS -generic_feature_set -score_type pep -out $f/merge_index_percolator_pep.idXML

                FalseDiscoveryRate -in $f/merge_index_percolator_pep.idXML -out $f/merge_index_percolator_pep_$FDR.idXML -protein false -threads $NTHREADS -FDR:PSM $FDR -algorithm:add_decoy_peptides -algorithm:add_decoy_proteins
                IDScoreSwitcher -in $f/merge_index_percolator_pep_$FDR.idXML -out $f/merge_index_percolator_pep_switched_$FDR.idXML -new_score 'Posterior Error Probability_score' -new_score_orientation lower_better -new_score_type pep -threads $NTHREADS
                rm $f/merge_index_percolator_pep_$FDR.idXML
        fi
	done
done
