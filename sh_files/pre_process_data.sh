#!/bin/bash
MAIN_PATH=$(pwd)

for data in 'jurkat' 'wtc11'
do
        echo '--- Converting mzID files to idXML with OpenMS toolkit ('$data' dataset) ---'
        cd Data/$data
        singularity exec ../../Containers/OpenMS.img sh ../../sh_files/get_openMS_file.sh > ../../Log_files/conversion_idXML_$data.txt 2>&1
        cd ../../

        if [ "$data" = "wtc11" ]
        then
                echo '--- Building mrna_isoform.tsv for wtc11 dataset ---'
                cd Data
                singularity exec ../Containers/python.img python3 accession_mapping.py --gencode_fasta wtc11/gencode_protein.fasta --pacbio_fasta wtc11/wtc11_protein_refined.fasta > ../Log_files/building_mrna_iso_wtc11.txt
                singularity exec ../Containers/IsoBayes.img Rscript get_mrna_isoform_wtc11.R >> ../Log_files/building_mrna_iso_wtc11.txt
                rm accession_map_gencode_uniprot_pacbio.tsv
                cd ../
        fi

        echo '--- Getting '$data' genes ---'
        singularity exec Containers/IsoBayes.img Rscript utils_function/get_gene_$data.R $MAIN_PATH

        echo '--- Building validation set ('$data') ---'
        singularity exec Containers/IsoBayes.img Rscript utils_function/build_validation_set.R $MAIN_PATH $data 
done
