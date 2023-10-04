#!/bin/bash
PATH_WD=/home/jbollon/prot_iso_mrna_dev/IsoBayes_paper

# Chapter 1
cp $PATH_WD/Benchmark_results/jurkat/ROC_main_result.png 1/
cp $PATH_WD/Benchmark_results/jurkat/SumTab_main_result.csv 1/

cp $PATH_WD/Abundance_correlation/jurkat/scatterplot_OpenMS*.png 1/
cp $PATH_WD/Abundance_correlation/jurkat/correlation_OpenMS*.csv 1/

cp $PATH_WD/Change_protein_mRNA_isoform/jurkat/main_result_*OpenMS.png 1/
cp $PATH_WD/Change_protein_mRNA_isoform/jurkat/main_result_*OpenMS_mRNA.png 1/
cp $PATH_WD/Change_protein_mRNA_isoform/jurkat/scatterplot_log2FC_OpenMS.png 1/
cp $PATH_WD/Change_protein_mRNA_isoform/jurkat/scatterplot_log2FC_OpenMS_mRNA.png 1/

# Chapter 2
cp $PATH_WD/Benchmark_results/jurkat/ROC_main_result_no*.png 2/
cp $PATH_WD/Benchmark_results/jurkat/SumTab_main_result_no*.csv 2/

cp $PATH_WD/Abundance_correlation/jurkat/scatterplot_no_unique_OpenMS*.png 2/
cp $PATH_WD/Abundance_correlation/jurkat/correlation_no_unique_OpenMS*.csv 2/

cp $PATH_WD/Change_protein_mRNA_isoform/jurkat/main_result_*OpenMS*no_unique.png 2/
cp $PATH_WD/Change_protein_mRNA_isoform/jurkat/scatterplot_log2FC_OpenMS*no_unique.png 2/

# Chapter 3
cp $PATH_WD/Robustness/jurkat/ROC_MM_vs_OpenMS.png 3/
cp $PATH_WD/Robustness/jurkat/ROC_MM_psm_vs_MM_intensities.png 3/

cp $PATH_WD/Robustness/jurkat/SumTab_MM_vs_OpenMS.csv 3/
cp $PATH_WD/Robustness/jurkat/SumTab_MM_psm_vs_MM_intensities.csv 3/

for model in 'MM_psm' 'MM_intensities' 'OpenMS'
do
	cp $PATH_WD/Abundance_correlation/jurkat/scatterplot_$model''*.png 3/
	cp $PATH_WD/Abundance_correlation/jurkat/correlation_$model''*.csv 3/
done

for mrna in '_mRNA' ''
do
	cp $PATH_WD/Change_protein_mRNA_isoform/jurkat/OpenMS_vs_MM''$mrna.png 3/
	cp $PATH_WD/Change_protein_mRNA_isoform/jurkat/MM_psm_vs_MM_intensities''$mrna.png 3/
	for model in 'MM_psm' 'MM_intensities' 'OpenMS'
	do
		cp $PATH_WD/Change_protein_mRNA_isoform/jurkat/main_result_extreme_$model''$mrna.png 3/
		cp $PATH_WD/Change_protein_mRNA_isoform/jurkat/scatterplot_log2FC_$model''$mrna.png 3/
	done
done

# Chapter 4
for model in 'MM_psm' 'MM_intensities'
do
	cp $PATH_WD/Robustness/jurkat/ROC_PEP_vs_no_PEP_mRNA_vs_no_mRNA_$model.png 4/
	cp $PATH_WD/Robustness/jurkat/SumTab_PEP_vs_no_PEP_mRNA_vs_no_mRNA_$model.csv 4/
	for pep in '_PEP' ''
	do
		for mrna in '_mRNA' ''
		do
			cp $PATH_WD/Abundance_correlation/jurkat/scatterplot_$model''$mrna''$pep.png 4/
			cp $PATH_WD/Abundance_correlation/jurkat/correlation_$model''$mrna''$pep.csv 4/
			cp $PATH_WD/Change_protein_mRNA_isoform/jurkat/PEP_no_PEP_$model''$mrna.png 4/
			cp $PATH_WD/Change_protein_mRNA_isoform/jurkat/main_result_extreme_$model''$mrna''$pep.png 4/
			cp $PATH_WD/Change_protein_mRNA_isoform/jurkat/scatterplot_log2FC_$model''$mrna''$pep.png 4/
		done
	done
done
