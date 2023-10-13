#!/bin/bash
PATH_WD=/home/jbollon/prot_iso_mrna_dev/IsoBayes_paper

# Chapter 1
cp $PATH_WD/Benchmark_results/ROC_main_result.png 1/
cp $PATH_WD/Benchmark_results/Average_Memory_usage.png 1/
cp $PATH_WD/Benchmark_results/Average_Run-Time.png 1/
cp $PATH_WD/Benchmark_results/scatterplot_abundance_*.png 1/
cp $PATH_WD/Benchmark_results/change_mrna_prot*.png 1/
cp $PATH_WD/Benchmark_results/scatterplot_log2fc*.png 1/

# Chapter 2
cp $PATH_WD/Benchmark_results/no_UP_ROC_main_result.png 2/
cp $PATH_WD/Benchmark_results/no_UP_scatterplot_abundance_*.png 2/
cp $PATH_WD/Benchmark_results/no_UP_change_mrna_prot*.png 2/
cp $PATH_WD/Benchmark_results/no_UP_scatterplot_log2fc*.png 2/

for data in 'jurkat' 'wtc11'
do
	# Chapter 1
	#cp $PATH_WD/Benchmark_results/$data/ROC_main_result.png 1/$data
	cp $PATH_WD/Benchmark_results/$data/SumTab_main_result.csv 1/$data

	#cp $PATH_WD/Abundance_correlation/$data/scatterplot_OpenMS*.png 1/$data

	#cp $PATH_WD/Change_protein_mRNA_isoform/$data/main_result_*OpenMS.png 1/$data
	#cp $PATH_WD/Change_protein_mRNA_isoform/$data/main_result_*OpenMS_mRNA.png 1/$data
	#cp $PATH_WD/Change_protein_mRNA_isoform/$data/scatterplot_log2FC_OpenMS.png 1/$data
	#cp $PATH_WD/Change_protein_mRNA_isoform/$data/scatterplot_log2FC_OpenMS_mRNA.png 1/$data

	# Chapter 2
	#cp $PATH_WD/Benchmark_results/$data/ROC_main_result_no*.png 2/$data
	cp $PATH_WD/Benchmark_results/$data/no_UP_SumTab_main_result.csv 2/$data

	#cp $PATH_WD/Abundance_correlation/$data/scatterplot_no_unique_OpenMS*.png 2/$data
	#cp $PATH_WD/Abundance_correlation/$data/correlation_no_unique_OpenMS*.csv 2/$data

	#cp $PATH_WD/Change_protein_mRNA_isoform/$data/main_result_*OpenMS*no_unique.png 2/$data
	#cp $PATH_WD/Change_protein_mRNA_isoform/$data/scatterplot_log2FC_OpenMS*no_unique.png 2/$data

	# Chapter 3
	cp $PATH_WD/Robustness/$data/ROC_MM_vs_OpenMS.png 3/$data
	cp $PATH_WD/Robustness/$data/ROC_MM_psm_vs_MM_intensities.png 3/$data

	cp $PATH_WD/Robustness/$data/SumTab_MM_vs_OpenMS.csv 3/$data
	cp $PATH_WD/Robustness/$data/SumTab_MM_psm_vs_MM_intensities.csv 3/$data

	for model in 'MM_psm' 'MM_intensities' 'OpenMS'
	do
		cp $PATH_WD/Abundance_correlation/$data/scatterplot_$model''*.png 3/$data
		cp $PATH_WD/Abundance_correlation/$data/correlation_$model''*.csv 3/$data
	done

	for mrna in '_mRNA' ''
	do
		cp $PATH_WD/Change_protein_mRNA_isoform/$data/OpenMS_vs_MM''$mrna.png 3/$data
		cp $PATH_WD/Change_protein_mRNA_isoform/$data/MM_psm_vs_MM_intensities''$mrna.png 3/$data
		for model in 'MM_psm' 'MM_intensities' 'OpenMS'
		do
			cp $PATH_WD/Change_protein_mRNA_isoform/$data/main_result_extreme_$model''$mrna.png 3/$data
			cp $PATH_WD/Change_protein_mRNA_isoform/$data/scatterplot_log2FC_$model''$mrna.png 3/$data
		done
	done

	# Chapter 4
	for model in 'MM_psm' 'MM_intensities'
	do
		cp $PATH_WD/Robustness/$data/ROC_PEP_vs_no_PEP_mRNA_vs_no_mRNA_$model.png 4/$data
		cp $PATH_WD/Robustness/$data/SumTab_PEP_vs_no_PEP_mRNA_vs_no_mRNA_$model.csv 4/$data
		for pep in '_PEP' ''
		do
			for mrna in '_mRNA' ''
			do
				cp $PATH_WD/Abundance_correlation/$data/scatterplot_$model''$mrna''$pep.png 4/$data
				cp $PATH_WD/Abundance_correlation/$data/correlation_$model''$mrna''$pep.csv 4/$data
				cp $PATH_WD/Change_protein_mRNA_isoform/$data/PEP_no_PEP_$model''$mrna.png 4/$data
				cp $PATH_WD/Change_protein_mRNA_isoform/$data/main_result_extreme_$model''$mrna''$pep.png 4/$data
				cp $PATH_WD/Change_protein_mRNA_isoform/$data/scatterplot_log2FC_$model''$mrna''$pep.png 4/$data
			done
		done
	done
done
