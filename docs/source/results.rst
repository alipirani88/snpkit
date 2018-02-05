.. _results:

Results
=======

All the final results will be saved into date_time_core_results directory inside the output folder. 

::

	2018_01_13_13_18_03_core_results
	├── core_snp_consensus
	└── data_matrix

core_snp_consensus folder contains the core consensus fasta and vcf files, whereas data_matrix folder contains all matrices and reports generated during report step.  

As the name suggests, data_matrix will contain various matrices that can be queried or plotted for further diagnosing the variant call results. Alternatively, you can run a R script provided inside the data_matrix folder to generate the plots. 

Require: ggplot2 and heatmap.3

::

	module load R

	Rscript generate_diagnostics_plots.R 



Variant Calling
---------------

Each sample folder under output directory contains standard variant calling outputs such as clean reads, aligned BAM files, filtered non-core vcf and various stats file.

Pending...

| Extension | Description |
| --------- | ----------- |
| . |  |
| . |  |
| . |  |
| . |  |
| . |  |
| . |  |
| . |  |




Core Variants
-------------

The final core step results will be moved to \*_core_results directory under output directory. There are two results folder, core_snp_consensus and data_matrix. The core SNP vcf and consensus fasta files are stored in core_snp_consensus while the data matrices for useful for variant diagnostics can be found in data_matrix.

Pending...

| Extension | Description |
| --------- | ----------- |
| . |  |
| . |  |
| . |  |
| . |  |
| . |  |
| . |  |
| . |  |


Diagnostics plots
-----------------

+------------+------------+
| Extension  | Description|
+============+============+
| barplot.pdf | Distribution of filter-pass variant positions(variants observed in all the samples) in each sample. colors represents the filter criteria that caused them to get filtered out in that particular sample.|
+------------+------------+
| barplot_DP.pdf | Distribution of filter-pass variant positions in each sample. color represents the read-depth range that they fall in.|
+------------+------------+
| temp_Only_filtered_positions_for_closely_matrix_FQ.pdf | Heatmap spanning reference genome and shows positions that were filtered out due to low FQ values |
+------------+------------+
| DP_position_analysis.pdf | same information as in barplot_DP.pdf but shown in heatmap format |
+------------+------------+
| temp_Only_filtered_positions_for_closely_matrix_DP.pdf | Heatmap spanning reference genome and shows positions that were filtered out due to low DP values |
+------------+------------+


Example Plots
-------------


.. figure::  image/barplot.png
   :align:   center

   Distribution of all the filter-pass variant positions in samples


.. figure::  image/barplot_DP.png
   :align:   center

   Distribution of depth for all the filter-pass variant positions in samples
