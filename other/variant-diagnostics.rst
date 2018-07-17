.. _variant-diagnostics:

Variant Diagnostics
===================

Types of Plot results
---------------------

+--------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Extension  						 | Description                                            							    |
+============+===========================================+==========================================================================================================================================================================================================+
| barplot.pdf 						 | Distribution of filter-pass variant positions(variants observed in all the samples) in each sample. colors represents the filter criteria that caused them to get filtered out in that particular sample.|
+------------+-------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| barplot_DP.pdf 					 | Distribution of filter-pass variant positions in each sample. color represents the read-depth range that they fall in                                                                                 |
+--------------------------------------------------------+---------------------------------------------------------------------------------------------------+
| temp_Only_filtered_positions_for_closely_matrix_FQ.pdf | Heatmap spanning reference genome and shows positions that were filtered out due to low FQ values |
+--------------------------------------------------------+------------------------------------------------------------------+
| DP_position_analysis.pdf 				 | same information as in barplot_DP.pdf but shown in heatmap format|
+--------------------------------------------------------+-------------------------------------------+
| temp_Only_filtered_positions_for_closely_matrix_DP.pdf | Heatmap spanning reference genome and shows positions that were filtered out due to low DP values |
+--------------------------------------------------------+------------+




.. figure::  image/barplot.png
   :align:   center

   Distribution of all the filter-pass variant positions in samples



.. figure::  image/barplot_DP.png
   :align:   center

   Distribution of depth for all the filter-pass variant positions in samples



.. figure::  image/DP_position_analysis.png
   :align:   center

   same information as in barplot_DP.pdf but shown in heatmap format



.. figure::  image/temp_Only_filtered_positions_for_closely_matrix_DP.png
   :align:   center

   Heatmap spanning reference genome and shows positions that were filtered out due to low DP values



.. figure::  image/temp_Only_filtered_positions_for_closely_matrix_FQ.png
   :align:   center

   Heatmap spanning reference genome and shows positions that were filtered out due to low FQ values
