#' ---
#' title: "Non-Core Variant Matrix QC"
#' params:
#'   mat: NULL
#' ---

#+ echo=F, warnings=T, message=T
# Input to script: SNP_matrix_code.csv

# Load libraries
library(pheatmap)
library(htmlTable)
# Load functions
source('/nfs/esnitkin/bin_group/pipeline/Github/scripts/variant_parser_functions.R')

# Get command line arguments
#args = commandArgs(trailingOnly = TRUE)

#params = list()
#params$mat = '/scratch/esnitkin_fluxod/apirani/Project_Penn_KPC/Analysis/Regional_KPC_transmission/2018_01_22_Penn_ST258_variant_calling/2018_12_06_14_33_46_core_results/data_matrix/matrices/SNP_matrix_code.csv'

# Read in and parse variant matrix
alt_mat = parse_snps(params$mat)

#' Total number of variants:

#+ echo=F, warnings=F, message=T
nrow(alt_mat$mat)

#' Total number of non-phage variants:
#+ echo=F, warnings=F, message=T
sum(rowSums(alt_mat$mat == -2) == 0)
# htmlTable(matrix(c(-4:3,'filtered (low MQ)','filtered (low FQ)','filtered (phage/repeat region)',
#                  'unmapped','reference allele','core variant','filtered (DP/QUAL)','non-core variant'),
#                  ncol=2),
#           align='c|l',header=c('Number','Description'))

#' What numbers mean: 
#+ echo=F, warnings=F, message=T

#' | Number | Description   |          
#' |:------:|:-------------:|
#' | -4     | filtered (low MQ) | 
#' | -3     | filtered (low FQ) | 
#' | -2     | filtered (phage/repeat region) | 
#' | -1     | unmapped | 
#' | 0      | reference allele |
#' | 1      | core variant  | 
#' | 2      | filtered (DP/QUAL) | 
#' | 3      | non-core variant |

#+ echo=F, warnings=F, message=T
var_type_counts = (sapply(as.data.frame(t(alt_mat$mat)), function(x) table(factor(x,levels=c(-1,0,1,2,3,-2,-3,-4)))))

#' Masked variant positions:
#+ echo=F, warnings=F, message=T
paste('Phage or repeat region:', sum(var_type_counts[6,] != 0), (sum(var_type_counts[6,] != 0)/length(var_type_counts[1,])))
paste('Low FQ:', sum(var_type_counts[7,] != 0),(sum(var_type_counts[7,] != 0)/length(var_type_counts[1,])))
paste('Low MQ:', sum(var_type_counts[8,] != 0),(sum(var_type_counts[8,] != 0)/length(var_type_counts[1,])))
#' Individual variant positions:
#+ echo=F, warnings=F, message=T
paste('Unmapped, filtered, and true:', sum(var_type_counts[1,] != 0 & var_type_counts[4,] != 0), (sum(var_type_counts[1,] != 0 & var_type_counts[4,] != 0)/length(var_type_counts[1,])))
paste('Only unmapped or true:', sum(var_type_counts[1,] != 0 & var_type_counts[4,] == 0),(sum(var_type_counts[1,] != 0 & var_type_counts[4,] == 0)/length(var_type_counts[1,])))
paste('Only filtered or true:', sum(var_type_counts[1,] == 0 & var_type_counts[4,] != 0),(sum(var_type_counts[1,] == 0 & var_type_counts[4,] != 0)/length(var_type_counts[1,])))
paste('Only true:', sum(var_type_counts[3,] > 0),(sum(var_type_counts[3,] > 0)/length(var_type_counts[1,])))

# Heatmap function
heatmap_colors = c('#b35806','#e08214','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788')
allele_heatmap = function(mat,col=heatmap_colors){
  par(mfrow=c(1,1))
  pos_not_filt_frac = colSums(mat == 3 | mat == 1)/(colSums(mat == 2) + colSums(mat == 3 | mat == 1))
  names(pos_not_filt_frac) = colnames(mat)
  genome_filt_count = rowSums(mat == 2)
  names(genome_filt_count) = rownames(mat)
  pheatmap(mat,color=col,
           show_rownames = F, show_colnames = F,
           annotation_row = as.data.frame(genome_filt_count), annotation_col = as.data.frame(pos_not_filt_frac),
           cluster_rows = F, cluster_cols = F)
}


#' Heatmap before gubbins recombination filtering. Phage regions masked. x axis is variants (ordered by position in genome), y axis is genomes. Rows are genomes and columns are variants.
#+ echo=F, warnings=F, message=T
mat = t(alt_mat$mat)
colnames(mat) = rownames(alt_mat$mat)
allele_heatmap(mat)
allele_heatmap(mat[,colSums(mat == -2) == 0],col=heatmap_colors[c(1,2,4:length(heatmap_colors))])

#' Histogram of number genomes with each variant. Phage regions masked.
#+ echo=F, warnings=F, message=T
hist(colSums(mat == 3 | mat == 1),1000,main='',xlab='Number of genomes with variant')

#' Heatmap before gubbins recombination filtering. Phage regions and variants positions filtered because of low FQ masked.
#+ echo=F, warnings=T, message=T
# Mask variants with low FQ (-3) 
mat_maskFQ = mat
mat_maskFQ[,colSums(mat == -3) != 0] = -3
#mat_maskFQ[mat == -4] = -2.5
allele_heatmap(mat_maskFQ)
allele_heatmap(mat_maskFQ[,colSums(mat_maskFQ == -2) == 0 & colSums(mat_maskFQ == -3) == 0],col=heatmap_colors[c(1,4:length(heatmap_colors))])

#' Histogram of number genomes with each variant. Phage regions and variants positions filtered because of low FQ masked.
#+ echo=F, warnings=T, message=T
hist(alt_mat$pos[colSums(mat == -3) != 0],1000,main='',xlab='Low FQ positions across genome')
hist(colSums(mat_maskFQ == 3 | mat_maskFQ == 1),1000,main='',xlab='Number of genomes with variant')

#' Heatmap before gubbins recombination filtering. Phage regions and variants positions filtered because of low FQ and low MQ masked.
#+ echo=F, warnings=T, message=T
# Mask variants with low MQ (-4) 
mat_maskMQ = mat_maskFQ
mat_maskMQ[,colSums(mat == -4) != 0] = -4
allele_heatmap(mat_maskMQ)
allele_heatmap(mat_maskMQ[,colSums(mat_maskMQ == -2) == 0 & colSums(mat_maskMQ == -3) == 0 & colSums(mat_maskMQ == -4) == 0],col=heatmap_colors[4:length(heatmap_colors)])

#' Histogram of number genomes with each variant. Phage regions and variants positions filtered because of low FQ and low MQ masked.
#+ echo=F, warnings=T, message=T
hist(alt_mat$pos[colSums(mat == -4) != 0],1000,main='',xlab='Low MQ positions across genome')
hist(colSums(mat_maskMQ == 3 | mat_maskMQ == 1),1000,main='',xlab='Number of genomes with variant')




