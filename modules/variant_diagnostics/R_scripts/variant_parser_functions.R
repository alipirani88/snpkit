# Forked from: /nfs/esnitkin/Project_MRSA/Analysis/2016-MRSA_colonization_vs_infection/2018-07-17_rewrite_parser_for_SNP_mat_with_newest_matrix_format_v2/lib
##################
# REQUIRE LIBRARIES
###################
suppressMessages(library(seqinr)) # Necessary for converting 3 letter amino acid code to 1 letter amino acid code
suppressMessages(library(Biostrings)) # SNT: Necessary for loading in BLOSUM matrix 
data(BLOSUM80) # Necessary for BLOSUM prediction
suppressMessages(library(magrittr)) # For piping commands 
suppressMessages(library(stringr)) # for counting characters, pipes 

########################
# LIBRARY 
########################

split_any_annotations <- function(variant_matrix, num_of_row_with_multiple_annotations){
  # This function takes in a variant matrix and a row number. 
  # This function outputs a variant matrix. 
  # This function takes in a matrix with a row.name with multiple annotations and: 
  #   1. Grabs all of the annotations within that row.name and formats them into separate annotations
  #   2. Appends rows to the matrix which are duplicates of the row in question.
  #   3. Updates the row.names of the appended rows and the original row in question. 
  
  # Check inputs
  if (class(variant_matrix) != "matrix"){
    stop("Input a variant matrix")
  }
  if (is.na(num_of_row_with_multiple_annotations)){
    return(variant_matrix)
    stop("Called split_any_annotations with no rows with multiple annotations")
  }
  if (class(num_of_row_with_multiple_annotations) != "integer" | num_of_row_with_multiple_annotations < 1 | num_of_row_with_multiple_annotations > nrow(variant_matrix)){
    stop("Input a valid row number")
  } 
  
  
  # Split the current row.name which has multiple annotations
  split_annotations <- strsplit(row.names(variant_matrix)[num_of_row_with_multiple_annotations], ";")
  # Count how many annotations are in this row.name
  num_annotations <- (length(split_annotations[[1]])-1) 
  # Initialize a list, each entry will correspond to a single annotation
  names_for_split_annotations <- list(NA, num_annotations)
  for (k in 1:(num_annotations)){
    start  <- split_annotations[[1]][1] # Grab the name of the variant ("Coding SNP at 60 > G;C)
    end <- split_annotations[[1]][k+1]
    full   <- paste(start, end, sep = ";") # concatenate the start, middle, and end
    names_for_split_annotations[[k]] <- full
  }
  
  # Initialize the output matrix
  variant_matrix_with_annotation_split <- matrix(NA, nrow = (nrow(variant_matrix) + num_annotations - 1), ncol = ncol(variant_matrix))
  # Populate the part of the matrix that will stay the same as the input matrix
  variant_matrix_with_annotation_split[1:nrow(variant_matrix), ] <- variant_matrix
  # Populate the new rows of the matrix (1 new row for each extra annotation)
  variant_matrix_with_annotation_split[(nrow(variant_matrix) + 1):nrow(variant_matrix_with_annotation_split), ] <- matrix(rep(variant_matrix[num_of_row_with_multiple_annotations, ], (num_annotations - 1)), (num_annotations - 1), ncol(variant_matrix), byrow = TRUE)
  # Attach the original matrix row.names for the unchanged part of the matrix and add the annotations you just split
  row.names(variant_matrix_with_annotation_split) <- c(row.names(variant_matrix), unlist(names_for_split_annotations)[2:length(names_for_split_annotations)])
  # Shorten the annotation for the row in question
  row.names(variant_matrix_with_annotation_split)[num_of_row_with_multiple_annotations] <- names_for_split_annotations[[1]]
  return(variant_matrix_with_annotation_split)
} # end split_any_annotations()

#### SNP Parser Function ####

# function to parse snp matrix
# input:
# 1) snpmat - path to snp matrix or snp matrix in matrix form
# output:
# list containing cleaned snpmat and parsed variables
parse_snps = function(snpmat){
  
  # if a path to the snp matrix is given
  if(is.character(snpmat)){
    #using the allele matrix, because eventually we will need to deal with duplicate alleles 
    snpmat <-   read.table(snpmat,
                           header = TRUE,
                           stringsAsFactors = FALSE,
                           sep = "\t",
                           quote = "", 
                           row.names = 1)
    
    # save snpmat as RData file
    #save.image(paste0(format(Sys.time(), "%Y-%m-%d"),'_snpmat.RData'))
  }
  snpmat = snpmat[!is.na(row.names(snpmat)),] #remove blank lines
  
  
  #### SPLIT UP THE MATRIX #### 
  
  # GET ROWS WITH MULTIPLE ALLELES:
  #gets the variant nucleotide, for duplicate alleles there will be a comma separating the dup alleles 
  alleles = strsplit(row.names(snpmat), ';') %>% sapply(., function(x){x[1]}) %>% gsub('functional.*$', '', .) %>% gsub('^.*>', '', .) %>% gsub(' ', '', .)
  
  rows_with_duplicate_alleles = as.integer(grep(',', alleles))
  
  # temporary change - need to fix when we figure out what's going on with these 
  if(length(rows_with_duplicate_alleles) != 0){
    snpmat_less <- snpmat[-rows_with_duplicate_alleles,]
  }else{
    snpmat_less <- snpmat
  }
  
  # KS ADDED 2 LINES: drop rows with "None". Temporary fix. What's with these?
  rows_with_none <- as.integer(grep("None", row.names(snpmat_less)))
  
  # ZL ADDED conditional
  if(length(rows_with_none) > 0){
    snpmat_less <- snpmat_less[-rows_with_none, ]
  }
  
  # SNT added 2018-09-14
  # If there are no pipes in the row, remove that row! 
  rows_with_empty_annots = which(str_count(row.names(snpmat_less), '\\|')==0)
  if(length(rows_with_empty_annots) > 0){
    snpmat_less <- snpmat_less[-rows_with_empty_annots, ]
  }
  
  # SNT added 2018-09-14
  # Might want to change this at some point - removing mutation in the chromosome end
  # because there are no locus tags that match that 
  rows_with_chr_end = grep('CHR_END', row.names(snpmat_less))
  if(length(rows_with_chr_end) > 0){
    snpmat_less <- snpmat_less[-rows_with_chr_end, ]
  }
  
  num_dividers <- sapply(1:nrow(snpmat_less), function(x) lengths(regmatches(row.names(snpmat_less)[x], gregexpr(";", row.names(snpmat_less)[x]))))
  
  # rows_with_multiple_annotations <- c(1:nrow(snpmat_less))[num_dividers > 2]
  # SNT - commented out (Above) and replaced with below

  
  # changed 2018-09-14
  # if there are multipe annotations, there must be more than 9 pipes 
  rows_with_multiple_annotations <- c(1:nrow(snpmat_less))[num_dividers > 2 & str_count(row.names(snpmat_less), '\\|') > 9]
  
  annotations_fixed_less <- as.matrix(snpmat_less)
  #print(head(rownames(annotations_fixed_less)))
  if (length(rows_with_multiple_annotations) != 0){
    for (j in 1:length(rows_with_multiple_annotations)){
      annotations_fixed_less <- split_any_annotations(annotations_fixed_less, rows_with_multiple_annotations[j])
    }
  }
  
  # SNT added 2019-09-14
  # Some annotations have  a semicolon in them, since we are splitting by semicolons, this
  # messes the parsing up - remove the lines that messed up 
  pipe_counts = sapply(row.names(annotations_fixed_less), function(i){str_count(i, '\\|')})
  rows_with_no_pipes = which(pipe_counts ==0)
  if(length(rows_with_no_pipes) > 0){
    annotations_fixed_less <- annotations_fixed_less[-rows_with_no_pipes, ]
  }
  
  # GET FUNCTIONAL ANNOTATION - PHAGE, REPEATS, MASKED 
  flag = strsplit(row.names(annotations_fixed_less), ';') %>% sapply(., function(x){x[1]}) %>% gsub('^.*functional=', '', .)
  
  phage = sapply(strsplit(flag, '_'), function(x){x[1]}) =='PHAGE'
  repeats = sapply(strsplit(flag, '_'), function(x){x[2]}) =='REPEATS'
  # note: MASK might not be the right word but I don't have any in my data set so need to ask Ali 
  masked = sapply(strsplit(flag, '_'), function(x){x[3]}) == 'MASK'
  
  # GET LOCUS TAG
  locus_tag = strsplit(row.names(annotations_fixed_less), ';') %>% sapply(., function(x){x[1]}) %>% gsub('^.*locus_tag=', '', .)
  locus_tag_ig_gene1 = sapply(strsplit(locus_tag, '-'), function(lt){lt[1]})
  locus_tag_ig_gene2 = sapply(strsplit(locus_tag, '-'), function(lt){lt[2]})
  
  # SPLIT UP COMPONENTS 
  annotation_components <- strsplit(row.names(annotations_fixed_less), "\\|")
  
  # GRAB PREDICTION OF FUNCTIONAL IMPACT OF EACH SNP 
  snpeff_prediction <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][3])
  snpeff_low <- snpeff_moderate <- snpeff_high <- snpeff_modifier <- snpeff_prediction
  
  snpeff_low[snpeff_low != "LOW"] = FALSE 
  snpeff_low[snpeff_low == "LOW"] = TRUE
  snpeff_low = as.logical(snpeff_low)
  
  snpeff_moderate[snpeff_moderate != "MODERATE"] <- FALSE
  snpeff_moderate[snpeff_moderate == "MODERATE"] <- TRUE
  snpeff_moderate = as.logical(snpeff_moderate)
  
  snpeff_high[snpeff_high != "HIGH"] <- FALSE
  snpeff_high[snpeff_high == "HIGH"] <- TRUE
  snpeff_high = as.logical(snpeff_high)
  
  snpeff_modifier[snpeff_modifier != "MODIFIER"] <- FALSE
  snpeff_modifier[snpeff_modifier == "MODIFIER"] <- TRUE
  snpeff_modifier = as.logical(snpeff_modifier)
  
  
  # GET SNP's GENOMIC POSITION
  pos <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][1])
  pos <- as.numeric(gsub(".* at (\\d+) >.*", "\\1", pos))
  
  # GET GENE ID 
  # sometimes it's the gene symbol and sometimes its USA300HOU_#### - i need them all in the latter format 
  genes <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][4]) # gene ID
  
  
  
  # GET NUCLEOTIDE CHANGE 
  nuc <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][5])
  nuc <- gsub( "[cn].\\d+", "" , nuc) #eg A>C
  nuc <- gsub(">", "->", nuc) # A->C
  
  
  
  # GRAB SnpEff's DESCRIPTION OF THE VARIANT
  #[1] "synonymous_variant"                                "missense_variant"                                 
  #[3] "intergenic_region"                                 "splice_region_variant&stop_retained_variant"      
  #[5] "intragenic_variant"                                "non_coding_transcript_variant"                    
  #[7] "stop_lost&splice_region_variant"                   "stop_gained"                                      
  #[9] "start_lost"                                        "initiator_codon_variant"                          
  #[11] "initiator_codon_variant&non_canonical_start_codon"
  
  
  var_type <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][2])
  
  s_mut = var_type %in% c('synonymous_variant', 'splice_region_variant&stop_retained_variant')
  
  ns_mut = var_type == 'missense_variant'
  
  intergenic = var_type == 'intergenic_region'
  
  intragenic = var_type == 'intragenic_variant'
  
  stop_lost = var_type == 'stop_lost&splice_region_variant'
  
  stop = var_type == 'stop_gained'
  
  start_lost = var_type == 'start_lost'
  
  # KS ADDED LINE: 
  initiator <- var_type == "initiator_codon_variant"
  
  
  
  # GET GENE ID OF THE GENES FLANKING AN INTERGENIC SNP
  ig_gene1 <- ig_gene2 <- genes
  ig_gene1[intergenic] <- gsub("[-].*", "", genes[intergenic])
  ig_gene2[intergenic] <- gsub(".*[-]", "", genes[intergenic])
  
  
  # GET REF_AA AND VAR_AA FOR S_MUT AND NS_MUT AND STOPS 
  
  aa <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][6])
  
  ref_aa <- aa %>% gsub("p[.]", "",.) %>% gsub("[0-9]+.*", "", .)  
  var_aa <- aa %>% gsub("p[.]", "",.) %>% gsub(".*[0-9]+", "", .) 
  
  ref_aa[ns_mut] <- a(ref_aa[ns_mut])
  var_aa[ns_mut] <- a(var_aa[ns_mut])
  
  ref_aa[s_mut] <- gsub('Ter', 'Stp', ref_aa[s_mut]) %>% a(.)
  var_aa[s_mut] <- gsub('Ter', 'Stp', var_aa[s_mut]) %>% a(.)
  
  ref_aa[stop] <- gsub('Ter', 'Stp', ref_aa[stop]) %>% a(.)
  #var_aa[stop] is already in single amino acid code format 
  
  
  
  # INITIALIZE SIFT_DEL AND PROVEAN_DEL FOR USE OUTSIDE OF SCRIPT
  sift_del <- provean_del <- rep(FALSE, length(ns_mut))
  
  aa_loc <- gsub("[A-z*?.]", "", aa)
  aa_var <- paste(ref_aa, aa_loc, var_aa, sep = "")
  aa <- paste(ref_aa, "->", var_aa, sep = "")
  
  # DEFINE BLOSUM DELETERIOUS MUTATIONS
  del_mut = rep(FALSE, length(ns_mut))
  for (i in 1:length(ref_aa[ns_mut])){
    del_mut[ns_mut][i] <- BLOSUM80[ref_aa[ns_mut][i], var_aa[ns_mut][i]] < 0
  }  
  
  
  # CALCULATE GENE LENGTH IN NUCLEOTIDES 
  gene_length <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][7])
  gene_length <- gsub("[0-9]+/", "", gene_length)
  gene_length <- as.numeric(gene_length)
  
  # GENE SYMBOL 
  gene_symbol <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][9])
  ig_gene1_symbol <- ig_gene2_symbol <- rep(NA, length(gene_symbol))
  ig_gene1_symbol[intergenic] = strsplit(gene_symbol[intergenic], ',') %>% sapply(., function(x){x[1]})
  ig_gene2_symbol[intergenic] = strsplit(gene_symbol[intergenic], ',') %>% sapply(., function(x){x[2]})
  
  # SAVE ANNOTATIONS AS A SEPARATE OBJECT
  annots <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][10])
  
  # label column names of snp matrix
  colnames(annotations_fixed_less) = colnames(snpmat)
  
  # SAVE MATRIX
  # ZL: changed to annotations_fixed_less?
  var_mat_bin <- annotations_fixed_less
  
  
  parsed = list(mat=annotations_fixed_less,
                     phage=phage,
                     repeats=repeats,
                     masked=masked,
                     locus_tag = locus_tag,
                     locus_tag_ig_gene1 = locus_tag_ig_gene1,
                     locus_tag_ig_gene2 = locus_tag_ig_gene2,
                     snpeff_prediction=snpeff_prediction,
                     snpeff_low=snpeff_low,
                     snpeff_moderate=snpeff_moderate,
                     snpeff_high=snpeff_high,
                     snpeff_modifier=snpeff_modifier,
                     pos=pos,
                     genes=genes,
                     nuc=nuc,
                     var_type=var_type,
                     s_mut=s_mut,
                     ns_mut=ns_mut,
                     intergenic=intergenic,
                     intragenic=intragenic,
                     stop_lost=stop_lost,
                     stop=stop,
                     start_lost=start_lost,
                     initiator=initiator,
                     ig_gene1=ig_gene1,
                     ig_gene2=ig_gene2,
                     aa=aa,
                     ref_aa=ref_aa,
                     var_aa=var_aa,
                     sift_del=sift_del,
                     aa_loc=aa_loc,
                     aa_var=aa_var,
                     del_mut=del_mut,
                     gene_length=gene_length,
                     gene_symbol=gene_symbol,
                     ig_gene1_symbol=ig_gene1_symbol,
                     ig_gene2_symbol=ig_gene2_symbol,
                     annots=annots)
  
  save(parsed, file = paste0(format(Sys.time(), "%Y-%m-%d"),'_snpmat_parsed.RData'))
  
  return(parsed)
  
} #end parse_snps

# ZL - not sure what's going on here, so commented it out until we can put it in a function
# # READ IN GENE_LOC FILE 
# gene_info = read.table('./KPNIH1_gene_loc.txt', sep = "\t", header = TRUE, row.names = 1)
# gene_loc_raw = cbind(gene_info[,1] - 27084, gene_info[,2]-27084)
# row.names(gene_loc_raw) = as.matrix(gsub("^USA300_TCH1516_genome:", "", row.names(gene_info), perl=TRUE))
# gene_loc = gene_loc_raw[grep("pUSA", row.names(gene_loc_raw) ,invert = TRUE),];
# 
# # RE-ASSIGN GENE NAME BASED ON FORMAT IN REF GENOME (in my case USA300HOU_####) AND NOT GENE SYMBOL
# 
# 
# # SENSE OF GENES
# # will assign genes with sense once I update the gene name above 
# sense = rep(NA, length(gene_loc[,1]))
# sense[gene_loc[,1] < gene_loc[,2]] = '+'
# sense[gene_loc[,1] > gene_loc[,2]] = '-'


#### Indel Parser Function ####

# 2018-07-30
# KS forked from ST'd indel_parser.R
# ------------------------------------------------------------------------------

# KATIE - can you add a logical to indicate where the split rows are (i.e. the duplicated annotations) -- STEPH: see new variable: second_split_of_annotation

########################
# LIBRARY 
########################

# parse indel matrix function
# input:
# 1) indel matrix
# output:
# list including cleaned indel matrix and parsed information
parse_indels = function(indelmat){
  
  indelmat <- read.table(indelmat,
                         header = TRUE, 
                         stringsAsFactors = FALSE, 
                         sep = "\t",
                         quote = "", 
                         row.names = 1
  )
  
  indelmat = indelmat[!is.na(row.names(indelmat)),] #remove blank lines 
  
  # GET ROWS WITH MULTIPLE ALLELES:
  #gets the variant nucleotide, for duplicate alleles there will be a comma separating the dup alleles 
  alleles = strsplit(row.names(indelmat), ';') %>% sapply(., function(x){x[1]}) %>% gsub('functional.*$', '', .) %>% gsub('^.*>', '', .) %>% gsub(' ', '', .)
  
  rows_with_duplicate_alleles = as.integer(grep(',', alleles))
  
  # temporary change - need to fix when we figure out what's going on with these 
  indelmat_less <- indelmat[-rows_with_duplicate_alleles,]
  
  rows_with_none <- as.integer(grep("None", row.names(indelmat_less)))
  if (length(rows_with_none) > 0){
    stop("There are rows in the matrix with the word 'None'. That means there is a bug.")
  }
  
  #if(length(rows_with_none) > 0){
  #  indelmat_less <- indelmat_less[-rows_with_none, ]    
  #}

  num_dividers <- sapply(1:nrow(indelmat_less), function(x) lengths(regmatches(row.names(indelmat_less)[x], gregexpr(";", row.names(indelmat_less)[x]))))
  rows_with_multiple_annotations <- c(1:nrow(indelmat_less))[num_dividers > 2]
  
  annotations_fixed_less <- as.matrix(indelmat_less)
  original_nrow <- nrow(annotations_fixed_less)
  for (j in 1:length(rows_with_multiple_annotations)){
    annotations_fixed_less <- split_any_annotations(annotations_fixed_less, rows_with_multiple_annotations[j])
  }
  colnames(annotations_fixed_less) = colnames(indelmat)
  
  # STEPH: this is the logical I added to capture where the split annotations end up. 
  # Each time through the loop a new row gets appended to the end of the matrix, 
  # which means that all of the rows between the original nrow and the new nrow
  # are these "duplicates"/second half of the split annotations
  second_split_of_annotation <- c(rep(FALSE, original_nrow), rep(TRUE, (nrow(annotations_fixed_less) - original_nrow)))
  
  # SPLIT UP COMPONENTS 
  annotation_components <- strsplit(row.names(annotations_fixed_less), "\\|")
  
  # GET FUNCTIONAL ANNOTATION - PHAGE, REPEATS, MASKED 
  flag = strsplit(row.names(annotations_fixed_less), ';') %>% sapply(., function(x){x[1]}) %>% gsub('^.*functional=', '', .)
  
  phage = sapply(strsplit(flag, '_'), function(x){x[1]}) =='PHAGE'
  repeats = sapply(strsplit(flag, '_'), function(x){x[2]}) =='REPEATS'
  # note: MASK might not be the right word but I don't have any in my data set so need to ask Ali 
  masked = sapply(strsplit(flag, '_'), function(x){x[3]}) == 'MASK'
  
  # GET LOCUS TAG
  locus_tag = strsplit(row.names(annotations_fixed_less), ';') %>% sapply(., function(x){x[1]}) %>% gsub('^.*locus_tag=', '', .)
  locus_tag_ig_gene1 = sapply(strsplit(locus_tag, '-'), function(lt){lt[1]})
  locus_tag_ig_gene2 = sapply(strsplit(locus_tag, '-'), function(lt){lt[2]})
  
  # GRAB PREDICTION OF FUNCTIONAL IMPACT OF EACH INDEL 
  snpeff_prediction <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][3])
  
  snpeff_low <- snpeff_moderate <- snpeff_high <- snpeff_modifier <- snpeff_prediction
  
  snpeff_low[snpeff_low != "LOW"] = FALSE 
  snpeff_low[snpeff_low == "LOW"] = TRUE
  snpeff_low = as.logical(snpeff_low)
  
  snpeff_moderate[snpeff_moderate != "MODERATE"] <- FALSE
  snpeff_moderate[snpeff_moderate == "MODERATE"] <- TRUE
  snpeff_moderate = as.logical(snpeff_moderate)
  
  snpeff_high[snpeff_high != "HIGH"] <- FALSE
  snpeff_high[snpeff_high == "HIGH"] <- TRUE
  snpeff_high = as.logical(snpeff_high)
  
  snpeff_modifier[snpeff_modifier != "MODIFIER"] <- FALSE
  snpeff_modifier[snpeff_modifier == "MODIFIER"] <- TRUE
  snpeff_modifier = as.logical(snpeff_modifier)
  
  # GET INDEL's GENOMIC POSITION
  pos <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][1])
  pos <- as.numeric(gsub(".* at (\\d+) >.*", "\\1", pos))
  
  # GRAB SnpEff's DESCRIPTION OF THE INDEL
  # [1] "conservative_inframe_deletion"                     
  # [2] "conservative_inframe_insertion"                    
  # [3] "disruptive_inframe_deletion"                       
  # [4] "disruptive_inframe_insertion"                      
  # [5] "frameshift_variant"                                
  # [6] "frameshift_variant&splice_region_variant"          
  # [7] "frameshift_variant&start_lost"                     
  # [8] "frameshift_variant&stop_gained"                    
  # [9] "frameshift_variant&stop_lost&splice_region_variant"
  # [10] "intergenic_region"                                 
  # [11] "intragenic_variant"                                
  # [12] "start_lost&conservative_inframe_deletion"
  
  
  var_type <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][2])
  
  intergenic = var_type == 'intergenic_region'
  
  intragenic = var_type == 'intragenic_variant'
  
  conservative_inframe_deletion = var_type == 'conservative_inframe_deletion'
  
  conservative_inframe_insertion_plus = var_type %in% c('conservative_inframe_insertion', "conservative_inframe_insertion&splice_region_variant")
  
  disruptive_inframe_deletion = var_type == 'disruptive_inframe_deletion'
  
  frameshift_variant = var_type %in% c('frameshift_variant', 
                                       'frameshift_variant&splice_region_variant', 
                                       'frameshift_variant&start_lost', 
                                       'frameshift_variant&stop_gained', 
                                       'frameshift_variant&stop_lost&splice_region_variant', 
                                       "frameshift_variant&splice_region_variant")
  
  start_lost_plus = var_type %in% c("start_lost&conservative_inframe_deletion",
                                    "start_lost&disruptive_inframe_insertion")
  
  gene_fusion_plus = var_type %in% c('bidirectional_gene_fusion', "gene_fusion")
  
  stop_gained_plus = var_type %in% c("stop_gained&conservative_inframe_insertion", 
                                     "stop_gained&disruptive_inframe_insertion")
  
  stop_lost_plus = var_type %in% c("stop_lost&disruptive_inframe_insertion&splice_region_variant", 
                                   "stop_lost&conservative_inframe_deletion&splice_region_variant",
                                   "stop_lost&splice_region_variant")
  
  synonymous = var_type == "synonymous_variant"
  
  disruptive_inframe_insertion = var_type == "disruptive_inframe_insertion"
  
  missense_variant = var_type == "missense_variant"
  
  # GET GENE ID 
  # sometimes it's the gene symbol and sometimes its USA300HOU_#### - i need them all in the latter format 
  genes <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][4]) # gene ID
  
  # GET GENE ID OF THE GENES FLANKING AN INTERGENIC SNP
  ig_gene1 <- ig_gene2 <- genes
  # FIRST BUG
  ig_gene1[intergenic] <- gsub("[-].*", "", genes[intergenic])
  ig_gene2[intergenic] <- gsub(".*[-]", "", genes[intergenic])
  
  
  
  #### GET NUCLEOTIDE CHANGE 
  #"n.62545A>C" - nuc's in this format should be a SNP, it is incorrectly in the indel table 
  # it occurs when there is a SNP and indel at that site, i think?
  nuc <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][5])
  
  #   LOGICALS INDICATING IF A INDEL IS A DELEITION (DEL), DUPLICATION (DUP), OR INSERTION (INS)
  del_dup_ins = gsub('^.*[0-9]','', nuc) %>% gsub('[A-Z].*$','', .)
  del = del_dup_ins == 'del'
  dup = del_dup_ins == 'dup'
  ins = del_dup_ins == 'ins'
  
  # NUCLEOTIDE SEQ OF DEL, DUP, OR INS
  nuc_del_dup_ins =  gsub('^.*[a-z]','',nuc)
  
  #n_or_c = gsub('[.].*$','', nuc)
  
  #nuc[del_dup_ins == 'del']
  #nuc[del_dup_ins == 'ins']
  #nuc[del_dup_ins == 'dup']
  
  # POSITION WHERE THE DELETION, INSERTION, OR DUPLICATION STARTS AND ENDS 
  # NOTE: If the del, ins, or dup is only 1 nuc. long, pos2 will be NA 
  pos_del_dup_ins = gsub('^.*[.]','',nuc) %>% gsub('[a-z].*$','',.) %>% strsplit(., split = '_')
  
  temp1 <- gsub('^.*[.]','',nuc) 
  temp2 <- gsub('[A-Za-z].*$','',temp1) 
  temp3 <- strsplit(temp2, split = '_')
  
  
  
  pos1 = unlist(lapply(pos_del_dup_ins, function(i) i[1]))
  pos2 = unlist(lapply(pos_del_dup_ins, function(i) i[2]))
  
  # LENGTH OF THE DELETION, INSERTION, OR DUPLICATION 
  length_del_dup_ins = nchar(nuc_del_dup_ins)
  
  # CALCULATE GENE LENGTH WHERE INDEL OCCURS IN NUCLEOTIDES 
  gene_length <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][7])
  gene_length <- gsub("[0-9]+/", "", gene_length)
  gene_length <- as.numeric(gene_length)
  
  # GENE SYMBOL 
  gene_symbol <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][9])
  ig_gene1_symbol <- ig_gene2_symbol <- rep(NA, length(gene_symbol))
  ig_gene1_symbol[intergenic] = strsplit(gene_symbol[intergenic], ',') %>% sapply(., function(x){x[1]})
  ig_gene2_symbol[intergenic] = strsplit(gene_symbol[intergenic], ',') %>% sapply(., function(x){x[2]})
  
  # SAVE ANNOTATIONS AS A SEPARATE OBJECT
  annots <- sapply(1:length(annotation_components), function(x) annotation_components[[x]][10])
  
  # SENSE OF GENE (is the gene on the + or - strand?)
  # orientation of the gene (+ or - strand) to be added later, once the gene locus is provided 
  
  parsed = list(mat=annotations_fixed_less,
                phage=phage,
                repeats=repeats,
                masked=masked,
                locus_tag = locus_tag,
                locus_tag_ig_gene1 = locus_tag_ig_gene1,
                locus_tag_ig_gene2 = locus_tag_ig_gene2,
                snpeff_prediction=snpeff_prediction,
                snpeff_low=snpeff_low,
                snpeff_moderate=snpeff_moderate,
                snpeff_high=snpeff_high,
                snpeff_modifier=snpeff_modifier,
                pos=pos,
                genes=genes,
                nuc=nuc,
                var_type=var_type,
                intergenic=intergenic,
                intragenic=intragenic,
                conservative_inframe_deletion=conservative_inframe_deletion,
                conservative_inframe_insertion_plus=conservative_inframe_insertion_plus,
                disruptive_inframe_deletion=disruptive_inframe_deletion,
                frameshift_variant=frameshift_variant,
                start_lost_plus=start_lost_plus,
                gene_fusion_plus=gene_fusion_plus,
                gene_fusion_plus=gene_fusion_plus,
                synonymous=synonymous,
                disruptive_inframe_insertion=disruptive_inframe_insertion,
                missense_variant=missense_variant,
                ig_gene1=ig_gene1,
                ig_gene2=ig_gene2,
                del_dup_ins=del_dup_ins,
                del=del,
                dup=dup,
                ins=ins,
                nuc_del_dup_ins=nuc_del_dup_ins,
                pos_del_dup_ins=pos_del_dup_ins,
                pos1=pos1,
                pos2=pos2,
                length_del_dup_ins=length_del_dup_ins,
                gene_length=gene_length,
                gene_symbol=gene_symbol,
                ig_gene1_symbol=ig_gene1_symbol,
                ig_gene2_symbol=ig_gene2_symbol,
                annots=annots)
  
  save(parsed, file = paste0(format(Sys.time(), "%Y-%m-%d"),'_indelmat_parsed.RData'))
  
  return(parsed)
  
  
  
} #end parse_indels


simplify_snp_code <- function(snp_matrix, keepMQ = FALSE){
  simplified_code_snpmat <- snp_matrix
  simplified_code_snpmat[simplified_code_snpmat == 3] <- 1 # true variant
  simplified_code_snpmat[simplified_code_snpmat == 2] <- 0 # filtered variant
  simplified_code_snpmat[simplified_code_snpmat == -1] <- 0 # unmapped
  simplified_code_snpmat[simplified_code_snpmat == -2] <- 0 # phage
  simplified_code_snpmat[simplified_code_snpmat == -3] <- 0 # lowFQ
  if (keepMQ){
    simplified_code_snpmat[simplified_code_snpmat == -4] <- 1 # keep lowMQ
  } else {
    simplified_code_snpmat[simplified_code_snpmat == -4] <- 0 # remove lowMQ
  }
  # remove columns with no variants (phage)
  #simplified_code_snpmat <- simplified_code_snpmat[rowSums(simplified_code_snpmat) != 0,]
  
  if (sum(sum(simplified_code_snpmat == 1) + sum(simplified_code_snpmat == 0)) != (ncol(simplified_code_snpmat) * nrow(simplified_code_snpmat))){
    stop("snpmat encoded incorrectly.")
  }
  return(simplified_code_snpmat)
} # end simplify_snp_code()
