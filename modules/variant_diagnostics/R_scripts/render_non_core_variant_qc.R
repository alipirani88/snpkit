# Get command line arguments
args = commandArgs(trailingOnly = TRUE)

if(length(args) == 2){
  pref=paste0(args[[2]],'_')
}else{
  pref=''
}

rmarkdown::render('non_core_variant_qc.R', 
                  params=list(mat=args[[1]]),
                  output_file=paste0(pref,"non_core_variant_qc.html"))
