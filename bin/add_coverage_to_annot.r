args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  if(length(unlist(res))==1) res[[1]][2]=""
  return(res)
}
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$table) | is.null(args$coverage) ) {
  cat("
      Mandatory arguments:
      --table=file_name           - annotated table file to add coverage (.txt)
      --coverage=path             - coverage file (output of mpileup-nf pipeline)
      --nthreads=1                - number of threads used in parallel
      --help                      - print this text
      Optional arguments:
      --remove_bc                 - remove sequencing barcodes in the annotated table (sample name in the form barcode-SM)
      --plot_coverage             - output a PDF for the distribution of coverage of variant in the annotated table (coverages.pdf)

      example: add_coverage_to_annot.r --table=test_annotated.txt --coverage=test_coverage.txt \n\n")
  q(save="no")
}

if(is.null(args$remove_bc)) {remove_bc=FALSE} else {remove_bc=TRUE}
if(is.null(args$plot_coverage)) {plot_coverage=FALSE} else {plot_coverage=TRUE}

var_table = read.table(args$table, quote="\"", stringsAsFactors=F, sep="\t", header=T)

if(remove_bc){
  var_table$old_SM = var_table$SM
  var_table$SM = unlist(lapply(var_table$SM, function(x) unlist(strsplit(x,"-"))[2]))
  var_table$SM = gsub("Qc_","",var_table$SM)
} 

coverage = read.table(args$coverage, quote="\"", stringsAsFactors=F, sep="\t", header=T)

rownames(coverage) = paste(coverage$chr, coverage$pos, sep="-")

var_table$coverage = unlist(lapply(1:nrow(var_table), function(i){
  dat = var_table[i,]
  bc = paste(dat$Chr,dat$Start,sep="-")
  cov=coverage[bc, match(dat$SM, colnames(coverage))]
  if(is.null(cov)) return(NA)
  cov
}))

if(plot_coverage) {
  pdf("coverages.pdf", 6.5, 5)
  hist(var_table$coverage, br=40, col=adjustcolor("cornflowerblue",0.85), xlab="coverage", ylab="counts", main="coverages of annotated table")
  dev.off()
}
  
write.table(var_table, file=paste(gsub(".txt","",args$table),"_with_coverage.txt",sep=""), quote=F, sep="\t", row.names = F)
