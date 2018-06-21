library(VariantAnnotation)
args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  if(length(unlist(res))==1) res[[1]][2]=""
  return(res)
}
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$VCFToAddCoverage) | is.null(args$coverage) ) {
  cat("
      Mandatory arguments:
      --VCFToAddCoverage=file     - VCF file to add coverage
      --coverage=path             - coverage file (output of mpileup-nf pipeline)
      --help                      - print this text
      Optional arguments:
      --remove_bc                 - remove sequencing barcodes in the annotated table (sample name in the form barcode-SM)
      --output_vcf                - output VCF name
      example: add_coverage_to_annot.r --VCFToAddCoverage=myfile.vcf --coverage=mycoverage.txt \n\n")
  q(save="no")
}

if(is.null(args$remove_bc)) {remove_bc=FALSE} else {remove_bc=TRUE}
vcf = args$VCFToAddCoverage
if(is.null(args$output_vcf)) {output_vcf=gsub(".vcf.gz","_addCoverage.vcf",vcf)} else {output_vcf=args$output_vcf}

coverage_dat = read.table(args$coverage, quote="\"", stringsAsFactors=F, sep="\t", header=T)
rownames(coverage_dat) = paste(coverage_dat$chr, coverage_dat$pos, sep = "-")

#initiate the first chunk
vcf_con <- open(VcfFile(vcf,  yieldSize=10000))
vcf_chunk = readVcf(vcf_con, "hg19")

#and continue
while(dim(vcf_chunk)[1] != 0) {
  
  #annotate the header of the chunk in first to don't have warnings after
  geno(header(vcf_chunk))["addCov",]=list("1","Integer","Added coverage from external table")
  
  for(i in 1:nrow(vcf_chunk)){
    print(i)
    dat = vcf_chunk[i,]
    bc = paste(as.character(seqnames(rowRanges(dat))), start(ranges(rowRanges(dat))), sep="-")
    samples = colnames(geno(dat)[[1]])
    if(remove_bc) samples = unlist(lapply(samples, function(x) unlist(strsplit(x,"-"))[2]))
    covs = matrix(unlist(lapply(samples, function(s){
      c = coverage_dat[bc,s]
      if(is.null(c)){return(NA)} else {return(c)}
    })), nrow=1)
    rownames(covs) = rownames(geno(dat)[[names(geno(dat))[1]]])
    colnames(covs) = colnames(geno(dat)[[names(geno(dat))[1]]])
    geno(vcf_chunk[i,])[["addCov"]] = covs
    names(geno(vcf_chunk[i])[[length(geno(vcf_chunk[i]))]]) = "addCov"
  }
  
  #write out the annotated VCF
  con = file(output_vcf, open = "a")
  writeVcf(vcf_chunk, con)
  vcf_chunk = readVcf(vcf_con, "hg19")
  close(con)
}
