args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  if(length(unlist(res))==1) res[[1]][2]=""
  return(res)
}
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$table) | (is.null(args$split_info) & is.null(args$split_geno)) ) {
  cat("
      Mandatory arguments:
      --table=file_name           - annotated table file to add coverage (.txt)

      --split_info                - split the INFO field
      AND/OR
      --split_geno                - split the GENOTYPE field

      Optional argumets:
      --plots                     - plot an histogram for each INFO/GENO fields (INFO_fields.pdf and/or GENO_fields.pdf)
      --help                      - print this text

      example: split_INFO_GENOTYPE.r --table=test_annotated.txt --split_info \n\n")
  q(save="no")
}

if(is.null(args$split_info)) {split_info=FALSE} else {split_info=TRUE}
if(is.null(args$split_geno)) {split_geno=FALSE} else {split_geno=TRUE}
if(is.null(args$plots)) {plots=FALSE} else {plots=TRUE}

table = read.table(args$table, quote="\"", stringsAsFactors=F, sep="\t", header=T)

if(split_info){
  field_value = lapply(table$INFO, function(x) unlist(strsplit(x,";")))
  # get for each line the fields (TYPE, AF, WARN, ...), then take the unique set 
  all_fields = unique(unlist(lapply(field_value, function(l) unlist(lapply(l, function(f) unlist(strsplit(f,"="))[1])))))
  # transform field_value: for each variant line, get each field(name) - value
  # when a field is present in a line and not in another one, this would give a NA
  field_value = lapply(field_value, function(x) {
    res = unlist(lapply(x, function(y) unlist(strsplit(y,"="))[2]))
    names(res) = unlist(lapply(x, function(y) unlist(strsplit(y,"="))[1]))
    res
  })
  for(i in 1:length(all_fields)){
    field = all_fields[i]
    table[,paste(field,"INFO",sep="_")] = unlist(lapply(field_value, function(x) x[field]))
  }
  
  if(plots){
    pdf("INFO_fields.pdf", 6.5, 5)
    for(f in all_fields){
      f = paste(field,"INFO",sep="_")
      if(!suppressWarnings(is.na(as.numeric(table[1,f])))){
        hist(as.numeric(table[,f]), br=20, main=f, ylab="counts", xlab=f, col=adjustcolor("tomato",0.75), border="tomato")
      }
    }
    dev.off()
  }
  table$VF = table$AF # in the case of needlestack, rename AF (allelic frequency) into VF (variant frequency) (cause AF is also in GENO)
}

if(split_geno){
  # get for each line the fields (GT, QVAL, QVAL_INV, ...), then take the unique set 
  all_fields = unique(unlist(lapply(table$FORMAT, function(x) unlist(strsplit(x,":")))))
  # create field_value: for each variant line, get each field(name) - value
  # when a field is present in a line and not in another one, this would give a NA
  field_value = lapply(1:nrow(table), function(i) {
    res = unlist(lapply(table$GENOTYPE[i], function(y) unlist(strsplit(y,":"))))
    names(res) = unlist(lapply(table$FORMAT[i], function(y) unlist(strsplit(y,":"))))
    res
  })
  for(i in 1:length(all_fields)){
    field = all_fields[i]
    table[,all_fields[i]] = unlist(lapply(field_value, function(x) x[field]))
  }
  
  if(plots){
    pdf("GENO_fields.pdf", 6.5, 5)
    for(f in all_fields){
      if(!suppressWarnings(is.na(as.numeric(table[which(!is.na(table[,f]))[1],f])))){ #return false if all NA or if the first non NA is character
        hist(as.numeric(table[,f]), br=20, main=f, ylab="counts", xlab=f, col=adjustcolor("purple3",0.75), border="purple3")
      }
    }
    dev.off()
  }
}

write.table(table, file=paste(paste(gsub(".txt","",args$table), ifelse(split_info,"INFO",""), ifelse(split_geno,"GENO",""),sep="_"),
                              ".txt",sep=""), quote=F, sep="\t", row.names = F)

