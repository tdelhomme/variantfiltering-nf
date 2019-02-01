args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  if(length(unlist(res))==1) res[[1]][2]=""
  return(res)
}
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$help)) {help=FALSE} else {help=TRUE}

if(is.null(args$table) | help) {
  cat("

      Method: add features to each variant of a table to be used in the classifier 

      Mandatory arguments:
      --table=file_name           - annotated table file to add a status
      Optional arguments:
      --remove_bc                 - remove sequencing barcodes in the annotated table (sample name in the form barcode-SM)
      --help                      - print this text

      example: add_calling_features.r --table=test_annotated.txt \n\n")
  q(save="no")
}

if(is.null(args$remove_bc)) {remove_bc=FALSE} else {remove_bc=TRUE}

table = read.table(args$table, quote="\"", stringsAsFactors=F, sep="\t", header=T)

if(remove_bc){
  table$old_SM = table$SM
  table$SM = unlist(lapply(table$SM, function(x) unlist(strsplit(x,"-"))[2]))
  table$SM = gsub("Qc_","",table$SM)
} 

# minimum distance

min_distance <-function(all_my_data, row_number_variant, af_ratio=10){
  cur_data = all_my_data[-row_number_variant,] #remove current variant from distances computations
  cur_data = cur_data[which(cur_data$Chr==all_my_data[row_number_variant,"Chr"] & 
                              cur_data$old_SM==all_my_data[row_number_variant,"old_SM"] &
                              cur_data$AF>=af_ratio*all_my_data[row_number_variant,"AF"]) ,]
  indels = which(as.numeric(cur_data$End) > as.numeric(cur_data$Start))
  min( unlist(lapply(all_my_data[row_number_variant,"Start"]:all_my_data[row_number_variant,"End"], function(var_pos) {
    #create a vector of all position where we observed a variant
    all_pos_sort = sort(as.numeric(unique(c(cur_data$Start, unlist(lapply(indels, function(i) cur_data[i,"Start"]:cur_data[i,"End"]))))))
    min(var_pos - all_pos_sort[which(all_pos_sort<=var_pos)[length(which(all_pos_sort<=var_pos))]] ,
        all_pos_sort[which(all_pos_sort>=var_pos)[1]] - var_pos, na.rm = T) #na if only one variant compared
  })) ) #min to minimum distance over all bases of the variant 
}

table$MIN_DIST = unlist(lapply(1:nrow(table), function(i) log10(min_distance(table, i)) ))

# transitions/transversions ratio

get_TsTv <- function(table, old_SM){
  transitions = c("AG","GA","CT","TC")
  data = table[which(old_SM == table$old_SM),] # get SM data
  data = data[which(data$TYPE_INFO=="snv"),] # get only SNV data
  subst = paste(data$Ref, data$Alt, sep="")
  Ts = length(which(subst %in% transitions)) 
  TsTv = Ts / (length(subst) - Ts)
  return(TsTv)
}

TsTv = unlist(lapply(unique(table$old_SM), function(old_SM){ get_TsTv(table, old_SM) }))
names(TsTv) = unique(table$old_SM)
table$TsTv = as.numeric(TsTv[table$old_SM])

# #SNV/#indels ratio

get_snv_indel_ratio <- function(table, old_SM){
  data = table[which(old_SM == table$old_SM),] # get SM data
  res = length(which(data$TYPE_INFO=="snv")) /   length(which(data$TYPE_INFO=="del" | data$TYPE_INFO=="ins"))
  return(res)
}

snv_indel_ratio = unlist(lapply(unique(table$old_SM), function(old_SM){ get_snv_indel_ratio(table, old_SM) }))
names(snv_indel_ratio) = unique(table$old_SM)
table$snv_indel_ratio = as.numeric(snv_indel_ratio[table$old_SM])

# MaxRatioWin

to_MaxRatioWin <-function(all_my_data, row_number_variant, wind_length=100){
  cur_data = all_my_data[-row_number_variant,] #remove current variant from distances computations
  cur_data = cur_data[which(cur_data$Chr==all_my_data[row_number_variant,"Chr"] & 
                              cur_data$old_SM==all_my_data[row_number_variant,"old_SM"] &
                              cur_data$Start <= (all_my_data[row_number_variant,"Start"] +  wind_length/2) &
                              cur_data$Start >= (all_my_data[row_number_variant,"Start"] -  wind_length/2)),]
  max(cur_data$AF / all_my_data[row_number_variant,"AF"])
}

table$MaxRatioWin = unlist(lapply(1:nrow(table), function(i) to_MaxRatioWin(table, i) ))

# NbVarWin

to_NbVarWin <-function(all_my_data, row_number_variant, wind_length=100){
  cur_data = all_my_data[-row_number_variant,] #remove current variant from distances computations
  cur_data = cur_data[which(cur_data$Chr==all_my_data[row_number_variant,"Chr"] & 
                              cur_data$old_SM==all_my_data[row_number_variant,"old_SM"] &
                              cur_data$Start <= (all_my_data[row_number_variant,"Start"] +  wind_length/2) &
                              cur_data$Start >= (all_my_data[row_number_variant,"Start"] -  wind_length/2)),]
  nrow(cur_data)
}

table$NbVarWin = unlist(lapply(1:nrow(table), function(i) to_NbVarWin(table, i) ))

# HpLength
count_rep = function(region, ref_pos=4){ #in needlestack vcf cont is 3x3 so ref is in pos 4
  right = substr(region, ref_pos, nchar(region))
  left = substr(region, 1, ref_pos);  left = paste(rev(unlist(strsplit(left,NULL))),collapse="")
  count_rep_first_letter(right) + count_rep_first_letter(left) - 1
}

count_rep_first_letter = function(word){
  counts=1
  i=2
  while(i<=nchar(word)){
    if(substr(word,i,i) == substr(word,i-1,i-1)) {counts=counts+1;i=i+1} else {i=nchar(word)+1}
  }
  counts
}

to_HpLength <- function(all_my_data, row_number_variant){
  cont = all_my_data[row_number_variant,"CONT_INFO"]
  cont = gsub("x",all_my_data[row_number_variant,"Ref"],cont)
  count_rep(cont)
}

table$HpLength = unlist(lapply(1:nrow(table), function(i) to_HpLength(table, i) ))

# output

write.table(table, file=paste(paste(gsub(".txt","",args$table),"supp_features",sep="_"),
                              ".txt",sep=""), quote=F, sep="\t", row.names = F)


