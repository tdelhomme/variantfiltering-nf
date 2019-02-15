library(reprtree)
library(caret)
library(randomForest)

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

if(is.null(args$train_table) | help) {
  cat("

      Method: add features to each variant of a table to be used in the classifier

      Mandatory arguments:
      --train_table=file_name     - table output from add_calling_features.r
      --features=AO,DP            - list of features separated by commas
      Optional arguments:
      --output_folder             - folder to save the models (default: .)
      --help                      - print this text

      example: compute_model.r --train_table=table_for_training.txt \n\n")
  q(save="no")
}

if(is.null(args$output_folder)) {output_folder="."} else {output_folder=args$output_folder}
if(is.null(args$features)) {stop("Please provide --features parameter")}

set.seed(123) # here I fix the seed of the random generator to have the same random numbers if re-run

train_table_all = read.table(args$train_table, quote="\"", stringsAsFactors=F, sep="\t", header=T)

train_table_all$IoD = 1 + train_table_all$SIG * 10^(train_table_all$ERR)

# keep only SNVs
for(type in c("snv","indel")){
  if(type=="snv") train_table = train_table_all[which(train_table_all$TYPE_INFO=="snv"),]
  if(type!="snv") train_table = train_table_all[which(train_table_all$TYPE_INFO!="snv"),]

  train_table[which(is.infinite(train_table$MIN_DIST)),"MIN_DIST"] = 1000000000
  train_table[which(is.infinite(train_table$MaxRatioWin)),"MaxRatioWin"] = 1000000000

  propTP = as.numeric(table(train_table$status)["TP"] / nrow(train_table))
  propFP = as.numeric(table(train_table$status)["FP"] / nrow(train_table))

  my_features = unlist(strsplit(args$features,","))

  rf = randomForest(as.factor(status) ~ .,
                    data = train_table[,my_features],
                    importance = TRUE, # to allow us to inspect variable importance
                    ntree = 50000, sampsize = as.numeric(table(train_table$status)["TP"], probs=T) #classwt = c(propTP, propFP), # weighted FP by propTP and TP by propFP
                    # ,maxnodes=10, nodesize=20
  )
  assign(paste("rf_",type,sep=""), rf)
  save(list=paste("rf_",type,sep=""), file = paste(output_folder,"/RF_model_",type,".Rdata",sep=""))
}
