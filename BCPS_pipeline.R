rm(list=ls())
graphics.off()

######################### UPLOAD LIBRARIES ##########################   
library(singscore)
library(stats)
library(Biobase)

######################### UPLOAD DATA ###############################
# load("../total_signatures_list.RData") 	#It uploads a list of signatures ("l.signature" object)
load("../expression_set.RData") 	#Upload the ExpressionSet object of your dataset ("eset" object)

l.signature <- list(BCPS=list(BCPS_UP=c("AP1M2", "CDK5", "PAFAH1B3", "SLC25A10", "SMG5"), BCPS_DOWN=c("CXCL12", "IFFO1", "MFAP4",                                                                                                      "TGFBR2")))

######################### Singscore function ########################   
# The function applyes singscore identifying up-regulated and down-regulated genes.
# It considers a list of signatures and for each of them looks for sublists of 
# up-regulated and down-regulated lists of genes. If they are so it performs 
# singscore considering the directions. If they are not inside the list, the 
# directions are ignored by singscore. Signatures with sublists of up- and down-
# -regolated genes smaller or equal to 3, are ignored.

mySimpleScore<-function(rankData,mysetlist,knownDir=TRUE){
  require(singscore)
  score<-matrix(0, nrow=length(mysetlist),ncol=ncol(rankData))
  rownames(score)<-names(mysetlist)
  colnames(score)<-colnames(rankData)
  for(i in 1:length(mysetlist)){
    if(class(mysetlist[[i]])=="character"){
      sl<-mysetlist[[i]]
      sl<-sl[sl%in%rownames(rankData)]
      if(length(sl)>=3){
        scoretemp<-simpleScore(rankData,upSet = sl, knownDirection = knownDir)
        score[rownames(score)==names(mysetlist)[i],]<-scoretemp$TotalScore
      }
    } else {
      sl.up<-mysetlist[[i]][[grep("UP$",names(mysetlist[[i]]), ignore.case = TRUE)]]
      sl.up<-sl.up[sl.up%in%rownames(rankData)]
      sl.dn<-mysetlist[[i]][[grep("DN$|DOWN$",names(mysetlist[[i]]), ignore.case = TRUE)]]
      sl.dn<-sl.dn[sl.dn%in%rownames(rankData)]
      if(length(sl.up)>=3 & length(sl.dn)>=3){
        scoretemp<-simpleScore(rankData,upSet = sl.up, downSet = sl.dn)
        score[rownames(score)==names(mysetlist)[i],]<-scoretemp$TotalScore
      }
    }
  }
  score<-score[rowSums(score!=0)>0,]
  return(score)
}

######################### Singscore Application ##############################
rankData <- rankGenes(exprs(eset))
result <- mySimpleScore(rankData, l.signature)
# result <- t(result)
result <- as.data.frame(result)
result$sample <- row.names(result)
colnames(result)[1] <- "BCPS"
result <- result[,c(2,1)]

# in "result" object you can find the BCPS evaluated for each sample of your dataset

######################### BCPS correction ######################################   
expression_matrix <- exprs(eset)
patients <- colnames(expression_matrix)
corrected_matrix <- expression_matrix[1:100,]

for (i in row.names(expression_matrix)[1:100]) {
  l<-lm(expression_matrix[i,]~result[patients,"BCPS"])
  corrected_matrix[i,]<-l$residuals 
  rm("l")
}

# "corrected_matrix" contained a gene expression matrix in which each gene was 
# corrected by its correlation with BCPS.