# BCPS
The BCPS (Breast Cancer Purity Score) is a transcription-based score with a high correlation with the tumour content in Breast Cancer bulk samples. 
The score showed promising results in different applications:
- In the emulation of tumour purity;
- As a prognostic biomarker;
- In sampling bias correction.

We present a pipeline to evaluate the BCPS and to perform BCPS correction.

# Applications of BCPS
For all the BCPS applications developed so far, you need the language program ```R``` and the packages ```singscore``` (version 1.14.0), ```stats``` (version 4.1.3), and ```Biobase``` (version 2.54.0). 

```
library(singscore)
library(stats)
library(Biobase)
```

We also recommend using the ExpressionSet object to manage the expression profiles of your samples and uploading it on ```R```. The 9-gene set is made of 5 tumour-associated genes and 4 stroma-associated genes and it is managed as a list of lists.
You can also add other gene sets to your list and evaluate them contemporaneously.

```
load("../expression_set.RData") 	#Upload the ExpressionSet object of your dataset

l.signature <- list(BCPS=list(BCPS_UP=c("AP1M2", "CDK5", "PAFAH1B3", "SLC25A10", "SMG5"),
                    BCPS_DOWN=c("CXCL12", "IFFO1", "MFAP4","TGFBR2")))

```
## Evaluation of BCPS
To evaluate BCPS we used a function able to manage the list of genesets and recognize the direction of subsets ("BCPS_UP" and "BCPS_DOWN").


```
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
```

Once the function has been defined, we proceed with the actual calculation:

```
rankData <- rankGenes(exprs(eset))
result <- mySimpleScore(rankData, l.signature)
result <- as.data.frame(result)
result$sample <- row.names(result)
colnames(result)[1] <- "BCPS"
result <- result[,c(2,1)]
```

In the "result" object you can find the BCPS (and other possible genesets) evaluated for each sample of your dataset.

## Sampling bias correction through BCPS
To correct the sampling bias you can use the BCPS evaluated in the previous paragraph to correct the expression of each gene in your expression matrix.

```
expression_matrix <- exprs(eset)
patients <- colnames(expression_matrix)
corrected_matrix <- expression_matrix[1:100,]

for (i in row.names(expression_matrix)[1:100]) {
  l<-lm(expression_matrix[i,]~result[patients,"BCPS"])
  corrected_matrix[i,]<-l$residuals 
  rm("l")
}
```

"corrected_matrix" contained a gene expression matrix in which each gene was corrected by its correlation with BCPS.

# Publication
A manuscript is currently under consideration for publication, to cite currently please refer to the bioRxiv preprint:

###### Add link

# Further information
Created by Dr Marco Barreca at the University of Milano-Bicocca in collaboration with Dr Matteo Dugo at the IRCCS San Raffaele Hospital, under the supervision of Dr Maurizio Callari at Fondazione Michelangelo.


