# BCPS
The BCPS (Breast Cancer Purity Score) is a transcription-based score that estimates the tumour content in bulk transcriptomic data from clinical breast cancer samples. 
Here we present a pipeline to compute the BCPS and to adjust gene expression for tumour purity using the BCPS.

## Applications of BCPS
For all the BCPS applications developed so far, the programming language ```R``` and the packages ```singscore``` (version 1.14.0), ```stats``` (version 4.1.3), and ```Biobase``` (version 2.54.0). 

```
library(singscore)
library(stats)
library(Biobase)
```

We also recommend using the ExpressionSet object to handle your dataset of expression profiles and uploading it in ```R```. 
The 9-gene BCPS is made of 5 tumour-associated genes and 4 stroma-associated genes and it is handled as a list of lists.
You can also add other gene sets to your list and evaluate them at the same time.

```
load("../expression_set.RData") 	#Upload the ExpressionSet object of your dataset

l.signature <- list(BCPS=list(BCPS_UP=c("AP1M2", "CDK5", "PAFAH1B3", "SLC25A10", "SMG5"),
                              BCPS_DOWN=c("CXCL12", "IFFO1", "MFAP4","TGFBR2")))

```
### Evaluation of BCPS
To evaluate BCPS we used a function able to manage the list of genesets and recognize the expected direction of expression in case of high tumour purity (i.e. "BCPS_UP" and "BCPS_DOWN").

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

Once the function has been defined, the BCPS can be computed:

```
rankData <- rankGenes(exprs(eset))
result <- mySimpleScore(rankData, l.signature)
result <- as.data.frame(result)
result$sample <- row.names(result)
colnames(result)[1] <- "BCPS"
result <- result[,c(2,1)]
```

In the "result" object you can find the BCPS (and other possible genesets) evaluated for each sample of your dataset.

### Adjusting for sampling bias using the BCPS
To adjust for the sampling bias, you can use the BCPS evaluated in the previous paragraph to adjust the expression of each gene in your expression matrix.

```
expression_matrix <- exprs(eset)
patients <- colnames(expression_matrix)

adjusted_matrix <- apply(expression_matrix, 1, function(x) {
  l <- lm(x~result[patients,"BCPS"])
  x <- l$residuals})
adjusted_matrix <- as.matrix(t(adjusted_matrix))
```

"adjusted_matrix" contains a gene expression matrix in which each gene was adjusted based on its linear relationship with BCPS.

# Publication
A manuscript is currently under consideration for publication, to cite currently please refer to the bioRxiv preprint:

###### Add link

# Further information
Created by Dr Marco Barreca at the University of Milano-Bicocca in collaboration with Dr Matteo Dugo at the IRCCS San Raffaele Hospital, under the supervision of Dr Maurizio Callari at Fondazione Michelangelo and Daniela Besozzi at the University of Milano-Bicocca.

If you need further information, please write an e-mail at: m.barreca@campus.unimib.it.

