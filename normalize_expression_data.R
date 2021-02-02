### Maud Fagny
### 2021-02-01
### normalize_expression_data.R
### Normalize expression data with YARN
###_______________________________

### load libraries
library(RColorBrewer)
library(yarn)
library(corrplot)

### Set variables
args=commandArgs(trailingOnly=TRUE)
folder=args[1]
threshold=args[2]

### Load counts data and sample annotation
counts=read.table(paste0(folder, "expr_raw.txt"),
                  header=T, row.names=1, stringsAsFactors=F)
samples.annot=read.table(paste0(folder, "samples_annotation.txt"),
                         header=T, row.names=1, stringsAsFactors=F)
samples.annot=samples.annot[colnames(counts),]

###Filter counts

#To filter out genes not expressed in at least one tissue
filter=apply(counts, 1, function(x, annot, threshold){
  any(unlist(tapply(x, annot$tissue, function(y){sum(y>0)}))>threshold)
}, annot=samples.annot[colnames(counts),], threshold=threshold)

counts.nonull=counts[filter,]

### Compute Reads per million
RPM=apply(counts.nonull, 2, function(x){x/sum(x)})

### Normalize expression with YARN 
# prepare Esets to use YARN

metadata <- data.frame(labelDescription=c("Tissue", "Replicates", "Batch", "Tissue Color", "Batch Color"), 
                       row.names=c("tissue", "rep", "batch", "tissue.color", "batch.color"))
phenoData <- new("AnnotatedDataFrame", data=samples.annot[colnames(counts.nonull),], varMetadata=metadata)
annotation <- "RNAseq 100bp single-end and 150bp paired-end"

experimentData <- new("MIAME",  name="Maud Fagny", lab="GQE - Le Moulon",
                      contact="maud.fagny@inra.fr", title="Tissue-specific enhancer Experiment")
ESet2 <- ExpressionSet(assayData=as.matrix(counts.nonull), phenoData=phenoData,
                       experimentData=experimentData, annotation="hgu95av2")


# Normalize expression and correct for batch effect
ESet2 <- normalizeTissueAware(ESet2, "tissue", "qsmooth")
norm.ESet2 <- yarn:::extractMatrix(ESet2, normalized = TRUE, log = FALSE)

batch=samples.annot[colnames(counts.nonull),"batch"]
expr=removeBatchEffect(norm.ESet2, batch=factor(batch))

#For counts equal to zero, give minimum value for each individual.
for(i in 1:ncol(expr)){
  expr[counts.nonull[,i]==0,i] <- min(expr[,i])
}

save(expr, file=paste0(folder, "normalized_exprs.Rdata"))
write.table(expr, file=paste0(folder, "exprs.txt"), col.names=F, quote=F, sep="\t")

### Plot PCA of gene expression data
batch.color=unlist(tapply(samples.annot$batch.color, samples.annot$batch, unique))
tissue.color=unlist(tapply(samples.annot$tissue.color, samples.annot$tissue, unique))

#Function
plot.pca <- function(pca, s.pca, i, j, col, pch=16, legend, main=""){
  m <- matrix(c(1,2),nrow = 2,ncol = 1,byrow = TRUE)
  layout(mat = m,heights = c(0.3,0.7))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  par(xpd=TRUE)
  legend("top", legend=names(legend), pch=16, 
         col=legend, bg='white', box.col='white', box.lw=0, box.lty = 1,
         ncol=3)
  par(xpd=NA)
  plot(pca$x[,i], pca$x[,j], 
       xlab=paste0("PC", i, " (", round(s.pca$importance[2,i]*100, digits = 2), "%)"),
       ylab=paste0("PC", j, " (", round(s.pca$importance[2,j]*100, digits = 2), "%)"),
       col=col, pch=pch, main=main)
  lines(c(0,0), par("usr")[3:4], lty=3)
  lines( par("usr")[1:2], c(0,0),lty=3)
  
}

wrap.PCA=function(count.tab, folder, tag, samples.annot, tissue.color, batch.color){
  pca=prcomp(t(count.tab), center=T, scale.=TRUE, retx=T)
  s.pca=summary(pca)
  
  tiff(paste0(folder, "Figures/PCA/", "PCA_", tag, "_variance.tiff"))
  par(mar=c(5,4,3,1)+.1)
  b=barplot(s.pca$importance[2,], col='lightseagreen', 
            ylab="Proportion of variance explained", xlab="Principal component", 
            names=1:ncol(count.tab), main="Eigenvectors")
  lines(b, s.pca$importance[2,], col='blue', lwd=2)
  dev.off()
  
  for(i in c(1,3)){
    j=i+1
    tiff(paste0(folder, "Figures/PCA/", "PCA_", tag, "_", i, "-", j, "_tissue.tiff"))
    par(mar=c(5,4,3,1)+.1)
    plot.pca(pca=pca, s.pca=s.pca, i=i, j=j, 
             col=tissue.color[tissue.corres[samples.annot$tissue]], 
             pch=16, main=paste0("PC", i, " vs PC", j, " - Tissue"), legend=tissue.color)
    dev.off()
    
    tiff(paste0(folder, "Figures/PCA/", "PCA_", tag, "_", i, "-", j, "_center.tiff"))
    par(mar=c(5,4,3,1)+.1)
    plot.pca(pca=pca, s.pca=s.pca, i=i, j=j, col=batch.color[samples.annot$batch], 
             pch=16, main=paste0("PC", i, " vs PC", j, " - Laboratory"), legend=batch.color)
    dev.off()
  }
  return(list("pca"=pca, "summary"=s.pca))
}

dir.create(paste0(folder, "Figures/PCA/"), recursive=T)

tpm.res=wrap.PCA(count.tab=RPM, folder=folder), 
                 tag="rpm", samples.annot=samples.annot, 
                 tissue.color=tissue.col, batch.color=batch.color)
yarn.res=wrap.PCA(count.tab=nmatrix.libscale, folder=folder, 
                 tag="YARN", samples.annot=samples.annot, 
                 tissue.color=tissue.col, batch.color=univbatch.color)
save(tpm.res, yarn.res, file=paste0(folder, "PCA_exprs.Rdata"))

### Plot correlations between counts 
calc.cor <- function(count.tab){
  tmp=cor(count.tab)
  tmp.cs=(tmp-mean(tmp[upper.tri(tmp)]))/sd(tmp[upper.tri(tmp)])
  tmp.cs[tmp.cs>1]=1
  tmp.cs[tmp.cs<(-1)]=-1
  return(tmp.cs)
}
pdf(paste0(folder,  "Figures/correlation_exprs.pdf"))
colnames(tpm)=colnames(expr)=paste(samples.annot$tissue, samples.annot$rep, sep="_")
corrplot(calc.cor(tpm), title = "RPM", tl.cex=0.5, mar=c(0,0,5,0)+.1,
         tl.col=tissue.col[tissue.corres[samples.annot$tissue]])
corrplot(calc.cor(nmatrix.nolibscale), title = "YARN", tl.cex=0.5, mar=c(0,0,5,0)+.1, 
         tl.col=tissue.col[tissue.corres[samples.annot$tissue]])

