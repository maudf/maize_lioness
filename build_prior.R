### Maud Fagny
### 2021-02-01
### build_prior.R
### Build prior
###_______________________________

### Load libraries
library(data.table)
library(tidyr)
library(purrr)

### Set variables
args=commandArgs(trailingOnly=TRUE)
folder=args[1]
window=args[2]
threshold=args[3]

### Create output folders
dir.create(file.path(folder, "Networks/"), showWarnings = F)
### Load data

##Enhancers
enh=read.table(paste0(folder, "enhancers.bed"), 
                    header=F, stringsAsFactors = F, sep="\t", quote="")

##TFBS
tfbs=read.table(paste0(folder, "FIMO/fimo.gff"),
                     header=F, stringsAsFactors = F, sep="\t", quote="")
tfbs=data.table(tfbs)
tfbs[, V1:=gsub(".*::", "", V1)]
setkey(tfbs, V1, V2, V3, V4, V5, V6, V7, V8,V9)

## Genes
B73.genes=read.table(paste0(folder, "genes.gff3"), header=F, stringsAsFactors = F, quote="", sep="\t")
load(paste0(folder, "normalized_exprs.Rdata"))

genes=B73.genes[grep("ID=gene", B73.genes$V9 ),]
genes=genes[genes$V1 %in% as.character(1:10),]
genenames=unlist(lapply(genes$V9, function(x){strsplit(x, ";")[[1]][1]}))
genes=genes[genenames %in% rownames(expr),]
genenames=genenames[genenames %in% rownames(expr)]

tss=data.frame("chr"=as.numeric(genes$V1), "tss"=ifelse(genes$V7=="+", genes$V4, genes$V5), 
               "geneID"=gsub("ID=gene:", "", genenames), stringsAsFactors = F)
tss=tss[tss$geneID %in% rownames(counts),]

tss$start=tss$tss
tss$end=tss$tss+1
tss=data.table(tss)
setkey(tss, chr, start, end)

### Find common and unique enhancers
enh=data.table(enh)
colnames(enh)[1:3]=c("chr", "start", "end")

### Find limits of enhancers +/- window

enhancers=data.table("chr"=h$chr, 
                          "start"=sapply(h$start, function(x) max(1, x-window)), 
                          "end"=h$end+window)
setkey(enhancers, chr, start, end)


corres.coord.enhancers=paste(enh$chr, paste(enh$start, enh$end, sep="-"), sep=":")
names(corres.coord.enhancers)=paste(enhancers$chr, paste(enhancers$start, enhancers$end, sep="-"), sep=":")


### Extract significant TFBS from TFBS results
find.signif = function(df, threshold){
  df=separate(df,col = "V9", 
              into = c("Name", "Alias", "ID", "pvalue", 
                       "qvalue", "sequence", "empty"), sep = ";")
  df[, empty:=NULL]
  df[, pvalue:=as.numeric(gsub("pvalue=", "", pvalue))]
  df[, qvalue:=as.numeric(gsub("qvalue= ", "", qvalue))]
  df=df[qvalue <= threshold,]
  df[, Alias:=gsub("Alias=", "", Alias)]
  return(df)
}

tfbs.signif=find.signif(tfbs)

### Find neighboring genes (+/- window)

genes.next.enhancers=foverlaps(enhancers, tss, nomatch = 0)
colnames(genes.next.enhancers)=c(colnames(genes.next.enhancers)[1:5], "enhancer.start", "enhancer.end")


genes.next.enhancers[, enhancer.name:=corres.coord.enhancers[paste(chr, paste(enhancer.start, enhancer.end, sep="-"), sep=":")]]

genes.by.enhancers=tapply(genes.next.enhancers$geneID, genes.next.enhancers$enhancer.name, function(x) x)

### Find prior

find.prior <- function(signif, genes.next){
  prior1=unique(unlist(tapply(1:nrow(signif), signif$V1, function(x, enh, gen){
    e=unique(enh$V1[x])
    a=unique(enh$Alias[enh$V1==e])
    b=unique(gen$geneID[gen$enhancer.name==e])
    tmp=expand.grid(a, b, stringsAsFactors=F)
    return(paste(tmp$Var1, tmp$Var2, sep="_"))
  }, enh=signif, gen=genes.next)))
  tmp=matrix(unlist(lapply(prior1, function(x) strsplit(x, "_")[[1]])), ncol=2, byrow=T)
  prior=data.frame("Var1"=tmp[,1], "Var2"=tmp[,2], "edge"=1)
  return(prior)
}
priors=find.prior(signif=tfbs.signif, genes.next=genes.next.enhancers)
save(priors, file=paste0(folder, "Networks/priors.Rdata"))
write.table(priors, file=paste0(folder, "Networks/priors.txt"), col.names=F, row.names=F, quote=F, sep="\t")
