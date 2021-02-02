### Maud Fagny
### 2021-02-01
### make_enhancers_tfbs_gene_te_table.R
### Cross all available data in terms of TE annotation, enhancers, TFBS annotation and cadidate target genes.
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
### Load data

##Enhancers
enhancers=read.table(paste0(folder, "enhancers.bed"), 
                    header=F, stringsAsFactors = F, sep="\t", quote="")
colnames(enhancers)=c("chr","start", "end", "name", "tissue")
enhancers=data.table(enhancers)
setkey(enhancers, chr, start, end)

##TFBS
tfbs=read.table(paste0(folder, "FIMO/fimo.gff"),
                     header=F, stringsAsFactors = F, sep="\t", quote="")
tfbs=tfbs[, c(1,4,5,9)]
colnames(tfbs)=c("enh.pos", "tfbs.start", "tfbs.end", "tfbs.annot")
tfbs=data.table(tfbs)
tfbs[, enh.pos:=gsub(".*::", "", enh.pos)]
setkey(tfbs, enh.pos, tfbs.start, tfbs.end, tfbs.annot)

## Genes
B73.genes=read.table(paste0(folder, "genes.gff3"), header=F, stringsAsFactors = F, quote="", sep="\t")
load(paste0(folder, "normalized_exprs.Rdata"))

B73.genes=B73.genes[, c(1,4,5,7,9)]
colnames(B73.genes)=c("chr", "gene.start", "gene.end", "gene.strand", "gene.annot")
genes=B73.genes[grep("ID=gene", B73.genes$gene.annot ),]
genes=genes[genes$chr %in% as.character(1:10),]

genenames=unlist(lapply(genes$gene.annot, function(x){strsplit(x, ";")[[1]][1]}))

tss=data.frame("chr"=as.numeric(genes$chr), "start.gene"=ifelse(genes$gene.strand=="+", genes$gene.start, genes$gene.end), 
               "geneID"=gsub("ID=gene:", "", genenames), stringsAsFactors = F)

tss$end.gene=tss$start.gene+1
tss=data.table(tss)
setkey(tss, chr, start.gene, end.gene)

## TE
TE=read.table(paste0(folder, "TransposableElements.gff"), 
              header=F, stringsAsFactors = F, sep="\t", quote="")
colnames(TE)=c("chr", "soft", "data", "TE.start", "TE.end", "score", "strand", "note", "TE.anno")
TE[,c("soft", "data", "score", "strand", "note")]=NULL
TE=data.table(TE)
setkey(TE, chr, TE.start, TE.end)
TE[, TE:=gsub("Motif:", "", unlist(lapply(TE$TE.anno, 
                                          function(x) strsplit(x, '"')[[1]][2])))]
annot.TE=read.table(paste0(folder, "data/annot_TE_names.txt"), 
                    header=F, stringsAsFactors = F, sep="\t", quote="")
colnames(annot.TE)=c("TE", "Superfamily", "consensus")
annot.TE.byconsensus=annot.TE[,2:3]
rownames(annot.TE.byconsensus)=annot.TE$TE
colnames(annot.TE.byconsensus)=c("superfamily", "consensus")
annot.TE.byconsensus$superfamily[annot.TE.byconsensus$superfamily=="DNA"]="TIR"
annot.TE.byconsensus$superfamily[annot.TE.byconsensus$superfamily=="TIR" & 
                                   annot.TE.byconsensus$consensus=="Helitron"]="Helitron"
annot.TE.byconsensus$consensus[annot.TE.byconsensus$superfamily=="Helitron" & 
                                   annot.TE.byconsensus$consensus=="Helitron"]="DHH"

TE[, superfamily:=annot.TE.byconsensus[TE, "superfamily"]]
TE[, family:=annot.TE.byconsensus[TE, "consensus"]]
TE=TE[!is.na(superfamily)]

### Compute proportion of each TE
a=by(data = TE, INDICES = TE$superfamily, FUN = function(x){
  res=as.numeric(by(data =x, INDICES = x$chr,  FUN = function(y){
    ir=IRanges(start=y$start, end=y$end)
    sum(width(IRanges::reduce(ir)))
  }))
  names(res)=sort(unique(x$chr))
  return(res)
})

TE.length=data.frame(matrix(0, nrow=length(unique(TE$chr)), ncol=length(a)))
colnames(TE.length)=names(a)
rownames(TE.length)=unique(TE$chr)
for(n in names(a)){
  TE.length[names(a[[`n`]]),n]=a[[`n`]]
}
tmp=read.table(paste0(folder,"genome.fa.fai"), header=F, stringsAsFactors=F)
rownames(tmp)=tmp$V1
TE.length$genome=tmp[rownames(TE.length), "V2"]
TE.prop.all=round(colSums(TE.length[1:10,c("TIR", "LINE", "LTR", "MITE", "Helitron")])/sum(TE.length$genome), digits = 3)
TE.prop.all=c("all"=sum(TE.prop.all), TE.prop.all)

### Functions
## Extract significant TFBS from TFBS results
find.signif = function(df, threshold){
  df=separate(df,col = "tfbs.annot", 
              into = c("Name", "Alias", "ID", "pvalue", 
                       "qvalue", "sequence", "empty"), sep = ";")
  df[, empty:=NULL]
  df[, pvalue:=as.numeric(gsub("pvalue=", "", pvalue))]
  df[, qvalue:=as.numeric(gsub("qvalue= ", "", qvalue))]
  df=df[qvalue <= threshold,]
  df[, Alias:=gsub("Alias=", "", Alias)]
  df[, enh:=gsub(":.*", "", gsub(".*_", "", df$Name))]
  df[, enh.pos:=gsub("[+-]$", "", gsub(".*::", "", Name))]
  df[,Name:=NULL]
  df[,ID:=NULL]
  df[,pvalue:=NULL]
  df[,qvalue:=NULL]
  df[,sequence:=NULL]
  return(df)
}

### Cross tables

## Merge tfbs
tfbs.signif=find.signif(tfbs)
setkey(tfbs.signif, enh.pos, tfbs.start, tfbs.end, Alias)
tfbs.signif[,tissue:=ifelse((enh.pos %in% enhancers[tissue=="Both"]$enh.pos), "Both", 
                     ifelse((enh.pos %in% enhancers[tissue=="HUSK"]$enh.pos), "Husk", "V2-IST"))]

### Merge tfbs loci
tfbs.unique.locus=foverlaps(tfbs.signif, tfbs.signif, mult ="first", nomatch=NA)
tfbs.unique.locus[,Name:=apply(data.frame(tfbs.unique.locus[,c("enh.pos", "tfbs.start", "tfbs.end", "Alias")], stringsAsFactors = F), 
                               1, function(x){gsub(" ", "", paste(x, collapse="_"))})]
tmp.res=tapply(1:nrow(tfbs.unique.locus), tfbs.unique.locus$Name, function(n, tab){
  tmp=tab[n,c(1,7:ncol(tab))]
  enh.pos=tmp$enh.pos[1]
  tfbs.start=min(tmp$i.tfbs.start)
  tfbs.end=max(tmp$i.tfbs.end)
  Alias=paste(unique(sort(tmp$i.Alias)), collapse=";")
  annot=unique(tmp[, c("i.enh","i.tissue")])
  colnames(annot)=c("enh", "tissue")
  return(data.frame(enh.pos, tfbs.start, tfbs.end, Alias, annot, stringsAsFactors = F))
}, tab=data.frame(tfbs.unique.locus, stringsAsFactors = F))

tfbs.all=NULL
for(i in 1:length(tmp.res)){
  tfbs.all=rbind(tfbs.all, tmp.res[[i]])
}
tfbs.all=data.table(tfbs.all)
tfbs.all[, chr:=as.numeric(gsub(":.*", "", enh.pos))]
tfbs.all[, enh.start:=as.numeric(gsub(".*:", "", gsub("-.*", "", enh.pos)))]
tfbs.all[, enh.end:=as.numeric(gsub(".*-", "", enh.pos))]
tfbs.all[, tfbs.genome.start:=enh.start+tfbs.start-1]
tfbs.all[, tfbs.genome.end:=enh.start+tfbs.end-1]
setkey(tfbs.all, chr, tfbs.genome.start, tfbs.genome.end)

##Overlap TFBS - TE
tfbs.TE=foverlaps(tfbs.all, TE,  nomatch = NA)
#Annotate TE
tfbs.TE[is.na(superfamily) & !is.na(TE), superfamily:=TE]
#Search start and end of TE in TFBS
tfbs.TE[,start.te.tfbs.overlap:=apply(tfbs.TE[,c("tfbs.genome.start","TE.start")], 1, max)]
tfbs.TE[, end.te.tfbs.overlap:=apply(tfbs.TE[,c("tfbs.genome.end","TE.end")], 1, min)]
#If no TE, replace coordinates by 0 and family and superfamily by None
tfbs.TE[is.na(start.te.tfbs.overlap), start.te.tfbs.overlap:=0]
tfbs.TE[is.na(end.te.tfbs.overlap), end.te.tfbs.overlap:=0]
tfbs.TE[is.na(TE.start), TE.start:=0]
tfbs.TE[is.na(TE.end), TE.end:=0]
tfbs.TE[is.na(TE.anno), TE.anno:="None"]
tfbs.TE[is.na(superfamily), superfamily:="None"]
tfbs.TE[is.na(family), family:="None"]
tfbs.TE[is.na(TE), TE:=""]
tfbs.TE[, tfbs.Name:=apply(data.frame(tfbs.TE[,c("chr", "tfbs.genome.start", "tfbs.genome.end", "Alias")], stringsAsFactors = F), 
                      1, function(x){gsub(" ", "", paste(x, collapse="_"))})]
tfbs.TE[,prop.TE:=ifelse(TE != "", ((end.te.tfbs.overlap-start.te.tfbs.overlap)+1)/((tfbs.end-tfbs.start)+1),0)]
tfbs.TE[, length.TE:=ifelse(TE=="", 0, ((end.te.tfbs.overlap-start.te.tfbs.overlap)+1))]
tfbs.TE=unique(tfbs.TE)
multiple.TEs=names(table(tfbs.TE$tfbs.Name)[table(tfbs.TE$tfbs.Name)>1])

tfbs.TE[, start.te.tfbs.overlap:=as.character(start.te.tfbs.overlap)]
tfbs.TE[, end.te.tfbs.overlap:=as.character(end.te.tfbs.overlap)]
tfbs.TE[, TE.start:=as.character(TE.start)]
tfbs.TE[, TE.end:=as.character(TE.end)]
tab=tfbs.TE[tfbs.Name %in% multiple.TEs]

for(x in multiple.TEs){
  tmp=data.frame(tab[tfbs.Name==x], stringsAsFactors = F)
  ir=IRanges::reduce(IRanges(start=as.numeric(tmp$start.te.tfbs.overlap), end=as.numeric(tmp$end.te.tfbs.overlap)))
  tfbs.TE[tfbs.Name==x, TE.start:=gsub(" ", "", 
                                  paste(unique(sort(tmp$TE.start)), collapse=";"))]
  tfbs.TE[tfbs.Name==x, TE.end:=gsub(" ", "", 
                                paste(unique(sort(tmp$TE.end)), collapse=";"))]
  tfbs.TE[tfbs.Name==x, TE.anno:=paste(unique(sort(tmp$TE.anno)), collapse=";")] 
  tfbs.TE[tfbs.Name==x, TE:=gsub(" ", "", paste(unique(sort(tmp$TE)), collapse=";"))]
  tfbs.TE[tfbs.Name==x, superfamily:=gsub(" ", "", 
                                     paste(unique(sort(tmp$superfamily)), collapse=";"))]
  tfbs.TE[tfbs.Name==x, family:=gsub(" ", "", paste(unique(sort(tmp$family)), collapse=";"))]
  tfbs.TE[tfbs.Name==x, start.te.tfbs.overlap:=gsub(" ", "", 
                                       paste(unique(sort(tmp$start.te.tfbs.overlap)), collapse=";"))]
  tfbs.TE[tfbs.Name==x, end.te.tfbs.overlap:=gsub(" ", "", 
                                     paste(unique(sort(tmp$end.te.tfbs.overlap)), collapse=";"))]
  tfbs.TE[tfbs.Name==x, prop.TE:=sum(width(ir))/((tmp$tfbs.end[1]-tmp$tfbs.start[1])+1)]
  tfbs.TE[tfbs.Name==x, length.TE:=ifelse(length(width(ir))>1, 
                                     gsub(" ", "", paste(as.character(width(ir)), collaspe=";")), 
                                     width(ir))]
}
tfbs.TE=unique(tfbs.TE)


### Cross TFBS with genes


tfbs.TE[, enh.window.start:=sapply(enh.start, function(x) max(1, x-window))]
tfbs.TE[, enh.window.end:=enh.end+window]
setkey(tfbs.TE, chr, enh.window.start, enh.window.end)
tfbs.TE.genes=foverlaps(tfbs.TE, tss, nomatch = NA)

enh.genes=unique(tfbs.TE.genes[, c("chr", "geneID", "start.gene", 
                                   "end.gene", "enh.pos", "enh.start", "enh.end",
                                   "enh", "tissue")])   
tfbs.genes=unique(tfbs.TE.genes[, c("chr", "geneID", "start.gene", 
                                   "end.gene", "enh.pos", "enh.start", "enh.end",
                                   "enh", "tissue", 
                                   "Alias", "tfbs.genome.start", "tfbs.genome.end")])   


### Cross enhancers with TE
### Cross enhancers with TE and annotate
annotate.TE.enh=function(enh, TE){
  enh.TE=foverlaps(enh, TE,  nomatch = NA)
  enh.TE[, TE:=gsub("Motif:", "", unlist(lapply(TE.anno, 
                                                   function(x) strsplit(x, '"')[[1]][2])))]
  enh.TE[,TE.enh.start:=apply(enh.TE[, c("TE.start", "enh.start")], 1, max)]
  enh.TE[, TE.enh.end:=apply(enh.TE[,c("TE.end", "enh.end")], 1, min)]
  enh.TE[is.na(TE.enh.start), TE.enh.start:=0]
  enh.TE[is.na(TE.enh.end), TE.enh.end:=0]
  enh.TE[is.na(superfamily), superfamily:="None"]
  enh.TE[is.na(family), family:="None"]
  enh.TE[is.na(TE), TE:=""]
  annotated.enh.TE=tapply(1:nrow(enh.TE), enh.TE$enh.pos, function(x, tab){
    te.length=unlist(tapply(x, tab$superfamily[x], function(y, tmp){
      ir=IRanges(start=tab$TE.enh.start[y][tab$TE.enh.start[y] >0], end=tab$TE.enh.end[y][tab$TE.enh.end[y] >0])
      return(sum(width(IRanges::reduce(ir))))
      
    }, tmp=tab))
    te.family=unlist(tapply(tab$family[x], tab$superfamily[x], 
                            function(x){paste(unique(sort(x)), collapse="/")}))
    te.superfamily=unlist(tapply(tab$superfamily[x], tab$superfamily[x], unique))
    return(data.frame(te.length, te.family, te.superfamily, 
                      stringsAsFactors = FALSE))
  }, tab=data.frame(enh.TE, stringsAsFactors = F))
  res=NULL
  for(i in 1:length(annotated.enh.TE)){
    res=rbind(res, cbind("enh.pos"=rep(names(annotated.enh.TE)[i], nrow(annotated.enh.TE[[i]])), 
                         annotated.enh.TE[[i]]))
  }
  enh$enh.pos %in% 
  for(e in enh$enh.pos)
  return(list(enh.TE, res))
}
TE.enh=annotate.TE.enh(enh=enhancers, TE=TE)
enh.TE=TE.enh[[1]]
annotated.enh.TE=TE.enh[[2]]
save(TE.length, TE.prop.all, tfbs.TE.genes, enh.genes, tfbs.genes, enh.TE, annotated.enh.TE, file=paste0(folder, "crosstables_TE_gene_enhancers_tfbs.Rdata"))
