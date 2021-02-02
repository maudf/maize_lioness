### Maud Fagny
### 2021-02-01
### find_enhancer_topgenes_pairs.R
### Identify modst likely target gene for each enhancer.
###_______________________________

### Load libraries
library(data.table)

### Set variables
args=commandArgs(trailingOnly=TRUE)
folder=args[1]

### Load data

load(paste0(folder,"output_lioness.Rdata"))
load(paste0(folder,"Networks/priors.Rdata"))
samples.annot=read.table(paste0(folder, "samples_annotation.txt"),
                         header=T, row.names=1, stringsAsFactors=F)
load(paste0(folder, "logtransformed_lioness_results.Rdata"))
load( paste0(folder,"crosstables_TE_gene_enhancers_tfbs.Rdata"))

tis=samples.annot$tissue

### Extract enhancers with transposable elements
enh=sort(unique(c(enh.genes$enh)))

### Compute targeting for all enhancer-gene pairs and pick highest target

enh.gene.pairs=t(apply(data.frame(enh, stringsAsFactors=F), 1, function(x, pairs, edges){
  print(x)
  tmp=pairs[pairs$enh==x,]
  g=unique(tmp$geneID)
  if( length(g)>0 & !all(is.na(g)) ){
      score.genes=t(apply(data.frame(g, stringsAsFactors=F), 1, function(y, tfbs, edges, gene){
          tfs=tfbs$Alias[tfbs$geneID==y]
          res=NULL
          for(tf in tfs){
              alltfs=unlist(strsplit(tf, ";"))
              e=edges[paste(y, alltfs, sep="_"), ]
              res=rbind(res, colSums(e,na.rm=T)/nrow(e))
          }
          if(is.null(res)){return(rep(0,12))} else {return(colSums(res)/nrow(res))}
      }, tfbs=tmp, edges=edges))
      rownames(score.genes)=g
      if(nrow(score.genes)==1){
        avg.husk=mean(score.genes[,1:6])
        avg.V2=mean(score.genes[,7:12])
        names(avg.husk)=names(avg.V2)=g
      } else {
        avg.husk=rowSums(score.genes[,1:6])/ncol(score.genes[,1:6])
        avg.V2=rowSums(score.genes[,7:12])/ncol(score.genes[,7:12])
      }
      out=data.frame("enh"=rep(x, length(avg.husk)), "gene"=names(avg.husk), 
                       "edge.husk"=avg.husk, "top.husk"=ifelse(avg.husk==max(avg.husk),1,0), 
                       "edge.V2"=avg.V2, "top.V2"=ifelse(avg.V2==max(avg.V2),1,0),
                       stringsAsFactors = F)
   

    } else {
      out=data.frame("enh"=x, "gene"=NA, "edge.husk"=NA, "top.husk"=NA, 
                     "edge.V2"=NA, "top.V2"=NA, stringsAsFactors = F)
   }
  return(out)
  },pairs=tfbs.genes, edges=edges[, colnames(exprs)]))

##Make result table
res=NULL
for(i in 1:length(enh.gene.pairs)){res=rbind(res, enh.gene.pairs[[i]][[1]]); print(i)}
### Extract enhancer-genes pairs 
## Annotate enhancers with position and if MITE, whether in Pif or not
pos.enh=unique(enh.genes[, c("chr", "enh.start", "enh.end", "enh")])
pos.enh[,enh:=ifelse(!is.na(enh.V2), enh.V2,enh.husk)]
setkey(pos.enh,enh)
setkey(pos.enh,enh)
res=data.table(res)
setkey(res, enh)
enh.gene.pairs=merge(res, pos.enh)
save(enh.gene.pairs, file=paste0(folder,"list_enhancers_topgenes.Rdata"))
write.table(enh.gene.pairs, file=paste0(folder,"Tables/list_enhancers_topgenes.txt"), quote=F, row.names=F, col.names=T, sep="\t")
