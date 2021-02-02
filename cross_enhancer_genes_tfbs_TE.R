### Maud Fagny
### 2021-02-01
### cross_enhancers_genes_tfbs_TE.R
### Perform enrichment analysis for TE superfamilies and families within enhancers and TFBS.
###_______________________________

### Load Libraries
library(data.table)
library(tidyr)
library(purrr)
library(IRanges)
library(topGO)
library(Rgraphviz)
library(ggplot2)


### Set variables
args=commandArgs(trailingOnly=TRUE)
folder=args[1]

### Load data
load(paste0(folder, "crosstables_TE_gene_enhancers_tfbs.Rdata"))


### Stats by tissues
tfbs.TE=tfbs.TE.genes
tfbs.TE[, start.gene:=NULL]
tfbs.TE[, end.gene:=NULL]
tfbs.TE[, geneID:=NULL]
tfbs.TE=unique(tfbs.TE)

pdf(paste0(folder, "Figures/overlap_TFBS_TE_propseq_distrib.pdf"))
hist(tfbs.TE$prop.TE[tfbs.TE$prop.TE>0], main="TFBS overlapping TE -Any tissue", 
     xlab="Proportion of overlapping TFBS sequence", ylab="Counts")
dev.off()

count.tfbs.te=function(tab){
  c("Any"=sum(tab[superfamily %in% c("LTR;TIR", "TIR", "MITE", "LTR", "LINE", "Helitron")]$prop.TE>0)/nrow(tab), 
    "Multiple"=sum(tab[superfamily=="LTR;TIR"]$prop.TE>0)/nrow(tab),
    "TIR"=sum(tab[superfamily=="TIR"]$prop.TE>0)/nrow(tab),
    "MITE"=sum(tab[superfamily=="MITE"]$prop.TE>0)/nrow(tab),
    "LTR"=sum(tab[superfamily=="LTR"]$prop.TE>0)/nrow(tab),
    "LINE"=sum(tab[superfamily=="LINE"]$prop.TE>0)/nrow(tab),
    "Helitron"=sum(tab[superfamily=="Helitron"]$prop.TE>0)/nrow(tab))
    }
summary.te.tfbs=data.frame("Genome"=c("Any"=TE.prop.all[1], "Multiple"=0, TE.prop.all[c("TIR", "MITE", "LTR", "LINE", "Helitron")]),
                           "All"=  c(count.tfbs.te(tfbs.TE)),
           "Both"=count.tfbs.te(tfbs.TE[tissue=="Both"]),
           "Husk"=count.tfbs.te(tfbs.TE[tissue=="HUSK"]),
           "V2-IST"=count.tfbs.te(tfbs.TE[tissue=="V2-IST"]))           
summary.te.tfbs=rbind("None"=1-summary.te.tfbs[1,], summary.te.tfbs)

pdf(paste0(folder, "Figures/prop_TFBS_withTE_byorder.pdf"))
par(mar=c(4,5,4,1)+.1)
barplot(as.matrix(summary.te.tfbs[c(1,3:nrow(summary.te.tfbs)),2:ncol(summary.te.tfbs)]), 
        xlab="Proportion of TFBS with at least 1 TE", 
        col=c("white", "gold", "forestgreen", "green3", "blue", "cyan", "purple"), 
        horiz=T)

par(xpd=T)
legend(0, par("usr")[4]+(par("usr")[4]-par("usr")[3])/7,
rownames(summary.te.tfbs)[c(1,3:nrow(summary.te.tfbs))], bty='n',
       fill=c("white", "gold", "forestgreen", "green3", "blue", "cyan", "purple"), ncol=3)
dev.off()

### Test for TE order enrichment among TFBS
res.enrich.TFBS.order=NULL
tissue=tfbs.TE$tissue
for(te in c("TIR", "MITE", "LTR", "LINE")){
  real.tmp=(tfbs.TE$superfamily== te)
  fh=fisher.test(table(real.tmp, (tissue=="HUSK")))
  fv=fisher.test(table(real.tmp, (tissue=="V2-IST")))
  res.enrich.TFBS.order=rbind(res.enrich.TFBS.order, 
                             c("p.husk"=fh$p.value, "OR.husk"=fh$estimate,
                               "p.v2"=fv$p.value, "OR.v2"=fv$estimate))
  rownames(res.enrich.TFBS.order)[nrow(res.enrich.TFBS.order)]=te
}

### Test for TE superfamily enrichment among TFBS
res.enrich.TFBS.sf=NULL
TE.order.sf=unique(paste(tfbs.TE$superfamily, tfbs.TE$family, sep="_"))
tissue=tfbs.TE$tissue
for(te in TE.order.sf[TE.order.sf != "None_None"]){
  real.tmp=(paste(tfbs.TE$superfamily, tfbs.TE$family, sep="_")== te)
  if(sum(real.tmp)>4 & length(grep(";", te))==0){
  fh=fisher.test(table(real.tmp, (tissue=="HUSK")))
  fv=fisher.test(table(real.tmp, (tissue=="V2-IST")))
  res.enrich.TFBS.sf=rbind(res.enrich.TFBS.sf, 
                              c("p.husk"=fh$p.value, "OR.husk"=fh$estimate, "OR.min.husk"=fh$conf.int[1], "OR.max.husk"=fh$conf.int[2],
                                "p.v2"=fv$p.value, "OR.v2"=fv$estimate, "OR.min.v2"=fv$conf.int[1], "OR.max.v2"=fv$conf.int[2]))
  rownames(res.enrich.TFBS.sf)[nrow(res.enrich.TFBS.sf)]=te
  }
}

res.enrich.TFBS.sf=data.frame(res.enrich.TFBS.sf)
res.enrich.TFBS.sf=res.enrich.TFBS.sf[order(rownames(res.enrich.TFBS.sf)),]
rownames(res.enrich.TFBS.sf)=c("RLC", "RLG", "RLX","DTA", "DTC", "DTM")
res.enrich.TFBS.sf=res.enrich.TFBS.sf[c(3,2,1,6,5,4),]
color=c(rep("blue",3), rep("forestgreen",3) )
pdf(paste0(folder, "Figures/OR_TFBS_withTE_bysf.pdf"), height=5, width=8)
par(las=1, mar=c(4,5,0,0)+.1)
plot(log2(res.enrich.TFBS.sf$OR.v2.odds.ratio), c(1:nrow(res.enrich.TFBS.sf)), 
     col=color, pch=18, cex=3, xlab="Odds ratio", ylab="", yaxt="n", ylim=c(0.5,6.5),
     xlim=c(min(log2(res.enrich.TFBS.sf$OR.min.v2)), max(log2(res.enrich.TFBS.sf$OR.max.v2))), cex.axis=2, cex.lab=2)
for(i in 1:6){
  lines(c(log2(res.enrich.TFBS.sf$OR.min.v2[i]), log2(res.enrich.TFBS.sf$OR.max.v2[i])), c(i,i), 
      col=color[i], lwd=3)
  lines(rep(log2(res.enrich.TFBS.sf$OR.min.v2[i]),2), c((i-0.1), (i+0.1)), 
        col=color[i], lwd=3)
  lines(rep(log2(res.enrich.TFBS.sf$OR.max.v2[i]),2), c((i-0.1), (i+0.1)) , 
        col=color[i], lwd=3)
}
lines(c(0,0), c(-1,8), lwd=3)
axis(side=2, labels=gsub("_", " ", rownames(res.enrich.TFBS.sf)), at=1:6,cex.axis=2, cex.lab=2)
dev.off()

### Summary table of TE order and family by enhancer

summarize.annot=function(tab.annot, tab.enh, levels){
  res=lapply(unique(tab.annot$enh.pos), function(n, annot, levels, enh){
    annot=data.table(annot)
    tot.length=enh[enh.pos==n, enh.end-enh.start+1]
    annotated.sf=summary(factor(annot[enh.pos == n]$te.superfamily, level=levels))
    annotated.f=tapply(annot[enh.pos == n]$te.family, 
                       factor(annot[enh.pos == n]$te.superfamily, level=levels), 
                      unique)
    length.annotated.te=tapply(annot[enh.pos == n]$te.length, 
                               factor(annot[enh.pos == n]$te.superfamily, level=levels), 
                      sum)
    length.annotated.te[is.na(length.annotated.te)]=0
    annotated.f[is.na(annotated.f)]=""
    return(c(as.character(n), tot.length, annotated.sf, annotated.f, 
             length.annotated.te, length.annotated.te/tot.length,
             sum(annotated.sf[levels %in% c("TIR", "LINE", "LTR", "MITE", "Helitron")]),
             sum(annotated.f[levels %in% c("TIR", "LINE", "LTR", "MITE", "Helitron")] !=""), 
             sum(length.annotated.te[levels %in% c("TIR", "LINE", "LTR", "MITE", "Helitron")]), 
             sum(length.annotated.te[levels %in% c("TIR", "LINE", "LTR", "MITE", "Helitron")])/tot.length))
   }, annot=tab.annot, levels=levels, enh=tab.enh)
  res2=data.frame(matrix(unlist(res), ncol=length(res[[1]]), byrow=TRUE), stringsAsFactors = FALSE)
  colnames(res2)=c("enhancer", "tot.length", levels, 
                  paste0("families.", levels), 
                  paste0("length.annotated", levels), paste0("prop.length.annotated", levels),
                  "nb.annotated.superfamilies","nb.annotated.families",
                  "tot.annotated.te", "prop.annotated.te")
  return(res2)
}
merged.enh=unique(enh.TE[,c("chr", "enh.start", "enh.end", "enh", "enh.pos", "tissue")])
tab.annot.all=summarize.annot(tab.annot=annotated.enh.TE, tab.enh=merged.enh, 
                          levels=unique(c(annot.TE.byconsensus$superfamily,
                                          enh.TE$superfamily))) 
tis=merged.enh$tissue
names(tis)=merged.enh$enh.pos
tab.annot.all$tissue=tis[tab.annot.all$enhancer]
summary.table=data.frame(rbind(both=c("total"=sum(merged.enh$tissue=="Both") ,
                                      "total.TE"=sum((tab.annot.all$tissue=="Both") & as.numeric(tab.annot.all$nb.annotated.superfamilies)>0), 
                                      "total.length.TE"=sum( as.numeric(tab.annot.all$tot.length[(tab.annot.all$tissue=="Both")])),
                                      "LTR"=sum((tab.annot.all$tissue=="Both") & tab.annot.all$LTR=="1"),
                                      "LINE"=sum((tab.annot.all$tissue=="Both") & tab.annot.all$LINE=="1"),
                                      "TIR"=sum( (tab.annot.all$tissue=="Both") & tab.annot.all$TIR=="1"),
                                      "MITE"=sum((tab.annot.all$tissue=="Both") & tab.annot.all$MITE=="1"),
                                      "Helitron"=sum((tab.annot.all$tissue=="Both") & tab.annot.all$Helitron=="1"),
                                      "Any.TE"=sum((tab.annot.all$tissue=="Both") & (tab.annot.all$LTR=="1" | tab.annot.all$LINE=="1") & (tab.annot.all$TIR=="1"| tab.annot.all$MITE=="1"| tab.annot.all$Helitron=="1")),
                                       "tot.cov.LTR"=sum(as.numeric(tab.annot.all$length.annotatedLTR[(tab.annot.all$tissue=="Both")])),
                                      "tot.cov.LINE"=sum( as.numeric(tab.annot.all$length.annotatedLINE[(tab.annot.all$tissue=="Both")])),
                                      "tot.cov.TIR"=sum( as.numeric(tab.annot.all$length.annotatedTIR[(tab.annot.all$tissue=="Both")])),
                                      "tot.cov.MITE"=sum(as.numeric(tab.annot.all$length.annotatedMITE[(tab.annot.all$tissue=="Both")])),
                                      "tot.cov.Helitron"=sum(as.numeric(tab.annot.all$length.annotatedHelitron[(tab.annot.all$tissue=="Both")]))
),
    husk=c("total"=sum(merged.enh$tissue=="HUSK") ,
           "total.TE"=sum((tab.annot.all$tissue=="HUSK") & as.numeric(tab.annot.all$nb.annotated.superfamilies)>0), 
           "total.length.TE"=sum( as.numeric(tab.annot.all$tot.length[(tab.annot.all$tissue=="Husk")])),
           "LTR"=sum((tab.annot.all$tissue=="HUSK") & tab.annot.all$LTR=="1"),
           "LINE"=sum((tab.annot.all$tissue=="HUSK") & tab.annot.all$LINE=="1"),
           "TIR"=sum( (tab.annot.all$tissue=="HUSK") & tab.annot.all$TIR=="1"),
           "MITE"=sum((tab.annot.all$tissue=="HUSK") & tab.annot.all$MITE=="1"),
           "Helitron"=sum((tab.annot.all$tissue=="HUSK") & tab.annot.all$Helitron=="1"),
           "Any.TE"=sum((tab.annot.all$tissue=="HUSK") & (tab.annot.all$LTR=="1" | tab.annot.all$LINE=="1") & (tab.annot.all$TIR=="1"| tab.annot.all$MITE=="1"| tab.annot.all$Helitron=="1")),
           "tot.cov.LTR"=sum(as.numeric(tab.annot.all$length.annotatedLTR[(tab.annot.all$tissue=="HUSK")])),
           "tot.cov.LINE"=sum( as.numeric(tab.annot.all$length.annotatedLINE[(tab.annot.all$tissue=="HUSK")])),
           "tot.cov.TIR"=sum( as.numeric(tab.annot.all$length.annotatedTIR[(tab.annot.all$tissue=="HUSK")])),
           "tot.cov.MITE"=sum(as.numeric(tab.annot.all$length.annotatedMITE[(tab.annot.all$tissue=="HUSK")])),
           "tot.cov.Helitron"=sum(as.numeric(tab.annot.all$length.annotatedHelitron[(tab.annot.all$tissue=="HUSK")]))
    ),
      V2=c("total"=sum(merged.enh$tissue=="V2-IST") ,
           "total.TE"=sum((tab.annot.all$tissue=="V2-IST") & as.numeric(tab.annot.all$nb.annotated.superfamilies)>0), 
           "total.length.TE"=sum( as.numeric(tab.annot.all$tot.length[(tab.annot.all$tissue=="V2-IST")])),
           "LTR"=sum((tab.annot.all$tissue=="V2-IST") & tab.annot.all$LTR=="1"),
           "LINE"=sum((tab.annot.all$tissue=="V2-IST") & tab.annot.all$LINE=="1"),
           "TIR"=sum( (tab.annot.all$tissue=="V2-IST") & tab.annot.all$TIR=="1"),
           "MITE"=sum((tab.annot.all$tissue=="V2-IST") & tab.annot.all$MITE=="1"),
           "Helitron"=sum((tab.annot.all$tissue=="V2-IST") & tab.annot.all$Helitron=="1"),
           "Any.TE"=sum((tab.annot.all$tissue=="V2-IST") & (tab.annot.all$LTR=="1" | tab.annot.all$LINE=="1") & (tab.annot.all$TIR=="1"| tab.annot.all$MITE=="1"| tab.annot.all$Helitron=="1")),
           "tot.cov.LTR"=sum(as.numeric(tab.annot.all$length.annotatedLTR[(tab.annot.all$tissue=="V2-IST")])),
           "tot.cov.LINE"=sum( as.numeric(tab.annot.all$length.annotatedLINE[(tab.annot.all$tissue=="V2-IST")])),
           "tot.cov.TIR"=sum( as.numeric(tab.annot.all$length.annotatedTIR[(tab.annot.all$tissue=="V2-IST")])),
           "tot.cov.MITE"=sum(as.numeric(tab.annot.all$length.annotatedMITE[(tab.annot.all$tissue=="V2-IST")])),
           "tot.cov.Helitron"=sum(as.numeric(tab.annot.all$length.annotatedHelitron[(tab.annot.all$tissue=="V2-IST")]))
      )))

### Plot distribution of proportion of enhancer sequence covered by each order of TE
pdf(paste0(folder, "Figures/TE_prop_seq_enhancers.pdf"), height=11.5, width=8)
par(mfrow=c(5,3), mar=c(5,4,3,0)+.1)
hist(as.numeric(tab.annot.all$prop.length.annotatedTIR)[tab.annot.all$prop.length.annotatedTIR !="0" & 
                                                          (tab.annot.all$tissue=="V2-IST") ], 
     breaks=seq(0,1,0.05), main="V2-IST (TIR)", xlab="% sequence covered by TE")
hist(as.numeric(tab.annot.all$prop.length.annotatedTIR)[tab.annot.all$prop.length.annotatedTIR !="0" & 
                                                          (tab.annot.all$tissue=="HUSK")], 
     breaks=seq(0,1,0.05), main="Husk (TIR)", xlab="% sequence covered by TE")
hist(as.numeric(tab.annot.all$prop.length.annotatedTIR)[tab.annot.all$prop.length.annotatedTIR !="0" & 
                                                          (tab.annot.all$tissue=="Both")], 
     breaks=seq(0,1,0.05), main="Both (TIR)", xlab="% sequence covered by TE")
hist(as.numeric(tab.annot.all$prop.length.annotatedMITE)[tab.annot.all$prop.length.annotatedMITE !="0" & 
                                                           (tab.annot.all$tissue=="V2-IST")], 
     breaks=seq(0,1,0.05), main="V2-IST (MITE)", xlab="% sequence covered by TE")
hist(as.numeric(tab.annot.all$prop.length.annotatedMITE)[tab.annot.all$prop.length.annotatedMITE !="0" & 
                                                           (tab.annot.all$tissue=="HUSK")], 
     breaks=seq(0,1,0.05), main="Husk (MITE)", xlab="% sequence covered by TE")
hist(as.numeric(tab.annot.all$prop.length.annotatedMITE)[tab.annot.all$prop.length.annotatedMITE !="0" & 
                                                           (tab.annot.all$tissue=="Both")], 
     breaks=seq(0,1,0.05), main="Both (MITE)", xlab="% sequence covered by TE")
hist(as.numeric(tab.annot.all$prop.length.annotatedLTR)[tab.annot.all$prop.length.annotatedLTR !="0" & 
                                                          (tab.annot.all$tissue=="V2-IST")], 
     breaks=seq(0,1,0.05), main="V2-IST (LTR)", xlab="% sequence covered by TE")
hist(as.numeric(tab.annot.all$prop.length.annotatedLTR)[tab.annot.all$prop.length.annotatedLTR !="0" & 
                                                          (tab.annot.all$tissue=="HUSK")], 
     breaks=seq(0,1,0.05), main="Husk (LTR)", xlab="% sequence covered by TE")
hist(as.numeric(tab.annot.all$prop.length.annotatedLTR)[tab.annot.all$prop.length.annotatedLTR !="0" & 
                                                          (tab.annot.all$tissue=="Both")], 
     breaks=seq(0,1,0.05), main="Both (LTR)", xlab="% sequence covered by TE")
hist(as.numeric(tab.annot.all$prop.length.annotatedLINE)[tab.annot.all$prop.length.annotatedLINE !="0" & 
                                                           (tab.annot.all$tissue=="V2-IST")], 
     breaks=seq(0,1,0.05), main="V2-IST (LINE)", xlab="% sequence covered by TE")
hist(as.numeric(tab.annot.all$prop.length.annotatedLINE)[tab.annot.all$prop.length.annotatedLINE !="0" & 
                                                           (tab.annot.all$tissue=="HUSK")], 
     breaks=seq(0,1,0.05), main="Husk (LINE)", xlab="% sequence covered by TE")
hist(as.numeric(tab.annot.all$prop.length.annotatedLINE)[tab.annot.all$prop.length.annotatedLINE !="0" & 
                                                           (tab.annot.all$tissue=="Both")], 
     breaks=seq(0,1,0.05), main="Both (LINE)", xlab="% sequence covered by TE")
hist(as.numeric(tab.annot.all$prop.length.annotatedHelitron)[tab.annot.all$prop.length.annotatedHelitron !="0" & 
                                                               (tab.annot.all$tissue=="V2-IST")], 
     breaks=seq(0,1,0.05), main="V2-IST (Helitron)", xlab="% sequence covered by TE")
hist(as.numeric(tab.annot.all$prop.length.annotatedHelitron)[tab.annot.all$prop.length.annotatedHelitron !="0" & 
                                                               (tab.annot.all$tissue=="HUSK")], 
     breaks=seq(0,1,0.05), main="Husk (Helitron)", xlab="% sequence covered by TE")
hist(as.numeric(tab.annot.all$prop.length.annotatedHelitron)[tab.annot.all$prop.length.annotatedHelitron !="0" & 
                                                               (tab.annot.all$tissue=="Both")], 
     breaks=seq(0,1,0.05), main="Both (Helitron)", xlab="% sequence covered by TE")

dev.off()

t(apply(summary.table, 1, function(x) c(x[1:2], round(x[4:8]/x[2], digits = 3))))

enh.seq.cov.te=cbind("Both"=as.numeric(summary.table[1,10:14]/summary.table[1,3]),
                     "Husk"=as.numeric(summary.table[2,10:14]/summary.table[2,3]),
                     "V2-IST"=as.numeric(summary.table[1,10:14]/summary.table[3,3]),
                               "Genome"=as.numeric(TE.prop.all[c(4,3,2,5,6)]))
rownames(enh.seq.cov.te)=c("LTR", "LINE", "TIR", "MITE", "Helitron")
enh.seq.cov.te=rbind("None"=1-colSums(enh.seq.cov.te), enh.seq.cov.te)

pdf(paste0(folder, "Figures/TE_coverage_enhancers.pdf"))
par( mar=c(4,5,4,1)+.1)
barplot(enh.seq.cov.te,  horiz=T, , 
col=c("white", "forestgreen", "green3", "blue", "cyan", "purple"),
xlab="Proportion of sequence covered with TE")
par(xpd=TRUE)
legend(0, par("usr")[4]+(par("usr")[4]-par("usr")[3])/7, legend=c("None", "LTR", "LINE", "TIR", "MITE", "Helitron"), 
       fill=c("white", "forestgreen", "green3", "blue", "cyan", "purple"), bty='n', ncol=3)
dev.off()


### Enrichment in TE among enhancers
enh.genes=enh.genes[!is.na(geneID) & (geneID %in% rownames(expr)) & (enh.pos %in% enh.TE$enh.pos) ]
enh.TE=enh.TE[ (enh.pos %in% unique(enh.genes$enh.pos))  & (enh.pos %in% unique(tfbs.TE.genes$enh.pos)) ]

##Number of enhancers with TE
real.TE=tapply(enh.TE$superfamily, enh.TE$enh.pos, function(x){
  u=unique(sort(x))
  if((length(u)>1) & ("None" %in% u)){u=u[u != "None"]}
  if(length(u)>1){u="Multiple"}
  u
  })
real.TE[real.TE=="None"]=""
tissue=tapply(enh.TE$tissue, enh.TE$enh.pos, unique)
tissue=tissue[names(real.TE)]
table(real.TE)/sum( table(real.TE))
tab=data.frame("All"=summary(as.factor(real.TE))/ length(real.TE),
               "Both"=summary(as.factor(real.TE)[tissue=="Both"])/sum(tissue=="Both"),
               "Husk"=summary(as.factor(real.TE)[tissue=="HUSK"])/sum(tissue=="HUSK"),
               "V2-IST"=summary(as.factor(real.TE)[tissue=="V2-IST"])/sum(tissue=="V2-IST"))

pdf(paste0(folder, "Figures/proportion_enhancer_with_TE.pdf"), height=5, width=8)
par(mar=c(3,7,7,2)+.1, las=1)
barplot(as.matrix(tab)[c(1,6,4,3,7,5,2),], horiz=T, 
        col=c("white", "gold", "forestgreen", "green", "blue", "cyan", "purple"),cex.axis=2, cex.lab=2,cex.names=2, names=c(colnames(tab)[1:3], "V2-IST"))
#barplot(as.matrix(tab)[c(1,6,4,3,7,5,2),3:4], horiz=T, 
#        col=c("white", "gold", "forestgreen", "green3", "blue", "cyan", "purple"), cex.names=1.5, cex.axis = 1.5)
par("xpd"=TRUE)
legend(par("usr")[1], ((par("usr")[4]-par("usr")[3])/2+par("usr")[4]), 
       legend=c("None",rownames(tab)[c(0,6,4,3,7,5,2)]), 
       fill=c("white", "gold", "forestgreen", "green", "blue", "cyan", "purple"),
       bty="n", ncol=3, cex=2)
dev.off()

## Enrichment of orders
res.enrich.TE.counts=NULL
for(te in c("TIR", "MITE", "LTR", "LINE", "Helitron")){
  real.tmp=tapply(enh.TE$superfamily, enh.TE$enh.pos, function(x){any(x == te)})
  fh=fisher.test(table(real.tmp, (tissue[names(real.tmp)]=="HUSK")))
  fv=fisher.test(table(real.tmp, (tissue[names(real.tmp)]=="V2-IST")))
  res.enrich.TE.counts=rbind(res.enrich.TE.counts, 
                             c("p.husk"=chisq.test(table(real.tmp, (tissue[names(real.tmp)]=="Husk")))$p.value, "OR.husk"=fh$estimate,
                               "p.v2"=chisq.test(table(real.tmp, (tissue[names(real.tmp)]=="V2-IST")))$p.value, "OR.v2"=fv$estimate))
}
rownames(res.enrich.TE.counts)=c("TIR", "MITE", "LTR", "LINE", "Helitron")

## Enrichment of superfamilies
res.enrich.TE.sf=NULL
TE.order.sf=unique(paste(enh.TE$superfamily, enh.TE$family, sep="_"))
tissue=enh.TE$tissue
for(te in TE.order.sf[TE.order.sf != "None_None"]){
  real.tmp=(paste(enh.TE$superfamily, enh.TE$family, sep="_")== te)
  fh=fisher.test(table(real.tmp, (tissue=="HUSK")))
  fv=fisher.test(table(real.tmp, (tissue=="V2-IST")))
  res.enrich.TE.sf=rbind(res.enrich.TE.sf, 
                           c("p.husk"=fh$p.value, "OR.husk"=fh$estimate,
                             "p.v2"=fv$p.value, "OR.v2"=fv$estimate))
  rownames(res.enrich.TE.sf)[nrow(res.enrich.TE.sf)]=te
}

write.table(res.enrich.TE.counts)

##Proportion of sequence covered
summarize.annot=function(tab.annot, tab.enh, levels){
  res=lapply(names(tab.annot), function(n, annot, levels, enh){
    tot.length=enh[names==n, end-start+1]
    annotated.sf=summary(factor(annot[[`n`]]$te.superfamily, level=levels))
    annotated.f=tapply(annot[[`n`]]$te.family, 
                       factor(annot[[`n`]]$te.superfamily, level=levels), 
                       unique)
    length.annotated.te=tapply(annot[[`n`]]$te.length, 
                               factor(annot[[`n`]]$te.superfamily, level=levels), 
                               sum)
    length.annotated.te[is.na(length.annotated.te)]=0
    annotated.f[is.na(annotated.f)]=""
    return(c(n, tot.length, annotated.sf, annotated.f, 
             length.annotated.te, length.annotated.te/tot.length,
             sum(annotated.sf[levels %in% c("DNA", "LINE", "LTR", "MITE")]),
             sum(annotated.f[levels %in% c("DNA", "LINE", "LTR", "MITE")] !=""), 
             sum(length.annotated.te[levels %in% c("DNA", "LINE", "LTR", "MITE")]), 
             sum(length.annotated.te[levels %in% c("DNA", "LINE", "LTR", "MITE")])/tot.length))
  }, annot=tab.annot, levels=levels, enh=tab.enh)
  res2=data.frame(matrix(unlist(res), ncol=length(res[[1]]), byrow=TRUE), stringsAsFactors = FALSE)
  colnames(res2)=c("enhancer", "tot.length", levels, 
                   paste0("families.", levels), 
                   paste0("length.annotated", levels), paste0("prop.length.annotated", levels),
                   "nb.annotated.superfamilies","nb.annotated.families",
                   "tot.annotated.te", "prop.annotated.te")
  return(res2)
}


## Get genes start and end
B73.genes=read.table(paste0(folder, "genes.gff3"), header=F, stringsAsFactors = F, quote="", sep="\t")
B73.genes=B73.genes[, c(1,4,5,7,9)]
colnames(B73.genes)=c("chr", "gene.start", "gene.end", "gene.strand", "gene.annot")
genes=B73.genes[grep("ID=gene", B73.genes$gene.annot ),]
genes=genes[genes$chr %in% as.character(1:10),]
genenames=unlist(lapply(genes$gene.annot, function(x){strsplit(x, ";")[[1]][1]}))
tss=data.frame("chr"=as.numeric(genes$chr), "start.gene"=ifelse(genes$gene.strand=="+", genes$gene.start, genes$gene.end), 
               "geneID"=gsub("ID=gene:", "", genenames), "strand"=genes$gene.strand, stringsAsFactors = F)

tss$end.gene=tss$start.gene+1
tss=data.table(tss)
setkey(tss, chr, start.gene, end.gene)
strand=tss$strand
names(strand)=tss$geneID

## Enrichment in husk-specific  short-distance enhancers compared to V2-IST-specific or shared enhancers
enh.genes=enh.genes[!is.na(geneID) & (geneID %in% rownames(expr)) & (enh.pos %in% enh.TE$enh.pos) ]
enh.genes[,gene.strand:=strand[enh.genes$geneID]]

enh.genes[,dist.start:=ifelse(strand[enh.genes$geneID]=="+", enh.genes$enh.start-enh.genes$start.gene, enh.genes$start.gene-enh.genes$enh.start)]
enh.genes[,dist.end:=ifelse(strand[enh.genes$geneID]=="+", enh.genes$enh.end-enh.genes$start.gene, enh.genes$start.gene-enh.genes$enh.end)]
enh.genes[, dist.gene:=apply(enh.genes[, c("dist.start", "dist.end")], 1, function(x){x[abs(x)==min(abs(x))][1]})]
min.dist=tapply(enh.genes$dist.gene, enh.genes$enh.pos, function(x){x[abs(x)==min(abs(x))][1]})
tissue=tapply(enh.genes$tissue, enh.genes$enh.pos, unique)
p=wilcox.test(abs(min.dist[tissue=="Husk"]), abs(min.dist[tissue !="HUSK"]), alternative="less")$p.value

hh=hist(abs(min.dist[tissue=="HUSK"]), breaks = seq(0, 250000, 5000), plot=F)
hv=hist(abs(min.dist[tissue=="V2-IST"]), breaks = seq(0, 250000, 5000), plot=F)
hb=hist(abs(min.dist[tissue=="Both"]), breaks = seq(0, 250000, 5000), plot=F)
hh$counts=hh$counts/sum(hh$counts)
hv$counts=hv$counts/sum(hv$counts)
hb$counts=hb$counts/sum(hb$counts)


pdf(paste0(folder, "Figures/distrib_dist_gene_enhancers.pdf"))
par(mar=c(4,5,1,1)+.1, las=0, xpd=FALSE)
plot(hb, border="black", main="", xlab="Distance Enhancer-Nearest gene", 
     ylim=c(0, max(hh$counts, hb$counts, hv$counts)))
plot(hh, border="olivedrab1",  add=T)
plot(hv, border="blue", add=T)
legend("topright", bty='n', legend=c("Both", "Husk", "V2-IST"), 
       col=c("black", "olivedrab1", "blue"), lty=1, pch=-1)

plot(hb, border="black", main="", xlab="Distance Enhancer-Nearest gene", 
     ylim=c(0, max(hh$counts, hb$counts, hv$counts)), xlim=c(0, 100000))
plot(hh, border="olivedrab1",  add=T)
plot(hv, border="blue", add=T)
legend("topright", bty='n', legend=c("Both", "Husk", "V2-IST"), 
       col=c("black", "olivedrab1", "blue"), lty=1, pch=-1)

bw=bw.nrd(sample(abs(min.dist), 50))
plot(density(abs(min.dist[tissue=="HUSK"]), bw = bw), col="olivedrab1", main="", xlab="Distance Enhancer-Nearest gene")
lines(density(abs(min.dist[tissue=="Both"]), bw = bw))
lines(density(abs(min.dist[tissue=="V2-IST"]), bw = bw), col="blue")
legend("topright", bty='n', legend=c("Both", "Husk", "V2-IST"), 
       col=c("black", "olivedrab1", "blue"), lty=1, pch=-1)
plot(density(abs(min.dist[tissue=="HUSK"]), bw = bw), col="olivedrab1", main="", 
     xlab="Distance Enhancer-Nearest gene", xlim=c(0, 100000))
lines(density(abs(min.dist[tissue=="Both"]), bw = bw))
lines(density(abs(min.dist[tissue=="V2-IST"]), bw = bw), col="blue")
legend("topright", bty='n', legend=c("Both", "Husk", "V2-IST"), 
       col=c("black", "olivedrab1", "blue"), lty=1, pch=-1)
dev.off()

### Resampling test of husk enhancer MITE content

## Get minimum distance from gene and enhancer length for each gene
enh.genes[,enh.length:=(enh.genes$enh.end-enh.genes$enh.start+1)]
enh.chr=tapply(enh.genes$chr, enh.genes$enh.pos, unique)
enh.length=tapply(enh.genes$enh.length, enh.genes$enh.pos, unique)

## Resample a dataset with similar size as the real dataset 
## with similar genomic characteristics (x1000)

## Functions

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
  return(list(enh.TE, res))
}
resample.genomicloc=function(nsam, chr.all, tss, min.d, len, real.TE, TE, tissue, tis){
  sampled.p.chr=NULL
  sampled.s.chr=NULL
  sampled.t.chr=NULL
  sampled.pcomp.chr=NULL
  sampled.tcomp.chr=NULL
  for(k in 1:nsam){
    if(k/10==floor(k/10)){print(k)}
    sampled=NULL
    for(i in 1:10){
      s=sample.int(nrow(tss[chr == i]), sum(chr.all==i))
      enh.tmp=tss[chr == i]$start.gene[s]
      str.tmp=tss[chr == i]$strand[s]
      sampled.enh.start=ifelse(str.tmp=="+", ifelse(min.d[chr.all==i]>0, enh.tmp+min.d[chr.all==i],
                                                    enh.tmp+min.d[chr.all==i]-len[chr.all==i]), 
                               ifelse(min.d[chr.all==i]>0,enh.tmp-min.d[chr.all==i]-len[chr.all==i],
                                      enh.tmp-min.d[chr.all==i]))
      sampled.enh.end=sampled.enh.start+abs(len[chr.all==i])
      sampled.enh.chr=chr.all[chr.all==i]
      sampled=rbind(sampled, data.frame("chr"=as.integer(sampled.enh.chr), "start"=as.integer(sampled.enh.start), 
                                        "end"=as.integer(sampled.enh.end),
                                        stringsAsFactors = F))
    }
    sampled$name=paste0("enh", 1:nrow(sampled))
    sampled=data.table(sampled)
    setkey(sampled, chr, start, end)
    sampled.TE=foverlaps(sampled, TE,  nomatch = 0)
    data=c(as.numeric(real.TE[tissue!=tis]), ifelse(sampled$name %in% sampled.TE$name, 1, 0))
    fac=c(rep(0, length(real.TE[tissue!=tis])), rep(1, length(sampled$name)))
    sampled.p.chr=c(sampled.p.chr,chisq.test(table(data, fac))$p.value)
    sampled.s.chr=c(sampled.s.chr,chisq.test(table(data, fac))$stat)
    sampled.t.chr=c(sampled.t.chr,fisher.test(table(data, fac))$estimate)
    data=c(as.numeric(real.TIR[tissue==tis]), ifelse(sampled$name %in% sampled.TE$name, 1, 0))
    fac=c(rep(1, length(real.TIR[tissue==tis])), rep(0, length(sampled$name)))
    sampled.pcomp.chr=c(sampled.pcomp.chr,chisq.test(table(data, fac))$p.value)
    sampled.tcomp.chr=c(sampled.tcomp.chr,fisher.test(table(data, fac))$estimate)
  }
  return(list("null.distrib.stat"=sampled.s.chr, "null.distrib.pval"=sampled.p.chr, "null.distrib.or"=sampled.t.chr, 
              "real.null.pval"=sampled.pcomp.chr, "real.null.or"=sampled.tcomp.chr))
}

## Run resamplings
tss=tss[start.gene>=250000]
load("Documents/INRA/GeneNetworks/Enhancers/data/TE_annot.Rdata")
enh.TE[, sf.family:=paste(enh.TE$superfamily, enh.TE$family, sep="_")]
real.MITE=tapply(enh.TE$superfamily, enh.TE$enh.pos, function(x){any(x == "MITE")})
real.MITE_DTT=tapply(enh.TE$sf.family, enh.TE$enh.pos, function(x){any(x == "MITE_DTT")})
real.MITE_DTM=tapply(enh.TE$sf.family, enh.TE$enh.pos, function(x){any(x == "MITE_DTM")})
real.TIR_DTM=tapply(enh.TE$sf.family, enh.TE$enh.pos, function(x){any(x == "TIR_DTM")})

resample.MITE.husk=resample.genomicloc(min.d=min.dist[tissue=="Husk"], len=enh.length[tissue=="Husk"], chr.all=enh.chr[tissue=="Husk"], 
                    real.TE=real.MITE, tissue=tissue, tis="Husk", TE=TE[!is.na(superfamily) & superfamily=="MITE"], nsam=1000, tss=tss)
real.MITE.c=chisq.test(table(real.MITE, ifelse(tissue=="Husk", TRUE, FALSE)))$statistic
sum((resample.MITE.husk$null.distrib.pval <= 0.05) & (resample.MITE.husk$null.distrib.or>1))/1000
sum((resample.MITE.husk$real.null.pval <= 0.05) & (resample.MITE.husk$real.null.or>1))/1000
sum((resample.MITE.husk$null.distrib.stat >= real.MITE.c) & (resample.MITE.husk$null.distrib.or>1))/1000

resample.MITE_DTT.husk=resample.genomicloc(min.d=min.dist[tissue=="Husk"], len=enh.length[tissue=="Husk"], chr.all=enh.chr[tissue=="Husk"], 
                                       real.TE=real.MITE_DTT, tissue=tissue, tis="Husk", TE=TE[!is.na(superfamily) & (superfamily=="MITE") & (family=="DTT")],
                                       nsam=1000, tss=tss)
real.MITE_DTT.c=chisq.test(table(real.MITE_DTT, ifelse(tissue=="Husk", TRUE, FALSE)))$statistic
sum((resample.MITE_DTT.husk$null.distrib.pval <= 0.05) & (resample.MITE_DTT.husk$null.distrib.or>1))/1000
sum((resample.MITE_DTT.husk$real.null.pval <= 0.05) & (resample.MITE_DTT.husk$real.null.or>1))/1000
sum((resample.MITE_DTT.husk$null.distrib.stat >= real.MITE_DTT.c) & (resample.MITE_DTT.husk$null.distrib.or>1))/1000

resample.MITE_DTM.husk=resample.genomicloc(min.d=min.dist[tissue=="Husk"], len=enh.length[tissue=="Husk"], chr.all=enh.chr[tissue=="Husk"], 
                                           real.TE=real.MITE_DTT, tissue=tissue, tis="Husk", TE=TE[!is.na(superfamily) & (superfamily=="MITE") & (family=="DTM")], 
                                           nsam=1000, tss=tss)
real.MITE_DTM.c=chisq.test(table(real.MITE_DTM, ifelse(tissue=="Husk", TRUE, FALSE)))$statistic
sum((resample.MITE_DTM.husk$null.distrib.pval <= 0.05) & (resample.MITE_DTM.husk$null.distrib.or>1))/1000
sum((resample.MITE_DTM.husk$real.null.pval <= 0.05) & (resample.MITE_DTM.husk$real.null.or>1))/1000
sum((resample.MITE_DTM.husk$null.distrib.stat >= real.MITE_DTM.c) & (resample.MITE_DTM.husk$null.distrib.or>1))/1000

resample.TIR_DTM.v2=resample.genomicloc(min.d=min.dist[tissue=="V2-IST"], len=enh.length[tissue=="V2-IST"], chr.all=enh.chr[tissue=="V2-IST"], 
                                    real.TE=real.TIR_DTM, tissue=tissue, tis="V2-IST", TE=TE[!is.na(superfamily) & (superfamily=="TIR") & (family=="DTM")], 
                                    nsam=1000, tss=tss)
real.TIR_DTM.c=chisq.test(table(real.TIR, ifelse(tissue=="V2-IST", TRUE, FALSE)))$statistic
sum((resample.TIR_DTM.v2$null.distrib.pval <= 0.01) & (resample.TIR_DTM.v2$null.distrib.or>1))/1000
sum((resample.TIR_DTM.v2$real.null.pval <= 0.01) & (resample.TIR_DTM.v2$real.null.or>1))/1000
sum((resample.TIR_DTM.v2$null.distrib.stat >= real.TIR_DTM.c) & (resample.TIR_DTM.v2$null.distrib.or>1))/1000

### Resampling TFBS
tfbs.genes[, tfbs.Name:=apply(data.frame(tfbs.genes[,c("chr", "tfbs.genome.start", "tfbs.genome.end", "Alias")], stringsAsFactors = F), 
                              1, function(x){gsub(" ", "", paste(x, collapse="_"))})]
tfbs.TE=tfbs.TE[tfbs.Name %in% tfbs.genes$tfbs.Name]
tfbs.genes=tfbs.genes[!is.na(geneID) & (geneID %in% rownames(expr)) & (tfbs.Name %in% tfbs.TE$tfbs.Name) ]
tfbs.genes[,gene.strand:=strand[tfbs.genes$geneID]]

tfbs.genes[,dist.start:=ifelse(strand[tfbs.genes$geneID]=="+", tfbs.genes$tfbs.genome.start-tfbs.genes$start.gene, tfbs.genes$start.gene-tfbs.genes$tfbs.genome.start)]
tfbs.genes[,dist.end:=ifelse(strand[tfbs.genes$geneID]=="+", tfbs.genes$tfbs.genome.end-tfbs.genes$start.gene, tfbs.genes$start.gene-tfbs.genes$tfbs.genome.end)]
tfbs.genes[, dist.gene:=apply(tfbs.genes[, c("dist.start", "dist.end")], 1, function(x){x[abs(x)==min(abs(x))][1]})]

tfbs.genes[,tfbs.length:=(tfbs.genes$tfbs.genome.end-tfbs.genes$tfbs.genome.start+1)]
tfbs.chr=tapply(tfbs.genes$chr, tfbs.genes$tfbs.Name, unique)
tfbs.length=tapply(tfbs.genes$tfbs.length, tfbs.genes$tfbs.Name, unique)

tfbs.TE[, sf.family:=paste(tfbs.TE$superfamily, tfbs.TE$family, sep="_")]

min.dist.tfbs=tapply(tfbs.genes$dist.gene, tfbs.genes$tfbs.Name, function(x){x[abs(x)==min(abs(x))][1]})
tissue=tapply(tfbs.TE$tissue, tfbs.TE$tfbs.Name, unique)
p=wilcox.test(abs(min.dist.tfbs[tissue=="Husk"]), abs(min.dist.tfbs[tissue !="Husk"]), alternative="less")$p.value

real.TIR.tfbs=tapply(tfbs.TE$superfamily, tfbs.TE$tfbs.Name, function(x){any(x == "TIR")})
real.TIR_DTM.tfbs=tapply(tfbs.TE$sf.family, tfbs.TE$tfbs.Name, function(x){any(x == "TIR_DTM")})
real.LTR_Copia.tfbs=tapply(tfbs.TE$sf.family, tfbs.TE$tfbs.Name, function(x){any(x == "LTR_Copia")})

resample.TIR.v2.tfbs=resample.genomicloc(min.d=min.dist.tfbs[tissue=="V2-IST"], len=tfbs.length[tissue=="V2-IST"], chr.all=tfbs.chr[tissue=="V2-IST"], 
                                       real.TE=real.TIR.tfbs, tissue=tissue, tis="V2-IST", TE=TE[!is.na(superfamily) & superfamily=="TIR"], nsam=1000, tss=tss)
real.TIR.TFBS.c=chisq.test(table(real.TIR.tfbs, ifelse(tissue=="V2-IST", TRUE, FALSE)))$statistic
sum((resample.TIR.v2.tfbs$null.distrib.pval <= 0.01) & (resample.TIR.v2.tfbs$null.distrib.or>1))/1000
sum((resample.TIR.v2.tfbs$real.null.pval <= 0.01) & (resample.TIR.v2.tfbs$real.null.or>1))/1000
sum((resample.TIR.v2.tfbs$null.distrib.stat >= real.TIR.TFBS.c) & (resample.TIR.v2.tfbs$null.distrib.or>1))/1000
resample.TIR_DTM.v2.tfbs=resample.genomicloc(min.d=min.dist.tfbs[tissue=="V2-IST"], len=tfbs.length[tissue=="V2-IST"], chr.all=tfbs.chr[tissue=="V2-IST"], 
                                         real.TE=real.TIR_DTM.tfbs, tissue=tissue, tis="V2-IST", 
                                         TE=TE[!is.na(superfamily) & superfamily=="TIR" & family=="DTM"], nsam=5000, tss=tss)
real.TIR_DTM.TFBS.c=chisq.test(table(real.TIR_DTM.tfbs, ifelse(tissue=="V2-IST", TRUE, FALSE)))$statistic
sum((resample.TIR_DTM.v2.tfbs$null.distrib.pval <= 0.01) & (resample.TIR_DTM.v2.tfbs$null.distrib.or>1))/5000
sum((resample.TIR_DTM.v2.tfbs$real.null.pval <= 0.01) & (resample.TIR_DTM.v2.tfbs$real.null.or>1))/5000
sum((resample.TIR_DTM.v2.tfbs$null.distrib.stat >= real.TIR_DTM.TFBS.c) & (resample.TIR_DTM.v2.tfbs$null.distrib.or>1))/5000
resample.LTR_Copia.v2.tfbs=resample.genomicloc(min.d=min.dist.tfbs[tissue=="V2-IST"], len=tfbs.length[tissue=="V2-IST"], chr.all=tfbs.chr[tissue=="V2-IST"], 
                                         real.TE=real.LTR_Copia.tfbs, tissue=tissue, tis="V2-IST", 
                                         TE=TE[!is.na(superfamily) & superfamily=="LTR" & family=="Copia"], nsam=1000, tss=tss)
real.LTR_Copia.TFBS.c=chisq.test(table(real.LTR_Copia.tfbs, ifelse(tissue=="V2-IST", TRUE, FALSE)))$statistic
sum((resample.LTR_Copia.v2.tfbs$null.distrib.pval <= 0.01) & (resample.LTR_Copia.v2.tfbs$null.distrib.or>1))/1000
sum((resample.LTR_Copia.v2.tfbs$real.null.pval <= 0.01) & (resample.LTR_Copia.v2.tfbs$real.null.or>1))/1000
sum((resample.LTR_Copia.v2.tfbs$null.distrib.stat >= real.LTR_Copia.TFBS.c) & (resample.LTR_Copia.v2.tfbs$null.distrib.or>1))/1000


### enhancers with mite GO analysis

husk.enh.mite=enh.TE[superfamily=="MITE" & tissue=="HUSK"]$enh.pos
genes.mite=unique(deg$geneID[deg$enh.pos %in% husk.enh.mite & deg$Husk>1 & deg$topgen.husk==1])
grp.mite=lapply(husk.genes, function(x, g){g[g %in% x]}, g=genes.mite )

library(topGO)
library(Rgraphviz)
all.genes=unique(unlist(husk.genes))
ZmB73_5a_xref.GO.topGO.FG <- read.delim(paste(folder, "GeneOntologyAnnotation.txt"), header=FALSE, stringsAsFactors = F)
geneNames<-ZmB73_5a_xref.GO.topGO.FG$V1
geneID2GO <- readMappings(file = paste0(folder,"GeneOntologyAnnotation.txt"))

allRes=list()
GOtermsGenes=NULL
allRes=NULL
myInterestingGenes=unique(genes.mite)
geneList <- factor(as.integer(geneNames[geneNames %in% all.genes] %in% myInterestingGenes))
names(geneList) <- geneNames[geneNames %in% all.genes]

GOtermsGenes=data.frame("Ontology"=character(0), "GOTerm"=character(0), "Genes"=character(0))
for(ont in c("BP", "MF")){
  GOdata <- new("topGOdata", ontology = ont, allGenes = geneList, 
                annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
  ## Compute stats
  test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultElimFisher <- getSigGroups(GOdata, test.stat)
   allRes[[`ont`]] <- GenTable(GOdata, 
                              elim = resultElimFisher, orderBy = "elim", 
                              ranksOf = "elim", topNodes = 100)
   allRes[[`ont`]]$ontology=rep(ont, nrow(allRes[[`ont`]]))
  
  pdf(paste0(folder, "Figures/GeneOntology/Husk_MITE_TFBS_targets_", ont, ".pdf"), width=11.5, height=8)
  showSigOfNodes(GOdata, score(resultElimFisher), firstSigNodes = 10, useInfo ='all')
  dev.off()
  
  myterms <- allRes[[`ont`]]$GO.ID
  mygenes  <- genesInTerm(GOdata, myterms)
  for (j in 1:length(myterms)) {
    myterm <- myterms[j]
    mygenesforterm <- mygenes[myterm][[1]]
    myfactor <- mygenesforterm %in% myInterestingGenes # find the genes that are in the list of genes of interest
    mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
    mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
    GOtermsGenes=rbind(GOtermsGenes, 
                       data.frame("Ontology"=ont, "GOTerm"=myterm,"Genes"=mygenesforterm2))
  }
}
allRes=rbind(allRes$BP, allRes$MF)
write.table(allRes[allRes$Significant >=3 &  
                     (allRes$elim <=0.01 ),], 
            file=paste0(folder, "Tables/GeneOntology//Gene_ontology_enrichment_Husk_MITE_TFBS_targets.txt"), 
            quote=FALSE,row.names=FALSE, sep="\t")
write.table(GOtermsGenes, file=paste0(folder, "Tables/GeneOntology/Gene_ontology_term_Husk_MITE_DTM_TFBS_targets.txt"), 
            quote=FALSE,row.names=FALSE, sep="\t")


### TFBS with DTM GO analysis

v2.tir.dtm.ERF=tfbs.TE[superfamily=="TIR" & family=="DTM"  & tissue  == "V2-IST" ][c(1,2,4,6,8),]
genes.v2=unique(tfbs.TE.genes[ tissue  == "V2-IST" & 
                                 geneID %in% deg$geneID[  deg$`V2-IST_Oka`>1 & deg$topgen.v2==1 ]]$geneID)
all.genes=unique(unlist(genes.v2))

ZmB73_5a_xref.GO.topGO.FG <- read.delim(paste(folder, "GeneOntologyAnnotation.txt"), header=FALSE, stringsAsFactors = F)
geneNames<-ZmB73_5a_xref.GO.topGO.FG$V1
geneID2GO <- readMappings(file = paste0(folder,"GeneOntologyAnnotation.txt"))

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, 
              annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)

allRes=list()
GOtermsGenes=NULL
allRes=NULL
myInterestingGenes=unique(unlist(grp))
geneList <- factor(as.integer(geneNames[geneNames %in% all.genes] %in% myInterestingGenes))
names(geneList) <- geneNames[geneNames %in% all.genes]

GOtermsGenes=data.frame("Ontology"=character(0), "GOTerm"=character(0), "Genes"=character(0))
for(ont in c("BP", "MF")){
  GOdata <- new("topGOdata", ontology = ont, allGenes = geneList, 
                annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
  ## Compute stats
  test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultElimFisher <- getSigGroups(GOdata, test.stat)
  allRes[[`ont`]] <- GenTable(GOdata, 
                                               elim = resultElimFisher, orderBy = "elim", 
                                               ranksOf = "elim", topNodes = 100)
  allRes[[`ont`]]$ontology=rep(ont, nrow(allRes[[`ont`]]))
      
  pdf(paste0(folder, "Figures/GeneOntology/V2-IST_TIR_DTM_TFBS_targets_", ont, ".pdf"), width=11.5, height=8)
  showSigOfNodes(GOdata, score(resultClassicFisher), firstSigNodes = 10, useInfo ='all')
  showSigOfNodes(GOdata, score(resultPCFisher), firstSigNodes = 10, useInfo ='all')
  showSigOfNodes(GOdata, score(resultElimFisher), firstSigNodes = 10, useInfo ='all')
  dev.off()
      
  myterms <- allRes[[`ont`]]$GO.ID
  mygenes  <- genesInTerm(GOdata, myterms)
  for (j in 1:length(myterms)) {
  myterm <- myterms[j]
  mygenesforterm <- mygenes[myterm][[1]]
  myfactor <- mygenesforterm %in% myInterestingGenes # find the genes that are in the list of genes of interest
  mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
  mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
  GOtermsGenes=rbind(GOtermsGenes, 
                     data.frame("Ontology"=ont, "GOTerm"=myterm,"Genes"=mygenesforterm2))
  }
}
allRes=rbind(allRes$BP, allRes$MF)
write.table(allRes[allRes$Significant >=3 &  
                     (allRes$elim <=0.01 ),], 
            file=paste0(folder, "Tables/GeneOntology//Gene_ontology_enrichment_V2-IST_TIR_DTM_TFBS_targets.txt"), 
            quote=FALSE,row.names=FALSE, sep="\t")
write.table(GOtermsGenes, file=paste0(folder, "Tables/GeneOntology/Gene_ontology_term_V2-IST_TIR_DTM_TFBS_targets.txt"), 
            quote=FALSE,row.names=FALSE, sep="\t")

pdf(paste0(folder, "Figures/GeneOntology//Gene_ontology_enrichment_V2-IST_TIR_DTM_TFBS_targets.pdf"), h=11, w=8)
allRes <- allRes[allRes$Significant >=5 &  allRes$adjPclassic <=0.01,]
allRes <- allRes[order(allRes$adjPclassic),]
allRes <- allRes[1:min(c(10, nrow(allRes))),]
allRes$OddsRatio=allRes$Significant/allRes$Expected
ggplot(allRes, aes(x=-log10(allRes$adjPclassic), y=allRes$Term, size=allRes$OddsRatio))+
  geom_point( shape=21, colour="blue", fill="dodgerblue3")+
  guides(fill=guide_legend(title.theme = element_text(size = 15)))+
  scale_size(range = c(min(allRes$OddsRatio), max(allRes$OddsRatio)*5))+
  labs(
          x = "-log10(FDR qvalue)",
          y = "",
          size = "Odds Ratio"
        )+
  theme_bw()+
  xlim(0, max(-log10(allRes$adjPclassic)+.5))+
  theme(axis.text.y = element_text(colour="black", size = rel(1.6)))+
  theme(plot.margin=unit(c(1, 0, 2, 3), unit="pt"), legend.position=c(0.25,0.85),
              legend.key.width=unit(2, unit="pt"),
              legend.key.height=unit(8, unit="pt"),
              legend.key = element_blank(),
        legend.text=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2))) # test

dev.off()
