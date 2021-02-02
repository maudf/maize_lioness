### Load libraries
library(umap)
library(RColorBrewer)
library(limma)


### Set variables
args=commandArgs(trailingOnly=TRUE)
folder=args[1]
thres=args[2]
GO.file=args[3]

### Create output folders
dir.create(file.path(folder, "Tables"), showWarnings = F)
dir.create(file.path(paste0(folder, "Figures/"), "Correlations"), showWarnings = F)
dir.create(file.path(paste0(folder, "Figures/"), "GeneOntology"), showWarnings = F)
dir.create(file.path(paste0(folder, "Tables/"), "GeneOntology"), showWarnings = F)
dir.create(file.path(paste0(folder, "Figures/GeneOntology"), "diff_targeted/"), showWarnings = F)
dir.create(file.path(paste0(folder, "Tables/GeneOntology"), "diff_targeted/"), showWarnings = F)

### Load data
load(paste0(folder,"output_lioness.Rdata"))
samples.annot=read.table(paste0(folder, "samples_annotation.txt"),
                         header=T, row.names=1, stringsAsFactors=F)
load(paste0(folder, "Networks/priors.Rdata"))

### Log transform edges weight for lioness
logtrans=function(x){
  y=log(exp(x)+1)
  y[y==Inf] <- x[y==Inf]
  return(y)
}

ind.net.log=lapply(lioness, logtrans)

### Compare edges between prior and non-prior

pdf(paste0(folder, "Figures/edges_by_prior.pdf"), height=11.5, width=8)
par(mar=c(4,5,3,0)+.1, mfrow=c(5,3))
plot(density(panda$edge[panda$prior==0]), col="grey", 
     lwd=2, cex.lab=1.5, cex.axis=1.5, xlim=c(0,40), 
     main="All", xlab="Edge weight")
lines(density(panda$edge[panda$prior==1]), col="red", lwd=2)
for(i in 1:length(ind.net.log)){
  cat("Plotting edges weight distribution for ", names(ind.net.log)[i], "\n")
  plot(density(as.numeric(t(as.matrix(lioness[[i]])))[panda$prior==0]), col="grey", 
       lwd=2, cex.lab=1.5, cex.axis=1.5, 
       main=names(ind.net.log)[i], xlab="Edge weight")
  lines(density(t(as.numeric(as.matrix(lioness[[i]])))[panda$prior==1]), col="red", lwd=2)
}
dev.off()

### filter by prior status
lioness.tab <- NULL
ind.net.log.tab <- NULL
for(i in 1:length(lioness)){
  cat(names(lioness)[i], "\n")
  lioness.tab=cbind(lioness.tab, as.numeric(t(as.matrix(lioness[[i]]))[panda$prior==1]))
  ind.net.log.tab=cbind(ind.net.log.tab, as.numeric(t(as.matrix(ind.net.log[[i]]))[panda$prior==1]))
}
lioness.tab=data.frame(t(lioness.tab))
ind.net.log.tab=data.frame(t(ind.net.log.tab))
colnames(lioness.tab)=colnames(ind.net.log.tab)=paste(panda$Gene, panda$TF, sep="_")[panda$prior==1]
rownames(lioness.tab)=rownames(ind.net.log.tab)=names(lioness)
ind.net.log.tab[,apply(ind.net.log.tab, 2, function(x) length(unique(x)))<=2] <- NULL

### Graphic representation of edges UMAP
custom.config = umap.defaults
if(nrow(ind.net.log.tab)<=15){
  custom.config$n_neighbors=nrow(ind.net.log.tab)-1
}
custom.config$n_neighbors=6
ind.net.log.umap=umap(ind.net.log.tab[,apply(ind.net.log.tab, 2, max)<=2], custom.config)

tissue.color=unlist(tapply(samples.annot$tissue.color, samples.annot$tissue, unique))

pdf(paste0(folder, "Figures/UMAP_lioness_logtransformed.pdf"))
par(mar=c(3,3,5,0)+.1, xpd=TRUE)
plot(ind.net.log.umap$layout[,1], ind.net.log.umap$layout[,2], 
     col=sample.annot[rownames(ind.net.log.umap$layout),"tissue.color"],
     pch=ifelse(samples.annot[rownames(ind.net.log.umap$layout),"tissue"] %in% c("HUSK", "V2-IST"),
                17, 16), xlab="", ylab="",cex.lab=2, cex.axis=2)
legend(par("usr")[1],par("usr")[4]+(par("usr")[4]-par("usr")[3])/7, legend=names(tissue.color), col=tissue.color, 
       pch=ifelse(names(tissue.color) %in% c("HUSK", "V2-IST"), 17, 16), ncol=5, bty='n')
dev.off()

save(lioness.umap, file=paste0(folder, "lioness_UMAP_results.Rdata"))


### Compute matrix of indegree, outdegree, and edges by sample
indegree=matrix(unlist(lapply(ind.net.log, rowSums)), ncol=length(ind.net.log))
colnames(indegree)=names(ind.net.log)
rownames(indegree)=rownames(ind.net.log[[1]])
outdegree=matrix(unlist(lapply(ind.net.log, colSums)), ncol=length(ind.net.log))
colnames(outdegree)=names(ind.net.log)
rownames(outdegree)=colnames(ind.net.log[[1]])
edges.names=unlist(lapply(rownames(ind.net.log[[1]]), paste, 
                          colnames(ind.net.log[[1]]), sep="_" ))
edges=data.frame(t(ind.net.log.tab))

### Compute matrix of average indegree, outdegree, and edge by tissue
tmp=by(ind.net.log.tab, samples.annot[rownames(ind.net.log.tab), "tissue"], 
                                function(x){colSums(x)/nrow(x)})
edges.bytissue=data.frame(matrix(unlist(tmp), 
                                 ncol=length(unique(samples.annot[rownames(ind.net.log.tab), "tissue"]))))
colnames(edges.bytissue)=names(tmp)
rownames(edges.bytissue)=names(tmp[[1]])
rm(tmp)

tmp=indegree.bytissue=by(t(indegree), 
                         tissue.corres[samples.annot[colnames(indegree), "tissue"]],
                         function(x){apply(x, 2, mean)})
indegree.bytissue=data.frame(matrix(unlist(tmp), 
                                    ncol=length(unique(samples.annot[colnames(indegree), "tissue"]))))
colnames(indegree.bytissue)=names(tmp)
rownames(indegree.bytissue)=names(tmp[[1]])
rm(tmp)

tmp=outdegree.bytissue=by(t(outdegree), 
                         tissue.corres[samples.annot[colnames(outdegree), "tissue"]],
                         function(x){apply(x, 2, mean)})
outdegree.bytissue=data.frame(matrix(unlist(tmp), 
                                    ncol=length(unique(samples.annot[colnames(outdegree), "tissue"]))))
colnames(outdegree.bytissue)=names(tmp)
rownames(outdegree.bytissue)=names(tmp[[1]])
rm(tmp)

save(indegree, indegree.bytissue, outdegree, outdegree.bytissue, edges, edges.bytissue, ind.net.log.tab, ind.net.log, file=paste0(folder, "logtransformed_lioness_results.Rdata"))

### Differential edges targeting analysis with limma between Husk and V2-IST

# Prepare data
tissues=samples.annot[rownames(ind.net.log.tab), "tissue"]
limmadiffexpr <- list()
limmadiffexprsignif <- list()
tab=ind.net.log.tab[tissues %in% c("HUSK", "V2-IST"),]
tissues.list=tissues[tissues %in% c("HUSK", "V2-IST")]
tissue=rep(0, nrow(tab))
tissue[tissues.list=="V2-IST"]=1

# Run limma
design <- model.matrix(~ tissue)
fit <- lmFit(t(tab), design=design)
fit <- eBayes(fit)

# Extract significant results
limmadiffexpr[["V2-IST"]]=topTable(fit, number=ncol(tab))
limmadiffexpr[["HUSK"]]=topTable(fit, number=ncol(tab))
limmadiffexprsignif[["V2-IST"]]=limmadiffexpr[["V2-IST"]][limmadiffexpr[["V2-IST"]]$adj.P.Val<=thres & 
                                                        limmadiffexpr[["V2-IST"]]$logFC>0,]
limmadiffexprsignif[["HUSK"]]=limmadiffexpr[["HUSK"]][limmadiffexpr[["HUSK"]]$adj.P.Val<=thres & 
                                                            limmadiffexpr[["HUSK"]]$logFC<0,]

a=lapply(limmadiffexprsignif, rownames)
diff.genes=lapply(a, function(x){unique(unlist(lapply(x, function(y) strsplit(y, "_")[[1]][1])))})
nb.diff.genes=unlist(lapply(genes, function(x) length(unique(x))))
diff.TFs=lapply(a, function(x){unlist(lapply(x, function(y) strsplit(y, "_")[[1]][2]))})
nb.diff.TFs=unlist(lapply(TFs, function(x) length(unique(x))))
save(diff.genes, diff.TFs, limmadiffexprsignif, limmadiffexpr, 
     file=paste0(folder, "limma_differential_targeting_results_thres_", thres, ".Rdata"))

### Gene Ontology analysis of differential targeting:
library(topGO)
library(Rgraphviz)
all.genes=unique(unlist(lapply(colnames(ind.net.log.tab), function(x){strsplit(x, "_")[[1]][1]})))
ZmB73_5a_xref.GO.topGO.FG <- read.delim(GO.file, header=FALSE, stringsAsFactors = F)
geneNames<-ZmB73_5a_xref.GO.topGO.FG$V1
geneID2GO <- readMappings(file = GO.file)

allRes=list()
GOtermsGenes=list()
for(i in 1:length(diff.genes)){
  tis=names(diff.genes)[[i]]
  myInterestingGenes=diff.genes[[i]]
  geneList <- factor(as.integer(geneNames[geneNames %in% all.genes] %in% myInterestingGenes))
  names(geneList) <- geneNames[geneNames %in% all.genes]
  #str(geneList)
  GOtermsGenes[[`tis`]]=data.frame("Ontology"=character(0), "GOTerm"=character(0), "Genes"=character(0))
  for(ont in c("BP", "MF", "CC")){
    GOdata <- new("topGOdata", ontology = ont, allGenes = geneList, 
                  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
    ## Compute stats
    test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
    resultElimFisher <- getSigGroups(GOdata, test.stat)
    
    allRes[[`tis`]][[`ont`]] <- GenTable(GOdata, 
                      elim = resultElimFisher, orderBy = "elim", 
                      ranksOf = "elim", topNodes = 100)
    allRes[[`tis`]][[`ont`]]$ontology=rep(ont, nrow(allRes[[`tis`]][[`ont`]]))
      
    pdf(paste0(folder, "Figures/GeneOntology/topdifftargetededges/topGO_B73_up_", ont, "_", tis, ".pdf"), width=11.5, height=8)
    showSigOfNodes(GOdata, score(resultElimFisher), firstSigNodes = 5, useInfo ='all')
    dev.off()

    myterms <- allRes[[`tis`]][[`ont`]]$GO.ID
    mygenes  <- genesInTerm(GOdata, myterms)
    for (i in 1:length(myterms)) {
      myterm <- myterms[i]
      mygenesforterm <- mygenes[myterm][[1]]
      myfactor <- mygenesforterm %in% myInterestingGenes # find the genes that are in the list of genes of interest
      mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
      mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
      GOtermsGenes[[`tis`]]=rbind(GOtermsGenes[[`tis`]], 
                                  data.frame("Ontology"=ont, "GOTerm"=myterm,"Genes"=mygenesforterm2))
    }
  }
  allRes[[`tis`]]=rbind(allRes[[`tis`]]$BP, allRes[[`tis`]]$MF, allRes[[`tis`]]$CC)
  write.table(allRes[[`tis`]], file=paste0(folder, "Tables/GeneOntology/topdifftargetededges/Gene_ontology_enrichment_", tis, "2.txt"), 
              quote=FALSE,row.names=FALSE, sep="\t")
  write.table(GOtermsGenes[[`tis`]], file=paste0(folder, "Tables/GeneOntology/topdifftargetededges/Gene_ontology_term_genes_", tis, "2.txt"), 
              quote=FALSE,row.names=FALSE, sep="\t")
}
save(allRes, GOtermsGenes, file=paste0(folder, "Gene_ontology_results_topdifftargetededges.Rdata"))

allRes=list()
GOtermsGenes=list()
for(i in 1:length(genes)){
  tis=names(genes)[i]
  myInterestingGenes=unique(genes[[i]])
  geneList <- factor(as.integer(geneNames[geneNames %in% all.genes] %in% myInterestingGenes))
  names(geneList) <- geneNames[geneNames %in% all.genes]
  #str(geneList)
  GOtermsGenes[[`tis`]]=data.frame("Ontology"=character(0), "GOTerm"=character(0), "Genes"=character(0))
  for(ont in c("BP", "MF", "CC")){
    GOdata <- new("topGOdata", ontology = ont, allGenes = geneList, 
                  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
    ## Compute stats
    test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
    resultClassicFisher <- getSigGroups(GOdata, test.stat)
    allRes[[`tis`]][[`ont`]] <- GenTable(GOdata, 
                                         classic = resultClassicFisher, orderBy = "elim", 
                                         ranksOf = "elim", pvalCutOff = 1)
    allRes[[`tis`]][[`ont`]]$adjPclassic=p.adjust(allRes[[`tis`]][[`ont`]]$classic, method="fdr")
    allRes[[`tis`]][[`ont`]]$ontology=rep(ont, nrow(allRes[[`tis`]][[`ont`]]))
    
    pdf(paste0(folder, "Figures/GeneOntology/limma/top10GO_B73_up_", ont, "_", tis, "_thres_",thres, ".pdf"), width=11.5, height=8)
    showSigOfNodes(GOdata, score(resultClassicFisher), firstSigNodes = 10, useInfo ='all')
     dev.off()
    
    myterms <- allRes[[`tis`]][[`ont`]]$GO.ID
    mygenes  <- genesInTerm(GOdata, myterms)
    for (j in 1:length(myterms)) {
      myterm <- myterms[j]
      mygenesforterm <- mygenes[myterm][[1]]
      myfactor <- mygenesforterm %in% myInterestingGenes # find the genes that are in the list of genes of interest
      mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
      mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
      GOtermsGenes[[`tis`]]=rbind(GOtermsGenes[[`tis`]], 
                                  data.frame("Ontology"=ont, "GOTerm"=myterm,"Genes"=mygenesforterm2))
    }
  }
  allRes[[`tis`]]=rbind(allRes[[`tis`]]$BP, allRes[[`tis`]]$MF, allRes[[`tis`]]$CC)
  write.table(allRes[[`tis`]][allRes[[`tis`]]$Significant >=5 & 
                                ( allRes[[`tis`]]$adjPclassic <=0.05 ),], 
              , file=paste0(folder, "Tables/GeneOntology/limma/Gene_ontology_enrichment_", tis, "_thres_",thres, ".txt"), 
              quote=FALSE,row.names=FALSE, sep="\t")
  write.table(GOtermsGenes[[`tis`]], file=paste0(folder, "Tables/GeneOntology/limma/Gene_ontology_term_genes_", tis, "_thres_",thres, ".txt"), 
              quote=FALSE,row.names=FALSE, sep="\t")
}
save(allRes, GOtermsGenes, file=paste0(folder, "Gene_ontology_results_limma_thres_",thres, ".Rdata"))


allRes=list()
GOtermsGenes=list()
for(i in 1:nrow(diff.target.genes)){
  tis=rownames(diff.target.genes)[i]
  myInterestingGenes=colnames(diff.target.genes)[as.numeric(diff.target.genes[i,])>1]
  geneList <- factor(as.integer(geneNames[geneNames %in% all.genes] %in% myInterestingGenes))
  names(geneList) <- geneNames[geneNames %in% all.genes]
  #str(geneList)
  GOtermsGenes[[`tis`]]=data.frame("Ontology"=character(0), "GOTerm"=character(0), "Genes"=character(0))
  for(ont in c("BP", "MF", "CC")){
    GOdata <- new("topGOdata", ontology = ont, allGenes = geneList, 
                  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
    ## Compute stats
    test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
    resultClassicFisher <- getSigGroups(GOdata, test.stat)
    test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
    resultElimFisher <- getSigGroups(GOdata, test.stat)
    test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test")
    resultWeightFisher <- getSigGroups(GOdata, test.stat)
    allRes[[`tis`]][[`ont`]] <- GenTable(GOdata, 
                                         classic = resultClassicFisher, weight = resultWeightFisher, 
                                         elim = resultElimFisher, orderBy = "weight", 
                                         ranksOf = "elim", topNodes = 100)
    allRes[[`tis`]][[`ont`]]$ontology=rep(ont, nrow(allRes[[`tis`]][[`ont`]]))
    
    pdf(paste0(folder, "Figures/GeneOntology/diff_targeted/topGO_B73_up_", ont, "_", tis, ".pdf"), width=11.5, height=8)
    showSigOfNodes(GOdata, score(resultClassicFisher), firstSigNodes = 5, useInfo ='all')
    showSigOfNodes(GOdata, score(resultWeightFisher), firstSigNodes = 5, useInfo ='all')
    showSigOfNodes(GOdata, score(resultElimFisher), firstSigNodes = 5, useInfo ='all')
    dev.off()
    
    myterms <- allRes[[`tis`]][[`ont`]]$GO.ID
    mygenes  <- genesInTerm(GOdata, myterms)
    for (i in 1:length(myterms)) {
      myterm <- myterms[i]
      mygenesforterm <- mygenes[myterm][[1]]
      myfactor <- mygenesforterm %in% myInterestingGenes # find the genes that are in the list of genes of interest
      mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
      mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
      GOtermsGenes[[`tis`]]=rbind(GOtermsGenes[[`tis`]], 
                                  data.frame("Ontology"=ont, "GOTerm"=myterm,"Genes"=mygenesforterm2))
    }
  }
  allRes[[`tis`]]=rbind(allRes[[`tis`]]$BP, allRes[[`tis`]]$MF, allRes[[`tis`]]$CC)
  write.table(allRes[[`tis`]], file=paste0(folder, "Tables/GeneOntology/diff_targeted/Gene_ontology_enrichment_", tis, ".txt"), 
              quote=FALSE,row.names=FALSE, sep="\t")
  write.table(GOtermsGenes[[`tis`]], file=paste0(folder, "Tables/GeneOntology/diff_targeted/Gene_ontology_term_genes_", tis, ".txt"), 
              quote=FALSE,row.names=FALSE, sep="\t")
}
save(allRes, GOtermsGenes, file=paste0(folder, "Gene_ontology_results_topdifftargetedgenes.Rdata"))

