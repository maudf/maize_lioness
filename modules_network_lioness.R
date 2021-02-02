### Maud Fagny
### 2021-02-01
### modules_network_lioness.R
### Identify Tissue-specific modules with ALPACA and perform gene ontology analysis
###_______________________________

###Load libraries
library(ALPACA)
library(condor)
library(topGO)
library(Rgraphviz)

### Set variables
args=commandArgs(trailingOnly=TRUE)
folder=args[1]
max.nb.genes=args[2]
GO.file=args[3]

### Load data
load(paste0(folder, "logtransformed_lioness_results.Rdata"))
load(paste0(folder,"output_lioness.Rdata"))
samples.annot=read.table(paste0(folder, "samples_annotation.txt"),
                         header=T, row.names=1, stringsAsFactors=F)
load(paste0(folder, "../../data/crosstables_TE_gene_enhancers_tfbs.Rdata"))
load(paste0(folder, "normalized_exprs.Rdata"))

### FUnction
alpaca <- function (net.table, file.stem, verbose = T) {
  net.table[, 1] <- paste(net.table[, 1], "A", sep = "_")
  net.table[, 2] <- paste(net.table[, 2], "B", sep = "_")
  print("Detecting communities in control network...")
  ctrl.pos <- net.table[net.table[, 3] >= 0, 1:3]
  ctrl.elist <- data.frame(red = ctrl.pos[, 2], blue = ctrl.pos[, 
                                                                1], weights = ctrl.pos[, 3])
  ctrl.condor <- create.condor.object(ctrl.elist)
  ctrl.condor <- condor.cluster(ctrl.condor, project = F)
  ctrl.memb <- c(ctrl.condor$red.memb[, 2], ctrl.condor$blue.memb[, 
                                                                  2])
  names(ctrl.memb) <- c(as.character(ctrl.condor$red.memb[, 
                                                          1]), as.character(ctrl.condor$blue.memb[, 1]))
  if (!(is.null(file.stem))) 
    write.table(ctrl.memb, paste(c(file.stem, "_ALPACA_ctrl_memb.txt"), 
                                 collapse = ""), row.names = T, col.names = F, quote = F, 
                sep = "\t")
  pos.table <- net.table[intersect(which(net.table[, 3] >= 
                                           0), which(net.table[, 4] >= 0)), ]
  pos.graph <- graph.edgelist(as.matrix(pos.table[, 1:2]), 
                              directed = T)
  if (length(setdiff(V(pos.graph)$name, names(ctrl.memb))) > 
      0) {
    uncounted <- setdiff(V(pos.graph)$name, names(ctrl.memb))
    unc.memb <- sample(1:max(ctrl.memb), length(uncounted), 
                       replace = T)
    names(unc.memb) <- uncounted
    ctrl.memb <- c(ctrl.memb, unc.memb)
  }
  print("Computing differential modularity matrix...")
  dwbm <- computeDWBMmat.mscale(pos.table, ctrl.memb[V(pos.graph)$name])
  if (verbose) {
    write.table(dwbm, paste(c(file.stem, "_DWBM.txt"), collapse = ""), 
                row.names = T, col.names = T, quote = F, sep = "\t")
    write.table(rownames(dwbm), paste(c(file.stem, "_DWBM_rownames.txt"), 
                                      collapse = ""), row.names = T, col.names = T, quote = F, 
                sep = "\t")
    write.table(colnames(dwbm), paste(c(file.stem, "_DWBM_colnames.txt"), 
                                      collapse = ""), row.names = T, col.names = T, quote = F, 
                sep = "\t")
  }
  ntfs <- nrow(dwbm)
  ngenes <- ncol(dwbm)
  this.B <- array(0, dim = c(ntfs + ngenes, ntfs + ngenes))
  this.B[1:ntfs, (ntfs + 1):(ntfs + ngenes)] <- dwbm
  this.B[(ntfs + 1):(ntfs + ngenes), 1:ntfs] <- t(dwbm)
  print("Computing differential modules...")
  louv.memb <- genlouvain(this.B)
  names(louv.memb) <- c(rownames(dwbm), colnames(dwbm))
  print("Computing node scores...")
  louv.Ascores <- NULL
  louv.Bscores <- NULL
  for (i in 1:max(louv.memb)) {
    this.comm <- names(louv.memb)[louv.memb == i]
    this.tfs <- this.comm[grep("_A", this.comm)]
    this.genes <- this.comm[grep("_B", this.comm)]
    if ((length(this.tfs) > 1 ) & (length(this.genes) > 1)) {
      tf.sums <- apply(dwbm[this.tfs, this.genes], 1, sum)
      gene.sums <- apply(dwbm[this.tfs, this.genes], 2, 
                         sum)
    }        else {
      tf.sums <- sum(dwbm[this.tfs, this.genes])
      gene.sums <- dwbm[this.tfs, this.genes]
    }
    this.denom <- sum(dwbm[this.tfs, this.genes])
    louv.Ascores <- c(louv.Ascores, tf.sums/this.denom)
    louv.Bscores <- c(louv.Bscores, gene.sums/this.denom)
  }
  louv.scores <- c(louv.Ascores, louv.Bscores)
  if (!is.null(file.stem)) {
    write.table(cbind(names(louv.memb), as.vector(louv.memb)), 
                paste(c(file.stem, "_ALPACA_final_memb.txt"), collapse = ""), 
                col.names = F, row.names = F, quote = F, sep = "\t")
    write.table(cbind(names(louv.scores), louv.scores), paste(c(file.stem, 
                                                                "_ALPACA_scores.txt"), collapse = ""), row.names = F, 
                col.names = F, quote = F, sep = "\t")
  }
  list(louv.memb, louv.scores)
}

### Compute tissue-specific condor networks
condor.network <- list()
core.tab <- lapply(rownames(edges.bytissue), function(x){return(unlist(strsplit(x,"_")[[1]]))})
core.tab <- data.frame(matrix(unlist(core.tab), ncol=2, byrow=TRUE),
                       stringsAsFactors = FALSE)
colnames(core.tab) <- c("Gene", "TF")
core.tab=core.tab[core.tab$TF %in% expr.tfs,]

for (tis in colnames(edges.bytissue)){
  cat("Finding community structure for tissue ", tis, "\n")
  tab <- cbind(core.tab, "edge"=edges.bytissue[paste(core.tab$Gene, core.tab$TF, sep="_"),tis])
  tab <- tab[tab$edge != Inf,]
  condor.network[[`tis`]] <- create.condor.object(tab)
  condor.network[[`tis`]] <- condor.cluster(condor.network[[`tis`]], project=F)
}

###Stats on networks
nb.com=sapply(condor.network, function(x) nrow(x$Qcoms))
modularity=sapply(condor.network, function(x) max(x$modularity))   

### ALPACA
dir.create(file.path(folder, "Alpaca"), showWarnings = F)
core.tab.alpaca=core.tab[,2:1]
TFs= paste0("tf",1:length(unique(core.tab$TF)))
names(TFs)=unique(core.tab$TF)
genes= paste0("gene", 1:length(unique(core.tab$Gene)))
names(genes)=unique(core.tab$Gene)

core.tab.alpaca=core.tab
core.tab.alpaca$Gene=genes[core.tab.alpaca$Gene]
core.tab.alpaca$TF=TFs[core.tab.alpaca$TF]

## Run ALPACA for V2-IST
husk.V2.edges=data.frame(core.tab[,2:1], "husk"=edges.bytissue[paste(core.tab$Gene, core.tab$TF, sep="_"),"HUSK"], 
                         "V2"=edges.bytissue[paste(core.tab$Gene, core.tab$TF, sep="_"),"V2-IST"], stringsAsFactors = FALSE)
husk.V2=alpaca(husk.V2.edges, file.stem=paste0(folder, "Alpaca/alpaca.results.huskref"), verbose=TRUE)

V2.genes=tapply(names(husk.V2[[1]]), husk.V2[[1]], function(x){gsub("_B", "", x[grep("_B",x)])})
V2.TFs=tapply(names(husk.V2[[1]]), husk.V2[[1]], function(x){gsub("_A", "", x[grep("_A",x)])})
V2.genes.scores=lapply(V2.genes, function(x, scores){
  res=scores[paste0(x, "_B")]
  names(res)=x
  return(sort(res))
}, scores=husk.V2[[2]])
V2.TFs.scores=lapply(V2.TFs, function(x, scores){
  res=scores[paste0(x, "_A")]
  names(res)=x
  return(sort(res))
}, scores=husk.V2[[2]])
V2.topgenes <- lapply(V2.genes.scores, function(x){
  n=names(sort(x, decreasing=T))[1:max.nb.genes]
  return(n[!is.na(n)])})

## Run ALPACA for Husk
V2.husk.edges=data.frame(core.tab[,2:1], "V2"=edges.bytissue[paste(core.tab$Gene, core.tab$TF, sep="_"),"V2-IST_Oka"], 
                         "husk"=edges.bytissue[paste(core.tab$Gene, core.tab$TF, sep="_"),"Husk"], stringsAsFactors = FALSE)
V2.husk=alpaca(V2.husk.edges, file.stem=paste0(folder, "Alpaca/alpaca.results.V2ref"), verbose=TRUE)

husk.genes=tapply(names(V2.husk[[1]]), V2.husk[[1]], function(x){gsub("_B", "", x[grep("_B",x)])})
husk.TFs=tapply(names(V2.husk[[1]]), V2.husk[[1]], function(x){gsub("_A", "", x[grep("_A",x)])})
husk.genes.scores=lapply(husk.genes, function(x, scores){
  res=scores[paste0(x, "_B")]
  names(res)=x
  return(sort(res))
}, scores=V2.husk[[2]])
husk.TFs.scores=lapply(husk.TFs, function(x, scores){
  res=scores[paste0(x, "_A")]
  names(res)=x
  return(sort(res))
}, scores=V2.husk[[2]])
husk.topgenes <- lapply(husk.genes.scores, function(x){
  n=names(sort(x, decreasing=T))[1:max.nb.genes]
  return(n[!is.na(n)])})

save(husk.V2, V2.husk, file=paste0(folder,"Alpaca/results_alpaca.Rdata"))

## Compare ALPACA results between husk and V2-IST
# Compare all genes from each modules
tab=NULL
for(i in 1:length(V2.genes)){
  tmp=c()
  for(j in 1:length(husk.genes)){
    tmp=c(tmp, (sum(V2.genes[[i]] %in% husk.genes[[j]])))  
  }
  tab=cbind(tab, tmp)
}

# Compare max.nb.genestop genes from each modules (max.nb.genes genes that participate the most to tissue-specific modularity)

tab.topgenes=NULL
for(i in 1:length(V2.topgenes)){
  tmp=c()
  for(j in 1:length(husk.topgenes)){
    tmp=c(tmp, (sum(V2.topgenes[[i]] %in% husk.topgenes[[j]])))  
  }
  tab.topgenes=cbind(tab.topgenes, tmp)
}

## Compute Jaccard index 
jaccard=NULL
for(i in 1:length(V2.genes)){
  tmp=c()
  for(j in 1:length(husk.genes)){
    tmp=c(tmp, (sum(V2.genes[[i]] %in% husk.genes[[j]])/length(union(V2.genes[[i]], husk.genes[[j]]))))  
  }
  jaccard=cbind(jaccard, tmp)
}
rownames(jaccard)=paste("husk", 1:67, sep="_")
colnames(jaccard)=paste("V2", 1:71, sep="_")

### Extract Genes in tissue-specific modules and genes in shared modules
husk.unique=apply(jaccard, 1, max)[unlist(lapply(husk.TFs, length))>=1 & unlist(lapply(husk.genes, length))>=10]
V2.unique=apply(jaccard, 2, max)[unlist(lapply(V2.TFs, length))>=1 & unlist(lapply(V2.genes, length))>=10]

find.shared=function(i, j, h=husk.genes.scores, v=V2.genes.scores){
  names(h[[i]][names(h[[i]]) %in% names(v[[j]]))
}
find.unique.husk=function(i, h=husk.genes.scores){
  names(head(sort(h[[i]], decreasing = T),max.nb.genes))
}
find.unique.V2=function(i, v=V2.genes.scores){
    names(head(sort(v[[i]], decreasing = T),max.nb.genes))
}
list.topgenes.shared=list("1.1"=find.shared(i=1,j=1), "2.2"=find.shared(i=2,j=2), "3.3"=find.shared(i=3,j=3), 
                   "4.4"=find.shared(i=4,j=4), "5.6"=find.shared(i=5,j=6), "6.7"=find.shared(i=6,j=7), 
                   "8.10"=find.shared(i=8,j=10), "12.12"=find.shared(i=12,j=12), "13.13"=find.shared(i=13,j=13))
list.topgenes.husk=list("7"=find.unique.husk(i=7),"9"=find.unique.husk(i=9))
list.topgenes.v2=list("5"=find.unique.V2(i=5),"8"=find.unique.V2(i=8),"9"=find.unique.V2(i=9))

### Print files to draw networks with cytoscape
for(i in 14:length(husk.TFs)){husk.TFs[[i]]=NULL; husk.genes[[i]]=NULL}
for(i in 14:length(V2.TFs)){V2.TFs[[i]]=NULL; V2.genes[[i]]=NULL}

V2.net=V2.husk.edges[V2.husk.edges$TF %in% unlist(V2.TFs) & 
                    V2.husk.edges$Gene  %in% unlist(V2.genes), 
                    c("TF", "Gene", "V2")]
husk.net=V2.husk.edges[V2.husk.edges$TF %in% unlist(husk.TFs) & 
                       V2.husk.edges$Gene  %in% unlist(husk.genes), 
                     c("TF", "Gene", "husk")]
gv=lapply(V2.genes, function(x){if(length(x)>max.nb.genes){sample(x, floor(length(x)/10))} else {x}})
gh=lapply(husk.genes, function(x){if(length(x)>max.nb.genes){sample(x, floor(length(x)/10))} else {x}})
gh.com=unlist(lapply(1:length(gh), function(x, g){rep(x, length(g[[x]]))}, g=gh))
names(gh.com)=unlist(gh)
gv.com=unlist(lapply(1:length(gv), function(x, g){rep(x, length(g[[x]]))}, g=gv))
names(gv.com)=unlist(gv)
th.com=unlist(lapply(1:length(husk.TFs), function(x, g){rep(x, length(g[[x]]))}, g=husk.TFs))
names(th.com)=unlist(husk.TFs)
tv.com=unlist(lapply(1:length(V2.TFs), function(x, g){rep(x, length(g[[x]]))}, g=V2.TFs))
names(tv.com)=unlist(V2.TFs)

V2.net=V2.net[V2.net$Gene %in% unlist(gv),]
husk.net=husk.net[husk.net$Gene %in% unlist(gh),]
V2.nodes=data.frame("Node"=c(unlist(V2.TFs) , unlist(gv)), stringsAsFactors = F)
husk.nodes=data.frame("Node"=c(unlist(husk.TFs), unlist(gh)), stringsAsFactors = F)
V2.nodes$type=ifelse(V2.nodes$Node %in% V2.net$TF, "TF", "Gene")
husk.nodes$type=ifelse(husk.nodes$Node %in% husk.net$TF, "TF", "Gene")
V2.nodes$com=c(tv.com, gv.com)
husk.nodes$com=c(th.com, gh.com)

write.table(V2.net, file=paste0(folder, "Network_V2_edges.txt"), row.names=F, sep="\t", quote=F)
write.table(husk.net, file=paste0(folder, "Network_Husk_edges.txt"), row.names=F, sep="\t", quote=F)
write.table(V2.nodes, file=paste0(folder, "Network_V2_nodes.txt"), row.names=F, sep="\t", quote=F)
write.table(husk.nodes, file=paste0(folder, "Network_Husk_nodes.txt"), row.names=F, sep="\t", quote=F)

### Enrichment in Gene Ontology terms for V2-specific

GOenrich=function(geneNames, geneID2GO, list.topgenes, all.genes, tag){
         allRes=list()
allRes[[tag]]=list()
GOtermsGenes=list()
GOtermsGenes[[tag]]=list()
for(i in names(list.topgenes)){
  allRes[[tag]][[i]] <- list()
  if(length(list.topgenes[[i]])>=5){
    myInterestingGenes=unique(list.topgenes[[i]])
    geneList <- factor(as.integer(geneNames[geneNames %in% all.genes] %in% myInterestingGenes))
    names(geneList) <- geneNames[geneNames %in% all.genes]
    #str(geneList)
    GOtermsGenes[["V2"]][[i]]=data.frame("Ontology"=character(0), "GOTerm"=character(0), "Genes"=character(0))
    for(ont in c("BP", "MF")){
      GOdata <- new("topGOdata", ontology = ont, allGenes = geneList, 
                    annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
      ## Compute stats
      test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
      resultElimFisher <- getSigGroups(GOdata, test.stat)
      allRes[[tag]][[i]][[`ont`]] <- GenTable(GOdata, 
                                           elim = resultElimFisher, orderBy = "elim", 
                                           ranksOf = "elim", topNodes = max.nb.genes)      
      allRes[[tag]][[i]][[`ont`]]$ontology=rep(ont, nrow(allRes[[tag]][[i]][[`ont`]]))
      
      pdf(paste0(folder, "Figures/GeneOntology/alpaca/top10GO_B73_up_",tag,"_", ont, "_", i, ".pdf"), width=11.5, height=8)
      showSigOfNodes(GOdata, score(resultElimFisher), firstSigNodes = 10, useInfo ='all')
      dev.off()
    
      myterms <- allRes[[tag]][[i]][[`ont`]]$GO.ID
      mygenes  <- genesInTerm(GOdata, myterms)
      for (j in 1:length(myterms)) {
        myterm <- myterms[j]
        mygenesforterm <- mygenes[myterm][[1]]
        myfactor <- mygenesforterm %in% myInterestingGenes # find the genes that are in the list of genes of interest
        mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
        mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
        GOtermsGenes[[tag]][[i]]=rbind(GOtermsGenes[[tag]][[i]], 
                                    data.frame("Ontology"=ont, "GOTerm"=myterm,"Genes"=mygenesforterm2))
      }
    }
  allRes[[tag]][[i]]=rbind(allRes[[tag]][[i]]$BP, allRes[[tag]][[i]]$MF)
  write.table(allRes[[tag]][[i]][allRes[[tag]][[i]]$Significant >=3 & 
                                    (allRes[[tag]][[i]]$elim <=0.01 | allRes[[tag]][[i]]$adjPclassic <=0.05 | allRes[[tag]][[i]]$adjPparentchild <=0.05),], 
              file=paste0(folder, "Tables/GeneOntology/alpaca/Gene_ontology_enrichment_", tag, "_", i, ".txt"), 
              quote=FALSE,row.names=FALSE, sep="\t")
  write.table(GOtermsGenes[[tag]][[i]], file=paste0(folder, "Tables/GeneOntology/alpaca/Gene_ontology_term_genes_",tag, "_", i, ".txt"), 
              quote=FALSE,row.names=FALSE, sep="\t")
  }
}
         return(list(allRes, GOtermsGenes))
}

  
ZmB73_5a_xref.GO.topGO.FG <- read.delim(GO.file, header=FALSE, stringsAsFactors = F)
geneNames<-ZmB73_5a_xref.GO.topGO.FG$V1
geneID2GO <- readMappings(file = GO.file)

## V2-specific
V2.GO.res=GOenrich(geneNames=geneNames, geneID2GO=geneID2GO, 
         list.topgenes=list.topgenes.v2, all.genes=unique(unlist(V2.genes)), tag="V2")
## Husk-specific
Husk.GO.res=GOenrich(geneNames=geneNames, geneID2GO=geneID2GO, 
         list.topgenes=list.topgenes.husk, all.genes=unique(unlist(husk.genes)), tag="Husk")
## Shared
Shared.GO.res=GOenrich(geneNames=geneNames, geneID2GO=geneID2GO, 
         list.topgenes=list.topgenes.shared, all.genes=unique(unlist(husk.genes)), tag="Shared")

  save(V2.GO.res, Husk.GO.res, Shared.GO.res, file=paste0(folder, "Alpaca/Gene_ontology_results_alpaca.Rdata"))


