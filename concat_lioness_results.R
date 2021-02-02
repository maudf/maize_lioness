### Maud Fagny
### 2021-02-01
### concat_lioness_results.R
### Concatenate lioness results
###_______________________________

library(RcppCNPy)

### Set Variables
args=commandArgs(trailingOnly=TRUE)
folder=args[1]
networkfolder="Networks/"

### Load files

panda=read.table(paste0(folder, networkfolder, "output_panda.txt"), header=F, sep="\t", stringsAsFactors=F)
colnames(panda)=c("TF", "Gene", "prior", "edge")
TF=panda$TF[panda$Gene==panda$Gene[1]]
genes=unique(panda$Gene)

### Read
files=list.files(path=paste0(folder, networkfolder, "lioness_output/"), pattern="*.npy")
lioness=list()
for (f in 1:length(files)){
  cat("Loading file", files[f], "...\n")
  n=gsub(".npy", '', files[f])
  n=as.numeric(gsub('lioness.', '', n))
  name=rownames(samples.annot)[n]
  lioness[[`name`]]=as.data.frame(t(npyLoad(paste0(folder,"lioness_output/", files[f]))))
  colnames(lioness[[`name`]])=TF
  rownames(lioness[[`name`]])=genes
}

save(lioness, panda, samples.annot, file=paste0(folder,"output_lioness.Rdata"))
