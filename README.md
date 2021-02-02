# Building and analyzing enhancer-driven TF-gene regulatory networks
Updated on 2021-02-01

## Introduction
These scripts were used to build TF-gene regulatory network from maize enhancers from Oka et al., Genome Biology, 2018.

The steps 1-4 are applicable to any datasets. The steps 5-8 have been
tailored to analyse husk vs. V2-IST networks. To apply it to your
dataset, you will need to modify the scripts.

To run this pipeline:
* install all required softwares and packages (see
**List of required softwares and packages** below).
* Create a main folder (eg.: ~/Documents/Networks/) to work from
* Place in a folder (eg.: ~/Documents/Networks/) the files listed in
**List of Files:**.



### List of Files:
* genome.fasta : FASTA file containing the genome sequence
* genes.gff3 : GFF3 file containing genes annotation in the genome
* GeneOntologyAnnotation.txt : TXT file with Gene Ontology annotation of the genes
* TransposableElements.gff : GFF file with transposable elements annotation for the genome
* TransposableElementsAnnotation.txt:  tab-separated TXT file containing transposable elements annotations, including 3 columns: TE name (TE), superfamily (Superfamily), and family (consensus)
* enhancers.bed : BED file containing the enhancers coordinates (chr, start, end), the enhancer name (name) and an optional column corresponding to tissue (here "HUSK", "V2-IST or "Both") that is necessary for step 7.
* motifs.meme : MEME file containing the transcription factor binding sites motifs
* expr\_raw.txt : tab-separated TXT file containing raw expression counts (genes in row and samples in columns)
* samples\_annotation.txt :  tab-separated TXT file containing samples annotation (in my case, tissue ID (tissue), replicate ID (rep), sequencing batch ID (batch), color to plot the tissue (tissue.color), color to plot the batch (batch.color))

### List of required softwares and packages
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [FIMO](http://gensoft.pasteur.fr/docs/meme/5.0.4/fimo.html) from the MEME suite
* [PyPanda](https://github.com/aless80/pypanda)
* [R](https://cran.r-project.org/) with [bioconductor](https://www.bioconductor.org/) and the following packages:
    * [purrr](https://cran.r-project.org/web/packages/purrr/index.html)
    * [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html)
    * [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
    * [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html)
    * [YARN](https://bioconductor.org/packages/release/bioc/html/yarn.html)
    * [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)
    * [corrplot](https://cran.r-project.org/web/packages/corrplot/index.html)
    * [condor](https://rdrr.io/github/jplatig/condor/)
    * [igraph](https://igraph.org/r/)
    * [umap](https://cran.r-project.org/web/packages/umap/index.html)
    * [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
    * [ALPACA](https://github.com/meghapadi/ALPACA/blob/master/DESCRIPTION)
    * [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html)
    * [Rgraphviz](https://bioconductor.org/packages/release/bioc/html/Rgraphviz.html)

## Step 1 -Map Transcription Factor Binding Sites (TFBS) to enhancers
To obtain enhancer sequences and map TFBS, run:

```
$ cd ~/Documents/Networks/
$ bedtools  getfasta -fo genome.fasta  -fi enhancers.fasta -bed enhancers.bed
$ fimo --o FIMO/ motifs.meme enhancers.fasta 1>FIMO/fimo.out 2>FIMO/fimo.err
```

## Step 2 - Normalize gene expression with YARN and plot principal component analyses results.
To normalize expression data and plot results for PC1 and PC2, run:
```
Rscript normalize_expression_data.R ./ 3
```
This R script takes 2 arguments as input:
* the path to the folder containing data
* the threshold for gene expression filtering. A gene is considered expressed if it has a raw counts >0 in at least *threshold* number of samples.

Results will be stored in the main folder and PCA will be plotted in Figures/

## Step 3 - Build prior from TFBS results, enhancers and gene coordinates and map Transposable elements.
To build a prior network that include all possible regulatory
relationships between TFBS within enhancers and candidate target
genes, run:
```
Rscript build_prior.R ./ 250000 0.01
```
This R script takes 3 arguments as input:
* the path to the folder containing data.
* the size of the window to search for potential target genes.
* the qvalue threshold to filter TFBS results. Only TFBS with a Benjamini-Hochberg corrected p-value (q-value in FIMO results) under this threshold will be kept.

To build an annotation table with enhancer, TFBS, all candidate target
genes, and Transposable elements overlapping TE, run:
```
Rscript make_enhancers_tfbs_gene_te_table.R ./ 250000 0.01
```
This R script takes 3 arguments as input:
* the path to the folder containing data.
* the size of the window to search for potential target genes.
* the qvalue threshold to filter TFBS results. Only TFBS with a Benjamini-Hochberg corrected p-value (q-value in FIMO results) under this threshold will be kept.

## Step 4 - Build networks from prior using PANDA/LIONESS.
```
ln -s expr.txt Network/expr.txt
cd Networks/
python run_panda.py -e ./expr.txt -m ./priors.txt -o output_panda.txt -q output_lioness.txt
cd ../
Rscript concat_lioness_results.R ./ 
```
The PANDA/LIONESS software takes 2 arguments as input:
* path to the file with expression data
* path to the file with network prior

The R scripts take 1 argument as input:
* path to the main folder

Results will be stored in the Network/ folder

## Step 5 - PANDA/LIONESS networks: differential targeting analysis.
To obtain the list of differentially targeted edges between two
tissues (here Husk and V2-IST), run:
```
Rscript differential_targeting.R ./ 0.01 GeneOntologyAnnotation.txt
```
This R script takes 3 arguments in input:
* the path to the main folder,
* the Benjamini-Hocheberg corrected p-value threshold for differential targeting analysis
* The name of the file containing the Gene Ontology annotations for
all genes.

## Step 6 - PANDA/LIONESS networks: tissue-specific module identification 
To identify the modules in each tissue-specific network and analyse
gene content of tissue-specific and shared modules, run:
``` 
Rscript modules_network_lioness.R ./ 100 GeneOntologyAnnotation.txt
```

This R script takes 3 arguments in input:
* the path to the main folder,
* the number of top genes to retain in each module for the Gene
Ontology enrichment analysis
* The name of the file containing the Gene Ontology annotations for
all genes.

## Step 7 - Identification of each enhancer top target

To identify enhancers most likely target genes, run:
``` 
Rscript find_enhancer_topgenes_pairs.R ./
```

This R script takes 1 argument in input:
* the path to the main folder.


## Step 8 - TE superfamilies and families enrichment analysis
To perform transposable elements superfamilies and families enrichment
analyses in enhancer and TFBS sequences in a tissue-specific function,
and obtain Gene Ontlogy enrichment analysis for candidate target genes, run:
``` 
Rscript cross_enhancers_genes_tfbs_TE.R ./
```

This R script takes 1 argument in input:
* the path to the main folder.

