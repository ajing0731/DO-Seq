# install enONE v0.0.1——enRUVseq
devtools::install_github("thereallda/enRUVseq")
library(enRUVseq)
library(tidyverse)
library(edgeR)
library(paintingr)
library(patchwork)
library(SummarizedExperiment)

#load sample information
meta <- read.csv('metadata.csv', comment.char = '#')
meta[1:3,]
#>  id    condition replicate
#>1 S1 embryo.Input         1
#>2 S3 embryo.Input         2
#>3 S4  larva.Input         1

#load sequencing data
counts.df <- read.csv('Counts.csv', row.names = 1)
counts.df[1:3,]
#>            embryo.Input_1 embryo.Input_2 larva.Input_1 larva.Input_2 pupa.Input_1 pupa.Input_2
#>FBgn0031208              0              0             7             0           48           72
#>FBgn0263584              0              0             0             0            0            0
#>FBgn0267987              1              0             0             0            1            0
#>            D3.Input_1 D3.Input_2 embryo.Enrich_1 embryo.Enrich_2 larva.Enrich_1 larva.Enrich_2
#>FBgn0031208         30         31               0               3              2              5
#>FBgn0263584          4          0               0               0              0              0
#>FBgn0267987          5          1               0               4              0              0
#>            pupa.Enrich_1 pupa.Enrich_2 D3.Enrich_1 D3.Enrich_2
#>FBgn0031208           109           168          70          81
#>FBgn0263584             0             0           0           0
#>FBgn0267987             0             2           4           3

keep <- filterByExpr(counts.df, group = meta$condition)
counts_keep <- counts.df[keep,]
enrich_group <- gsub(".*\\.", "", meta$condition)
spikeInPrefix <- '^ENS'

# create Enone
Enone <- createEnone(data = counts_keep,
                     bio.group = meta$condition,
                     enrich.group = enrich_group,
                     batch.group = NULL,
                     spike.in.prefix = spikeInPrefix,
                     input.id = "Input",
                     enrich.id = "Enrich")
#>class: Enone 
#>dim: 20616 16 
#>metadata(0):
#>  assays(1): ''
#>rownames(20616): FBgn0031208 FBgn0067779 ... ENSMUSG00000095041 ENSMUSG00000095742
#>rowData names(3): GeneID SpikeIn Synthetic
#>colnames(16): embryo.Input_1 embryo.Input_2 ... D3.Enrich_1 D3.Enrich_2
#>colData names(5): id condition enrich replicate batch

Enone <- enONE(Enone, 
               ruv.norm = TRUE, ruv.k = 3,
               pam.krange = 2:6, pc.k = 3,
               n.neg.eval = 1000, n.pos.eval = 1000)
#>The number of negative control genes for RUV: 1000 
#>The number of positive evaluation genes: 1000 
#>The number of negative evaluation genes: 1000 
#>Apply normalization...
#>Perform assessment...

enScore <- getScore(Enone)
pca.eval <- prcomp(enScore[,-c(3, 9)], scale = TRUE)
ggPCA_Biplot(pca.eval, score = enScore$SCORE)

# select normalization method
norm.method <- rownames(enScore[5,])
Enone <- UseNormalization(Enone, slot = 'sample', method = norm.method)
norm.data <- Counts(Enone, slot = 'sample', method = norm.method)
norm.factors <- getFactor(Enone, slot = 'sample', method = norm.method)
norm.method
#> [1] "DESeq_RUVs_k2"

#draw PCA plot
samples_name <- paste(meta$condition, meta$replicate, sep='.')
p1 <- ggPCA(log1p(Counts(Enone, slot='sample', 'Raw')), 
            color = meta$condition,
            label = samples_name, vst.norm = FALSE) + ggtitle('Before normalization')
p2 <- ggPCA(log1p(norm.data), 
            color = meta$condition,
            label = samples_name, vst.norm = FALSE) + ggtitle('After normalization')
p1 + p2

#set logFC>1 and padj<0.05 as cutoff
Enone <- FindEnrichment(Enone, slot='sample', method = norm.method, 
                        logfc.cutoff = 1, p.cutoff = 0.05)
res.best.ls <- getEnrichment(Enone, slot='sample', filter=TRUE)
unlist(lapply(res.best.ls, nrow))
#>embryo.Enrich_embryo.Input   larva.Enrich_larva.Input     pupa.Enrich_pupa.Input 
#>                      1611                       1754                       3224 
#>D3.Enrich_D3.Input 
#>              2752 

#get the list of identified NAD-RNAs
nad_df1 <- reduceRes(res.best.ls, fc.col = 'logFC')
nad_df1$Group <- gsub('\\..*', '', nad_df1$Group)
nad_df1$Group <- factor(nad_df1$Group, levels = unique(nad_df1$Group))
nad_df1[1:3,]
#>       GeneID    logFC  Group
#>1 FBgn0013688 8.697634 embryo
#>2 FBgn0026666 4.872506 embryo
#>3 FBgn0038449 4.679704 embryo