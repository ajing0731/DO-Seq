devtools::install_github("thereallda/enRUVseq")
library(enRUVseq)
library(dplyr)
library(gtools)
library(dplyr)

##find the transcriptome signatures
#read normalized RNA-seq counts
counts.rnaseq <- read.csv("counts.rnaseq.csv",row.names = 1)
meta.rnaseq <- read.csv("meta.rnaseq.csv",row.names = 1)
counts.rnaseq[1:3,]
#>            Embryo_1 Embryo_2 Embryo_3 Larva_1 Larva_2 Larva_3 Pupa_1 Pupa_2 Pupa_3 D3_1 D3_2 D3_3
#>FBgn0031208        8        4        8      12      30      17    105    165    139   72   71   69
#>FBgn0263584        0        0        1       0       0       0      0      1      0    0    0    0
#>FBgn0267987        8       12        5       1       1       0      6      3      1    1    5    5
meta.rnaseq[1:3,]
#>   id condition replicate
#>1 S25    Embryo         1
#>2 S26    Embryo         2
#>3 S27    Embryo         3

#ANOVA-like DE to filter stable-expressed genes during fly development
keep.rnaseq <- filterByExpr(counts.rnaseq, group = meta.rnaseq$condition)
counts.rnaseq.keep <- counts.rnaseq[keep.rnaseq,]
dev.de <- edgeRDE(counts.rnaseq.keep,
                  group = meta.rnaseq$condition,
                  design.formula = as.formula("~condition"),
                  coef = 2:length(unique(meta.rnaseq$condition)))
counts.rnaseq.norm <- cpm(dev.de$de.obj)
res.ls <- dev.de$res.ls[[1]]
res.sig.ls <- subset(res.ls,res.ls$FDR <0.01)
counts.rnaseq.norm <- counts.rnaseq.norm[res.sig.ls$GeneID,]

#transform the expression level of each gene into z-score
rnaseq.mean <- matrix(NA,nrow = length(counts.rnaseq.norm[,1]),ncol = 4)%>%data.frame()
colnames(rnaseq.mean) <- c("Embryo","Larva","Pupa","D3")
for (i in 1:4) {
  j = 3*i-2
  rnaseq.mean[,i] <- apply(counts.rnaseq.norm[,j:(j+2)], 1, mean)
}

rnaseq.mean <- log2(rnaseq.mean+1)
rnaseq.scale <- t(apply(rnaseq.mean, 1, scale))%>%data.frame()
colnames(rnaseq.scale) <- colnames(rnaseq.mean)
rownames(rnaseq.scale) <- rownames(counts.rnaseq.norm)

#pick out transcriptome signatures by z-score > 1.3
Judge <- rep("no", length(rnaseq.scale[,1]))
rnaseq.scale <- cbind(rnaseq.scale, Judge)
for (i in 1:length(Judge)){
     tmp.scale <- rnaseq.scale[i,1:4]
  if (sum(tmp.scale>1.3) == 1) {
     sig.pos <- which(tmp.scale>1.3)
     rnaseq.scale[i,]$Judge <- colnames(rnaseq.scale)[sig.pos]
  }
}
table(rnaseq.scale$Judge)
#>  D3 Embryo  Larva     no   Pupa 
#>1225    554   1089   8188    750

##find NAD-RNA epitranscriptome signatures
#transform the enrichment level of each NAD-RNA into z-score
nad_df1 <- read.csv("nad_df1.csv",row.names = 1)
nad_df1_all <- read.csv("nad_df1_all.csv",row.names = 1)
nad_df1_all <- nad_df1_all[mixedorder(nad_df1_all$GeneID),]
n <- unique(nad_df1_all$GeneID) %>% length()
logfc.scale <- matrix(NA, nrow = n,ncol = 5) %>% data.frame()
colnames(logfc.scale) <- c("GeneID", unique(nad_df1$Group))
for (i in 1:n){
  tmp.fc <- nad_df1_all$logFC[(4*i-3):(4*i)]
  tmp.fc <- t(scale(tmp.fc))
  logfc.scale[i,] <- c(nad_df1_all$GeneID[4*i], tmp.fc)
}

#pick out NAD-RNA epitranscriptome signatures by z-score > 1.3
Judge <- rep("no", length(logfc.scale[,1]))
logfc.scale <- cbind(logfc.scale, Judge)
for (i in 1:length(Judge)){
  tmp.scale <- logfc.scale[i,2:5]
  if (sum(tmp.scale>1.3) == 1) {
    sig.pos <- which(tmp.scale>1.3)
    sig.pos <- colnames(logfc.scale)[sig.pos+1]
    gene.group <- nad_df1[nad_df1$GeneID == logfc.scale$GeneID[i],]$Group
    if (sig.pos %in% gene.group){
      logfc.scale[i,]$Judge <- sig.pos
    }
  }
}
table(logfc.scale$Judge)
#> D3 Embryo  Larva     no   Pupa 
#>266    304    139   3056    441 

write.csv(rnaseq.scale,"rnaseq.scale.csv")
write.csv(logfc.scale,"logfc.scale.csv")
