library(rtracklayer)
library(dplyr)
#load NAD-RNA information
nad_df1 <- read.csv("nad_df1.csv",row.names = 1)
head(nad_df1)
#>       GeneID    logFC  Group external_gene_name      chromosome_name
#>1 FBgn0013688 8.697634 Embryo           mt:srRNA mitochondrion_genome
#>2 FBgn0026666 4.872506 Embryo               MagR                    X
#>3 FBgn0038449 4.679704 Embryo            CG17562                   3R
#>4 FBgn0033139 4.326826 Embryo            Tsp42Er                   2R
#>5 FBgn0035904 4.741153 Embryo              GstO3                   3L
#>6 FBgn0036909 3.800327 Embryo             Ccdc58                   3L

#load gene information from database (Drosophila melanogaster genome r6.36)
gtf_data <- import("dmel-all-r6.36.gtf")
gtf_data <- data.frame(gtf_data)
GenesInfo <- gtf_data[gtf_data$type == "gene",c(1,10,2,3,5)]
colnames(GenesInfo) <- c("chr","ID","start","stop","strand")
head(GenesInfo)
#>    chr          ID    start     stop strand
#>1     X FBgn0031081 19961297 19969323      +
#>68    X FBgn0052826 20025099 20025170      +
#>71    X FBgn0031085 20051294 20052519      +
#>79    X FBgn0062565 20094398 20095767      +
#>89    X FBgn0031088 20133579 20138878      +
#>120   X FBgn0041626 20141819 20143188      -

#create a data frame for storing results
random.res <- data.frame(matrix(1, nrow = 10000, ncol = 20))
colnames(random.res) <- c(paste0(dev.point, ".1kb"),paste0(dev.point, ".2kb"),
                          paste0(dev.point, ".3kb"),paste0(dev.point, ".4kb"),
                          paste0(dev.point, ".5kb"))

max.dist <- c(1000,2000,3000,4000,5000)
dev.point <- c("Embryo","Larva","Pupa","D3")
chr <- c("X","Y","2L","2R","3L","3R","4")
for (m in 1:10000){
  for (l in 1:5){
    for (i in 1:4){
      gene.number <- table(nad_df1[nad_df1$Group == dev.point[i],]$chromosome_name)%>% data.frame()
      rownames(gene.number) <- gene.number$Var1
      gene.number <- gene.number[chr,]
      gene.ls <- c()
      for (j in 1:7){
        ran.gene <- subset(GenesInfo$ID, GenesInfo$chr == chr[j],)
        ran.gene <- ran.gene[sample(1:length(ran.gene),gene.number$Freq[j])]
        gene.ls <- c(gene.ls,ran.gene)
      }
      nad.rna <- GenesInfo[GenesInfo$ID %in% gene.ls,]
      nad.rna <- arrange(nad.rna, nad.rna[,1],nad.rna[,3])
      num = 1
      cls.nad <- c()
      for(k in 1:(length(nad.rna$ID)-1)){
        gene.info <- GenesInfo[GenesInfo$ID == nad.rna$ID[k],]
        dist.range <- gene.info$stop + max.dist[l]
        next.gene <- nad.rna[(k+1):min((k+15),length(nad.rna$ID)),]
        next.gene <- subset(next.gene,next.gene$start >= gene.info$start &
                              next.gene$start <= dist.range)
        next.gene <- next.gene[next.gene$chr == nad.rna$chr[k],]
        if (length(next.gene$ID) > 0){
          cls.nad <- c(cls.nad, gene.info$ID, next.gene$ID)
        } else if(length(cls.nad) > 2) {
          if (nad.rna$ID[k] != last(cls.nad)){
            next
          }
          cls.nad <- c()
          num = num +1
        } else{
          cls.nad <- c()
        }
      }
      position = 4*(l-1) + i
      random.res[m, position] <- num-1
    }
  }
  if (m %% 1000 == 0){
    write.csv(random.res, "tmp.csv")
  }
}


write.csv(random.res,"randomres.csv")
