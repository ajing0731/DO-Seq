library(biomaRt)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(gghalves)
library(UpSetR)
library(patchwork)
library(paintingr)

#load NAD-RNA list
nad_df1 <- read.csv("nad_df1.csv", row.names = 1)
nad_rna_all <- unique(nad_df1$GeneID)

#fetch drosophila ensembl information
mart <- useMart("ensembl","dmelanogaster_gene_ensembl")
fly_anno <- getBM(attributes = c("external_gene_name","flybase_gene_id","gene_biotype",
                                 "description","chromosome_name"),
                  filters = "flybase_gene_id",values = nad_rna_all, mart = mart)
colnames(fly_anno)[2] <- "GeneID"

## gene type 
nad_df1 <- left_join(nad_df1, fly_anno, by = "GeneID")
gene.type <- nad_df1 %>% dplyr::count(gene_biotype, Group) %>% group_by(gene_biotype)%>%data.frame()
gene.type$gene_biotype <- factor(gene.type$gene_biotype, 
                                 levels = c("protein_coding","ncRNA","pseudogene","snoRNA","tRNA","rRNA","snRNA"))
gene.type <- gene.type[order(gene.type$Group),]
#draw pie chart 
p <- ggpie(gene.type[gene.type$Group == "Embryo",], "n",fill = "gene_biotype", color = "NA")+
  scale_fill_npg()+labs(title = "Embryo")+
  theme(legend.position = "right")

##chromosome distribution
gene.chromo <- nad_df1%>% dplyr::count(Group,chromosome_name)%>%group_by(Group)
gene.chromo$chromosome_name <- gsub("mitochondrion_genome","MT",gene.chromo$chromosome_name)
gene.chromo <- gene.chromo[str_length(gene.chromo$chromosome_name)<3,]
gene.chromo$chromosome_name <- factor(gene.chromo$chromosome_name, 
                                      levels = c("2L","2R","3L","3R","4","X","Y","MT"))
#draw barplot
p <- ggbarplot(gene.chromo[gene.chromo$Group == "Embryo",], x = "chromosome_name", y = "n", 
                fill = "chromosome_name", label = TRUE,legend = "none",
                xlab = "chromosome", ylab = "Numbers of NAD-RNAs", 
                color = "NA",title = "Embryo") +
  scale_fill_npg() + theme(plot.title = element_text(size = 18))

## gene length by fold change group
# retrieve gene length from featrueCounts files
Translen <- read.table('counts.txt', sep='\t', header=T)
gene.length <- inner_join(nad_df1, Translen[,c('Geneid','Length')], by = c('GeneID'='Geneid')) %>% 
  group_by(Group) %>% 
  mutate(quartile=ntile(logFC, 10)) %>%
  ungroup()
gene.length$Group <- factor(gene.length$Group, levels = c("Embryo","Larva","Pupa","D3"))
# half box half point
ggplot(gene.length, aes(x=factor(quartile), y=Length/1000)) +
  geom_half_boxplot(aes(fill=factor(quartile)), color="grey30",
                    outlier.shape = NA) +
  geom_half_point_panel(aes(color=logFC), side = 'r',size=1.5,alpha=.6, stroke=0,
                        position=position_jitter(width=.15)) +
  see::theme_modern() +
  theme(axis.ticks = element_line(color='black'),
        axis.text = element_text(color='black'),
        legend.position='right') +
  facet_wrap(.~Group, nrow=2, scales = 'free_y') +
  theme(strip.text = element_text(size = 18))+
  scale_fill_manual(values=colorRampPalette(c('#e5e5be','#C8B98C',"#567794",'#38678F'))(10), guide='none') +
  scale_color_gradientn(colors=colorRampPalette(c('#e5e5be','#C8B98C',"#567794",'#38678F'))(30),
                        breaks=1:3, labels=1:3, limits=c(0,3.5)) +
  labs(x='Deciles', y='Gene Length(kb)',color='log2FC')

##transcription factor in NAD-RNA
#fetch transcription factor list from AnimalTFDB v4.0
TF.database <- read.table("Drosophila_melanogaster_TF_AnimalTFDBv4.0.txt", sep = "\t", header = T)
TF.nad.rna <- inner_join(TF.database, nad_df1 %>% mutate(Ensembl = GeneID), by = "Ensembl")

a <- subset(TF.nad.rna$GeneID, TF.nad.rna$Group == "Embryo")
b <- subset(TF.nad.rna$GeneID, TF.nad.rna$Group == "Larva")
c <- subset(TF.nad.rna$GeneID, TF.nad.rna$Group == "Pupa")
d <- subset(TF.nad.rna$GeneID, TF.nad.rna$Group == "D3")

TF.list <- list(a,b,c,d)
names(TF.list) <- c('Embryo','Larva','Pupa','D3')
upset(fromList(TF.list),order.by = "freq", 
      mainbar.y.label = "Gene Intersection",
      sets.x.label = "Transcription Factor Number",
      text.scale = c(1.5, 1.5, 1.3, 1.3, 1.5, 1.5),
      sets = c("D3","Pupa","Larva","Embryo"), keep.order = T)

##transcription cofactor in NAD-RNA
#fetch transcription cofactor list from AnimalTFDB v4.0
cofactor.database <- read.table("Drosophila_melanogaster_Cof_AnimalTFDBv4.0.txt",
                                sep = "\t", header = T)
cofactor.nad.rna <- inner_join(cofactor.database, nad_df1 %>% mutate(Ensembl = GeneID), by = "Ensembl")

a <- subset(cofactor.nad.rna$GeneID, cofactor.nad.rna$Group == "Embryo")
b <- subset(cofactor.nad.rna$GeneID, cofactor.nad.rna$Group == "Larva")
c <- subset(cofactor.nad.rna$GeneID, cofactor.nad.rna$Group == "Pupa")
d <- subset(cofactor.nad.rna$GeneID, cofactor.nad.rna$Group == "D3")

cofactor.list <- list(a,b,c,d)
names(cofactor.list) <- c('Embryo','Larva','Pupa','D3')
upset(fromList(cofactor.list),order.by = "freq", 
      mainbar.y.label = "Gene Intersection",
      sets.x.label = "Transcription Cofactor Number",
      text.scale = c(1.5, 1.5, 1.3, 1.3, 1.5, 1.5),
      sets = c("D3","Pupa","Larva","Embryo"), keep.order = T)

##5utr length
utr5 <- getBM(attributes = c( "5utr",
                              "flybase_gene_id"),
              filters = "flybase_gene_id",values = rownames(counts.nsp), mart = mart)
utr5 <- utr5[!utr5$`5utr` == "Sequence unavailable",]
utr5$length <- str_length(utr5$`5utr`)
utr5.length <-utr5 %>% group_by(flybase_gene_id)%>% summarise(length = max(length))
utr5.length$log2length <- log2(utr5.length$length)
Judge <- matrix("non-NAD-RNA", nrow = length(utr5.length$flybase_gene_id),ncol = 4)%>%data.frame()
colnames(Judge) <- unique(nad_df1$Group)
utr5.length <- cbind(utr5.length, Judge)

a <- subset(nad_df1$GeneID, nad_df1$Group == "Embryo")
b <- subset(nad_df1$GeneID, nad_df1$Group == "Larva")
c <- subset(nad_df1$GeneID, nad_df1$Group == "Pupa")
d <- subset(nad_df1$GeneID, nad_df1$Group == "D3")
utr5.length[ utr5.length$flybase_gene_id %in% a,]$Embryo <- "NAD-RNA"
utr5.length[ utr5.length$flybase_gene_id %in% b,]$Larva <- "NAD-RNA"
utr5.length[ utr5.length$flybase_gene_id %in% c,]$Pupa <- "NAD-RNA"
utr5.length[ utr5.length$flybase_gene_id %in% d,]$D3 <- "NAD-RNA"

p1 <- ggboxplot(utr5.length, x = "Embryo",y = "log2length",legend = "none", 
                ylab = "log2(5'UTR length)",fill = "Embryo", 
                palette = paint_palette("Twilight")[5:4]) +
  stat_compare_means(aes(label = ..p.signif..),
                     label.x = 1.5, label.y = 15) +
  stat_summary(fun.y = mean, geom = "point", shape = 16, size =5)

p2 <- ggboxplot(utr5.length, x = "Larva",y = "log2length",legend = "none", 
                ylab = "log2(5'UTR length)",fill = "Larva", 
                palette = paint_palette("Twilight")[5:4]) +
  stat_compare_means(aes(label = ..p.signif..),
                     label.x = 1.5, label.y = 15) +
  stat_summary(fun.y = mean, geom = "point", shape = 16, size =5)+
  rremove("y.axis")+rremove("y.text") +rremove("y.ticks")+rremove("ylab")

p3 <- ggboxplot(utr5.length, x = "Pupa",y = "log2length",legend = "none", 
                ylab = "log2(5'UTR length)",fill = "Pupa", 
                palette = paint_palette("Twilight")[5:4]) +
  stat_compare_means(aes(label = ..p.signif..),
                     label.x = 1.5, label.y = 15) +
  stat_summary(fun.y = mean, geom = "point", shape = 16, size =5)+
  rremove("y.axis")+rremove("y.text") +rremove("y.ticks")+rremove("ylab")

p4 <- ggboxplot(utr5.length, x = "D3",y = "log2length",legend = "none", 
                ylab = "log2(5'UTR length)",fill = "D3", 
                palette = paint_palette("Twilight")[5:4]) +
  stat_compare_means(aes(label = ..p.signif..),
                     label.x = 1.5, label.y = 15) +
  stat_summary(fun.y = mean, geom = "point", shape = 16, size =5)+
  rremove("y.axis")+rremove("y.text") +rremove("y.ticks")+rremove("ylab")
p1+p2+p3+p4+plot_layout(nrow = 1,ncol = 4)