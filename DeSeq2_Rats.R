setwd("~/OneDrive - KI.SE/Mac/Documents/bulk_RNA/Miguel")
library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(tidyverse)
library(circlize)
rawd <- read.table("counts/all_sg.txt", header= TRUE)
data <- rawd[,-1]
rownames(data) = sub("gene-","",rawd$ENSEMBL_GeneID)
data<-data[which(rowSums(data)!=0),]
data <- data[-c(1:3),]

samples<-factor(colnames(data))
stages<-factor(rep(c("P1","P4","P7","P10","P60"),c(3,3,3,3,2)))
meta<-as.data.frame(cbind(samples,stages))
meta$samples<-as.factor(samples)
meta$stages<-as.factor(stages)
dds <- DESeqDataSetFromMatrix(data, meta, ~ stages)
dds
deseq <- DESeq(dds)
db <- varianceStabilizingTransformation(deseq)
saveRDS(db,"deseq_rats_normalized.rds")

pdf("../Lili/plots/PCA.pdf",width = 10,height = 8)
plotPCA(db, intgroup =c("stages")) #re-ran plotPCA function at the end of the code to plot differnt PCs
dev.off()

#heatmap
db_table <- as.data.frame(assay(db))
db_table_sum <- transform(db_table,sum=rowSums(db_table))
colnames(db_table_sum)
select<- order(db_table_sum$sum,decreasing = TRUE)[1:100]
#heatmap(as.matrix(db_table[select,]),Rowv = F,Colv=F)

ha<-columnAnnotation(Stages=db$stages,gp=gpar(col="white"))
#version1
df<-db_table[select,]
pdf("Lili/heatmap_topgenes.pdf",width = 10,height = 14)
Heatmap(as.matrix(df),cluster_rows = F,
        col=colorRamp2(c(5,10,15,20,25),c("#2c7fb8","#41b6c4","#a1dab4","#ffffcc","#fe9929")),
        heatmap_legend_param=list(
          at=seq(5,25,5),
          labels=c("5","10","15","20","25"),
          title="Expression Level",
          legen_height=unit(4,"cm"),
          title_position="lefttop-rot"
        ),
        top_annotation = ha)

dev.off()

##DE analysis
p1vsp10 <-results(deseq,contrast = c("stages","P1","P10"),alpha=0.05) #example
depair<-list()
stages<-c("P1","P4","P7","P10","P60")
paired_stages<-as.data.frame(t(combn(stages,2)))

for (i in 1:dim(paired_stages)[1]) {
    depair[[paste0(paired_stages[i,1],"vs",paired_stages[i,2])]]<-results(deseq,contrast = c("stages",paired_stages[i,1],paired_stages[i,2]),
                                                                          alpha=0.05)
  }
#save DE markers



sig_DEgenes_id<-list()
sig_DEgenes_id<-lapply(depair, function(x){
  as.factor(abs(x$log2FoldChange) >= 2 & x$padj < 0.05)
})

#plot top features

for (i in 1:length(depair)) {
  g<-ggplot(data=as.data.frame(depair[[i]]), 
                         aes(x=log2FoldChange, y=-log10(padj), 
                             # Colour based on the threshold defined before
                             colour=sig_DEgenes_id[[i]])) +
    scale_x_continuous(name= "log2 fold change",limits=c(-10,10)) +
    # Define the look of the points
    geom_point(alpha=0.7, size=1) +
    # Hide the legend
    theme(legend.position = "none") +
    # Apply another theme
    theme_bw() + theme(legend.position="none") +
    # Add the lines separating the DEGs
    geom_vline(xintercept = 2) +
    geom_vline(xintercept = -2) +
    geom_hline(yintercept = -log10(0.05)) +
    xlab("log2 fold change") + ylab("-log10 padj") +
    ggtitle(paste0(names(depair)[i]," DESeq2\nwith padj <= 0.05 and log2FoldChange >= 2")) +
    geom_text_repel(aes(label=ifelse(padj < 0.1e-6 & abs(log2FoldChange) >= 2,
                                     row.names(depair[[i]]), '')))+
    scale_color_manual(values = c("grey","#31a354"))
  print(g)
  rm(g)
}
dev.off()


#summarize the up/down-regulated features 
sig_DEgenes<-lapply(depair, function(x){
  id<-c()
  id<-(abs(x$log2FoldChange) >= 2 & x$padj < 0.05)
  as.data.frame(x[which(id=="TRUE"),])
})
sig_DEgenes<-lapply(sig_DEgenes, function(x)x[order(x$log2FoldChange,decreasing = T),])

saveRDS(sig_DEgenes,"Lili/Sig_DEgenes_pairs.rds")

library(dplyr)
upregulated_top20<-lapply(sig_DEgenes, function(x){
  if(dim(x)[1]>20){
    x<-x[x$log2FoldChange>2,]
    x %>% top_n(20,log2FoldChange)
  }
  else{
    x 
  }
})

downregulated_bottom20<-lapply(sig_DEgenes, function(x){
  if(dim(x)[1]>20){
    x<-x[x$log2FoldChange<(-2),]
    x %>% top_n(-20,log2FoldChange)
  }
  else{
    x 
  }
})

up_genes<-unique(unlist(lapply(upregulated_top20, function(x){
  rownames(x)
})))

down_genes<-unique(unlist(lapply(downregulated_bottom20, function(x){
  rownames(x)
})))

#p1vsRest<-Reduce(intersect, lapply(sig_DEgenes[c(1:4)], function(x)rownames(x)))

df<-db_table[which(rownames(db_table) %in% c(up_genes,down_genes)),]
df<-df[order(rowSums(df),decreasing = T),]
pdf("../Lili/heatmap_top_bottom_20.pdf",width = 10,height = 14)

ha = HeatmapAnnotation(Stages = db$stages,
                       col = list(Stages = c("P1" = "#d7191c", 
                                             "P4" = "#fdae61", 
                                             "P7" = "#ffffbf",
                                             "P10" = "#abd9e9",
                                             "P60" = "#2c7bb6"
                                             )))
Heatmap(as.matrix(df),cluster_rows = T,
        cluster_columns = F,
        #col=colorRamp2(c(5,10,15,20,25),c("#2c7fb8","#41b6c4","#a1dab4","#ffffcc","#fe9929")),
        heatmap_legend_param=list(
          at=seq(3,15,3),
          labels=seq(3,15,3),
          title="Expression Level",
          legen_height=unit(6,"cm"),
          title_position="lefttop-rot"
        ),
        top_annotation = ha,
        row_names_gp = grid::gpar(fontsize = 4),
        column_names_gp = grid::gpar(fontsize=6))

dev.off()


#top50
upregulated_top50<-lapply(sig_DEgenes, function(x){
  if(dim(x)[1]>50){
    x<-x[x$log2FoldChange>2,]
    x %>% top_n(50,log2FoldChange)
  }
  else{
    x 
  }
})

up_genes_50<-unique(unlist(lapply(upregulated_top50, function(x){
  rownames(x)
})))

df<-db_table[which(rownames(db_table) %in% c(up_genes_50)),]
df<-df[order(rowSums(df),decreasing = T),]
pdf("../Lili/heatmap_top50.pdf",width = 10,height = 14)

ha = HeatmapAnnotation(Stages = db$stages,
                       col = list(Stages = c("P1" = "#d7191c", 
                                             "P4" = "#fdae61", 
                                             "P7" = "#ffffbf",
                                             "P10" = "#abd9e9",
                                             "P60" = "#2c7bb6"
                       )))
Heatmap(as.matrix(df),cluster_rows = T,
        cluster_columns = F,
        #col=colorRamp2(c(5,10,15,20,25),c("#2c7fb8","#41b6c4","#a1dab4","#ffffcc","#fe9929")),
        heatmap_legend_param=list(
          at=seq(3,15,3),
          labels=seq(3,15,3),
          title="Expression Level",
          legen_height=unit(6,"cm"),
          title_position="lefttop-rot"
        ),
        top_annotation = ha,
        row_names_gp = grid::gpar(fontsize = 4),
        column_names_gp = grid::gpar(fontsize=6))

dev.off()

#bottom50
downregulated_top50<-lapply(sig_DEgenes, function(x){
  if(dim(x)[1]>50){
    x<-x[x$log2FoldChange<(-2),]
    x %>% top_n(-50,log2FoldChange)
  }
  else{
    x 
  }
})

down_genes_50<-unique(unlist(lapply(downregulated_top50, function(x){
  rownames(x)
})))
#subset only genes were not included in up_genes_50
id<-intersect(up_genes_50,down_genes_50)
down_genes_50<-down_genes_50[-which(down_genes_50 %in% id)]

df<-db_table[which(rownames(db_table) %in% c(down_genes_50)),]
df<-df[order(rowSums(df),decreasing = T),]
pdf("../Lili/heatmap_bottom50.pdf",width = 10,height = 14)

ha = HeatmapAnnotation(Stages = db$stages,
                       col = list(Stages = c("P1" = "#d7191c", 
                                             "P4" = "#fdae61", 
                                             "P7" = "#ffffbf",
                                             "P10" = "#abd9e9",
                                             "P60" = "#2c7bb6"
                       )))
Heatmap(as.matrix(df),cluster_rows = T,
        cluster_columns = F,
        #col=colorRamp2(c(5,10,15,20,25),c("#2c7fb8","#41b6c4","#a1dab4","#ffffcc","#fe9929")),
        heatmap_legend_param=list(
          at=seq(3,15,3),
          labels=seq(3,15,3),
          title="Expression Level",
          legen_height=unit(6,"cm"),
          title_position="lefttop-rot"
        ),
        top_annotation = ha,
        row_names_gp = grid::gpar(fontsize = 4),
        column_names_gp = grid::gpar(fontsize=6))

dev.off()


#function to plot PCA based on different PCs
plotPCA <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:4]
    return(d)
  }
  #change PCs can plot PCA of differnt PCs
  ggplot(data = d, aes_string(x = "PC2", y = "PC4", color = "group")) + 
    geom_point(size = 3) + xlab(paste0("PC2: ", round(percentVar[2] * 
                                                        100), "% variance")) + ylab(paste0("PC4: ", round(percentVar[4] * 
                                                                                                            100), "% variance")) + coord_fixed()
}

de<-readRDS('Sig_DEgenes_pairs.rds')
write.table(de$P1vsP4,'bulk_P1vsP4_de.txt',sep = '\t',col.names = TRUE,row.names = TRUE)
write.table(de$P1vsP7,'bulk_P1vsP7_de.txt',sep = '\t',col.names = TRUE,row.names = TRUE)
write.table(de$P1vsP10,'bulk_P1vsP10_de.txt',sep = '\t',col.names = TRUE,row.names = TRUE)
write.table(de$P1vsP60,'bulk_P1vsP60_de.txt',sep = '\t',col.names = TRUE,row.names = TRUE)
write.table(de$P4vsP7,'bulk_P4vsP7_de.txt',sep = '\t',col.names = TRUE,row.names = TRUE)
write.table(de$P4vsP10,'bulk_P4vsP10_de.txt',col.names = TRUE,row.names = TRUE)
write.table(de$P4vsP60,'bulk_P4vsP60_de.txt',sep = '\t',col.names = TRUE,row.names = TRUE)
write.table(de$P7vsP10,'bulk_P7vsP10_de.txt',sep = '\t',col.names = TRUE,row.names = TRUE)
write.table(de$P7vsP60,'bulk_P7vsP60_de.txt',sep = '\t',col.names = TRUE,row.names = TRUE)
write.table(de$P10vsP60,'bulk_P10vsP60_de.txt',sep = '\t',col.names = TRUE,row.names = TRUE)


###comparison between bulk and scRNA

epith_bulk_sc=read_csv(file = '../../snRNA/MS/mouse/pipeline_withoutMT/new_analysis_20220309/subcell_analysis/results/Epith_DE_bulk_sc_overlap.csv')
epith_bulk_sc = as.data.frame(epith_bulk_sc)
#pdf("Lili/volcanos_stages.pdf",width = 10,height = 8)

label_name=epith_bulk_sc$...1
i=7
pdf('../Lili/Epith_bulk_scRNA.pdf',width = 8,height = 8)
g<-ggplot(data=as.data.frame(depair[[i]]), 
          aes(x=log2FoldChange, y=-log10(padj), 
              # Colour based on the threshold defined before
              colour=sig_DEgenes_id[[i]])) +
  scale_x_continuous(name= "log2 fold change",limits=c(-10,10)) +
  # Define the look of the points
  geom_point(alpha=0.7, size=1) +
  # Hide the legend
  theme(legend.position = "none") +
  # Apply another theme
  theme_bw() + theme(legend.position="none") +
  # Add the lines separating the DEGs
  geom_vline(xintercept = 2) +
  geom_vline(xintercept = -2) +
  geom_hline(yintercept = -log10(0.05)) +
  xlab("log2 fold change") + ylab("-log10 padj") +
  ggtitle(paste0(names(depair)[i]," DESeq2\nwith padj <= 0.05 and log2FoldChange >= 2")) +
  geom_text_repel(aes(label=ifelse(rownames(depair[[i]]) %in% label_name,row.names(depair[[i]]), ''),size=5))+
  scale_color_manual(values = c("grey","#31a354"))+
  theme( plot.title=element_text(size=16,face="bold"),
         axis.text=element_text(size=18),
         axis.title=element_text(size=20,face="bold"))

print(g)
rm(g)
dev.off()

#immune
immune_bulk_sc=read_csv(file = '../../snRNA/MS/mouse/pipeline_withoutMT/new_analysis_20220309/subcell_analysis/results/Immune_DE_bulk_sc_overlap.csv')
immune_bulk_sc = as.data.frame(immune_bulk_sc)

label_name=immune_bulk_sc$...1
i=7
pdf('../Lili/Immune_bulk_scRNA.pdf',width = 8,height = 8)
g<-ggplot(data=as.data.frame(depair[[i]]), 
          aes(x=log2FoldChange, y=-log10(padj), 
              # Colour based on the threshold defined before
              colour=sig_DEgenes_id[[i]])) +
  scale_x_continuous(name= "log2 fold change",limits=c(-10,10)) +
  # Define the look of the points
  geom_point(alpha=0.7, size=1) +
  # Hide the legend
  theme(legend.position = "none") +
  # Apply another theme
  theme_bw() + theme(legend.position="none") +
  # Add the lines separating the DEGs
  geom_vline(xintercept = 2) +
  geom_vline(xintercept = -2) +
  geom_hline(yintercept = -log10(0.05)) +
  xlab("log2 fold change") + ylab("-log10 padj") +
  ggtitle(paste0(names(depair)[i]," DESeq2\nwith padj <= 0.05 and log2FoldChange >= 2")) +
  geom_text_repel(aes(label=ifelse(rownames(depair[[i]]) %in% label_name,row.names(depair[[i]]), ''),
                      size=5))+
  #geom_label_repel()+
  scale_color_manual(values = c("grey","#31a354"))+
  theme( plot.title=element_text(size=16,face="bold"),
         axis.text=element_text(size=18),
         axis.title=element_text(size=20,face="bold"))+
  geom_label_repel(max.overlaps=20)

print(g)
rm(g)
dev.off()
#max.overlaps = getOption("ggrepel.max.overlaps", default = 10)

