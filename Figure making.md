# Visualization and Analysis of Hematopoietic Cell Populations in scRNA-seq Data

This R script processes and visualizes single-cell RNA sequencing (scRNA-seq) data to investigate hematopoietic cell populations. It includes cell annotation refinement, dimensional reduction plots (openTSNE, UMAP), gene expression analysis, stemness and cell cycle profiling, and differential gene expression analysis. The script leverages various bioinformatics packages such as Seurat, monocle, tricycle, and clusterProfiler to analyze hematopoietic stem and progenitor cells (HSPCs), myeloid, and lymphoid lineages under different conditions, including chemotherapy and CDK inhibition. The figures generated illustrate cellular distributions, gene expression dynamics, and functional implications of stemness and differentiation.

~~~R
suppressPackageStartupMessages({
    library(dplyr)
    library(Seurat)
    library(Matrix)
    library(proxy)
    library(gplots)
    library(Rtsne)
    library(densityClust)
    library(irlba)
    library(monocle)
    library(plyr)
    library(DOSE)
    library(clusterProfiler)
    library(topGO)
    library(pathview)
    library(AnnotationDbi)
    library(cowplot)
    library(ggplot2)
    library(velocyto.R)
    library(trqwe)
    library(Rsamtools)
    library(GenomicFeatures)
    library(GenomicAlignments)
    library(BiocParallel)
    library(pheatmap)
    library(RColorBrewer)
    library(PoiClaClu)
    library(org.Mm.eg.db)
    library(org.Hs.eg.db)
    library(DESeq2)
    library(data.table)
    library(stringr)
    library(iTALK)
    library(nichenetr)
    library(tidyr)
    library(GenomicRanges)
    library(viridis)
    library(chromVAR)
    library(ggpubr)
    library(corrplot)
    library(SingleCellExperiment)
    library(scater)
    library(flexmix)
    library(splines)
    library(biomaRt)
    library(miQC)
    library(scales)
    library(BuenColors)
    library(PCAtools)
    library(Seurat)
	library(SeuratData)
	library(SeuratWrappers)
	library(flexmix)
})
source("/mnt/d/xiangyu.ubuntu/code/log-summery/MyBestFunction_scRNA.R.v4.R")
library(future)
library(future.apply)
options(future.globals.maxSize = 200 * 1024^3)
plan("multicore", workers = 60)
plan()
require(rJava)
require(xlsx)
~~~



~~~r
pool2 <- mcreadRDS("./All_filter_miQC.pool2.merge.Sv4.rds",mc.cores=20)
pool2$v2_Cell_annotation <- as.character(pool2$Cell_annotation)
pool2$v2_Cell_annotation[pool2$Cell_annotation %in% c("Mmp8+mNeu","Retnlg+mNeu")] <- "mNeu"
pool2$v2_Cell_annotation[pool2$Cell_annotation %in% c("Camp+immNeu","Fcnb+immNeu","Fmnl2+immNeu","Hist1h3c+immNeu","Ifitm6+immNeu","Ltf+immNeu","Pclaf+immNeu","Top2a+immNeu","Ube2c+immNeu")] <- "immNeu"
pool2$v2_Cell_annotation[pool2$Cell_annotation %in% c("Apoe+Mono.P","S100a4+Mono.P")] <- "Mono.P"
pool2$v2_Cell_annotation[pool2$Cell_annotation %in% c("Elane+GMP","F13a1+GMP")] <- "GMP"
pool2$v2_Cell_annotation[pool2$Cell_annotation %in% c("Stfa3+PMNa","Ifitm1+PMNa","Isg15+PMNb")] <- "PMN"
pool2$v2_Cell_annotation[pool2$Cell_annotation %in% c("Vpreb1+Early.B","Vpreb3+Early.B")] <- "Early.B"
pool2$v2_Cell_annotation <- factor(pool2$v2_Cell_annotation,levels=c("MPP","MEP","Ery","CLP","GMP","immNeu","mNeu","PMN","Mono.P","DC","NK.T","Early.B","mature.B","Plasma"))

col <- jdb_palette("corona")[c(2,1,3,4:7,9:length(jdb_palette("corona")))]
cell_sub <- c("MPP","MEP","Ery","CLP","GMP","immNeu","mNeu","PMN","Mono.P","DC","NK.T","Early.B","mature.B","Plasma")
col <- col[1:length(cell_sub)]
names(col) <- cell_sub
p1 <- DimPlot(object = pool2, reduction = "openTSNE",group.by="v2_Cell_annotation",label=TRUE,raster=FALSE,cols= col,pt.size=.1) +labs(title="Cell_annotation")
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig1B.openTSNE_pos.svg", plot=p1,width = 6, height = 5,dpi=300)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig1B.openTSNE_pos.pdf", plot=p1,width = 6, height = 5,dpi=300)
~~~

![image-20250226144938832](./Figure%20making.assets/image-20250226144938832.png)

~~~r
p1 <- DimPlot(object = pool2, reduction = "openTSNE",group.by="Cell_annotation",label=FALSE,raster=FALSE,cols= jdb_palette("corona"),pt.size=.1,split.by="condition") +labs(title="Cell_annotation")
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/pool2.openTSNE_pos.svg", plot=p1,width = 35, height = 5,dpi=300)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/pool2.openTSNE_pos.png", plot=p1,width = 35, height = 5,dpi=300)
~~~

![image-20250226145043339](./Figure%20making.assets/image-20250226145043339.png)

~~~r
scRNA.hematopoiesis <- readRDS(file="./Fig_New_T00_obj.markers.rds")
filter_scRNA.hematopoiesis <- scRNA.hematopoiesis[scRNA.hematopoiesis$p_val_adj<=0.05,]
Stemness <- filter_scRNA.hematopoiesis[filter_scRNA.hematopoiesis$cluster=="HSPC",]
Differentiated <- filter_scRNA.hematopoiesis[filter_scRNA.hematopoiesis$cluster=="Neutro",]
Lineage_marker <- intersect(rownames(GetAssayData(object = pool2, slot = "data")),Stemness$gene)
Stemness1 <- Stemness[Stemness$gene %in% Lineage_marker,]
write.csv(Stemness1,"./Stemness.markers.csv")
speci_raw <- FetchData(object = pool2, vars = Lineage_marker,slot="data")
pool2[["Stemness"]] <- (rowSums(speci_raw))/length(Lineage_marker)
Sel_sig <- c("Stemness")
All_sum <- as.data.frame(FetchData(object = pool2, vars = c(Sel_sig,"condition","Cell_annotation"),slot="data"))
only_MPP <- All_sum[All_sum$Cell_annotation %in% c("MPP"),]
p1 <- ggboxplot(only_MPP[only_MPP$condition %in% c("Vav_vehicle","Vav_carbo","Vav_trila","Vav_combo"),], x = "condition", y = "Stemness", fill="condition",title=paste0("Stemness",".exp"), legend = "none",outlier.shape = NA,notch = FALSE) +rotate_x_text(angle = 45)+
stat_summary(fun.y = median, geom="point",colour="darkred", size=3) + stat_summary(fun = median, geom = "line",aes(group = 1),col = "red",size=1)+
stat_compare_means(comparisons =list(c("Vav_vehicle","Vav_carbo"),c("Vav_vehicle","Vav_trila"),c("Vav_vehicle","Vav_combo"),c("Vav_carbo","Vav_combo")),method = "wilcox.test")
p2 <- ggboxplot(only_MPP[only_MPP$condition %in% c("TP53_veh","TP53_carbo","TP53_trila","TP53_combo"),], x = "condition", y = "Stemness", fill="condition",title=paste0("Stemness",".exp"), legend = "none",outlier.shape = NA,notch = FALSE) +rotate_x_text(angle = 45)+
stat_summary(fun.y = median, geom="point",colour="darkred", size=3) + stat_summary(fun = median, geom = "line",aes(group = 1),col = "red",size=1)+
stat_compare_means(comparisons =list(c("TP53_veh","TP53_carbo"),c("TP53_veh","TP53_trila"),c("TP53_veh","TP53_combo"),c("TP53_carbo","TP53_combo")),method = "wilcox.test")
plot <- plot_grid(p1,p2)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig1C.boxplots.svg", plot=plot,width = 6, height = 5,dpi=300)
~~~

![image-20250226145119603](./Figure%20making.assets/image-20250226145119603.png)

~~~r
Idents(pool2) <- pool2$v2_Cell_annotation
only_MPP <- subset(pool2,idents=c("MPP"))
library(tricycle)
only_MPP <- Runtricycle(object = only_MPP, slot = "data", reduction.name = "tricycleEmbedding", 
    reduction.key = "tricycleEmbedding_", gname = NULL, gname.type = "SYMBOL", species = "mouse",center.pc1=0.5,center.pc2=0.5)
only_MPP$tricycleGroup <- "G01"
only_MPP$tricycleGroup[only_MPP$tricyclePosition > 0.5*pi & only_MPP$tricyclePosition <= 1*pi] <- "S"
only_MPP$tricycleGroup[only_MPP$tricyclePosition > 1*pi & only_MPP$tricyclePosition <= 1.5*pi] <- "G2M"
only_MPP$tricycleGroup[only_MPP$tricyclePosition > 1.5*pi & only_MPP$tricyclePosition <= 1.75*pi] <- "M"
aa <- as.data.frame(table(only_MPP$condition,only_MPP$tricycleGroup))
aa <- aa[order(aa$Var2),]
aa$Var1 <- as.character(aa$Var1)
aa_all <- c()
for (i in unique(aa$Var1)){
  group_sel <- subset(aa,Var1==i)   
  group_sel$sum_number <- sum(group_sel$Freq)
  group_sel$normal_ratio <- (group_sel$Freq/group_sel$sum_number)*100
  group_sel$normal_ratio <- round(group_sel$normal_ratio,2)
  aa_all <- rbind(aa_all,group_sel)
}
aa_all <- aa_all[aa_all$Var1 %in% c("TP53_veh","TP53_carbo","TP53_trila","TP53_combo"),]
aa_all$Var1 <- factor(aa_all$Var1,levels=c("TP53_veh","TP53_carbo","TP53_trila","TP53_combo"))
aa_all <- aa_all[order(aa_all$Var1),]
library(ggalluvial)
p1 <- ggplot(aa_all, aes(x = Var1, y = normal_ratio, fill = Var2, stratum = Var2, alluvium = Var2)) + geom_stratum(width = 0.75) + geom_flow(alpha = 0.5) + 
theme_classic() + theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1)) +labs(x = '', y = 'Relative Abundance(%)',title="MPP.CD45.2.scRNA")+scale_fill_manual(values = jdb_palette("corona"))+scale_color_manual(values = jdb_palette("corona"))
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig1D.cellcycle.svg", plot=p1,width = 4, height = 5,dpi=300)
~~~

![image-20250226145146087](./Figure%20making.assets/image-20250226145146087.png)

~~~r
col <- jdb_palette("corona")[c(2,1,3,4:7,9:length(jdb_palette("corona")))]
cell_sub <- c("MPP","MEP","Ery","CLP","GMP","immNeu","mNeu","PMN","Mono.P","DC","NK.T","Early.B","mature.B","Plasma")
col <- col[1:length(cell_sub)]
names(col) <- cell_sub
col <- c(col,"#e1dfdd")
names(col)[length(col)] <- "OTS"
pool2$new_anno5 <- as.character(pool2$v2_Cell_annotation)
pool2$new_anno5[!pool2$v2_Cell_annotation %in% c("MPP","GMP","immNeu","mNeu","PMN","Mono.P")] <- "OTS"
pool2$new_anno5 <- factor(pool2$new_anno5,levels=c("OTS","MPP","GMP","immNeu","mNeu","PMN","Mono.P"))
p1 <- DimPlot(pool2, reduction = 'openTSNE', label = FALSE,repel=FALSE, group.by="new_anno5",cols=col[levels(pool2$new_anno5)],pt.size=1)+labs(title="Myeolid")
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig1F.openTSNE_pos.png", plot=p1,width = 6, height = 5,dpi=300)
~~~

![image-20250226145218152](./Figure%20making.assets/image-20250226145218152.png)

~~~R
only_Myeolid_seurat_seurat <- mcreadRDS("./only_Myeolid.pool2.merge.Sv4.rds",mc.cores=20)
only_Myeolid_seurat_seurat$v2_Cell_annotation[only_Myeolid_seurat_seurat$Cell_annotation %in% c("Mmp8+mNeu","Retnlg+mNeu")] <- "mNeu"
only_Myeolid_seurat_seurat$v2_Cell_annotation[only_Myeolid_seurat_seurat$Cell_annotation %in% c("Camp+immNeu","Fcnb+immNeu","Fmnl2+immNeu","Hist1h3c+immNeu","Ifitm6+immNeu","Ltf+immNeu","Pclaf+immNeu","Top2a+immNeu","Ube2c+immNeu")] <- "immNeu"
only_Myeolid_seurat_seurat$v2_Cell_annotation[only_Myeolid_seurat_seurat$Cell_annotation %in% c("Apoe+Mono.P","S100a4+Mono.P")] <- "Mono.P"
only_Myeolid_seurat_seurat$v2_Cell_annotation[only_Myeolid_seurat_seurat$Cell_annotation %in% c("Elane+GMP","F13a1+GMP")] <- "GMP"
only_Myeolid_seurat_seurat$v2_Cell_annotation[only_Myeolid_seurat_seurat$Cell_annotation %in% c("Stfa3+PMNa","Ifitm1+PMNa","Isg15+PMNb")] <- "PMN"
only_Myeolid_seurat_seurat$v2_Cell_annotation <- factor(only_Myeolid_seurat_seurat$v2_Cell_annotation,levels=c("MPP","GMP","immNeu","mNeu","PMN","Mono.P"))

cell_sub <- c("MPP","Elane+GMP",
    "F13a1+GMP","Camp+immNeu","Fcnb+immNeu","Fmnl2+immNeu","Hist1h3c+immNeu","Ifitm6+immNeu","Ltf+immNeu","Pclaf+immNeu",
    "Top2a+immNeu","Ube2c+immNeu","Mmp8+mNeu","Retnlg+mNeu","Stfa3+PMNa","Ifitm1+PMNa","Isg15+PMNb","Apoe+Mono.P",
    "S100a4+Mono.P")
pal <- jdb_palette("corona")[c(2,1,3:length(cell_sub))]
names(pal) <- cell_sub
p1 <- DimPlot(object = only_Myeolid_seurat_seurat, reduction = "umap",group.by="Cell_annotation",label=TRUE,raster=FALSE,cols= pal,pt.size=0.1)+labs(title="Cell_annotation")
~~~

![image-20250226145317730](./Figure%20making.assets/image-20250226145317730.png)

~~~R
p1 <- FeaturePlot(object = only_Myeolid_seurat_seurat, features = c("Neu.dev","Mono.dev"),ncol=2,pt.size=.1,raster=FALSE,reduction="umap",label=FALSE,cols = CustomPalette(low ="#007BBF", mid = "#FFF485",high = "#FF0000")) +NoAxes()
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig1G.UMAP_pos.svg", plot=p1,width = 11, height = 5,dpi=300)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig1G.UMAP_pos.pdf", plot=p1,width = 11, height = 5,dpi=300)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig1G.UMAP_pos.png", plot=p1,width = 11, height = 5,dpi=300)
~~~

![image-20250226145352116](./Figure%20making.assets/image-20250226145352116.png)

~~~R
only_Myeolid_seurat_seurat$CDK4_6 <- colSums(GetAssayData(object = only_Myeolid_seurat_seurat, slot = "data")[c("Cdk4","Cdk6"),])
scRNA.hematopoiesis <- readRDS(file="./Fig_New_T00_obj.markers.rds")
filter_scRNA.hematopoiesis <- scRNA.hematopoiesis[scRNA.hematopoiesis$p_val_adj<=0.05,]
Stemness <- filter_scRNA.hematopoiesis[filter_scRNA.hematopoiesis$cluster=="HSPC",]
Differentiated <- filter_scRNA.hematopoiesis[filter_scRNA.hematopoiesis$cluster=="Neutro",]
Lineage_marker <- intersect(rownames(GetAssayData(object = only_Myeolid_seurat_seurat, slot = "data")),Stemness$gene)
speci_raw <- FetchData(object = only_Myeolid_seurat_seurat, vars = Lineage_marker,slot="data")
only_Myeolid_seurat_seurat[["Stemness"]] <- (rowSums(speci_raw))/length(Lineage_marker)
cc.genes <- readRDS("./mouse_cell_cycle_genes.rds")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
g2m.genes_list <- intersect(g2m.genes,rownames(only_Myeolid_seurat_seurat))
s.genes_list <- intersect(s.genes,rownames(only_Myeolid_seurat_seurat))
only_Myeolid_seurat_seurat <- CellCycleScoring(object = only_Myeolid_seurat_seurat,g2m.features =g2m.genes, s.features = s.genes,set.ident = TRUE,ctrl=2)
p1 <- FeaturePlot(object = only_Myeolid_seurat_seurat, features = c("G2M.Score","Stemness","CDK4_6"),raster=FALSE,
  ncol=3,pt.size=.1,reduction="umap",label=FALSE,cols = CustomPalette(low ="#007BBF", mid = "#FFF485",high = "#FF0000"),sort=TRUE)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig1H.UMAP_pos.svg", plot=p1,width = 16, height = 5,dpi=300)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig1H.UMAP_pos.pdf", plot=p1,width = 16, height = 5,dpi=300)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig1H.UMAP_pos.png", plot=p1,width = 16, height = 5,dpi=300)
~~~

![image-20250226145445272](./Figure%20making.assets/image-20250226145445272.png)

~~~R
Sel_sig <- c("G2M.Score","Stemness","CDK4_6")
All_sum <- as.data.frame(FetchData(object = only_Myeolid_seurat_seurat, vars = c(Sel_sig,"condition","v2_Cell_annotation"),slot="data"))
All_plots1 <- lapply(1:length(Sel_sig),function(x) {
  p1 <- ggboxplot(All_sum, x = "v2_Cell_annotation", y = Sel_sig[x], fill="v2_Cell_annotation",title=paste0(Sel_sig[x],""), legend = "none",outlier.shape = NA,notch = FALSE) +rotate_x_text(angle = 45)+
  stat_summary(fun.y = median, geom="point",colour="darkred", size=3) +scale_fill_manual(values = col)+scale_color_manual(values = col)+ 
  stat_summary(fun = median, geom = "line",aes(group = 1),col = "red",size=1)
  return(p1)
  })
p1 <- CombinePlots(c(All_plots1),nrow=1)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig1I.box.svg", plot=p1,width = 8, height = 5,dpi=300)
~~~

![image-20250226145527611](./Figure%20making.assets/image-20250226145527611.png)

~~~R
All.DEGs <- mcreadRDS("./pool2.v1_cell_anno.DEGs.rds",mc.cores=20)
unique(All.DEGs$CT)
table(All.DEGs$cluster)
CDK_KO <- subset(All.DEGs,cluster=="TP53_trila")
CDK_KO$group <- ifelse(CDK_KO$p_val < 0.05 & abs(CDK_KO$avg_log2FC) > 0.1, "sig","None")
CDK_KO$group[CDK_KO$p_val < 0.05 & CDK_KO$avg_log2FC > 0.1] <- "UP"
CDK_KO$group[CDK_KO$p_val < 0.05 & CDK_KO$avg_log2FC < -0.1] <- "DN"
CDK_KO <- subset(CDK_KO,group!="None")
Comb_KO <- subset(All.DEGs,cluster=="TP53_combo")
Comb_KO$group <- ifelse(Comb_KO$p_val < 0.05 & abs(Comb_KO$avg_log2FC) > 0.1, "sig","None")
Comb_KO$group[Comb_KO$p_val < 0.05 & Comb_KO$avg_log2FC > 0.1] <- "UP"
Comb_KO$group[Comb_KO$p_val < 0.05 & Comb_KO$avg_log2FC < -0.1] <- "DN"
Comb_KO <- subset(Comb_KO,group!="None")
CDK_KO.both.up <- CDK_KO[CDK_KO$gene %in% intersect(subset(CDK_KO,group=="UP")$gene,subset(Comb_KO,group=="UP")$gene),]
Comb_KO.both.up <- Comb_KO[Comb_KO$gene %in% intersect(subset(CDK_KO,group=="UP")$gene,subset(Comb_KO,group=="UP")$gene),]
CDK_KO.both.dn <- CDK_KO[CDK_KO$gene %in% intersect(subset(CDK_KO,group=="DN")$gene,subset(Comb_KO,group=="DN")$gene),]
Comb_KO.both.dn <- Comb_KO[Comb_KO$gene %in% intersect(subset(CDK_KO,group=="DN")$gene,subset(Comb_KO,group=="DN")$gene),]
CDK_KO <- rbind(CDK_KO.both.up,CDK_KO.both.dn)
CDK_KO.MPP <- subset(CDK_KO,CT=="MPP")
CDK_KO.MPP <- CDK_KO.MPP[order(CDK_KO.MPP$avg_log2FC),]
CDK_KO.MPP$group <- ifelse(CDK_KO.MPP$p_val < 0.05 & abs(CDK_KO.MPP$avg_log2FC) > 0.1, "sig","None")
CDK_KO.MPP$group[CDK_KO.MPP$p_val < 0.05 & CDK_KO.MPP$avg_log2FC > 0.1] <- "CDK_KO.UP"
CDK_KO.MPP$group[CDK_KO.MPP$p_val < 0.05 & CDK_KO.MPP$avg_log2FC < -0.1] <- "CDK_KO.DN"
CDK_KO.MPP <- subset(CDK_KO.MPP,group!="None")
CDK_KO.MPP$cluster <- paste0(CDK_KO.MPP$CT,".",CDK_KO.MPP$group)
info <- CDK_KO.MPP[,c("cluster","gene")]
top_marker <- c()
number.group <- length(unique(info$cluster))
for (i in c(1:number.group)){
  y <- info$cluster
  marker <- info[with(info,y==unique(info$cluster)[i]),]
  top_marker[[i]] <- marker
  names(top_marker)[i] <- unique(info$cluster)[i]
}
gcSampl <- c()
for (i in c(1:length(top_marker))){
t <- top_marker[[i]]
symbol <- as.character(t$gene)
DD <- symbol
t$entrez <- mapIds(x = org.Mm.eg.db, keys = DD,keytype ="SYMBOL",column ="ENTREZID",multiVals="first")
names <- na.omit(t)
entrez <- as.character(names$entrez)
gcSampl[[i]] <- entrez
names(gcSampl)[i] <- names(top_marker)[i]
print(paste(names(top_marker)[i],"is done",sep = " "))
}
HLH_T1_OFF_HIGH_GO <- mcreadRDS(file="./MPP.GO_v2.rds",mc.cores=20)
df <- as.data.frame(HLH_T1_OFF_HIGH_GO)
UP.df <- subset(df,Cluster=="MPP.CDK_KO.UP")
Myeolid <- UP.df[grep("mye",UP.df$Description,value=FALSE),]
UP.df <- rbind(subset(df,Cluster=="MPP.CDK_KO.UP")[1:5,],Myeolid)
DN.df <- subset(df,Cluster=="MPP.CDK_KO.DN")[1:10]
SE_all_clu1 <- rbind(UP.df,DN.df)
data_sel <- new("compareClusterResult", compareClusterResult = SE_all_clu1, geneClusters = gcSampl[unique(SE_all_clu1$Cluster)],fun = "enrichGO")
plot <- clusterProfiler::dotplot(data_sel,showCategory=10,includeAll=FALSE,label_format=100,color ="pvalue") + theme(axis.text.x  = element_text(angle=45, vjust=1,size=8,hjust = 1)) + labs(title = "MPP GO BP")
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig1J.GO.svg", plot=plot,width = 7, height = 5,dpi=300)
~~~

![image-20250226145613827](./Figure%20making.assets/image-20250226145613827.png)

~~~R
MPP <- c("Cd34","Hlf","Cd48","Gcnt2")
MEP <- c("Gata1","Itga2b","Pf4")
Ery <- c("Hbb-bt","Hba-a2","Hba-a1","Alas2")
CLP <- c("Flt3","Dntt")
GMP <- c("Mpo")
immNeu <- c("Ltf","Ngp","Camp")
mNeu <- c("Retnlg","Mmp8","Cxcl2")
PMNa <- c("Ccl6","Stfa2l1")
PMNb <- c("Isg15","Rsad2","Ifit3")
Mono.P <- c("Csf1r","Klf4")
DC <- c("Sirpa", "Itgam")
NK.T <- c("Cd3d","Cd3g","Cd3e","Lck","Nkg7", "Gzma","Cd160")
Early.B <- c("Vpreb2")
mature.B <- c("Cd79a","Cd74","Cd79b","Mzb1")
Plasma <- c("Slamf7","Sdc1","Prdm1","Jchain")
Sel.Mar <- c(MPP,MEP,Ery,CLP,GMP,"Elane",immNeu, mNeu,PMNa,PMNb,Mono.P,"Apoe","S100a4",DC,NK.T,Early.B,mature.B,Plasma)
plot <- DotPlot(pool2, features = Sel.Mar, cols=c("#ffffff", "#B30000"),scale = TRUE,group.by="v2_Cell_annotation",
    col.min = 0,col.max = 5) + RotatedAxis()+ labs(title="All cell annotation")
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S1A.dotplot.svg", plot=plot,width = 13, height = 4,dpi=300)
~~~

![image-20250226145645378](./Figure%20making.assets/image-20250226145645378.png)

~~~R
pool2$CDK4_6 <- colSums(GetAssayData(object = pool2, slot = "data")[c("Cdk4","Cdk6"),])
All_sum <- as.data.frame(FetchData(object = pool2, vars = c("CDK4_6","condition","Cell_annotation"),slot="data"))
only_MPP <- All_sum[All_sum$Cell_annotation %in% c("MPP"),]
p1 <- ggboxplot(only_MPP[only_MPP$condition %in% c("Vav_vehicle","TP53_veh"),], x = "condition", y = "CDK4_6", fill="condition",title=paste0("CDK4_6",".exp"), legend = "none",outlier.shape = NA,notch = FALSE) +rotate_x_text(angle = 45)+
stat_summary(fun.y = median, geom="point",colour="darkred", size=3) + stat_summary(fun = median, geom = "line",aes(group = 1),col = "red",size=1)+
stat_compare_means(comparisons =list(c("Vav_vehicle","TP53_veh")),method = "wilcox.test")
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S1B.CDK46.svg", plot=p1,width = 6, height = 5,dpi=300)
~~~

![image-20250226145800731](./Figure%20making.assets/image-20250226145800731.png)

~~~R
pool2$condition.v2 <- as.character(pool2$condition)
pool2$condition.v2[grep("veh",pool2$condition.v2,value=FALSE)] <- "veh"
pool2$condition.v2[grep("carb",pool2$condition.v2,value=FALSE)] <- "carbo"
pool2$condition.v2[grep("tril",pool2$condition.v2,value=FALSE)] <- "trila"
pool2$condition.v2[grep("combo",pool2$condition.v2,value=FALSE)] <- "comb"
pool2$condition.v2 <- factor(pool2$condition.v2,levels=c("veh","carbo","trila","comb"))
Idents(pool2) <- pool2$v2_Cell_annotation
library(tricycle)
pool2 <- Runtricycle(object = pool2, slot = "data", reduction.name = "tricycleEmbedding", 
    reduction.key = "tricycleEmbedding_", gname = NULL, gname.type = "SYMBOL", species = "mouse",center.pc1=0.5,center.pc2=0.5)
pool2$tricycleGroup <- "G01"
pool2$tricycleGroup[pool2$tricyclePosition > 0.5*pi & pool2$tricyclePosition <= 1*pi] <- "S"
pool2$tricycleGroup[pool2$tricyclePosition > 1*pi & pool2$tricyclePosition <= 1.5*pi] <- "G2M"
pool2$tricycleGroup[pool2$tricyclePosition > 1.5*pi & pool2$tricyclePosition <= 1.75*pi] <- "M"
aa <- as.data.frame(table(pool2$condition.v2,pool2$tricycleGroup))
aa <- aa[order(aa$Var2),]
aa$Var1 <- as.character(aa$Var1)
aa_all <- c()
for (i in unique(aa$Var1)){
  group_sel <- subset(aa,Var1==i)   
  group_sel$sum_number <- sum(group_sel$Freq)
  group_sel$normal_ratio <- (group_sel$Freq/group_sel$sum_number)*100
  group_sel$normal_ratio <- round(group_sel$normal_ratio,2)
  aa_all <- rbind(aa_all,group_sel)
}
aa_all$Var1 <- factor(aa_all$Var1,levels=c("veh","carbo","trila","comb"))
aa_all <- aa_all[order(aa_all$Var1),]
library(ggalluvial)
p1 <- ggplot(aa_all, aes(x = Var1, y = normal_ratio, fill = Var2, stratum = Var2, alluvium = Var2)) + geom_stratum(width = 0.75) + geom_flow(alpha = 0.5) + 
theme_classic() + theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1)) +labs(x = '', y = 'Relative Abundance(%)',title="CD45.2.scRNA")+scale_fill_manual(values = jdb_palette("corona"))+scale_color_manual(values = jdb_palette("corona"))
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S1C.cellcycle.svg", plot=p1,width = 4, height = 5,dpi=300)
~~~

![image-20250226145957161](./Figure%20making.assets/image-20250226145957161.png)

~~~R
col <- jdb_palette("corona")[c(2,1,3,4:7,9:length(jdb_palette("corona")))]
cell_sub <- c("MPP","MEP","Ery","CLP","GMP","immNeu","mNeu","PMN","Mono.P","DC","NK.T","Early.B","mature.B","Plasma")
col <- col[1:length(cell_sub)]
names(col) <- cell_sub
library(ggalluvial)
aa <- as.data.frame(table(pool2$condition,pool2$v2_Cell_annotation))
aa <- aa[order(aa$Var2),]
aa$Var1 <- as.character(aa$Var1)
aa_all <- c()
for (i in unique(aa$Var1)){
  group_sel <- subset(aa,Var1==i)   
  group_sel$sum_number <- sum(group_sel$Freq)
  group_sel$normal_ratio <- (group_sel$Freq/group_sel$sum_number)*100
  group_sel$normal_ratio <- round(group_sel$normal_ratio,2)
  aa_all <- rbind(aa_all,group_sel)
}
aa_all$Var1 <- factor(aa_all$Var1,levels=c("Vav_vehicle","Vav_carbo","Vav_trila","Vav_combo","TP53_veh","TP53_carbo","TP53_trila","TP53_combo"))
aa_all <- aa_all[order(aa_all$Var1),]
p1 <- ggplot(aa_all, aes(x = Var1, y = normal_ratio, fill = Var2, stratum = Var2, alluvium = Var2)) + geom_stratum(width = 0.75) + geom_flow(alpha = 0.5) + 
theme_classic() + theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1)) +labs(x = '', y = 'Relative Abundance(%)',title="pool2.scRNA")+scale_fill_manual(values = col)+scale_color_manual(values = col)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S1D.proporation_1.svg", plot=p1,width = 4, height = 5,dpi=300)
~~~

![image-20250226150031359](./Figure%20making.assets/image-20250226150031359.png)

~~~R
Idents(pool2) <- pool2$condition.v2
pool2.DS <- subset(pool2,downsample=8000)
cell_sub <- c("MPP","MEP","Ery","CLP","GMP","immNeu","mNeu","PMN","Mono.P","DC","NK.T","Early.B","mature.B","Plasma")
pal <- jdb_palette("corona")[c(2,1,3:7,9:14,16:length(jdb_palette("corona")))]
names(pal) <- cell_sub
library(tricycle)
pool2.DS <- Runtricycle(object = pool2.DS, slot = "data", reduction.name = "tricycleEmbedding", 
    reduction.key = "tricycleEmbedding_", gname = NULL, gname.type = "SYMBOL", species = "mouse")
pool2.DS$tricycleGroup <- "G01"
pool2.DS$tricycleGroup[pool2.DS$tricyclePosition > 0.5*pi & pool2.DS$tricyclePosition <= 1*pi] <- "S"
pool2.DS$tricycleGroup[pool2.DS$tricyclePosition > 1*pi & pool2.DS$tricyclePosition <= 1.5*pi] <- "G2M"
pool2.DS$tricycleGroup[pool2.DS$tricyclePosition > 1.5*pi & pool2.DS$tricyclePosition <= 1.75*pi] <- "M"
col <- jdb_palette("corona")[c(1,3,2,4:length(jdb_palette("corona")))]
col <- col[1:4]
names(col) <- c("G01","S","G2M","M")
col <- c(col,"#e1dfdd")
names(col)[length(col)] <- "OTS"
con.sel <- levels(pool2.DS$condition.v2)
All_plot <- lapply(1:length(con.sel),function(x) {
    pool2.DS$new_anno5 <- as.character(pool2.DS$tricycleGroup)
    pool2.DS$new_anno5[pool2.DS$condition.v2!=con.sel[x]] <- "OTS"
    pool2.DS$new_anno5 <- factor(pool2.DS$new_anno5,levels=c("OTS",c("G01","S","G2M","M")))
    plot <- DimPlot(pool2.DS, reduction = 'openTSNE', label = FALSE,repel=FALSE, pt.size =.1,group.by="new_anno5",cols=col[levels(pool2.DS$new_anno5)],order=TRUE) +labs(title=paste0("All.",con.sel[x]))
    return(plot)
    })
p1 <- CombinePlots(All_plot,ncol=4)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S1E.UMAP_pos.svg", plot=p1,width = 24, height = 5,dpi=300)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S1E.UMAP_pos.pdf", plot=p1,width = 24, height = 5,dpi=300)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S1E.UMAP_pos.png", plot=p1,width = 24, height = 5,dpi=300)
~~~

![image-20250226150123518](./Figure%20making.assets/image-20250226150123518.png)

~~~R
aa_all_merge <- mcreadRDS("./CD45.2.cells.propo.v2.tricycleGroup.v2.rds", mc.cores = 20)
cell_sub <-c("G01","S","G2M","M")
comb <- list(c("veh","carbo"),c("veh","trila"),c("veh","comb"),c("carbo","comb"))
All_plot_merge <- lapply(1:length(cell_sub),function(x) {
  plot <- ggviolin(aa_all_merge[aa_all_merge$Var2==cell_sub[x],], "group", "normal_ratio", fill = "group", add = c("boxplot","jitter"),outlier.shape = NA, 
    add.params = list(fill = "white"),legend = "none",width=1,title=paste0("%",cell_sub[x],".CD45.2"),) +rotate_x_text(angle = 45)+
  stat_summary(fun.y = median, geom="point",colour="darkred", size=3) + stat_summary(fun = median, geom = "line",aes(group = 1),col = "red",size=1)+
  stat_compare_means(comparisons =comb, method = "wilcox.test",paired=FALSE)
  return(plot)
})
p1 <- CombinePlots(All_plot_merge,nrow=1)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S1F.UMAP_pos.svg", plot=p1,width = 10, height = 4,dpi=300)
~~~

![image-20250226150211997](./Figure%20making.assets/image-20250226150211997.png)

~~~R
cell_sub <- c("MPP","Elane+GMP","F13a1+GMP","Camp+immNeu","Fcnb+immNeu","Fmnl2+immNeu","Hist1h3c+immNeu","Ifitm6+immNeu","Ltf+immNeu","Pclaf+immNeu","Top2a+immNeu","Ube2c+immNeu","Mmp8+mNeu","Retnlg+mNeu","Stfa3+PMNa","Ifitm1+PMNa","Isg15+PMNb","Apoe+Mono.P","S100a4+Mono.P")
pal <- jdb_palette("corona")[c(2,1,3:length(cell_sub))]
names(pal) <- cell_sub
Sel_sig <- c("CDK4_6","Stenmness")
All_sum <- as.data.frame(FetchData(object = only_Myeolid_seurat_seurat, vars = c("Neu.dev","Mono.dev",Sel_sig,"Cell_annotation"),slot="data"))
All_sum1 <- All_sum[!is.na(All_sum$Neu.dev),]
All_sum1 <- All_sum1[order(All_sum1$Neu.dev,decreasing=FALSE),]
All_sum1$order <- 1:nrow(All_sum1)
All_plots <- lapply(1:length(Sel_sig),function(x) {
    All_sum1.tmp <- All_sum1
    All_sum1.tmp[,Sel_sig[x]] <- (All_sum1.tmp[,Sel_sig[x]]-mean(All_sum1.tmp[,Sel_sig[x]]))/(max(All_sum1.tmp[,Sel_sig[x]])-min(All_sum1.tmp[,Sel_sig[x]]))
    plot <- ggplot(data = All_sum1.tmp,aes_string(x = "order", y = Sel_sig[x],color="Cell_annotation")) + geom_point(alpha = 0.1,size = 0.1)+ 
        xlab("Neu.dev (by Pesudo time)") + ylab(paste(Sel_sig[x]," (Expression)",sep=""))+
        geom_smooth(colour = "orange",se=TRUE)+ylim(c(min(All_sum1.tmp[,Sel_sig[x]]),max(All_sum1.tmp[,Sel_sig[x]])))+labs(title=paste0(Sel_sig[x]))+theme_classic()+NoLegend()+scale_fill_manual(values = pal)+scale_color_manual(values = pal)
    return(plot)
    })
p1 <- CombinePlots(c(All_plots),nrow=1)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S1G.exp.svg", plot=p1,width = 10, height = 5,dpi=300)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S1G.exp.pdf", plot=p1,width = 10, height = 5,dpi=300)
~~~

![image-20250226150250178](./Figure%20making.assets/image-20250226150250178.png)

~~~R
library(ggridges)
All_sum <- as.data.frame(FetchData(object = only_Myeolid_seurat_seurat, vars = c("Mono.dev","Mono.dev",Sel_sig,"Cell_annotation","condition"),slot="data"))
p1 <- ggplot(All_sum[All_sum$condition %in% c("Vav_vehicle","Vav_carbo","Vav_trila","Vav_combo"),], 
    aes(x = Mono.dev, y = condition, color = condition, fill = condition)) +
geom_density_ridges(alpha=0.8) +scale_y_discrete(expand = c(0, 0)) +coord_cartesian(clip = "off") +
ggtitle("Mono.dev.in.Vav") +theme_ridges(center = TRUE)+NoLegend()+scale_fill_manual(values = jdb_palette("corona"))+scale_color_manual(values = jdb_palette("corona"))
p2 <- ggplot(All_sum[All_sum$condition %in% c("TP53_veh","TP53_carbo","TP53_trila","TP53_combo"),], 
    aes(x = Mono.dev, y = condition, color = condition, fill = condition)) +
geom_density_ridges(alpha=0.8) +scale_y_discrete(expand = c(0, 0)) +coord_cartesian(clip = "off") +
ggtitle("Mono.dev.in.TP53") +theme_ridges(center = TRUE)+NoLegend()+scale_fill_manual(values = jdb_palette("corona"))+scale_color_manual(values = jdb_palette("corona"))
plot <- plot_grid(p1,p2)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S1H.Mono.svg", plot=plot,width = 8, height = 4,dpi=300)
~~~

![image-20250226150316310](./Figure%20making.assets/image-20250226150316310.png)

~~~R
col <- jdb_palette("corona")[c(2,1,3,4:7,9:length(jdb_palette("corona")))]
cell_sub <- c("MPP","MEP","Ery","CLP","GMP","immNeu","mNeu","PMN","Mono.P","DC","NK.T","Early.B","mature.B","Plasma")
col <- col[1:length(cell_sub)]
names(col) <- cell_sub
col <- c(col,"#e1dfdd")
names(col)[length(col)] <- "OTS"
pool2$new_anno5 <- as.character(pool2$v2_Cell_annotation)
pool2$new_anno5[!pool2$v2_Cell_annotation %in% c("MPP")] <- "OTS"
pool2$new_anno5 <- factor(pool2$new_anno5,levels=c("MPP","OTS"))
p1 <- DimPlot(pool2, reduction = 'openTSNE', label = TRUE,repel=TRUE, group.by="new_anno5",cols=col[levels(pool2$new_anno5)],pt.size=1,order=TRUE)+labs(title="Myeolid")
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S1I.openTSNE_pos.svg", plot=p1,width = 6, height = 5,dpi=300)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S1I.openTSNE_pos.pdf", plot=p1,width = 6, height = 5,dpi=300)
~~~

![image-20250226150400362](./Figure%20making.assets/image-20250226150400362-1740600240827-1.png)

~~~R
All.DEGs <- mcreadRDS("./pool2.v1_cell_anno.DEGs.rds",mc.cores=20)
unique(All.DEGs$CT)
table(All.DEGs$cluster)
CDK_KO <- subset(All.DEGs,cluster=="TP53_trila")
CDK_KO$group <- ifelse(CDK_KO$p_val < 0.05 & abs(CDK_KO$avg_log2FC) > 0.1, "sig","None")
CDK_KO$group[CDK_KO$p_val < 0.05 & CDK_KO$avg_log2FC > 0.1] <- "UP"
CDK_KO$group[CDK_KO$p_val < 0.05 & CDK_KO$avg_log2FC < -0.1] <- "DN"
CDK_KO <- subset(CDK_KO,group!="None")
Comb_KO <- subset(All.DEGs,cluster=="TP53_combo")
Comb_KO$group <- ifelse(Comb_KO$p_val < 0.05 & abs(Comb_KO$avg_log2FC) > 0.1, "sig","None")
Comb_KO$group[Comb_KO$p_val < 0.05 & Comb_KO$avg_log2FC > 0.1] <- "UP"
Comb_KO$group[Comb_KO$p_val < 0.05 & Comb_KO$avg_log2FC < -0.1] <- "DN"
Comb_KO <- subset(Comb_KO,group!="None")
CDK_KO.both.up <- CDK_KO[CDK_KO$gene %in% intersect(subset(CDK_KO,group=="UP")$gene,subset(Comb_KO,group=="UP")$gene),]
Comb_KO.both.up <- Comb_KO[Comb_KO$gene %in% intersect(subset(CDK_KO,group=="UP")$gene,subset(Comb_KO,group=="UP")$gene),]
CDK_KO.both.dn <- CDK_KO[CDK_KO$gene %in% intersect(subset(CDK_KO,group=="DN")$gene,subset(Comb_KO,group=="DN")$gene),]
Comb_KO.both.dn <- Comb_KO[Comb_KO$gene %in% intersect(subset(CDK_KO,group=="DN")$gene,subset(Comb_KO,group=="DN")$gene),]
CDK_KO <- rbind(CDK_KO.both.up,CDK_KO.both.dn)
CDK_KO.MPP <- subset(CDK_KO,CT=="MPP")

only_Myeolid_seurat_seurat <- mcreadRDS("./only_Myeolid.pool2.merge.Sv4.rds",mc.cores=20)
only_Myeolid_seurat_seurat$v2_Cell_annotation[only_Myeolid_seurat_seurat$Cell_annotation %in% c("Mmp8+mNeu","Retnlg+mNeu")] <- "mNeu"
only_Myeolid_seurat_seurat$v2_Cell_annotation[only_Myeolid_seurat_seurat$Cell_annotation %in% c("Camp+immNeu","Fcnb+immNeu","Fmnl2+immNeu","Hist1h3c+immNeu","Ifitm6+immNeu","Ltf+immNeu","Pclaf+immNeu","Top2a+immNeu","Ube2c+immNeu")] <- "immNeu"
only_Myeolid_seurat_seurat$v2_Cell_annotation[only_Myeolid_seurat_seurat$Cell_annotation %in% c("Apoe+Mono.P","S100a4+Mono.P")] <- "Mono.P"
only_Myeolid_seurat_seurat$v2_Cell_annotation[only_Myeolid_seurat_seurat$Cell_annotation %in% c("Elane+GMP","F13a1+GMP")] <- "GMP"
only_Myeolid_seurat_seurat$v2_Cell_annotation[only_Myeolid_seurat_seurat$Cell_annotation %in% c("Stfa3+PMNa","Ifitm1+PMNa","Isg15+PMNb")] <- "PMN"
only_Myeolid_seurat_seurat$v2_Cell_annotation <- factor(only_Myeolid_seurat_seurat$v2_Cell_annotation,levels=c("MPP","GMP","immNeu","mNeu","PMN","Mono.P"))

mouse_TF <- read.table("./Mus_musculus_TF.txt",sep="\t",header=TRUE)
ALL_GSEA_GMT <- read.gmt("./msigdb.v7.1.symbols.gmt")
ALL_GSEA_GMT1 = ALL_GSEA_GMT %>% mutate(mouse_gene = convert_human_to_mouse_symbols(gene))
ALL_GSEA_GMT2 <- ALL_GSEA_GMT1[,c("term","mouse_gene")]
colnames(ALL_GSEA_GMT2) <- c("ont","gene")
MYELOID <- intersect(mouse_TF$Symbol,ALL_GSEA_GMT2[ALL_GSEA_GMT2$ont %in% grep("MYELOID",unique(ALL_GSEA_GMT2$ont),value=TRUE),"gene"])
LYMPHOID <- intersect(mouse_TF$Symbol,ALL_GSEA_GMT2[ALL_GSEA_GMT2$ont %in% grep("LYMPHOID",unique(ALL_GSEA_GMT2$ont),value=TRUE),"gene"])

Idents(only_Myeolid_seurat_seurat) <- only_Myeolid_seurat_seurat$Cell_annotation
MPP <- subset(only_Myeolid_seurat_seurat,idents="MPP")
Idents(MPP) <- MPP$condition
All_gsva_seura_ <- future_lapply(1:length(levels(MPP$condition)),function(i) {
  sel_tmp <- subset(MPP,idents=levels(MPP$condition)[i])
  sel_tmp <- pseudo_bulk_seurat_mean_random(seurat_obj=sel_tmp,num_split=20,seed.use=1,prefix=levels(MPP$condition)[i],slot="data",assay="RNA")
  metadata <- data.frame(cell_type=c(rep(levels(MPP$condition)[i],20)),
    row.names=colnames(sel_tmp))
  sel_gsva_seurat <- CreateSeuratObject(counts = sel_tmp,assay = 'RNA',project = 'RNA',min.cells = 0,meta.data = metadata)
  message(levels(MPP$condition)[i], " is done")
  return(sel_gsva_seurat)
})
new_group_pseudo_bulk <- merge(x = All_gsva_seura_[[1]], y = All_gsva_seura_[c(2:length(All_gsva_seura_))])
Idents(new_group_pseudo_bulk) <- new_group_pseudo_bulk$cell_type
new_group_pseudo_bulk1 <- subset(new_group_pseudo_bulk,idents=c("TP53_veh","TP53_carbo","TP53_trila","TP53_combo"))
new_group_pseudo_bulk1$cell_type <- factor(new_group_pseudo_bulk1$cell_type,levels=c("TP53_veh","TP53_carbo","TP53_trila","TP53_combo"))
CDK_KO <- rbind(CDK_KO.both.up,CDK_KO.both.dn)
CDK_KO.MPP <- subset(CDK_KO,CT=="MPP")
CDK_KO.MPP <- CDK_KO.MPP[order(CDK_KO.MPP$avg_log2FC),]
mark_gene1 <- intersect(subset(CDK_KO.MPP,avg_log2FC < 0)$gene,MYELOID)
mark_gene2 <- intersect(subset(CDK_KO.MPP,avg_log2FC > 0)$gene,LYMPHOID)
pdf("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S1J.MPP.DEG.TF2.pdf",width =5, height =8)
XY_heatmap(seurat_obj=new_group_pseudo_bulk1,group="cell_type",genes=c("Gata1","Gata2",unique(CDK_KO.MPP$gene)),all_num=FALSE,new_names=NULL,labels_rot=90,
  assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),min_and_max_cut=2,show_row_names=FALSE,mark_gene=c("Gata1","Gata2",mark_gene1,mark_gene2),label_size=0,scale=TRUE)
dev.off()
~~~

![image-20250226150524785](./Figure%20making.assets/image-20250226150524785.png)

~~~R
only_Myeolid_seurat_seurat <- mcreadRDS("./only_Myeolid.pool2.merge.Sv4.rds",mc.cores=20)
library(CytoTRACE)
results <- mcreadRDS("./All.CytoTRACE.DS.pool2.rds",mc.cores=20)
only_Myeolid_seurat_seurat.DS <- subset(only_Myeolid_seurat_seurat,cells=intersect(colnames(only_Myeolid_seurat_seurat),names(results$CytoTRACE)))
only_Myeolid_seurat_seurat.DS$CytoTRACE <- results$CytoTRACE[rownames(only_Myeolid_seurat_seurat.DS[[]])]
only_Myeolid_seurat_seurat.DS$CytoTRACErank <- results$CytoTRACErank[rownames(only_Myeolid_seurat_seurat.DS[[]])]
ALL_GSEA_GMT <- read.gmt("./msigdb.v7.1.symbols.gmt")
ALL_GSEA_GMT1 = ALL_GSEA_GMT %>% mutate(mouse_gene = convert_human_to_mouse_symbols(gene))
ALL_GSEA_GMT2 <- ALL_GSEA_GMT1[,c("term","mouse_gene")]
colnames(ALL_GSEA_GMT2) <- c("ont","gene")
ALL_GSEA_GMT2$ont <- as.character(ALL_GSEA_GMT2$ont)
CELL_CYCLE <- ALL_GSEA_GMT2[ALL_GSEA_GMT2$ont=="KEGG_CELL_CYCLE",]
Lineage_marker <- intersect(rownames(GetAssayData(object = only_Myeolid_seurat_seurat.DS, slot = "data")),CELL_CYCLE$gene)
speci_raw <- FetchData(object = only_Myeolid_seurat_seurat.DS, vars = Lineage_marker,slot="data")
only_Myeolid_seurat_seurat.DS[["CELL_CYCLE"]] <- (rowSums(speci_raw))/length(Lineage_marker)
Sel_sig <- c("CytoTRACE","CELL_CYCLE","Cdk4","Cdk6")
All_sum <- as.data.frame(FetchData(object = only_Myeolid_seurat_seurat.DS, vars = c(Sel_sig,"condition","Cell_annotation","v2_Cell_annotation"),slot="data"))
Sel.sum <- All_sum[All_sum$Cell_annotation %in% c("MPP"),]
Sel.sum <- Sel.sum[Sel.sum$condition %in% c("TP53_veh","TP53_carbo","TP53_trila","TP53_combo"),]
Sel.sum$condition  <- factor(Sel.sum$condition,levels=c("TP53_veh","TP53_carbo","TP53_trila","TP53_combo"))
comb <- split(t(combn(levels(Sel.sum$condition), 2)), seq(nrow(t(combn(levels(Sel.sum$condition), 2)))))
p1 <- ggboxplot(Sel.sum, x = "condition", y = "CytoTRACE", fill="condition",title=paste0("CytoTRACE",".in.MPP"), legend = "none",outlier.shape = NA,notch = FALSE) +rotate_x_text(angle = 45)+
stat_summary(fun.y = median, geom="point",colour="darkred", size=3) + stat_summary(fun = median, geom = "line",aes(group = 1),col = "red",size=1)+ylim(0.85,1.2)+
stat_compare_means(comparisons =comb,method = "wilcox.test")
p2 <- ggboxplot(Sel.sum, x = "condition", y = "CELL_CYCLE", fill="condition",title=paste0("CELL_CYCLE",".in.MPP"), legend = "none",outlier.shape = NA,notch = FALSE) +rotate_x_text(angle = 45)+
stat_summary(fun.y = median, geom="point",colour="darkred", size=3) + stat_summary(fun = median, geom = "line",aes(group = 1),col = "red",size=1)+
stat_compare_means(comparisons =comb,method = "wilcox.test")
plot <- plot_grid(p1,p2,nrow=1)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig1E.svg", plot=plot,width = 6, height = 5,dpi=300)
~~~

![image-20250226150701250](./Figure%20making.assets/image-20250226150701250.png)

~~~R
pool1 <- mcreadRDS("./All_filter_miQC.pool1.merge.Sv4.rds",mc.cores=20)
pool1$v2_Cell_annotation <- as.character(pool1$Cell_annotation)
pool1$v2_Cell_annotation[pool1$Cell_annotation %in% c("Retnlg+mNeu","Il1b+mNeu")] <- "mNeu"
pool1$v2_Cell_annotation[pool1$Cell_annotation %in% c("Camp+immNeu","Fcnb+immNeu","Fmnl2+immNeu","Hist1h3c+immNeu","Ltf+immNeu","Ube2c+immNeu")] <- "immNeu"
pool1$v2_Cell_annotation[pool1$Cell_annotation %in% c("Apoe+Mono.P","F13a1+Mono.P","S100a4+Mono.P")] <- "Mono.P"
pool1$v2_Cell_annotation[pool1$Cell_annotation %in% c("Mpo+GMP","Elane+GMP")] <- "GMP"
pool1$v2_Cell_annotation[pool1$Cell_annotation %in% c("Stfa3+PMNa","Ifitm1+PMNa","Isg15+PMNb")] <- "PMN"
pool1$v2_Cell_annotation[pool1$Cell_annotation %in% c("Vpreb1+Early.B","Vpreb3+Early.B")] <- "Early.B"
pool1$v2_Cell_annotation <- factor(pool1$v2_Cell_annotation,levels=c("MPP","MEP","Ery","CLP","GMP","immNeu","mNeu","PMN","Mono.P","DC","NK.T","Early.B","mature.B","Plasma","ILC2p"))
col <- jdb_palette("corona")[c(2,1,3,4:7,9:length(jdb_palette("corona")))]
cell_sub <- c("MPP","MEP","Ery","CLP","GMP","immNeu","mNeu","PMN","Mono.P","DC","NK.T","Early.B","mature.B","Plasma","ILC2p")
col <- col[1:length(cell_sub)]
names(col) <- cell_sub
p1 <- DimPlot(object = pool1, reduction = "openTSNE",group.by="v2_Cell_annotation",label=TRUE,raster=FALSE,cols= col,pt.size=.1) +labs(title="Cell_annotation")
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig2A.openTSNE_pos.svg", plot=p1,width = 6, height = 5,dpi=300)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig2A.openTSNE_pos.pdf", plot=p1,width = 6, height = 5,dpi=300)
p1 <- DimPlot(object = pool1, reduction = "openTSNE",group.by="v2_Cell_annotation",label=FALSE,raster=FALSE,cols= col,pt.size=.1) +labs(title="Cell_annotation")
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig2A.openTSNE_pos.png", plot=p1,width = 6, height = 5,dpi=300)
~~~

![image-20250226151033496](./Figure%20making.assets/image-20250226151033496.png)

~~~R
CD45.1.only_T_cells <- mcreadRDS("./only_T_cells.pool1.rds", mc.cores = 20)
CD45.1.only_T_cells$v2_Cell_annotation <- as.character(CD45.1.only_T_cells$v2_Cell_annotation)
CD45.1.only_T_cells$v2_Cell_annotation[grep("Ifit3",CD45.1.only_T_cells$v2_Cell_annotation,value=FALSE)] <- "CD8Tm"
CD45.1.only_T_cells$v2_Cell_annotation[grep("Gzma",CD45.1.only_T_cells$v2_Cell_annotation,value=FALSE)] <- "CD8Tem"
CD45.1.only_T_cells$v2_Cell_annotation[grep("Tox",CD45.1.only_T_cells$v2_Cell_annotation,value=FALSE)] <- "CD8Tpex"
CD45.1.only_T_cells$v2_Cell_annotation[grep("Tnfrsf4",CD45.1.only_T_cells$v2_Cell_annotation,value=FALSE)] <- "Treg"
cell_sub <- c("Naive.T","Naive.CD8T","Naive.CD4T","Cxcr3+Th1","Cx3cr1+Th1","Treg","CD8Tm","CD8Tem","CD8Tpex")
CD45.1.only_T_cells$v2_Cell_annotation <- factor(CD45.1.only_T_cells$v2_Cell_annotation,levels=c(cell_sub))
pal <- jdb_palette("corona")[c(2,1,3:7,9:14,16:length(jdb_palette("corona")))]
names(pal) <- cell_sub
p1 <- DimPlot(object = CD45.1.only_T_cells, reduction = "openTSNE",label=TRUE,group.by="v2_Cell_annotation",cols=pal,pt.size=1.5) +labs(title="openTSNE")
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig2B.onlyT.openTSNE_pos.svg", plot=p1,width = 8, height = 5,dpi=300)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig2B.onlyT.openTSNE_pos.pdf", plot=p1,width = 8, height = 5,dpi=300)
~~~

![image-20250226151201501](./Figure%20making.assets/image-20250226151201501.png)

~~~R
aa_all_merge <- mcreadRDS("./only_T.cells.propo.v2.rds", mc.cores = 20)
comb <- list(c("TP53_veh","TP53_carbo"),c("TP53_veh","TP53_trila"),c("TP53_veh","TP53_combo"),c("TP53_carbo","TP53_combo"))
aa_all_merge <- aa_all_merge[!aa_all_merge$group %in% c("Vav_vehicle","Vav_carbo","Vav_trila","Vav_combo"),]
aa_all_merge$group <- factor(aa_all_merge$group,levels=c("TP53_veh","TP53_carbo","TP53_trila","TP53_combo"))
cell_sub <- c("CD8Tm","CD8Tem","CD8Tpex")
All_plot_merge <- lapply(1:length(cell_sub),function(x) {
  plot <- ggboxplot(aa_all_merge[aa_all_merge$Var2==cell_sub[x],], x = "group", y = "normal_ratio", fill="group",
    title=paste0("%",cell_sub[x],".CD45.1"), legend = "none",outlier.shape = NA,notch = FALSE,add = "jitter") +rotate_x_text(angle = 45)+
  stat_summary(fun.y = median, geom="point",colour="darkred", size=3) + stat_summary(fun = median, geom = "line",aes(group = 1),col = "red",size=1)+
  stat_compare_means(comparisons =comb, method = "wilcox.test",paired=FALSE)
  return(plot)
})
p1 <- CombinePlots(All_plot_merge,nrow=1)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig2C.onlyT.prop.svg", plot=p1,width = 8, height = 5,dpi=300)
~~~

![image-20250226151059586](./Figure%20making.assets/image-20250226151059586.png)

~~~R
openTSNE_pos <- as.data.frame(CD45.1.only_T_cells[["openTSNE"]]@cell.embeddings)
openTSNE_pos <- openTSNE_pos[rownames(CD45.1.only_T_cells[[]]),]
openTSNE_pos$v2_Cell_annotation <- CD45.1.only_T_cells$v2_Cell_annotation
openTSNE_pos$condition <- CD45.1.only_T_cells$condition
TP53_veh <- openTSNE_pos[openTSNE_pos$condition=="TP53_veh",]
TP53_carbo <- openTSNE_pos[openTSNE_pos$condition=="TP53_carbo",]
TP53_trila <- openTSNE_pos[openTSNE_pos$condition=="TP53_trila",]
TP53_combo <- openTSNE_pos[openTSNE_pos$condition=="TP53_combo",]
p1 <- ggplot(openTSNE_pos, aes(x =openTSNE_1, y = openTSNE_2,color=group.v3))+geom_point(alpha = 0.1, size = 2, aes(color = v2_Cell_annotation))+
geom_density_2d(data=TP53_veh,aes(alpha = ..nlevel..), size=1,bins = 10,color=c(jdb_palette("corona")[2]))+scale_alpha_continuous(range = c(-0.2, 1))+
guides(color = guide_legend(override.aes = list(size = 2,alpha=1)))+scale_fill_manual(values = pal)+scale_color_manual(values = pal)+
labs(title=paste0("TP53_veh"))+theme_void()
p2 <- ggplot(openTSNE_pos, aes(x =openTSNE_1, y = openTSNE_2,color=group.v3))+geom_point(alpha = 0.1, size = 2, aes(color = v2_Cell_annotation))+
geom_density_2d(data=TP53_carbo,aes(alpha = ..nlevel..), size=1,bins = 10,color=c(jdb_palette("corona")[2]))+scale_alpha_continuous(range = c(-0.2, 1))+
guides(color = guide_legend(override.aes = list(size = 2,alpha=1)))+scale_fill_manual(values = pal)+scale_color_manual(values = pal)+
labs(title=paste0("TP53_carbo"))+theme_void()
p3 <- ggplot(openTSNE_pos, aes(x =openTSNE_1, y = openTSNE_2,color=group.v3))+geom_point(alpha = 0.1, size = 2, aes(color = v2_Cell_annotation))+
geom_density_2d(data=TP53_trila,aes(alpha = ..nlevel..), size=1,bins = 10,color=c(jdb_palette("corona")[2]))+scale_alpha_continuous(range = c(0, 2))+
guides(color = guide_legend(override.aes = list(size = 2,alpha=1)))+scale_fill_manual(values = pal)+scale_color_manual(values = pal)+
labs(title=paste0("TP53_trila"))+theme_void()
p4 <- ggplot(openTSNE_pos, aes(x =openTSNE_1, y = openTSNE_2,color=group.v3))+geom_point(alpha = 0.1, size = 2, aes(color = v2_Cell_annotation))+
geom_density_2d(data=TP53_combo,aes(alpha = ..nlevel..), size=1,bins = 10,color=c(jdb_palette("corona")[2]))+scale_alpha_continuous(range = c(-0.2, 1))+
guides(color = guide_legend(override.aes = list(size = 2,alpha=1)))+scale_fill_manual(values = pal)+scale_color_manual(values = pal)+
labs(title=paste0("TP53_combo"))+theme_void()
plot <- plot_grid(p1,p2,p3,p4,nrow=1)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig2D.onlyT.dis.svg", plot=plot,width = 24, height = 5,dpi=300)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig2D.onlyT.dis.pdf", plot=plot,width = 24, height = 5,dpi=300)
~~~

![image-20250226151235849](./Figure%20making.assets/image-20250226151235849.png)

~~~R
Idents(CD45.1.only_T_cells) <- CD45.1.only_T_cells$condition
only_T_lin_seurat1 <- subset(CD45.1.only_T_cells,idents=c("TP53_veh","TP53_carbo","TP53_trila","TP53_combo"))
Idents(only_T_lin_seurat1) <- only_T_lin_seurat1$v2_Cell_annotation
only_T_lin_seurat1 <- subset(only_T_lin_seurat1,idents=c("Naive.T","Naive.CD8T","CD8Tm","CD8Tem","CD8Tpex"))
Idents(only_T_lin_seurat1) <- only_T_lin_seurat1$condition
group <- c("TP53_veh","TP53_carbo","TP53_trila","TP53_combo")
All_seurat <- lapply(1:length(group),function(x) {
    sel_tmp_sub <- subset(only_T_lin_seurat1,idents=group[x])
    if (ncol(sel_tmp_sub) < 10) {if (ncol(sel_tmp_sub) < 10) {nbin=ncol(sel_tmp_sub) } else {nbin=10}} else {nbin=10}
    sel_tmp_sub_DS <- pseudo_bulk_seurat_mean_random(seurat_obj=sel_tmp_sub,num_split=nbin,seed.use=1,slot="data",prefix=paste0(group[x]),assay="RNA")
    metadata <- data.frame(cell_type=rep(group[x],nbin),group=c(rep(group[x],nbin)),row.names=colnames(sel_tmp_sub_DS))
    sel_tmp_seurat <- CreateSeuratObject(counts = sel_tmp_sub_DS,assay = 'RNA',project = 'RNA',min.cells = 0,meta.data = metadata)
    message(group[x], " is done")
    return(sel_tmp_seurat)
    })
All_gsva_seura <- merge(x = All_seurat[[1]], y = All_seurat[c(2:length(All_seurat))])
Sel_G <- c("Kif2","Ccr7","Sell","Itgb7","Il7r","Slamf6","Tcf7","Lef1","S1pr1","Tox","Havcr2","Nr4a2","Lag3","Pdcd1","Id2","Tigit","Cxcr6")
Sel_G <- intersect(Sel_G,rownames(All_gsva_seura))
setdiff(Sel_G,rownames(All_gsva_seura))
All_gsva_seura$cell_type <- factor(All_gsva_seura$cell_type,levels=group)
pdf("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/Fig2E.heatmap.pdf",width =5, height =4)
XY_heatmap(seurat_obj=All_gsva_seura,group="cell_type",genes=Sel_G,all_num=FALSE,new_names=NULL,labels_rot=90,
  assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),min_and_max_cut=2,show_row_names=TRUE,mark_gene=NULL,label_size=0,scale=TRUE)
dev.off()
~~~

![image-20250226151316305](./Figure%20making.assets/image-20250226151316305.png)

~~~r
MPP <- c("Cd34","Hlf","Cd48","Gcnt2")
MEP <- c("Gata1","Itga2b","Pf4")
Ery <- c("Hbb-bt","Hba-a2","Hba-a1","Alas2")
CLP <- c("Flt3")
GMP <- c("Mpo")
immNeu <- c("Ngp","Camp")
mNeu <- c("Retnlg","Mmp8","Cxcl2")
PMNa <- c("Ccl6","Stfa2l1")
PMNb <- c("Isg15","Rsad2","Ifit3")
Mono.P <- c("Csf1r","Klf4")
NK.T <- c("Cd3d","Cd3g","Cd3e","Lck", "Gzma","Cd160")
Early.B <- c("Vpreb2")
mature.B <- c("Cd79a","Cd74","Cd79b","Mzb1","Ebf1")
Plasma <- c("Slamf7","Sdc1","Prdm1","Jchain")
Sel.Mar <- c(MPP,MEP,Ery,CLP,GMP,"Elane",mNeu,PMNa,PMNb,
  Mono.P,"Apoe","F13a1","S100a4",NK.T,Early.B,mature.B,Plasma,"Rora","Gata3","Il2r","Tcf7")
p1 <- DotPlot(pool1, features = Sel.Mar, cols=c("#ffffff", "#B30000"),scale = TRUE,group.by="v2_Cell_annotation",
  col.min = 0,col.max = 5) + RotatedAxis()+ labs(title="CD45.1 cell annotation")
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S2A.sotplot.svg", plot=p1,width = 12, height = 4,dpi=300)
~~~

![image-20250226151342430](./Figure%20making.assets/image-20250226151342430.png)

~~~R
aa_all_merge <- mcreadRDS("./CD45.1.and.CD45.2.cells.propo.v2.rds", mc.cores = 20)
aa_all_merge1 <- aa_all_merge[aa_all_merge$Var2 %in% c("Myeloid","Lymphoid"),]
aa_all_merge1$group <- factor(aa_all_merge1$group,levels=c("CD45.1","CD45.2"))
aa_all_merge1$Var2 <- factor(aa_all_merge1$Var2,levels=c("Myeloid","Lymphoid"))
CD45.2.sum <- aa_all_merge1[grep("CD45.2",aa_all_merge1$group,value=FALSE),]
sample <- unique(aa_all_merge1$Var1)
aa_all1_ <- lapply(1:length(sample),function(x) {
    tmp <- aa_all_merge1[aa_all_merge1$Var1==sample[x],]
    tmp$relative <- 0
    tmp[tmp$Var2=="Lymphoid",]$relative <- tmp[tmp$Var2=="Lymphoid",]$Freq/median(CD45.2.sum[CD45.2.sum$Var2=="Lymphoid",]$Freq)
    tmp[tmp$Var2=="Myeloid",]$relative <- tmp[tmp$Var2=="Myeloid",]$Freq/median(CD45.2.sum[CD45.2.sum$Var2=="Myeloid",]$Freq)
    return(tmp)
    })
aa_all1 <- do.call(rbind,aa_all1_)
aa_all1$relative <- log(aa_all1$relative+1,10)
aa_all1$mouse_tmp <- gsub("CD45.1.","",aa_all1$Var1)
aa_all1$mouse_tmp <- gsub("CD45.2.","",aa_all1$mouse_tmp)
aa_all1$mouse_id <- paste0(aa_all1$mouse_tmp,"_",aa_all1$mouse_id)
aa_all1 <- aa_all1[aa_all1$condition %in% c("Vav_trila","Vav_combo","TP53_trila","TP53_combo"),]
Lymphoid.paired.df <- reshape2::dcast(aa_all1[aa_all1$Var2=="Lymphoid",],mouse_id~group,value.var = "relative")
Myeloid.paired.df <- reshape2::dcast(aa_all1[aa_all1$Var2=="Myeloid",],mouse_id~group,value.var = "relative")
Lymphoid.paired.df$group <- "Lymphoid"
Myeloid.paired.df$group <- "Myeloid"
paired.df <- rbind(Myeloid.paired.df,Lymphoid.paired.df)
paired.df[is.na(paired.df)] <- 0
paired.df$Diff <- paired.df$CD45.1 - paired.df$CD45.2
p1 <- ggpaired(paired.df, cond1  = "CD45.2", cond2  = "CD45.1",fill = "condition", palette = "jco",legend = "none",
    outlier.shape = NA,notch = TRUE, line.color = "gray", line.size = 0.1 ,title=paste0("CDKi & Combined"))+ stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),method = "t.test",paired=FALSE)+
    geom_line(data = ~subset(., Diff > 0), aes(group = id), color = "#ff4757")+ geom_line(data = ~subset(., Diff < 0), aes(group = id), color = "#3498db")+facet_wrap(~group)
aa_all1 <- do.call(rbind,aa_all1_)
aa_all1$relative <- log(aa_all1$relative+1,10)
aa_all1$mouse_tmp <- gsub("CD45.1.","",aa_all1$Var1)
aa_all1$mouse_tmp <- gsub("CD45.2.","",aa_all1$mouse_tmp)
aa_all1$mouse_id <- paste0(aa_all1$mouse_tmp,"_",aa_all1$mouse_id)
aa_all1 <- aa_all1[aa_all1$condition %in% c("Vav_carbo","Vav_combo","TP53_carbo","TP53_combo"),]
Lymphoid.paired.df <- reshape2::dcast(aa_all1[aa_all1$Var2=="Lymphoid",],mouse_id~group,value.var = "relative")
Myeloid.paired.df <- reshape2::dcast(aa_all1[aa_all1$Var2=="Myeloid",],mouse_id~group,value.var = "relative")
Lymphoid.paired.df$group <- "Lymphoid"
Myeloid.paired.df$group <- "Myeloid"
paired.df <- rbind(Myeloid.paired.df,Lymphoid.paired.df)
paired.df[is.na(paired.df)] <- 0
paired.df$Diff <- paired.df$CD45.1 - paired.df$CD45.2
p2 <- ggpaired(paired.df, cond1  = "CD45.2", cond2  = "CD45.1",fill = "condition", palette = "jco",legend = "none",
    outlier.shape = NA,notch = TRUE, line.color = "gray", line.size = 0.1 ,title=paste0("carbo & Combined"))+ stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),method = "t.test",paired=FALSE)+
    geom_line(data = ~subset(., Diff > 0), aes(group = id), color = "#ff4757")+ geom_line(data = ~subset(., Diff < 0), aes(group = id), color = "#3498db")+facet_wrap(~group)
plot <- plot_grid(p1,p2)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S2B.paired.plots.svg", plot=plot,width = 8, height = 4,dpi=300)
~~~

![image-20250226151413692](./Figure%20making.assets/image-20250226151413692.png)

~~~R
Naive <- c("Ccr7","Tcf7","Lef1")
Th1 <- c("Tbx21","Cd40lg")
Tregs <- c("Foxp3","Il2ra","Ctla4","Ikzf2")
Memory <- c("Eomes")
Effector_Memory <- c("Klrd1")
Exhausted_activated <- c("Lag3","Havcr2","Pdcd1","Entpd1")
Sel.Mar <- c(Naive,Th1,"Cxcr6","Zeb2","Cxcr3","Cx3cr1",Tregs,"Tnfrsf4",Memory,"Ifit3",Effector_Memory,"Gzma",Exhausted_activated,"Tox")
Idents(CD45.1.only_T_cells) <- CD45.1.only_T_cells$v2_Cell_annotation
group <- levels(CD45.1.only_T_cells$v2_Cell_annotation)
All_seurat <- lapply(1:length(group),function(x) {
    sel_tmp_sub <- subset(CD45.1.only_T_cells,idents=group[x])
    if (ncol(sel_tmp_sub) < 20) {if (ncol(sel_tmp_sub) < 20) {nbin=ncol(sel_tmp_sub) } else {nbin=20}} else {nbin=20}
    sel_tmp_sub_DS <- pseudo_bulk_seurat_mean_random(seurat_obj=sel_tmp_sub,num_split=nbin,seed.use=1,slot="data",prefix=paste0(group[x]),assay="RNA")
    metadata <- data.frame(cell_type=rep(group[x],nbin),group=c(rep(group[x],nbin)),row.names=colnames(sel_tmp_sub_DS))
    sel_tmp_seurat <- CreateSeuratObject(counts = sel_tmp_sub_DS,assay = 'RNA',project = 'RNA',min.cells = 0,meta.data = metadata)
    message(group[x], " is done")
    return(sel_tmp_seurat)
    })
only_Zeb2.All_gsva_seura <- merge(x = All_seurat[[1]], y = All_seurat[c(2:length(All_seurat))])
cell_sub <- c("Naive.T","Naive.CD8T","Naive.CD4T","Cxcr3+Th1","Cx3cr1+Th1","Treg","CD8Tm","CD8Tem","CD8Tpex")
only_Zeb2.All_gsva_seura$cell_type <- factor(only_Zeb2.All_gsva_seura$cell_type,levels=cell_sub)
pdf("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S2C.heatmap.pdf",width =5, height =6)
XY_heatmap(seurat_obj=only_Zeb2.All_gsva_seura,group="cell_type",genes=Sel.Mar,all_num=FALSE,new_names=NULL,labels_rot=90,
  assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),min_and_max_cut=2,show_row_names=TRUE,mark_gene=NULL,label_size=0,scale=TRUE)
dev.off()
~~~

![image-20250226151539105](./Figure%20making.assets/image-20250226151539105.png)

~~~R
only_Myeolid_seurat_seurat <- mcreadRDS("./only_Myeolid.pool2.merge.Sv4.rds",mc.cores=20)
col <- jdb_palette("corona")[c(2,1,3,4:7,9:length(jdb_palette("corona")))]
cell_sub <- c("MPP","MEP","Ery","CLP","GMP","immNeu","mNeu","PMN","Mono.P","DC","NK.T","Early.B","mature.B","Plasma","ILC2p")
col <- col[1:length(cell_sub)]
names(col) <- cell_sub
Sel_sig <- c("Cdk4","Cdk6","Stenmness")
All_sum <- as.data.frame(FetchData(object = only_Myeolid_seurat_seurat, vars = c(Sel_sig,"Mono.dev","condition","v2_Cell_annotation"),slot="data"))
All_sum$CDK4_6 <- All_sum$Cdk4+All_sum$Cdk6
All_sum1 <- All_sum[!is.na(All_sum$Mono.dev),]
All_sum1 <- All_sum1[order(All_sum1$Mono.dev,decreasing=FALSE),]
All_sum1$order <- 1:nrow(All_sum1)
Sel_sig <- c("CDK4_6","Stenmness")
All_plots <- lapply(1:length(Sel_sig),function(x) {
    All_sum1.tmp <- All_sum1
    All_sum1.tmp[,Sel_sig[x]] <- (All_sum1.tmp[,Sel_sig[x]]-mean(All_sum1.tmp[,Sel_sig[x]]))/(max(All_sum1.tmp[,Sel_sig[x]])-min(All_sum1.tmp[,Sel_sig[x]]))
    plot <- ggplot(data = All_sum1.tmp,aes_string(x = "order", y = Sel_sig[x],color="v2_Cell_annotation")) + geom_point(alpha = 0.1,size = 0.1)+ 
        xlab("Mono.dev (by Pesudo time)") + ylab(paste(Sel_sig[x]," (Expression)",sep=""))+
        geom_smooth(colour = "orange",se=TRUE)+ylim(c(min(All_sum1.tmp[,Sel_sig[x]]),max(All_sum1.tmp[,Sel_sig[x]])))+labs(title=paste0(Sel_sig[x]))+theme_classic()+NoLegend()+scale_fill_manual(values = col)+scale_color_manual(values = col)
    return(plot)
    })
plot <- CombinePlots(c(All_plots),nrow=1)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/S2C.CDK4_6.and.Stenmness.plots.svg", plot=plot,width = 8, height = 4,dpi=300)
~~~

![image-20250226152222421](./Figure%20making.assets/image-20250226152222421.png)

~~~R
pool2 <- mcreadRDS("./All_filter_miQC.pool2.merge.Sv4.rds",mc.cores=20)
Idents(pool2) <- pool2$v2_Cell_annotation
only_MPP <- subset(pool2,idents=c("MPP"))
library(tricycle)
cc.genes <- readRDS("./mouse_cell_cycle_genes.rds")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
g2m.genes_list <- intersect(g2m.genes,rownames(only_MPP))
s.genes_list <- intersect(s.genes,rownames(only_MPP))
only_MPP <- CellCycleScoring(object = only_MPP,g2m.features =g2m.genes, s.features = s.genes,set.ident = TRUE,ctrl=2)
only_MPP <- Runtricycle(object = only_MPP, slot = "data", reduction.name = "tricycleEmbedding", 
    reduction.key = "tricycleEmbedding_", gname = NULL, gname.type = "SYMBOL", species = "mouse",center.pc1=0.5,center.pc2=0.5)
only_MPP$tricycleGroup <- "G01"
only_MPP$tricycleGroup[only_MPP$tricyclePosition > 0.5*pi & only_MPP$tricyclePosition <= 1*pi] <- "S"
only_MPP$tricycleGroup[only_MPP$tricyclePosition > 1*pi & only_MPP$tricyclePosition <= 1.5*pi] <- "G2M"
only_MPP$tricycleGroup[only_MPP$tricyclePosition > 1.5*pi & only_MPP$tricyclePosition <= 1.75*pi] <- "M"
tricycleEmbedding <- as.data.frame(only_MPP[["tricycleEmbedding"]]@cell.embeddings)
tricycleEmbedding <- tricycleEmbedding[rownames(only_MPP[[]]),]
tricycleEmbedding$condition <- as.character(only_MPP$condition)
tricycleEmbedding$tricycleGroup <- as.character(only_MPP$tricycleGroup)
tricycleEmbedding <- tricycleEmbedding[tricycleEmbedding$condition %in% c("TP53_veh","TP53_carbo","TP53_trila","TP53_combo"),]
tricycleEmbedding$condition <- factor(tricycleEmbedding$condition,levels=c("TP53_veh","TP53_carbo","TP53_trila","TP53_combo"))
tricycleEmbedding$tricycleEmbedding_1 <- (-1)*tricycleEmbedding$tricycleEmbedding_1
col <- jdb_palette("corona")[c(1,3,2,4:length(jdb_palette("corona")))]
col <- col[1:4]
names(col) <- c("G01","S","G2M","M")
col <- c(col,"#e1dfdd")
names(col)[length(col)] <- "OTS"


All_plots1 <- lapply(1:length(levels(tricycleEmbedding$condition)),function(x) {
    tmp.tricycleEmbedding <- tricycleEmbedding[tricycleEmbedding$condition==levels(tricycleEmbedding$condition)[x],]
    if (levels(tricycleEmbedding$condition)[x] %in% c("Vav_combo","TP53_combo")){min=0} else {min=-0}
    if (levels(tricycleEmbedding$condition)[x] %in% c("Vav_carbo")){min=0}
    sp <- ggplot(tmp.tricycleEmbedding, aes(x=tricycleEmbedding_1, y=tricycleEmbedding_2,color=tricycleGroup)) +
    geom_point(alpha=0.5,size=2)+ stat_ellipse(data=tricycleEmbedding,aes(fill = tricycleGroup), geom="polygon",level=0.9,alpha=0.1,lwd =0.5)+
    geom_density_2d(data=tmp.tricycleEmbedding,aes(alpha = ..nlevel..), size=1,bins = 8,color=c(jdb_palette("corona")[2]))+scale_alpha_continuous(range = c(min, 1))+
    theme_classic()+scale_fill_manual(values = col)+scale_color_manual(values = col) +  guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(title=paste0("HSC/MPP.",levels(tricycleEmbedding$condition)[x]))+xlim(range(tricycleEmbedding$tricycleEmbedding_1))+ylim(range(tricycleEmbedding$tricycleEmbedding_2))
    return(sp)
    })
All_plots2 <- lapply(1:length(levels(tricycleEmbedding$condition)),function(x) {
    tmp.tricycleEmbedding <- tricycleEmbedding[tricycleEmbedding$condition==levels(tricycleEmbedding$condition)[x],]
    if (levels(tricycleEmbedding$condition)[x] %in% c("Vav_combo","TP53_combo")){min=0} else {min=-0}
    if (levels(tricycleEmbedding$condition)[x] %in% c("Vav_carbo")){min=0}
    sp <- ggplot(tmp.tricycleEmbedding, aes(x=tricycleEmbedding_1, y=tricycleEmbedding_2,color=tricycleGroup)) +
    geom_point(alpha=0.5,size=2)+ stat_ellipse(data=tricycleEmbedding,aes(fill = tricycleGroup), geom="polygon",level=0.9,alpha=0.1,lwd =0.5)+
    theme_classic()+scale_fill_manual(values = col)+scale_color_manual(values = col) +  guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(title=paste0("HSC/MPP.",levels(tricycleEmbedding$condition)[x]))+xlim(range(tricycleEmbedding$tricycleEmbedding_1))+ylim(range(tricycleEmbedding$tricycleEmbedding_2))
    return(sp)
    })
All_plots3 <- lapply(1:length(levels(tricycleEmbedding$condition)),function(x) {
    tmp.tricycleEmbedding <- tricycleEmbedding[tricycleEmbedding$condition==levels(tricycleEmbedding$condition)[x],]
    if (levels(tricycleEmbedding$condition)[x] %in% c("Vav_combo","TP53_combo")){min=0} else {min=-0}
    if (levels(tricycleEmbedding$condition)[x] %in% c("Vav_carbo")){min=0}
    sp <- ggplot(tmp.tricycleEmbedding, aes(x=tricycleEmbedding_1, y=tricycleEmbedding_2,color=tricycleGroup)) +
    stat_ellipse(data=tricycleEmbedding,aes(fill = tricycleGroup), geom="polygon",level=0.9,alpha=0.1,lwd =0.5)+
    geom_density_2d(data=tmp.tricycleEmbedding,aes(alpha = ..nlevel..), size=1,bins = 8,color=c(jdb_palette("corona")[2]))+scale_alpha_continuous(range = c(min, 1))+
    theme_classic()+scale_fill_manual(values = col)+scale_color_manual(values = col) +  guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(title=paste0("HSC/MPP.",levels(tricycleEmbedding$condition)[x]))+xlim(range(tricycleEmbedding$tricycleEmbedding_1))+ylim(range(tricycleEmbedding$tricycleEmbedding_2))
    return(sp)
    })
plot <- CombinePlots(c(All_plots1,All_plots2,All_plots3),ncol=4)
ggsave("/mnt/d/xiangyu.ubuntu/projects/OTS/Figure_making.X.Pan.v1/cellCycle.all.svg", plot=plot,width = 20, height = 12,dpi=300)
~~~

![image-20250226152421333](./Figure%20making.assets/image-20250226152421333.png)





