library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggrepel) 
library(magrittr)
library(scRNAtoolVis)
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
set.seed(123)

marker <- readRDS("marker.RDS")
scRNA <- readRDS("03_FindClusters_3.RDS")

library(scRNAtoolVis)
colnames(scRNA@meta.data)
table(scRNA@meta.data$region.main)
Idents(scRNA) <- 'region.main'
list <- rownames(table(Idents(scRNA)))
list

marker.dataframe <- data.frame()
logfc.threshold=0.2 
for (i in list){ 
  marker <- FindMarkers(scRNA, ident.1 = "sti", group.by = 'type',subset.ident = i,raster=FALSE,logfc.threshold = logfc.threshold,)
  output_dir <- paste0('.')
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("Dir already exists!")
  }
  write.csv(marker,file = paste0(output_dir,'/DEG.bytype.',i,'.csv')) 
  marker.filter <- marker[marker$p_val_adj<0.05 & abs(marker$avg_log2FC)>logfc.threshold,] 
  marker.filter$cluster <- i
  marker.filter$gene <- rownames(marker.filter) #这一步为了防止因基因名重复自动带上尾号数字
  marker.dataframe <- rbind(marker.dataframe,marker.filter)

  write.csv(marker.filter,file = paste0(output_dir,'/DEG.filtered.bytype.',i,'.csv'))
  assign(paste0('marker.filter.',i),marker.filter) #分配循环中的变量到环境中

sce2.markers <- marker.dataframe
sce2.markers <- sce2.markers[-grep(pattern = "^Hb[^(p)]",x = rownames(sce2.markers)),]  #由于结果中较多红细胞基因，在此去除红细胞相关基因
}

sce2.markers <- marker.dataframe
table(sce2.markers$cluster)
sce2.markers$cluster[sce2.markers$cluster=='NAc'] <- '1_NAc'
sce2.markers$cluster[sce2.markers$cluster=='Cortex'] <- '2_Cortex'
sce2.markers$cluster[sce2.markers$cluster=='LSI'] <- '3_LSI'
sce2.markers$cluster[sce2.markers$cluster=='CC'] <- '4_CC'
sce2.markers$cluster[sce2.markers$cluster=='CPu'] <- '5_CPu'
sce2.markers$cluster[sce2.markers$cluster=='VDB'] <- '6_VDB'
sce2.markers <- sce2.markers[sce2.markers$cluster!='Meninges',]
table(sce2.markers$cluster)

global <- read.csv(file = 'DEG.global.filtered.csv',header = T)
global_gene <- global[abs(global$avg_log2FC) >= 0.5,'X']
global_gene
# -----filter genes
global_gene <- c(global_gene,'Gm3417','Gm45208','Gm21887','Gh','Gm15601','Rpl13a-ps1')
global_gene
dim(sce2.markers)
sce2.markers <- sce2.markers[!sce2.markers$gene %in% global_gene, ]
dim(sce2.markers)

sce2.markers <- sce2.markers[-grep(pattern = "^Hb[^(p)]",x = sce2.markers$gene),]  #由于结果中较多红细胞基因，在此去除红细胞相关基因

logfc.threshold <- 0.2
markerVocalno(markers = sce2.markers,
              pforce=10,
              nforce=10,
              log2FC=logfc.threshold,
              topn = 5,
              labelCol = ggsci::pal_npg()(9))

saveRDS(object = sce2.markers,file ='./sce2.markers.RDS' )
write.csv(x = sce2.markers,file = './sce2.markers.csv')
ggsave(filename = '3_Volcano.pdf',width =8.5 ,height =5 )


#-----GO analysis-----
library(ggpubr)
library(clusterProfiler) 

marker <- sce2.markers 

library("org.Mm.eg.db")
Entrez = mapIds(x = org.Mm.eg.db,
                keys = marker$gene, 
                keytype = "SYMBOL", 
                column = "ENTREZID") 
head(Entrez )
marker$Entrez <- Entrez

marker <- marker[abs(marker$avg_log2FC) >= logfc.threshold,]

marker$group <- 'nochange'
marker$group[marker$avg_log2FC >= logfc.threshold] <- "up"
marker$group[marker$avg_log2FC <= logfc.threshold] <- "down"
table(marker$group)

formula_res <- compareCluster(OrgDb=org.Mm.eg.db,Entrez~group+cluster, 
                              data=marker, 
                              fun='enrichGO')
dotplot(formula_res)
ggsave(filename = './compareCluster.go.pdf',width =16 ,height =26 )

formula_res@compareClusterResult

library(ggpubr)
compareClusterResult <- formula_res@compareClusterResult
head(formula_res)
ggplot(formula_res,
       aes(x=cluster,y=Description))+
  geom_point(aes(size=`GeneRatio`,
                 color=`p.adjust`))+
  scale_color_gradientn(colours=c("#25c1f2","#8adcb4","#bdeb84","#f9fb42"), 
                        guide=guide_colorbar(reverse=T, 
                                             order=1))+
  scale_size_continuous(range=c(2,10))+
  xlab("RichFactor")+ 
  ylab(NULL)+
  ggtitle("BiologicalProcesses")+
  theme_set(theme_linedraw())+ 
  facet_grid(~group)

write.csv(compareClusterResult,'./compareClusterResult.csv')
saveRDS(object = formula_res,file ='./formula_res.RDS' )
ggsave(filename = './compareCluster.2.pdf',width = 11,height = 7)