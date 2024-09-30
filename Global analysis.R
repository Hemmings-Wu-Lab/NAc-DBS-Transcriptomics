library(dplyr)
library(tidyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(clustree)
library(kableExtra)
library(harmony)
library(lisi)
library(Matrix)
set.seed(123)

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
# -----Data Arrangement-----
# Load the processed data (RDS file)
control<- readRDS(file = "/SS200000830BL_D1.bin100.RDS")
sti <- readRDS(file = "/B01320C3.bin100.RDS")

# Adding Group Information into metadata
sti$type = "sti"
control$type = "control"

# Merge all datasets into one single seurat object
alldata <- merge(sti, c(control), add.cell.ids = c("sti", "control"))
as.data.frame(alldata@assays$Spatial@counts[1:10, 1:2])
head(alldata@meta.data, 10)

# -----QC caculation-----
alldata$log10GenesPerUMI <- log10(alldata$nFeature_Spatial)/log10(alldata$nCount_Spatial)
grep('^mt',rownames(alldata),value=T)
grep('^Rp[sl]',rownames(alldata),value=T)
grep('^Hb[^(p)]',rownames(alldata),value=T)
alldata <- PercentageFeatureSet(alldata, "^mt-", col.name = "percent_mito")
alldata <- PercentageFeatureSet(alldata, "^Rp[sl]", col.name = "percent_ribo")
alldata <- PercentageFeatureSet(alldata, "^Hb[^(p)]", col.name = "percent_hb")

features <- c("nFeature_Spatial", "nCount_Spatial", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(alldata, features = features, pt.size = 0.1, ncol = 3,raster=FALSE) + NoLegend()
dir.create('QC')
ggsave(filename = paste0('./QC/QC_before','.pdf'),width=10,height=10)


# -----filter cells
filtered_seurat <- subset(
    x = alldata,
    nFeature_Spatial >= 200 &
    nFeature_Spatial <= 6200 & 
    percent_mito < 25
)

# -----filter genes
# Filter MALAT1 (Usually,based on paper, the high-level MALAT1 genes are judged as mainly
x <- c('Gm42418','AY036118',"Malat1","mt-Atp6","mt-Atp8","mt-Co1","mt-Co2","mt-Co3","mt-Cytb","mt-Nd1","mt-Nd2","mt-Nd3","mt-Nd4","mt-Nd4l","mt-Nd5","mt-Nd6")
dim(filtered_seurat)
filtered_seurat <- subset(filtered_seurat,features=setdiff(rownames(filtered_seurat),x))
dim(filtered_seurat)

features <- c("nFeature_Spatial", "nCount_Spatial", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(filtered_seurat, features = features, pt.size = 0.1, ncol = 3,raster=FALSE) + NoLegend()
dir.create('QC')
ggsave(filename = paste0('./QC/QC_after','.pdf'),width=10,height=10)

# -----filter clusters
colnames(filtered_seurat@meta.data)
table(filtered_seurat@meta.data$anno_leiden)
filtered_seurat@meta.data$region.main <- NA
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='CC'] <- 'CC'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='Cortex.1'] <- 'Meninges'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='Cortex.2'] <- 'Cortex'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='Cortex.3'] <- 'Cortex'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='Cortex.4'] <- 'Cortex'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='Cortex.5'] <- 'Cortex'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='Cortex.6'] <- 'Cortex'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='Cortex.7'] <- 'Cortex'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='Cortex.8'] <- 'Cortex'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='CPu'] <- 'CPu'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='Interval'] <- 'Interval'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='Lesion'] <- 'Lesion'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='LSI'] <- 'LSI'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='NAc.c'] <- 'NAc'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='NAc.s'] <- 'NAc'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden=='VDB'] <- 'VDB'

table(filtered_seurat@meta.data$region.main)
saveRDS(filtered_seurat,file = '01_filtered_seurat.RDS')


# -----Harmony-----
filtered_seurat <- readRDS("01_filtered_seurat.RDS")
seurat_harmony <- SCTransform(assay = "Spatial",filtered_seurat, vst.flavor = "v2", method = "glmGamPoi"
                              ) %>% 
  RunPCA(verbose = TRUE)

table(seurat_harmony@meta.data$batch)
seurat_harmony <- RunHarmony(seurat_harmony, group.by.vars = c("batch"), plot_convergence = TRUE,
                             theta = c(1), assay.use = "SCT")
seurat_harmony <- RunUMAP(seurat_harmony, reduction = "harmony", dims = 1:25)
colnames(seurat_harmony@meta.data)
table(seurat_harmony@meta.data$type)
DimPlot(seurat_harmony,reduction = "umap", group.by='type',raster=FALSE)
DimPlot(seurat_harmony, reduction = "umap", group.by='type',split.by = "type",raster=FALSE)
plot11 & theme(plot.title = element_text(hjust = 0.5))
head(seurat_harmony@meta.data, 10)
saveRDS(seurat_harmony,file = ('./02_seurat_harmonyed.RDS'))


# -----cluster(聚类)-----
scRNA <- load('./02_seurat_harmonyed.RDS')
resolution=3
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
ElbowPlot(scRNA, ndims = 20, reduction = "pca")
dir.create('cluster')
ggsave(filename = paste0('./cluster/ElbowPlot.pdf'),width=10,height=10)
pc.num <- 1:15
scRNA <- FindNeighbors(scRNA, dims = pc.num)
scRNA <- FindClusters(scRNA, resolution = resolution)
scRNA <- RunUMAP(scRNA, dims = pc.num)
saveRDS(scRNA,file =paste0('./03_FindClusters_',resolution,'.RDS'))



# -----DEGs analysis-----
scRNA <- load("./cluster/FindClusters_3.RDS")
marker <- FindMarkers(scRNA,features=rownames(scRNA), ident.1 = "sti", group.by = 'type',raster=FALSE,
                      logfc.threshold = 0,min.pct=0.1 
                      )
marker.filter <- marker[abs(marker$avg_log2FC) >=0.25 & marker$p_val_adj <0.05,]
dir.create('DEG')
write.csv(marker,file = paste0('./DEG/DEG.global.csv'))
write.csv(marker.filter,file = paste0('./DEG/DEG.global.filtered.csv'))
saveRDS(marker,'./DEG/marker.RDS')

# -----GO analysis-----
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

dir.create('s_fig2')

logfc.threshold <- 0.2 
marker <- readRDS("/dev16T/ccw/Depression/02_trial/04_finalFig/output/fig2_global及脑区DEG/DEG/marker.RDS")
marker$gene <- rownames(marker)
marker.filter <- marker[abs(marker$avg_log2FC) >=logfc.threshold & marker$p_val_adj <0.05,]

marker$group <- 'nochange'
marker$group[marker$avg_log2FC >= logfc.threshold] <- "up"
marker$group[marker$avg_log2FC <= -logfc.threshold] <- "down"
table(marker$group)

up_genes <- marker[marker$group=='up','gene']
down_genes <- marker[marker$group=='down','gene']

# -----BP
ego <- enrichGO(gene         = up_genes,
                   OrgDb        = org.Mm.eg.db,
                   keyType      = "SYMBOL",
                   ont          = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable     = TRUE)
ego_bp.up <- simplify(ego, measure="Wang", cutoff=0.7, by="p.adjust", select_fun=min)

ego <- enrichGO(gene         = down_genes,
                     OrgDb        = org.Mm.eg.db,
                     keyType      = "SYMBOL",
                     ont          = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05,
                     readable     = TRUE)
ego_bp.down <- simplify(ego, measure="Wang", cutoff=0.7, by="p.adjust", select_fun=min)

# -----CC
ego <- enrichGO(gene         = up_genes,
                      OrgDb        = org.Mm.eg.db,
                      keyType      = "SYMBOL",
                      ont          = "CC",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable     = TRUE)
ego_cc.up <- simplify(ego, measure="Wang", cutoff=0.7, by="p.adjust", select_fun=min)

ego <- enrichGO(gene         = down_genes,
                        OrgDb        = org.Mm.eg.db,
                        keyType      = "SYMBOL",
                        ont          = "CC",
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05,
                        readable     = TRUE)
ego_cc.down <- simplify(ego, measure="Wang", cutoff=0.7, by="p.adjust", select_fun=min)

# -----MF
ego <- enrichGO(gene         = up_genes,
                      OrgDb        = org.Mm.eg.db,
                      keyType      = "SYMBOL",
                      ont          = "MF",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable     = TRUE)
ego_mf.up <- simplify(ego, measure="Wang", cutoff=0.7, by="p.adjust", select_fun=min)

ego <- enrichGO(gene         = down_genes,
                        OrgDb        = org.Mm.eg.db,
                        keyType      = "SYMBOL",
                        ont          = "MF",
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05,
                        readable     = TRUE)
ego_mf.down <- simplify(ego, measure="Wang", cutoff=0.7, by="p.adjust", select_fun=min)



# ------barplot
# -----BP
barplot_ego.bp.up <- barplot(ego_bp.up, showCategory=15) +
  scale_fill_gradient(low = "#BC9FFA", high = "#EAE0FE") +
  theme_minimal()
print(barplot_ego.bp.up)
ggsave(filename = './s_fig2/ego_bp.up.pdf',width = 7,height = 4)

barplot_ego.bp.down <- barplot(ego_bp.down, showCategory=15) +
  scale_fill_gradient(low = "#BC9FFA", high = "#EAE0FE") +
  theme_minimal()
print(barplot_ego.bp.down)
ggsave(filename = './s_fig2/ego_bp.down.pdf',width = 7,height = 4)

# -----CC
barplot_ego.cc.up <- barplot(ego_cc.up, showCategory=15) +
  # scale_fill_gradient(low = "#4FBAF6", high = "#B6E3FB") +
  scale_fill_gradient(low = "#7B90EF", high = "#BFC9F7") +
  theme_minimal()
print(barplot_ego.cc.up)
ggsave(filename = './s_fig2/ego_cc.up.pdf',width = 7,height = 4)

barplot_ego.cc.down <- barplot(ego_cc.down, showCategory=15) +
  scale_fill_gradient(low = "#7B90EF", high = "#BFC9F7") +
  theme_minimal()
print(barplot_ego.cc.down)
ggsave(filename = './s_fig2/ego_cc.down.pdf',width = 7,height = 4)

# -----MF
barplot_ego.mf.up <- barplot(ego_mf.up, showCategory=15) +
  scale_fill_gradient(low = "#AAAAAA", high = "#E2E2E2") +
  theme_minimal()
print(barplot_ego.mf.up)
ggsave(filename = './s_fig2/ego_mf.up.pdf',width = 7,height = 4)

barplot_ego.mf.down <- barplot(ego_mf.down, showCategory=15) +
  scale_fill_gradient(low = "#AAAAAA", high = "#E2E2E2") +
  theme_minimal()
print(barplot_ego.mf.down)
ggsave(filename = './s_fig2/ego_mf.down.pdf',width = 7,height = 4)


# -----Volcano plot
library(tidyverse)
library(ggrepel)
library(ggplot2)
marker <- readRDS("./DEG/marker.RDS")
DEG=marker
head(DEG) 
avg_log2FC_cutoff <- 0.2
data <- 
  DEG %>% 
  mutate(change = as.factor(ifelse(p_val < 0.05 & abs(avg_log2FC) > avg_log2FC_cutoff,
                                   ifelse(avg_log2FC > avg_log2FC_cutoff ,'Up','Down'),'No change')))
data$gene <- rownames(data)
head(data)
table(data$change)

ggplot(data,aes(avg_log2FC, -log10(p_val_adj)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1.2,1.2), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 1.75, 
             alpha = 0.5) +
             # size = 0.2, 
             # alpha = 0.5) +
  theme_bw(base_size = 12)+
  ggsci::scale_color_jama() +
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  geom_text_repel(data = filter(data, abs(avg_log2FC) > avg_log2FC_cutoff & -log10(p_val_adj) > 38),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = gene, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(FDR q-value)")

core_gene <- data$gene[data$change!='No change']

ggplot(data,aes(avg_log2FC, -log10(p_val_adj)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-avg_log2FC_cutoff,avg_log2FC_cutoff), linetype = "dashed", color = "black")+
  geom_point(aes(size = 1.5,alpha=0.4,
                 color = -log10(p_val_adj)))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#f9fb42","#bdeb84","#8adcb4","#25c1f2"))+
  scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.2)) + # 设置x轴范围为0到10，间距为2
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        legend.justification = c(0,1))+

  guides(col = 
           guide_colorbar(title = "-Log10_q-value",
                          ticks.colour = NA,
                          reverse = T,
                          title.vjust = 0.8,
                          barheight = 8,
                          barwidth = 1),
         size = "none") +
  xlab("Log2FC")+
  ylab("-Log10(FDR q-value)")


ggsave(width = 9.5 ,height = 5.5 ,filename = './Volcano.pdf')

# -----GO analysis
head(data)
table(data$change)
gene_down <- as.data.frame(data$gene[data$change=='Down'])
gene_up   <- as.data.frame(data$gene[data$change=='Up'])
colnames(gene_down) <- 'gene'
colnames(gene_up) <- 'gene'
dir.create('enrich')
write_csv(gene_down,'./enrich/gene.0.2_down.csv')
write_csv(gene_up,'./enrich/gene.0.2_up.csv')

library(clusterProfiler)
data$ensemblid <- as.character(na.omit(bitr(data$gene, 
                                            fromType="SYMBOL", 
                                            toType="ENSEMBL", 
                                            OrgDb="org.Hs.eg.db")[,2])) 
gene_down_entrez <- as.character(na.omit(bitr(gene_down, 
                                              fromType="SYMBOL", 
                                              toType="ENTREZID", 
                                              OrgDb="org.Mm.eg.db")[,2])) 
gene_diff_entrez <- unique(c(gene_up_entrez ,gene_down_entrez ))