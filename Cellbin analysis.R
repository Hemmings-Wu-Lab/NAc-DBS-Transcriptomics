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
library(reshape2)
library(RColorBrewer)
set.seed(123)

options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
set.seed(123)

# -----Data Arrangement-----
# Load the processed data (RDS file)
control<- readRDS(file = "/dev16T/ccw/Depression/02_trial/04_finalFig/data/RDS/SS200000830BL_D1.cellbin.RDS")
sti <- readRDS(file = "/dev16T/ccw/Depression/02_trial/04_finalFig/data/RDS/B01320C3.cellbin.RDS")

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
ggsave(filename = paste0('QC_before','.pdf'),width=10,height=10)

# -----filter cells
filtered_seurat <- subset(
  x = alldata,
  nFeature_Spatial >= 200 &
    nFeature_Spatial <= 1500 & 
    percent_mito < 10
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
table(filtered_seurat@meta.data$anno_leiden_from_bins)
filtered_seurat@meta.data$region.main <- NA
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden_from_bins=='CC'] <- 'CC'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden_from_bins=='Cortex.1'] <- 'Meninges'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden_from_bins=='Cortex.2'] <- 'Cortex'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden_from_bins=='Cortex.3'] <- 'Cortex'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden_from_bins=='Cortex.4'] <- 'Cortex'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden_from_bins=='Cortex.5'] <- 'Cortex'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden_from_bins=='Cortex.6'] <- 'Cortex'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden_from_bins=='Cortex.7'] <- 'Cortex'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden_from_bins=='Cortex.8'] <- 'Cortex'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden_from_bins=='CPu'] <- 'CPu'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden_from_bins=='LSI'] <- 'LSI'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden_from_bins=='NAc.c'] <- 'NAc'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden_from_bins=='NAc.s'] <- 'NAc'
filtered_seurat@meta.data$region.main[filtered_seurat@meta.data$anno_leiden_from_bins=='VDB'] <- 'VDB'

table(filtered_seurat@meta.data$region.main)
saveRDS(filtered_seurat,paste0('seurat_filtered.RDS'))


# -----Harmony-----
resolution=3
filtered_seurat <- NormalizeData(filtered_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_seurat <- FindVariableFeatures(filtered_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(filtered_seurat)
filtered_seurat <- ScaleData(filtered_seurat, features = all.genes)
filtered_seurat <- RunPCA(filtered_seurat, features = VariableFeatures(object = filtered_seurat))
ElbowPlot(filtered_seurat, ndims = 20, reduction = "pca")
ggsave(filename = paste0('ElbowPlot.pdf'),width=10,height=10)
pc.num <- 1:10
filtered_seurat <- FindNeighbors(filtered_seurat, dims = pc.num)
filtered_seurat <- FindClusters(filtered_seurat, resolution = resolution)
filtered_seurat <- RunUMAP(filtered_seurat, dims = pc.num)
DimPlot(filtered_seurat,reduction = "umap", group.by='type',raster=FALSE,split.by = 'type')
saveRDS(filtered_seurat,paste0('BFharmony.reso_',resolution,'.RDS'))

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
saveRDS(seurat_harmony,paste0('seurat_harmonyed.RDS'))


# -----cluster-----
scRNA <- seurat_harmony
scRNA <- readRDS("seurat_harmonyed.RDS")
SpatialFeaturePlot(scRNA, features = c("Gfap"),
                   # cols = c("#fef4ef", "#d82422"), label = F, ncol = 3,raster=FALSE
) 


# -----subset NAc-----
Cells.sub <- subset(scRNA@meta.data, region.main=='NAc')
scRNAsub <- subset(scRNA, cells=row.names(Cells.sub))
rownames(scRNAsub)

scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2000)
scale.genes <-  rownames(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub <- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
ElbowPlot(scRNAsub, ndims=20, reduction="pca")
pc.num=1:10

scRNAsub <- FindNeighbors(scRNAsub, dims = pc.num) 
scRNAsub <- FindClusters(scRNAsub, resolution = 0.2)#--此处需要修改
table(scRNAsub@meta.data$seurat_clusters)
metadata <- scRNAsub@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cell_cluster.csv',row.names = F)

#tSNE
scRNAsub = RunTSNE(scRNAsub, dims = pc.num)
embed_tsne <- Embeddings(scRNAsub, 'tsne')
write.csv(embed_tsne,'embed_tsne.csv')
#group_by_cluster
plot1 = DimPlot(scRNAsub, reduction = "tsne",label = TRUE) 
ggsave("tSNE.pdf", plot = plot1, width = 8, height = 7)
ggsave("tSNE.png", plot = plot1, width = 8, height = 7)
#group_by_sample
plot2 = DimPlot(scRNAsub, reduction = "tsne", group.by='type',label = TRUE) 
ggsave("tSNE_sample.png", plot = plot2, width = 8, height = 7)
ggsave("tSNE_sample.pdf", plot = plot2, width = 8, height = 7)
#combinate
plotc <- plot1+plot2
ggsave("tSNE_cluster_sample.png", plot = plotc, width = 10, height = 5)
ggsave("tSNE_cluster_sample.pdf", plot = plotc, width = 10, height = 5)


#UMAP
scRNAsub <- RunUMAP(scRNAsub, dims = pc.num)
embed_umap <- Embeddings(scRNAsub, 'umap')
write.csv(embed_umap,'embed_umap.csv') 
#group_by_cluster
plot3 = DimPlot(scRNAsub, reduction = "umap", label=T) 
ggsave("UMAP.png", plot = plot3, width = 8, height = 7)
ggsave("UMAP.pdf", plot = plot3, width = 8, height = 7)
#group_by_sample
plot4 = DimPlot(scRNAsub, reduction = "umap", group.by='type',label = TRUE)
ggsave("UMAP_sample.png", plot = plot4, width = 8, height = 7)
ggsave("UMAP_sample.pdf", plot = plot4, width = 8, height = 7)
#combinate
plotc <- plot3+plot4
ggsave("UMAP_cluster_sample.png", plot = plotc, width = 10, height = 5)
ggsave("UMAP_cluster_sample.pdf", plot = plotc, width = 10, height = 5)

plotc <- plot2+plot4+ plot_layout(guides = 'collect')
ggsave("tSNE_UMAP.png", plot = plotc, width = 10, height = 5)
ggsave("tSNE_UMAP.pdf", plot = plotc, width = 10, height = 5)


# D1
scRNAsub_Drd1 <- subset(x = scRNAsub, subset = Drd2 == 0 & Drd1 > 0)
FeaturePlot(scRNAsub_Drd1,features = 'Drd1',reduction = 'tsne',pt.size = 0.5)
dim(scRNAsub_Drd1) #1031
cell_Drd1 <- colnames(scRNAsub_Drd1)

# D2
scRNAsub_Drd2 <- subset(x = scRNAsub, subset = Drd2 > 0 & Drd1==0)
FeaturePlot(scRNAsub_Drd2,features = 'Drd2',reduction = 'tsne',pt.size = 0.5)
dim(scRNAsub_Drd2) #757
cell_Drd2 <- colnames(scRNAsub_Drd2)

# D12
scRNAsub_Drd12 <- subset(x = scRNAsub, subset = Drd2 > 0 & Drd1 > 0)
FeaturePlot(scRNAsub_Drd12,features = c('Drd1','Drd2'),reduction = 'tsne',pt.size = 0.5)
dim(scRNAsub_Drd12) #50
cell_Drd12 <- colnames(scRNAsub_Drd12)

# neuron
scRNAsub_neuron <- subset(x = scRNAsub, subset = Snap25 > 0)
FeaturePlot(scRNAsub_neuron,features = c('Snap25'),reduction = 'tsne',pt.size = 0.5)
dim(scRNAsub_neuron) #7255
cell_neuron <- colnames(scRNAsub_neuron)

# IN                                                                                                      
scRNAsub_IN <- subset(x = scRNAsub, subset = Resp18 > 0)
FeaturePlot(scRNAsub_IN,features = c('Resp18'),reduction = 'tsne',pt.size = 0.5)
dim(scRNAsub_IN) #733
cell_IN <- colnames(scRNAsub_IN)

# astro
scRNAsub_astro <- subset(x = scRNAsub, subset = Gja1 > 0)
FeaturePlot(scRNAsub_astro,features = c('Gja1'),reduction = 'tsne',pt.size = 0.5)
dim(scRNAsub_astro) #905
cell_astro <- colnames(scRNAsub_astro)

# micro
scRNAsub_micro <- subset(x = scRNAsub, subset = C1qa > 0)
FeaturePlot(scRNAsub_micro,features = c('C1qa'),reduction = 'tsne',pt.size = 0.5)
dim(scRNAsub_micro) #776
cell_micro <- colnames(scRNAsub_micro)

# endo
scRNAsub_endo <- subset(x = scRNAsub, subset = Cldn5 > 0)
FeaturePlot(scRNAsub_endo,features = c('Cldn5'),reduction = 'tsne',pt.size = 0.5)
dim(scRNAsub_endo) #291
cell_endo <- colnames(scRNAsub_endo)

# oligo
scRNAsub_oligo <- subset(x = scRNAsub, subset = Mog > 0)
FeaturePlot(scRNAsub_oligo,features = c('Mog'),reduction = 'tsne',pt.size = 0.5)
dim(scRNAsub_oligo) #398
cell_oligo <- colnames(scRNAsub_oligo)

# OPC
scRNAsub_OPC <- subset(x = scRNAsub, subset = Pdgfra > 0)
FeaturePlot(scRNAsub_OPC,features = c('Pdgfra'),reduction = 'tsne',pt.size = 0.5)
dim(scRNAsub_OPC) #127
cell_OPC <- colnames(scRNAsub_OPC)

# NB
scRNAsub_NB <- subset(x = scRNAsub, subset = Top2a > 0)
FeaturePlot(scRNAsub_NB,features = c('Top2a'),reduction = 'tsne',pt.size = 0.5)
dim(scRNAsub_NB) #11
cell_NB <- colnames(scRNAsub_NB)

# neuron_NMDA
scRNAsub_NMDA <- subset(x = scRNAsub, subset = Grin1 > 0)
FeaturePlot(scRNAsub_NMDA,features = c('Grin1'),reduction = 'tsne',pt.size = 0.5)
dim(scRNAsub_NMDA) #1816
cell_NMDA <- colnames(scRNAsub_NMDA)

# neuron_AMPA
scRNAsub_AMPA <- subset(x = scRNAsub, subset = Gria2 > 0)
FeaturePlot(scRNAsub_AMPA,features = c('Gria2'),reduction = 'tsne',pt.size = 0.5)
dim(scRNAsub_AMPA) #3064
cell_AMPA <- colnames(scRNAsub_AMPA)

characters <- list(cell_Drd1,cell_Drd2,cell_IN,cell_astro,cell_micro,cell_endo,cell_oligo,cell_OPC,cell_NB,cell_NMDA,cell_AMPA)

unique_elements <- lapply(seq_along(characters), function(i) {
  other_chars <- unlist(characters[-i])
  characters[[i]][!characters[[i]] %in% other_chars]
})

names(unique_elements) <- c('cell_Drd1','cell_Drd2','cell_IN','cell_astro','cell_micro','cell_endo','cell_oligo','cell_OPC','cell_NB','cell_NMDA','cell_AMPA')
unique_elements

# -----add metadata
table(scRNAsub@meta.data$region.main)
table(scRNAsub@meta.data[cell_Drd1,'type'])
table(scRNAsub@meta.data[cell_Drd2,'type'])
table(scRNAsub@meta.data[cell_IN,'type'])
table(scRNAsub@meta.data[cell_astro,'type'])
table(scRNAsub@meta.data[cell_micro,'type'])
table(scRNAsub@meta.data[cell_endo,'type'])
table(scRNAsub@meta.data[cell_oligo,'type'])
table(scRNAsub@meta.data[cell_OPC,'type'])
table(scRNAsub@meta.data[cell_NB,'type'])

scRNAsub@meta.data$celltype_bymarker = 'UD'
scRNAsub@meta.data[cell_Drd1,'celltype_bymarker'] = 'cell_Drd1'
scRNAsub@meta.data[cell_Drd2,'celltype_bymarker']  = 'cell_Drd2'
scRNAsub@meta.data[cell_NMDA,'celltype_bymarker']  = 'cell_NMDA'
scRNAsub@meta.data[cell_AMPA,'celltype_bymarker']  = 'cell_AMPA'
scRNAsub@meta.data[cell_IN,'celltype_bymarker'] = 'cell_IN'
scRNAsub@meta.data[cell_astro,'celltype_bymarker']  = 'cell_astro'
scRNAsub@meta.data[cell_micro,'celltype_bymarker'] = 'cell_micro'
scRNAsub@meta.data[cell_endo,'celltype_bymarker']  = 'cell_endo'
scRNAsub@meta.data[cell_oligo,'celltype_bymarker'] = 'cell_oligo'
scRNAsub@meta.data[cell_OPC,'celltype_bymarker']  = 'cell_OPC'
scRNAsub@meta.data[cell_NB,'celltype_bymarker'] = 'cell_NB'
table(scRNAsub@meta.data$celltype_bymarker)

saveRDS(scRNAsub,file = 'scRNAsub.NAc.RDS')


# -------------cell count / NAc bin100 count--------------------
load("scRNAsub_final.Rdata")
table(scRNAsub@meta.data$celltype_bymarker)

# ---cellbin cell count
pbmc <- scRNAsub
data_plotC <- table(pbmc@meta.data$type, pbmc@meta.data$celltype_bymarker) %>% melt()
colnames(data_plotC) <- c("Sample", "CellType","Number")
data_plotC
# ---NAc bin100 count
scRNA <- readRDS("03_FindClusters_3.RDS")
data_NAc <- table(scRNA@meta.data$type, scRNA@meta.data$anno_leiden) %>% melt()
colnames(data_NAc) <- c("Sample", "CellType","Number")
data_NAc 
# ----cell count / NAc bin100 count
data_plotC$ratio <- 0
data_plotC[data_plotC$Sample=='control','ratio'] <- data_plotC[data_plotC$Sample=='control','Number']/2042
data_plotC[data_plotC$Sample=='sti','ratio'] <- data_plotC[data_plotC$Sample=='sti','Number']/1461
data_plotC$combinetype <- paste0(data_plotC$Sample,data_plotC$CellType)
data_plotD <- data_plotC[!data_plotC$CellType%in%c('cell_neuron_other','cell_undefined','cell_endo'),]

library(RColorBrewer)
color_count <- length(unique(data_plotD$combinetype))
color_palette <- brewer.pal(n = min(color_count, 12), name = "Set3")

if(color_count > 12) {
  color_palette <- color_palette[rep(1:12, length.out = color_count)]
}

ggplot(data_plotD, aes(x = combinetype, y = ratio, fill = combinetype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  labs(title = "Value by Cell Type", x = "Cell Type", y = "Value")
ggsave(filename = 'cell.count_NAc.bin100.count.pdf',width =10 ,height =5 )
write.csv(x = data_plotD,file = 'cell.count_NAc.bin100.count.csv')

write.csv(x = data_plotC,file = './data_plotC_before.csv')
data_plotC <- read.csv(file = './data_plotC_after.csv',header = T)

colourCount = length(unique(pbmc@meta.data$celltype_bymarker))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
celltype_colors <- getPalette(colourCount)
celltype_colors <- c('#E6E9A5','#AFE4DE','#3DA5B8','#A7E2B3','#256A6F','#E1B3E6','#71A2D0','#44C1A5','#2E8A62','#D18C75','#B85A3D')
celltype_colors <- c('#E6E9A5','#3DA5B8','#A7E2B3','#71A2D0','#256A6F','#E1B3E6','#AFE4DE','#44C1A5','#2E8A62','#D18C75','#B85A3D')
pC1 <- ggplot(data = data_plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="stack")+
  scale_fill_manual(values=celltype_colors) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Average number")+
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6)) 

pC2 <- ggplot(data = data_plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=celltype_colors) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6)) 

pC <- pC1 + pC2 + plot_layout(ncol = 2, widths = c(1,1),guides = 'collect')
pC

ggsave(filename = './proportion.pdf',width =10 ,height =10 )
save(scRNAsub,file = 'scRNAsub_final.Rdata')


# -----DEGs------
# filter cells
table(scRNAsub@meta.data$celltype_bymarker)
seurat_obj <-subset(x = seurat_obj,subset = celltype_bymarker !='cell_endo')  
table(seurat_obj@meta.data$celltype_bymarker)

Idents(seurat_obj)
table(Idents(seurat_obj))
Idents(seurat_obj) <- seurat_obj@meta.data$celltype_bymarker
Idents(seurat_obj)
table(Idents(seurat_obj))

cell_types <- unique(seurat_obj@meta.data$celltype_bymarker)
cell_types

diff_genes_list_cell <- list()
for (cell_type in cell_types) {
  diff_genes <- FindMarkers(seurat_obj, ident.1 = cell_type)
  diff_genes$gene <- rownames(diff_genes)
  diff_genes_list_cell[[cell_type]] <- diff_genes
}

# ---DEG1 sti vs control 
diff_genes_list <- list()
for (cell_type in cell_types) {
  diff_genes <- FindMarkers(seurat_obj, ident.1 = "sti", group.by = 'type', subset.ident = cell_type) # ---logfc.threshold默认为0.1
  diff_genes$gene <- rownames(diff_genes)
  diff_genes_list[[cell_type]] <- diff_genes
}

final_results <- data.frame()
p_value_threshold <- 0.05

for (cell_type in names(diff_genes_list)) {
  current_results <- diff_genes_list[[cell_type]]
  filtered_results <- current_results[current_results$p_val < p_value_threshold, ]
  sorted_results <- filtered_results[order(filtered_results$avg_log2FC), ]
  sorted_results$cell_type <- cell_type
  final_results <- rbind(final_results, sorted_results)
}

table(final_results$cell_type)

table(final_results[,'cell_type'])
table(final_results[final_results$p_val_adj<0.05,'cell_type'])
final_results.filtedbyadP <- final_results[final_results$p_val_adj<0.05,]

write.csv(x = final_results,file = 'final_results.csv')
write.csv(x = final_results.filtedbyadP,file = 'final_results.filtedbyadP.csv')
saveRDS(object = final_results,file = 'final_results.RDS')
saveRDS(object = final_results.filtedbyadP,file = 'final_results.filtedbyadP.RDS')

seurat_obj$combine_type <- paste0(seurat_obj$celltype_bymarker,"_",seurat_obj$type)
table(seurat_obj$combine_type)
saveRDS(seurat_obj,file = 'seurat_obj_fig4.dotplot.DEG.RDS')

seurat_obj <- readRDS("seurat_obj_fig4.dotplot.DEG.RDS")
final_results.filtedbyadP <- read.csv(file = './final_results.filtedbyadP.csv',row.names = 1)

table(final_results.filtedbyadP$cell_type)
final_results.filtedbyadP <- final_results.filtedbyadP[unique(final_results.filtedbyadP$gene),]
final_results.filtedbyadP <- na.omit(final_results.filtedbyadP)

genes_of_interest2 <- c('Drd1','Drd2','Gria2','Grin1','Resp18','Gja1','C1qa','Mog','Pdgfra')

genes_of_interest <- final_results.filtedbyadP$gene
genes_of_interest1 <- c('Cfl1','Fth1','Nkain2','Kcnq5','Gm15564','Nlgn1','Phactr1','Celf2','Pde10a','Syt1','Atp2b2','Snca','Hspa8','Rpl38','Gm48099','Stk31','Rps27') 

p1 <- DotPlot(seurat_obj,
        features = genes_of_interest1,group.by = 'combine_type',
        cols = c("#ffffff", "#448444")
  ) +
  RotatedAxis() + 
  theme(
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(1, "mm"),
    axis.title = element_blank(),
  )+
  scale_color_gradientn(colours = c('#F4BD0B','white','#25DAA6')) 
p1
ggsave(filename = 'fig4_dotplot_combine_p1_v3.pdf',width =11 ,height =5)

p2 <- DotPlot(seurat_obj,
        features = genes_of_interest2,group.by = 'combine_type',
        cols = c("#ffffff", "#448444")
) +
  RotatedAxis() + 
  theme(
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(1, "mm"),
    axis.title = element_blank(),
  )+
  scale_color_gradientn(colours = c('#F4BD0B','#25DAA6')) #颜色
p2

pC <- p1 + p2 + plot_layout(ncol = 2, widths = c(1,1),guides = 'collect')
pC
ggsave(filename = 'fig4_dotplot_combine_p1_v3.pdf',width =12 ,height =5)
ggsave(filename = 'fig4_dotplot_combine_p2_v3.pdf',width =8.75 ,height =5 )