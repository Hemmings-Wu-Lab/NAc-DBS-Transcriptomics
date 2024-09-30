library(dplyr)
library(rjson)
library(Seurat)
library(ggplot2)
library(argparser)
library(SeuratDisk)

data_all <- readRDS("/dev16T/ccw/Depression/02_trial/02_seurat/data/cellbin/seurat_filtered.RDS")
table(data_all@meta.data$type)

data_contol <- subset(data_all,type='control')
out_path <- '/dev16T/ccw/Depression/02_trial/02_seurat/data/cellbin/seurat_filtered.control.h5Seurat'
SaveH5Seurat(data_contol,filename=out_path)
Convert(out_path,dest='h5ad')

data_sti <- subset(data_all,type='sti')
out_path <- '/dev16T/ccw/Depression/02_trial/02_seurat/data/cellbin/seurat_filtered.sti.h5Seurat'
SaveH5Seurat(data_sti,filename=out_path)
Convert(out_path,dest='h5ad')

