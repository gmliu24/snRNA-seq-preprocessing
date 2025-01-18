#Annotate cells following preprocess20250111
#reference https://www.bilibili.com/video/BV16N411g7wA/?spm_id_from=333.788.videopod.sections&vd_source=2b730183f1e3aae2288d173e8c97f4c9


library(tidyverse) # dplyr and ggplot2
library(Seurat) # Seurat toolkit
library(hdf5r) # for data import
library(patchwork) # for plotting
library(reticulate)
library(clustree)
library(harmony)
set.seed(12345)

basedir <- "/projectnb/mccall/guangmeiliu/snrnaseq_gmliu/analysis/"
outputdir <- paste0(basedir,"preprocess20250111/")

# 1. Fetch data ===================================================
data_harmony <- readRDS('/projectnb/mccall/guangmeiliu/snrnaseq_gmliu/analysis/preprocess20250111/data_harmony.RDS')

DefaultAssay(data_harmony) <- 'soupX_SCT'
Idents(data_harmony) <- 'soupX_SCT_snn_res.0.1'

DimPlot(data_harmony, reduction = "umap", label =TRUE)
p2 <- DimPlot(data_harmony) & theme_classic(base_size = 16)
#ggsave(p2,filename=paste0(outputdir,"figures/","w1118_42d_soupX_SCT_snn_res.0.03_dimplot.pdf"),width=10,height=10)

# 2. Find differentially expressed features for broad annotation (cluster biomarkers) ===================================================
library(presto)
 # find all markers
allmarkers <- FindAllMarkers(data_harmony, test.use = 'wilcox', 
                             only.pos = TRUE, logfc.threshold = 0.25,
                             min.pct = 0.25)
 # filter markers genes based on p_val
markers = allmarkers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
 # filter top 50 genes
top50 = markers %>% group_by(cluster) %>% top_n(n=50, wt = avg_log2FC)
#view(top50)
 # Save top50 marker genes into csv
write.csv(top50, paste0(outputdir, 'data_harmony_soupX_SCT_snn_res.0.1_12c_top50genes.csv'), row.names = T)

# 3. Add broad annotation ===================================================
celltype=data.frame(clusterID=0:11, 
                    celltype='unknown')
celltype[celltype$clusterID %in% c(0, 1, 2, 4, 5, 6, 7, 9, 10, 11), 2] = c('neuron', 'epithelial cell', 'photoreceptor cell', 
                                                                          'cone cell', 
                                                                          'pigment cell', 'fat cell', 
                                                                          'muscle cell', 'lamina monopolar neuron',
                                                                          'columnar neuron T1', 'hemocyte')
celltype[celltype$clusterID %in% c(3, 8), 2] = 'glial cell'

table(celltype$celltype)
data_harmony@meta.data$celltype = 'NA'

for(i in 1:nrow(celltype)){
  data_harmony@meta.data[which(data_harmony@meta.data$soupX_SCT_snn_res.0.1 == celltype$clusterID[i]),
                           'celltype'] <- celltype$celltype[i]
}

view(data_harmony@meta.data)
tab.1=table(data_harmony@meta.data$orig.ident, data_harmony$celltype)
library(gplots)
# https://www.rdocumentation.org/packages/gplots/versions/3.2.0/topics/balloonplot
balloonplot(tab.1,
            dotsize=4/max(strwidth(19),strheight(19)),
            text.size = 1.5,
            label.size = 1.5,
            sorted = TRUE)

p2 = DimPlot(data_harmony, reduction = 'umap', group.by = 'celltype', label=T) + ggtitle('data_harmony_Broad Annotation') + theme_classic(base_size = 16) 
ggsave(p2,filename=paste0(outputdir,"figures/","data_harmony_Broad Annotation.pdf"),width=12,height=10)

p1 = DimPlot(data_harmony, reduction = 'umap', group.by = 'orig.ident', label=T) + ggtitle('data_harmony_umap') + theme_classic(base_size = 16) 
ggsave(p1,filename=paste0(outputdir,"figures/","data_harmony_umap.pdf"),width=12,height=10)

saveRDS(data_harmony, file=paste0(outputdir,"data_harmony.RDS"))

# 4. Subcluster Glial cells ===================================================
Idents(data_harmony) = 'celltype'
glialcell <- subset(data_harmony, idents = c('glial cell')) 
table(glialcell@meta.data$celltype)

 ##4.1 Set variable features before clustering ===================================================

library(glmGamPoi)
setvarfeature <- function(sobj){
  DefaultAssay(sobj)="soupXcounts"
  sobj <- SCTransform(sobj, vars.to.regress = "percent.mt", verbose = FALSE, assay="soupXcounts", new.assay.name="soupX_SCT")
  sobj <- RunPCA(sobj, npcs = 50, verbose = FALSE)
  print(DimPlot(sobj, reduction = 'pca', group.by = 'orig.ident') + ggtitle('Before Integration'))
  print(DimPlot(sobj, reduction = 'harmony', group.by = 'orig.ident') + ggtitle('Harmony Integration') )
  ElbowPlot(sobj, reduction='harmony')
}

setvarfeature(glialcell)

 ##4.2 Subclustering ===================================================
glialcell <- FindNeighbors(glialcell, reduction='harmony', dims = 1:20 )
glialcell <- FindClusters(glialcell, reduction='harmony', resolution = seq(0.01, to = 0.1, by = 0.01), verbose = FALSE)
library(clustree)
clustree(glialcell) #find optimal resolution

Idents(glialcell) <- 'soupX_SCT_snn_res.0.08'
glialcell <- RunUMAP(glialcell, reduction='harmony', dims = 1:20, verbose = FALSE)
glialcell <- RunTSNE(glialcell, reduction='harmony', dims = 1:20, verbose = FALSE)

glialcell@meta.data$soupX_SCT_snn_res.0.02 <- NULL

res = 'soupX_SCT_snn_res.0.08'
name = 'glialcell'

umap_tsne_plot <- function(sobj){
  print(DimPlot(sobj, reduction = 'umap', group.by = 'orig.ident') + ggtitle(paste0(name,'_',res)) )
  print(DimPlot(sobj, reduction = 'umap', label = T) + ggtitle(paste0(name,'_',res)))
  print(DimPlot(sobj, reduction = 'tsne', group.by = 'orig.ident') + ggtitle(paste0(name,'_',res)) )
  print(DimPlot(sobj, reduction = 'tsne', label = T) + ggtitle(paste0(name,'_',res))) 
  p1<-DimPlot(sobj, reduction = 'umap', label = T) + ggtitle(paste0(name,'_',res)) 
  ggsave(p1,filename=paste0(outputdir,"figures/","umap_", name, "_", res, ".pdf"),width=10,height=10)
  p2<-DimPlot(sobj, reduction = 'tsne', label = T)
  ggsave(p2,filename=paste0(outputdir,"figures/","tsne_", name, "_", res, ".pdf"),width=10,height=10)
}

umap_tsne_plot(glialcell)


##4.3 Find marker genes ===================================================
library(presto)

findtop50 <- function(sobj){
  # find all markers
  allmarkers <- FindAllMarkers(sobj, test.use = 'wilcox', 
                               only.pos = TRUE, logfc.threshold = 0.25,
                               min.pct = 0.25)
  # filter markers genes based on p_val
  markers = allmarkers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
  # filter top 50 genes
  top50 = markers %>% group_by(cluster) %>% top_n(n=50, wt = avg_log2FC)
  #view(top50)
  # Save top50 marker genes into csv
  write.csv(top50, paste0(outputdir, name, "_", res, '_5c_top50genes.csv'), row.names = T)
}

findtop50(glialcell)

markers <- c('repo',
              'wrapper', 'zyd', # cortex glial
             'Eaat1', 'Gat', 'alrm', 'e', 'repo', #astrocyte
             'CG4797', 'Indy', 'Vmat', 'repo', 'Tret1-1', #perineurial
             'alrm', 'Indy', 'moody', 'Gli', 'Mdr65', 'repo', #subperineurial
             'bi', #optic chiasma glial cell
             'Optix', 'bi', #lamina epithelial/marginal glial cell
  
)

DotPlot(data_list[[1]], features = markers) +
  scale_color_viridis_c(option = "magma",direction = -1) +
  theme_classic(base_size = 12) + RotatedAxis() +theme(
    panel.grid.major = element_line(color = "grey80", linetype = "dashed")
  )
##4.4. Add annotation ===================================================
celltype=data.frame(clusterID=0:11, 
                    celltype='glial cell')
celltype[celltype$clusterID %in% c(0, 1, 2, 3), 2] = c('neuron', 'epithelial cell', 'photoreceptor cell', 
                                                                           'cone cell', 
                                                                           'pigment cell', 'fat cell', 
                                                                           'muscle cell', 'lamina monopolar neuron',
                                                                           'columnar neuron T1', 'hemocyte')
celltype[celltype$clusterID %in% c(3, 8), 2] = 'glial cell'

table(celltype$celltype)
data_harmony@meta.data$celltype = 'NA'

for(i in 1:nrow(celltype)){
  data_harmony@meta.data[which(data_harmony@meta.data$soupX_SCT_snn_res.0.1 == celltype$clusterID[i]),
                         'celltype'] <- celltype$celltype[i]
}

view(data_harmony@meta.data)
tab.1=table(data_harmony@meta.data$orig.ident, data_harmony$celltype)
library(gplots)
balloonplot(tab.1)

p2 = DimPlot(data_harmony, reduction = 'umap', group.by = 'celltype', label=T) + ggtitle('data_harmony_Broad Annotation') + theme_classic(base_size = 16) 
ggsave(p2,filename=paste0(outputdir,"figures/","data_harmony_Broad Annotation.pdf"),width=12,height=10)

p1 = DimPlot(data_harmony, reduction = 'umap', group.by = 'orig.ident', label=T) + ggtitle('data_harmony_umap') + theme_classic(base_size = 16) 
ggsave(p1,filename=paste0(outputdir,"figures/","data_harmony_umap.pdf"),width=12,height=10)

saveRDS(data_harmony, file=paste0(outputdir,"data_harmony.RDS"))






