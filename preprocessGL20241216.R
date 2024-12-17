####################################################################################
# Run QC on drprnull42d and w111842d cells separately 
# QC using median absolute deviations (MAD) method, Ambient RNA correction using SoupX package, Normalization using log10 transformation
# Outputs AnnData .h5ad files for both samples separately
# Code adpated from https://biostatsquid.com/scrnaseq-preprocessing-workflow-seurat/#step3
####################################################################################
#Didn't remove doublets
## Load required packages

library(tidyverse) # dplyr and ggplot2
library(Seurat) # Seurat toolkit
library(hdf5r) # for data import
library(patchwork) # for plotting
library(SeuratDisk)

# 1. Import data ===================================================
 ## Set up I/O paths

 basedir <- "/projectnb/mccall/guangmeiliu/snrnaseq_gmliu/analysis/"
 inputdir <- paste0(basedir,"cellranger/")
 outputdir <- paste0(basedir,"preprocess20241216/")
 dir.create(outputdir,showWarnings = FALSE)
 dir.create(paste0(outputdir,"figures"),showWarnings = FALSE)

 ## Load raw counts
 read_raw_sobj <- function(sample_id){
  inputdir <- "/projectnb/mccall/guangmeiliu/snrnaseq_gmliu/analysis/cellranger/"
  path <- paste0(inputdir, sample_id, "/outs/filtered_feature_bc_matrix/")
  sobj <- Read10X(data.dir = path)
  sobj <- CreateSeuratObject(counts = sobj, min.cells = 3, min.features = 200) #Seurat Object
  sobj$sample_id <- sample_id
  #add QC metrics
  sobj$log1p_total_counts <- log1p(sobj@meta.data$nCount_RNA)
  sobj$log1p_n_genes_by_counts <- log1p(sobj@meta.data$nFeature_RNA)
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt:") #add mitochondrial percentage
  
  return(sobj)
 }
 
# 2. Create your Seurat object (raw counts) ===========================================
samples <- c('w1118_42D', 'drprnull_42D')
data_list <- sapply(samples, read_raw_sobj)

SaveH5Seurat(data_list$w1118_42D, filename = paste0(outputdir,"w1118_42d_raw.h5seurat"), overwrite = TRUE)
SaveH5Seurat(data_list$drprnull_42D, filename = paste0(outputdir,"drprnull_42d_raw.h5seurat"), overwrite = TRUE)

 ##Visualize features
 plot_features <- function(sobjmetadata){
  p0 <- ggplot(as.data.frame(sobjmetadata), aes(x=nCount_RNA,y=nFeature_RNA, color=percent.mt))+
    geom_point()+
    scale_color_viridis_c(option = "magma")+
    theme_classic(base_size = 16)
  
  p1 <- ggplot(as.data.frame(sobjmetadata),aes(x=nCount_RNA))+
    geom_histogram()+
    theme_classic(base_size = 16)
  
  p2 <- ggplot(as.data.frame(sobjmetadata),aes(x=nFeature_RNA))+
    geom_histogram()+
    theme_classic(base_size = 16)
  
  p3 <- ggplot(as.data.frame(sobjmetadata),aes(x=log1p_total_counts))+
    geom_histogram()+
    theme_classic(base_size = 16)
  
  p4 <- ggplot(as.data.frame(sobjmetadata),aes(x=log1p_n_genes_by_counts))+
    geom_histogram()+
    theme_classic(base_size = 16)
  
  p5 <- ggplot(as.data.frame(sobjmetadata), aes(x=percent.mt))+
    geom_histogram(bins=20)+
    theme_classic(base_size = 16)
  
  #featureplots <- list(p0,p1,p2,p3,p4,p5)
  p0+p1+p2+p3+p4+p5 
  
}
 ## Plot features before filtering
 plot_features(data_list[[1]]@meta.data)
 plot_features(data_list[[2]]@meta.data)

# 3. QC ===================================================
 #QC steps contain filtering low quality cells (based on nCounts, nFeatures, percent.mt), correction of ambient RNA, doublet detection. https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
 ## perc_mt -----------------------
 ### View w1118_42d
 View(data_list[[1]]@meta.data)
 VlnPlot(data_list[[1]], features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
 FeatureScatter(data_list[[1]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
   geom_smooth(method = 'lm')
 
 ### View drprnull_42d
 View(data_list[[2]]@meta.data)
 VlnPlot(data_list[[2]], features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
 FeatureScatter(data_list[[2]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
   geom_smooth(method = 'lm')
 
 ###QC metrics (stored as metadata); 
 ##### A value is considered an outlier if it is greater than 3 MADs away (both directions) from the median of these two metrics

 mad_outlier <- function(sobj, metric, nmads){
  M <- sobj@meta.data[[metric]]
  median_M <- median(M, na.rm = TRUE)
  mad_M <- mad(M, na.rm = TRUE)
  outlier <- (M < (median_M - nmads * mad_M)) | (M > (median_M + nmads * mad_M))
  return(outlier)
 }

 flag_outliers <- function(sobj){
  
  ncount_outlier <- mad_outlier(sobj, 'log1p_total_counts', 5) 
  names(ncount_outlier) <- colnames(sobj)
  sobj <- AddMetaData(object = sobj,metadata = ncount_outlier,col.name = 'nCount_outlier')
  
  nfeature_outlier <- mad_outlier(sobj, 'log1p_n_genes_by_counts', 5) 
  names(nfeature_outlier) <- colnames(sobj)
  sobj <- AddMetaData(object = sobj,metadata = nfeature_outlier,col.name = 'nFeature_outlier')
  
  mito_outlier <- mad_outlier(sobj, 'percent.mt', 3) 
  names(mito_outlier) <- colnames(sobj)
  sobj <- AddMetaData(object = sobj,metadata = mito_outlier,col.name = 'mt_outlier')
 }

 data_list <- sapply(data_list, flag_outliers)

## Calculate the percent of outliers
outlier_counts <- lapply(data_list, function(sobj) {
  list(
    nCount_outliers = sum(sobj@meta.data$nCount_outlier, na.rm = TRUE)/sum(sobj@assays$RNA$counts, na.rm = TRUE)*100,
    nFeature_outliers = sum(sobj@meta.data$nFeature_outlier, na.rm = TRUE)/sum(sobj@assays$RNA$counts, na.rm = TRUE)*100,
    mito_outliers = sum(sobj@meta.data$mt_outlier, na.rm = TRUE)/sum(sobj@assays$RNA$counts, na.rm = TRUE)*100
  )
})
print(outlier_counts)

## Plot outliers before filtering
plot_outliers <- function(sobjmetadata){
  ncount_outlier_plot <- ggplot(data = sobjmetadata, aes(x=sample_id,y=log1p_n_genes_by_counts,color=nFeature_outlier))+
    #geom_violin(alpha=0.5)+
    geom_jitter(alpha=0.5)+
    theme_classic(base_size = 16) +
    labs(title = "nCount Outliers")
  nfeature_outlier_plot <- ggplot(data = sobjmetadata, aes(x=sample_id,y=log1p_total_counts,color=nCount_outlier))+
    #geom_violin(alpha=0.5) +
    geom_jitter(alpha=0.5) +
    theme_classic(base_size = 16) +
    labs(title = "nFeature Outliers")
  mt_outlier_plot <- ggplot(data = sobjmetadata, aes(x=sample_id,y=percent.mt,color=mt_outlier)) +
    #geom_violin(alpha=0.5)+
    geom_jitter(alpha=0.5) +
    theme_classic(base_size = 16)+
    labs(title = "Mitochondrial Outliers")
  return(ncount_outlier_plot + nfeature_outlier_plot + mt_outlier_plot)
}

plot_outliers(data_list[[1]]@meta.data)
plot_outliers(data_list[[2]]@meta.data)

## Filtering (best-practice method MAD)-----------------------
 filter_by_counts <- function(sobj){
  ##find outliers and subset
  bool_vector <- !mad_outlier(sobj, 'log1p_total_counts', 5) & !mad_outlier(sobj, 'log1p_n_genes_by_counts', 5) & !mad_outlier(sobj, 'percent.mt', 3) 
  sobj <- subset(sobj, cells = which(bool_vector))
  sobj <- subset(sobj, subset = percent.mt < 5) ## Remove cells with percent.mt > 5% because mito counts in snRNAseq indicates incomplete cell lysis.
  return(sobj)
 }

 data_list <- sapply(data_list, filter_by_counts)

## Plot features after filtering
 plot_features(data_list[[1]]@meta.data)
 plot_features(data_list[[2]]@meta.data)

 data_list[1]$w1118_42D ##9306 cells x 11652 genes 
 data_list[2]$drprnull_42D ##7219 cells x 11903 genes 

 SaveH5Seurat(data_list$w1118_42D, filename = paste0(outputdir,"w1118_42d_filtered.h5Seurat"), overwrite = TRUE)
 SaveH5Seurat(data_list$drprnull_42D, filename = paste0(outputdir,"drprnull_42d_filtered.h5Seurat"), overwrite = TRUE)

## Correction of ambient RNA using SoupX -----------------------
library(SoupX)

get_soup_groups <- function(sobj) {
  sobj <- NormalizeData(sobj, verbose = FALSE)
  sobj <- FindVariableFeatures(object = sobj, verbose = FALSE, selection.method = 'vst')
  sobj <- ScaleData(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, npcs = 20, verbose = FALSE)
  sobj <- FindNeighbors(sobj, dims = 1:20, verbose = FALSE)
  sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)
  sobj <- RunUMAP(sobj, dims = 1:10)
  return(sobj)
}

add_soup_groups <- function(sobj) {
  sobj <- get_soup_groups(sobj)
  sobj$soup_group <- sobj@meta.data[['seurat_clusters']]
  return(sobj)
}
## Apply to each Seurat object in the list
data_list <- lapply(data_list, add_soup_groups)

## Run SoupX
get_soup_groups <- function(sobj){
  sobj <- NormalizeData(sobj, verbose = FALSE)
  sobj <- FindVariableFeatures(object = sobj, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')
  sobj <- ScaleData(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, npcs = 20, verbose = FALSE)
  sobj <- FindNeighbors(sobj, dims = 1:20, verbose = FALSE)
  sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)
  sobj <- RunUMAP(sobj, dims = 1:10)
  return(sobj)
}
 add_soup_groups <- function(sobj){
  sobj <- get_soup_groups(sobj)
  sobj$soup_group <- sobj@meta.data[['seurat_clusters']]
  return(sobj)
 }

 data_list <- sapply(data_list, add_soup_groups)

 make_soup <- function(sobj){
  inputdir <- "/projectnb/mccall/guangmeiliu/snrnaseq_gmliu/analysis/cellranger/"
  sample_id <- as.character(sobj$sample[1]) 
  path <- paste0(inputdir,sample_id, "/outs/raw_feature_bc_matrix/")
  raw <- Read10X(data.dir = path)
  
  sc = SoupChannel(raw,sobj@assays$RNA$counts) ## Create SoupChannel with raw and filtered data
  sc = setClusters(sc,sobj$soup_group) # Assign soup groups
  sc = autoEstCont(sc, doPlot=FALSE) # Estimate contamination
  out = adjustCounts(sc, roundToInt = TRUE) # Correct counts
  sobj[["soupXcounts"]] <- CreateAssayObject(counts = out)   #save SoupX corrected counts
  
  return(sobj)
 }

###  Check SoupX correction - % of RNA removed

sum(data_list[[1]]@assays$RNA$counts)
sum(data_list[[1]]@assays$soupXcounts$counts)/sum(data_list[[1]]@assays$RNA$counts)*100

sum(data_list[[2]]@assays$RNA$counts)
sum(data_list[[2]]@assays$soupXcounts$counts)/sum(data_list[[2]]@assays$RNA$counts)*100

detach("package:SoupX", unload = TRUE)

# Doublet Detection
library(scater)
library(scDblFinder)
library(BiocParallel)

set.seed(123)

doublets <- function(sobj) {
  # Extract counts from the RNA assay
  counts_matrix <- GetAssayData(sobj, slot = "counts", assay = "RNA")
  # Create SingleCellExperiment object
  sce <- SingleCellExperiment(list(counts = counts_matrix))
  # Run scDblFinder to identify doublets
  sce <- scDblFinder(sce)
  # Extract doublet scores and classifications
  doublet_score <- sce$scDblFinder.score
  doublet_class <- sce$scDblFinder.class
  
  # Add to Seurat object's metadata
  sobj <- AddMetaData(sobj, metadata = doublet_score, col.name = "scDblFinder_score")
  sobj <- AddMetaData(sobj, metadata = doublet_class, col.name = "scDblFinder_class")
  
  return(sobj)
}

# Apply to all objects in data_list
data_list <- lapply(data_list, doublets)

summary(data_list[[1]]@meta.data$scDblFinder_class) #portion of doublets vs. singlet
summary(data_list[[2]]@meta.data$scDblFinder_class) #portion of doublets vs. singlet

## Normalisation -----------------------
 #data_list[[1]] <- NormalizeData(data_list[[1]], normalization.method = 'LogNormalize', scale.factor = 10000)
 data_list[[1]] <- NormalizeData(data_list[[1]])
 data_list[[2]] <- NormalizeData(data_list[[2]])
 str(data_list[[1]])

## Identify highly-variable features ===========================================
data_list[[1]] <- FindVariableFeatures(data_list[[1]], selection.method =  'vst', nfeatures = 2000)
data_list[[2]] <- FindVariableFeatures(data_list[[2]], selection.method =  'vst', nfeatures = 2000)

# Identify the top 10 HVGs
w1118_42d_top10 <- head(VariableFeatures(data_list[[1]]), 10)
w1118_42d_top10_plot <- VariableFeaturePlot(data_list[[1]])
LabelPoints(plot = w1118_42d_top10_plot, points = w1118_42d_top10, repel = TRUE)

drprnull_42d_top10 <- head(VariableFeatures(data_list[[2]]), 10)
drprnull_42d_top10_plot <- VariableFeaturePlot(data_list[[2]])
LabelPoints(plot = drprnull_42d_top10_plot, points = drprnull_42d_top10, repel = TRUE)

# Scaling ==================================================
all_genes1 <- rownames(data_list[[1]])
data_list[[1]] <- ScaleData(data_list[[1]], features = all_genes1)
View(data_list[[1]]@assays$RNA)

all_genes2 <- rownames(data_list[[2]])
data_list[[2]] <- ScaleData(data_list[[2]], features = all_genes2)
View(data_list[[2]]@assays$RNA)

## Set the corrected assay as default for further analysis
DefaultAssay(data_list[[1]]) <- "RNA"
DefaultAssay(data_list[[2]]) <- "RNA"

SaveH5Seurat(data_list$w1118_42D, filename = paste0(outputdir,"w1118_42d_QCnorm.h5Seurat"), overwrite = TRUE)
SaveH5Seurat(data_list$drprnull_42D, filename = paste0(outputdir,"drprnull_42d_QCnorm.h5Seurat"), overwrite = TRUE)

# Clustering ===================================================
### Run clustering with different resolutions
run_clustering <- function(sobj){
  sobj <- RunPCA(sobj, verbose = FALSE)
  sobj <- RunUMAP(sobj, dims = 1:30, verbose = FALSE)
  
  sobj <- FindNeighbors(sobj, dims = 1:30, verbose = FALSE)
  sobj <- FindClusters(sobj, verbose = FALSE, resolution = c(0.8, 1, 1.5, 2.5, 3.5, 4.5, 8, 12))
  
  return(sobj)
}

### Run clustering 
data_list <- lapply(data_list, run_clustering)

### Plot UMAP and marker gene dotplots for each correction/clustering resolution combination
norm_methods <- c("RNA")
cluster_res <- c("snn_res.0.8","snn_res.1","snn_res.1.5","snn_res.2.5","snn_res.3.5","snn_res.4.5",
                 "snn_res.8","snn_res.12")

markers <- c("elav","lncRNA:noe","VAChT","VGlut","Gad1","Vmat","SerT","Tdc2","ple", # neurons
             "repo","lncRNA:CR34335","alrm","wrapper","Indy","moody",#glia 
             "ninaC",	"trp",	"trpl", #photoreceptors
             "Hml", #hemocytes
             "ppl",#fatbody
             "drpr")

for(nm in norm_methods){
  for(cr in cluster_res){
    nmcr <- paste0(nm,"_",cr)
    p1 <- DimPlot(data_list[1]$w1118_42D,group.by = nmcr) &
      theme_classic(base_size = 16)
    
    ggsave(p1,filename=paste0(outputdir,"figures/","w1118_42d_",nmcr,"_dimplot.pdf"),width=10,height=10)
    
    p2 <- DimPlot(data_list[2]$drprnull_42D,group.by = nmcr) &
      theme_classic(base_size = 16)
    
    ggsave(p2,filename=paste0(outputdir,"figures/","drprnull_42d_",nmcr,"_dimplot.pdf"),width=10,height=10)
    
    p3 <- DotPlot(data_list[1]$w1118_42D, features = markers, group.by = nmcr) +
      scale_color_viridis_c(option = "magma",direction = -1) +
      theme_classic(base_size = 16) + RotatedAxis()
    
    ggsave(p3,filename=paste0(outputdir,"figures/","w1118_42d_",nmcr,"_dotplot.pdf"), width=16,height=10)
    
    p4 <- DotPlot(data_list[2]$drprnull_42D, features = markers, group.by = nmcr)&
      scale_color_viridis_c(option = "magma",direction = -1)&
      theme_classic(base_size = 16)+ RotatedAxis()
    
    ggsave(p4,filename=paste0(outputdir,"figures/","drprnull_42d_",nmcr,"_dotplot.pdf"), width=16,height=10)
    #break
  }
}

## Save seurat objects together
saveRDS(data_list, file=paste0(outputdir,"data_list.RDS"))

## Save seurat objects as h5ad files separately
data_list[[1]][["RNA3"]] <- as(object = data_list[[1]][["RNA"]], Class = "Assay")
DefaultAssay(data_list[[1]]) <- "RNA3"
data_list[[1]][["RNA"]] <- NULL
data_list[[1]] <- RenameAssays(object = data_list[[1]], RNA3 = 'RNA')
DefaultAssay(data_list[[1]]) <- "RNA"

data_list[[2]][["RNA3"]] <- as(object = data_list[[2]][["RNA"]], Class = "Assay")
DefaultAssay(data_list[[2]]) <- "RNA3"
data_list[[2]][["RNA"]] <- NULL
data_list[[2]] <- RenameAssays(object = data_list[[2]], RNA3 = 'RNA')
DefaultAssay(data_list[[2]]) <- "RNA"

SaveH5Seurat(data_list[[1]], filename = paste0(outputdir,"w1118_42d_QCnorm.h5Seurat"), overwrite = TRUE)
SaveH5Seurat(data_list[2]$drprnull_42D, filename = paste0(outputdir,"drprnull_42d_QCnorm.h5Seurat"), overwrite = TRUE)


## Convert h5Seurat to h5ad
Convert(paste0(outputdir,"w1118_42d_QCnorm.h5Seurat"), dest = paste0(outputdir,"w1118_42d_QCnorm.h5ad"), overwrite = TRUE)
Convert(paste0(outputdir,"drprnull_42d_QCnorm.h5Seurat"), dest = paste0(outputdir,"drprnull_42d_QCnorm.h5ad"), overwrite = TRUE)

sessionInfo()














