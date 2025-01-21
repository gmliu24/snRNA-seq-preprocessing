# DEG analysis based on preprocessGL_20250111
# Output Folder structure: 
# - /projectnb/mccall/guangmeiliu/snrnaseq_gmliu/analysis
#   - preprocessDate (same name as the script’s)
#     - DEG
#       - Different method (filtering parameters)
#         - GO
#           - GO_up
#           - GO_down

library(tidyverse) # dplyr and ggplot2
library(Seurat) # Seurat toolkit
library(dplyr) 
library(patchwork) # for plotting
library(scales)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
library(ggplotify)
set.seed(12345)

basedir <- "/projectnb/mccall/guangmeiliu/snrnaseq_gmliu/analysis/preprocess20250111/"
inputdir <- basedir

# 1. Fetch data ===================================================
data_harmony <- readRDS(paste0(basedir,'data_harmony.RDS'))

DefaultAssay(data_harmony) <- 'soupX_SCT'
Idents(data_harmony) <- 'celltype'

DimPlot(data_harmony, reduction = "umap", group.by = "celltype", label =TRUE, repel = F, #shuffle aviods ovelapping
        shuffle = T,  pt.size = 1, label.size = 5)

# 2. DEG analysis in specific cell type ===================================================
outputdir <- paste0(basedir,"DEG/soupX_SCT_data/") #specify method in the outputdir
dir.create(outputdir,showWarnings = FALSE)

glialcell_degs <- FindMarkers(subset(data_harmony, celltype == 'glial cell'), logfc.threshold = 0.25,
                           only.pos = F,
                           test.use = 'wilcox',
                           slot = 'data',
                           ident.1 = "drprnull_42D", ident.2 = "w1118_42D",
                           group.by = 'orig.ident'
                           ) %>%
  filter(p_val_adj < 0.05) %>%
  mutate(gene = rownames(.)) #add gene names to another column for convenient downstream analysis
##Add a column to specify up- or down- regulation
up_down_spec <- function(data){
  data$significance <- "Non-significant"
  data$significance[data$avg_log2FC > 1 & data$p_val < 0.05] <- 'Up-regulated'
  data$significance[data$avg_log2FC < -1 & data$p_val < 0.05] <- 'Down-regulated'
  return(data)
}
glialcell_degs <- up_down_spec(glialcell_degs)

write.csv(glialcell_degs, paste0(outputdir, 'glial_cell_DEG.csv'), row.names = T)

##volcano plot

top10_labels <- glialcell_degs %>%
  filter(significance != "Non-significant") %>%
  group_by(significance) %>%
  slice_max(order_by = abs(avg_log2FC), n = 10)

p1 <- ggplot(data = glialcell_degs, aes(x = avg_log2FC, y = -log10(p_val_adj), col = significance)) +
  geom_point(size = 3) +
  geom_text_repel(data = top10_labels, 
             aes(label = gene), 
             max.overlaps = Inf,
             point.padding = 0.3,
             size = 5, show.legend = F) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  ggtitle ('DEGs of 42D drprnull vs w1118') +
  scale_color_manual(values = c("blue", "grey", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  theme_classic(base_size = 16) +
  labs(
    x='Log2 Fold Change',
    y = '-Log10 Adjusted P value',
    tag = 'Glial cell' #change according to tissue selection
  ) +
  theme(
    legend.key = element_blank(), # Simplify the legend
    legend.title = element_blank(), # Remove the legend title
    legend.position = 'top',
    legend.text = element_text(size = 14, face = 'bold'),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.line = element_line(linewidth = 1.5),   # Change axis line thickness
    axis.ticks = element_line(linewidth = 1.5),
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.1, 0.99)
  ) 

ggsave(p1, filename = paste0(outputdir,'DEGs_glialcell.png'), width = 12, height = 10)

# 3. Function to run DEG on each broad cluster ===================================================
outputdir <- paste0(basedir,"DEG/soupX_SCT_data2/") #specify method in the outputdir
dir.create(outputdir,showWarnings = FALSE)
figdir <- paste0(outputdir,"figures/")
dir.create(figdir,showWarnings = FALSE)

broadcell <- unique(data_harmony@meta.data$celltype)
for (cell in broadcell) {
  # Subset the data for the current cell type
  current_cell_data <- subset(data_harmony, celltype == cell)
  
  # Find differentially expressed genes
  degs <- FindMarkers(current_cell_data, 
                      logfc.threshold = 0.25,
                      only.pos = F,
                      test.use = 'wilcox',
                      slot = 'data',
                      ident.1 = "drprnull_42D", 
                      ident.2 = "w1118_42D",
                      group.by = 'orig.ident'
  ) %>%
    filter(-log10(p_val_adj) > 0) %>%
    mutate(gene = rownames(.)) # Add gene names
  
  # Add up-/down-regulation significance
  up_down_spec <- function(data){
    data$significance <- "Non-significant"
    data$significance[data$avg_log2FC > 1 & data$p_val_adj < 0.05] <- 'Up-regulated'
    data$significance[data$avg_log2FC < -1 & data$p_val_adj < 0.05] <- 'Down-regulated'
    return(data)
  }
  degs <- up_down_spec(degs)
 
  # Write the DEG results to a CSV file
  write.csv(degs, paste0(outputdir, cell, "_DEG.csv"), row.names = T)
  
  # Make a histogram to see the distribution of avg_log2FC
  p1 <- ggplot(subset(degs, significance == "Up-regulated"), aes(x=avg_log2FC)) + geom_histogram() +
    geom_density(alpha=0.2, fill="lightpink") +
    ggtitle(paste("Log2FC distribution of upgenes in",cell))
  
  p2 <- ggplot(subset(degs, significance == "Down-regulated"), aes(x=avg_log2FC)) + geom_histogram() +
    geom_density(alpha=0.2, fill="lightblue") +
    ggtitle(paste("Log2FC distribution of downgenes in",cell))
  
  ggsave(p1, filename = paste0(figdir, cell, "_up_hist.png"), width = 12, height = 10)
  ggsave(p2, filename = paste0(figdir, cell, "_down_hist.png"), width = 12, height = 10)
  # Categorize up- and down-regulated genes as lowly, moderately and higly 
  
  
  # Create a volcano plot
  top20_labels <- degs %>%
    filter(significance != "Non-significant") %>%
    group_by(significance) %>%
    slice_max(order_by = abs(avg_log2FC), n = 20)
  
  p <- ggplot(data = degs, aes(x = avg_log2FC, y = -log10(p_val_adj), 
                               col = significance)) +
    geom_point(size = 3) +
    geom_text_repel(data = top20_labels, 
                    aes(label = gene), 
                    max.overlaps = Inf,
                    point.padding = 0.3,
                    size = 5, show.legend = F) +
    geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    ggtitle(expression("DEGs of 42D " * italic("drprnull") * " vs " * italic('w1118'))) +
    scale_color_manual(values = c("#0072B2", "grey", "#CC0033"), 
                       labels = c("Downregulated", "Not significant", "Upregulated")) +
    theme_classic(base_size = 16) +
    labs(
      x = 'Log2 Fold Change',
      y = '-Log10 Adjusted P value',
      tag = cell
    ) +
    theme(
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.position = 'top',
      legend.text = element_text(size = 14, face = 'bold'),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.line = element_line(linewidth = 1.5),
      axis.ticks = element_line(linewidth = 1.5),
      plot.tag = element_text(size = 18, face = "bold"),
      plot.tag.position = c(0.1, 0.99)
    )
  
  # Save the volcano plot
  ggsave(p, filename = paste0(figdir, cell, "_DEG.png"), width = 12, height = 10)
}

## Merge all .csv files into one file ===================================================
file_list <- list.files(path = outputdir, pattern = "\\.csv$", full.names = T)
data_list = list()

for (file in file_list) {
  # Extract the filename without extension
  filename <- basename(file)
  filename_no_ext <- tools::file_path_sans_ext(filename)
  
  # Extract the string before "_" as the cell type
  celltype <- strsplit(filename_no_ext, "_")[[1]][1]
  
  # Read the CSV file
  data <- read.csv(file, stringsAsFactors = FALSE)
  
  # Add the celltype column
  data$celltype <- celltype
  
  # Append to the list
  data_list[[length(data_list) + 1]] <- data
}

# Merge all data frames in the list into one
all.degs <- bind_rows(data_list)

# Write the merged data to a new CSV file
write.csv(all.degs, file = paste0(outputdir, "All_celltypes_DEGs.csv"), row.names = FALSE)

## Organize genes into different sheets in xlxs file ===================================================
library(openxlsx)
data <- read_csv(paste0(outputdir, "All_celltypes_DEGs.csv"))
broadcell <- unique(data$celltype)

output_file <- paste0(outputdir, "Up_and_Down_Genes_by_celltype.xlsx")
wb <- createWorkbook()

for (cell in broadcell) {
  # Subset and order the data for upregulated genes
  upgene <- subset(data, celltype == cell & significance == "Up-regulated")
  upgene <- upgene[order(upgene$avg_log2FC, decreasing = TRUE), ]
  
  # Subset and order the data for downregulated genes
  downgene <- subset(data, celltype == cell & significance == "Down-regulated")
  downgene <- downgene[order(downgene$avg_log2FC, decreasing = T), ]
  # Add upregulated genes to a sheet
  if (nrow(upgene) > 0) {
    addWorksheet(wb, paste0(cell, "_up"))
    writeData(wb, sheet = paste0(cell, "_up"), upgene)
  }
  
  # Add downregulated genes to a sheet
  if (nrow(downgene) > 0) {
    addWorksheet(wb, paste0(cell, "_down"))
    writeData(wb, sheet = paste0(cell, "_down"), downgene)
  }
}
# Save the workbook
saveWorkbook(wb, file = output_file, overwrite = TRUE)

## Create a histogram to see the distribution of log2FC across all up- and down- regulated genes ===================================================
## fig size and font size are optimal for publication
p3 <- ggplot(subset(data, significance != "Non-significant"), aes(x=avg_log2FC, fill = significance, color = significance)) + geom_histogram() +
  geom_density(alpha=0.2, fill="grey") +
  ggtitle(paste("Log2FC distribution of all genes")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "top")
ggsave(p3, filename = paste0(figdir, "All_genes_hist.png"), width = 8, height = 8)

p4 <- ggplot(subset(data, significance == "Up-regulated"), aes(x=avg_log2FC)) + geom_histogram(aes(y=after_stat(density)), fill = "lightblue", color = "black") +
  geom_density(alpha=0.2, fill = "grey") +
  ggtitle(paste("Log2FC distribution of up-regulated genes")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank())
ggsave(p4, filename = paste0(figdir, "Upgenes_hist.png"), width = 8, height = 8)

p5 = ggplot(subset(data, significance == "Down-regulated"), aes(x=avg_log2FC)) + geom_histogram(aes(y=after_stat(density)), fill = "lightpink", color = "black") +
  geom_density(alpha=0.2, fill = "grey") +
  ggtitle(paste("Log2FC distribution of down-regulated genes")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank())
ggsave(p5, filename = paste0(figdir, "Downgenes_hist.png"), width = 8, height = 8)


# 4. GO analysis ===================================================
##https://www.melbournebioinformatics.org.au/tutorials/tutorials/seurat-go/seurat-go/
##Convert gene symbols to Entrez IDs
library(org.Dm.eg.db)
library(enrichplot)
library(clusterProfiler)

basedir <- "/projectnb/mccall/guangmeiliu/snrnaseq_gmliu/analysis/preprocess20250111/"
## GO analysis for upregulated genes
outputdir <- paste0(basedir,"DEG/soupX_SCT_data2/GO/GO_up/")
dir.create(outputdir,showWarnings = FALSE)

for (cell in broadcell){
  up.degs <- subset(all.degs, celltype == cell & significance == 'Up-regulated')
  cell.gene_ids <- bitr(up.degs$gene, fromType = "SYMBOL",
                        toType = "ENTREZID", OrgDb = org.Dm.eg.db)
  
  ##Perform GO enrichment analysis for Cluster 2
  cell.ego <- enrichGO(gene = cell.gene_ids$ENTREZID, 
                       OrgDb = org.Dm.eg.db, 
                       ont = "ALL", # biological process
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05, 
                       readable = TRUE)
  
  # Skip if enrichment result is empty
  if (is.null(cell.ego) || nrow(cell.ego) == 0) {
    message(paste("No significant GO terms for", cell, "- skipping."))
    next
  }
  ##Visualise the GO enrichment results
  upbar <- barplot(cell.ego, showCategory=10) + 
    aes(fill = ONTOLOGY) +
    scale_fill_manual(values = c("BP" = "dodgerblue1", "CC" = "orange", "MF" = "green3"),
                      labels = c("BP" = "Biological Process", "CC" = "Cellular Component", "MF" = "Molecular Function")) +
    ggtitle(str_wrap(paste("GO Term Analysis for", cell), width = 30)) +
    theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 25, face = "bold"),
      axis.title = element_text(size = 20, face = "bold")
    )
  
  ggsave(upbar, filename = paste0(outputdir, cell, "_up_GO_barplot.png"), width = 12, height = 10)
  
  updot <- dotplot(cell.ego, showCategory=10) +
    aes(color = ONTOLOGY) +
    scale_color_manual(values = c("BP" = "black", "CC" = "blue", "MF" = "green"),
                       labels = c("BP" = "Biological Process", "CC" = "Cellular Component", "MF" = "Molecular Function")) +
    guides(color = guide_legend(override.aes = list(size = 5, stroke = 1.5))) +
    theme(
      axis.text.y = element_text(size = 20, face = "bold"),
      axis.text.x = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 20, face = "bold"),
      plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
      axis.line = element_line(linewidth = 1.5),
      axis.ticks = element_line(linewidth = 1.5)
    ) + 
    ggtitle(str_wrap(paste("GO Term Analysis for", cell), width = 30))
  ggsave(updot, filename = paste0(outputdir, cell, "_up_GO_dotplot.png"), width = 12, height = 10)
  
}

#broadcell = "epithelial cell"
## GO analysis for down-regulated genes
outputdir <- paste0(basedir,"GO_down/")
#dir.create(outputdir,showWarnings = FALSE)
for (cell in broadcell){
  down.degs <- subset(all.degs, celltype == cell & significance == 'Down-regulated')
  cell.gene_ids <- bitr(down.degs$gene, fromType = "SYMBOL",
                        toType = "ENTREZID", OrgDb = org.Dm.eg.db)
  
  ##Perform GO enrichment analysis for Cluster 2
  cell.ego <- enrichGO(gene = cell.gene_ids$ENTREZID, 
                       OrgDb = org.Dm.eg.db, 
                       ont = "ALL", # biological process
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05, 
                       readable = TRUE)
  
  # Skip if enrichment result is empty
  if (is.null(cell.ego) || nrow(cell.ego) == 0) {
    message(paste("No significant GO terms for", cell, "- skipping."))
    next
  }
  ##Visualise the GO enrichment results
  downbar <- barplot(cell.ego, showCategory=10) + 
    aes(fill = ONTOLOGY) +
    scale_fill_manual(values = c("BP" = "dodgerblue1", "CC" = "orange", "MF" = "green3"),
                      labels = c("BP" = "Biological Process", "CC" = "Cellular Component", "MF" = "Molecular Function")) +
    ggtitle(str_wrap(paste("GO Term Analysis for", cell), width = 35)) +
    theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 25, face = "bold"),
      axis.title = element_text(size = 20, face = "bold")
    )
  
  ggsave(downbar, filename = paste0(outputdir, cell, "_down_GO_barplot.png"), width = 12, height = 10)
  
  downdot <- dotplot(cell.ego, showCategory=10) +
    aes(color = ONTOLOGY) +
    scale_color_manual(values = c("BP" = "black", "CC" = "blue", "MF" = "green"),
                      labels = c("BP" = "Biological Process", "CC" = "Cellular Component", "MF" = "Molecular Function")) +
    guides(color = guide_legend(override.aes = list(size = 5, stroke = 1.5))) +
    theme(
      axis.text.y = element_text(size = 20, face = "bold"),
      axis.text.x = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 20, face = "bold"),
      plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
      axis.line = element_line(linewidth = 1.5),
      axis.ticks = element_line(linewidth = 1.5)
    ) + 
    ggtitle(str_wrap(paste("GO Term Analysis for", cell), width = 35))
  ggsave(downdot, filename = paste0(outputdir, cell, "_down_GO_dotplot.png"), width = 12, height = 10)
  
}

# supp. Add meta.data column using join function ===================================================
###Load glialcell annotation
glialcell_metadata <- glialcell@meta.data
###Merge subcluster in glialcell with data_harmony
###Fetch celltype info
data_harmony_metadata <- FetchData(data_harmony, 'celltype')
data_harmony_metadata$cell_id <- rownames(data_harmony_metadata)
###merge using left_join
data_harmony_metadata <- left_join(x=data_harmony_metadata, y = glialcell_metadata, by='celltype')
### Readd rownames
rownames(data_harmony_metadata) <- data_harmony_metadata$cell_id
view(data_harmony_metadata)
###Add metadata
data_harmony <- AddMetaData(data_harmony, metadata = data_harmony_metadata)
table(data_harmony@meta.data$celltype)






