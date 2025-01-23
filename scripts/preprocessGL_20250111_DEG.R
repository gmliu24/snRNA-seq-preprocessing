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

# Calculate the 0.33 and 0.66 percentiles
Up_per <- quantile(subset(data, significance == "Up-regulated")$avg_log2FC, probs = c(0.33, 0.66))
Down_per <- quantile(subset(data, significance == "Down-regulated")$avg_log2FC, probs = c(0.33, 0.66))
Up_per
Down_per

p4 <- ggplot(subset(data, significance == "Up-regulated"), aes(x=avg_log2FC)) + geom_histogram(aes(y=after_stat(density)), fill = "lightblue", color = "black") +
  geom_density(alpha=0.2, fill = "grey") +
  ggtitle(paste("Log2FC distribution of up-regulated genes")) +
  geom_vline(xintercept = Up_per[1], color = "grey", linetype = "dashed", size = 1) + # 33rd percentile
  geom_vline(xintercept = Up_per[2], color = "grey", linetype = "dashed", size = 1) + # 66th percentile
  annotate("text", x = Up_per[1], y = 0.5, label = "33%", color = "black", size = 5) +
  annotate("text", x = Up_per[2], y = 0.5, label = "66%", color = "black", size = 5) +
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
  geom_vline(xintercept = Down_per[1], color = "grey", linetype = "dashed", size = 1) + # 33rd percentile
  geom_vline(xintercept = Down_per[2], color = "grey", linetype = "dashed", size = 1) + # 66th percentile
  annotate("text", x = Down_per[1], y = 0.7, label = "33%", color = "black", size = 5) +
  annotate("text", x = Down_per[2], y = 0.7, label = "66%", color = "black", size = 5) +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank())
ggsave(p5, filename = paste0(figdir, "Downgenes_hist.png"), width = 8, height = 8)

## Add category for lowly/moderately/highly regulated genes ===================================================
## Calculate percentiles

basedir <- "/projectnb/mccall/guangmeiliu/snrnaseq_gmliu/analysis/preprocess20250111/"
outputdir <- paste0(basedir,"DEG/soupX_SCT_data2/") #specify method in the outputdir
data <- read_csv(paste0(outputdir, "All_celltypes_DEGs.csv"))

Up_per <- quantile(subset(data, significance == "Up-regulated")$avg_log2FC, probs = c(0.33, 0.66))
Down_per <- quantile(subset(data, significance == "Down-regulated")$avg_log2FC, probs = c(0.33, 0.66))
Up_per
Down_per

data <- data %>%
  mutate(
    category = case_when(
      significance == "Up-regulated" & avg_log2FC >= 1 & avg_log2FC < Up_per[1] ~ "low",
      significance == "Up-regulated" & avg_log2FC >= Up_per[1] & avg_log2FC < Up_per[2] ~ "mid",
      significance == "Up-regulated" & avg_log2FC >= Up_per[2] ~ "high",
      significance == "Down-regulated" & avg_log2FC < Down_per[1] ~ "high",
      significance == "Down-regulated" & avg_log2FC >= Down_per[1] & avg_log2FC < Down_per[2] ~ "mid",
      significance == "Down-regulated" & avg_log2FC >= Down_per[2] ~ "low",
      TRUE ~ NA_character_
    )
  ) 

write.csv(data, file = paste0(outputdir, "All_celltypes_DEGs.csv"), row.names = FALSE)

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

# 4. GO analysis 1===================================================
##https://www.melbournebioinformatics.org.au/tutorials/tutorials/seurat-go/seurat-go/
##Convert gene symbols to Entrez IDs
library(org.Dm.eg.db)
library(enrichplot)
library(clusterProfiler)

basedir <- "/projectnb/mccall/guangmeiliu/snrnaseq_gmliu/analysis/preprocess20250111/"
## GO analysis for upregulated genes
outputdir <- paste0(basedir,"DEG/soupX_SCT_data2/GO/GO_up/")
dir.create(outputdir,showWarnings = FALSE)

 ##Map gene id ============================
##Add EntrezID column in the data frame
cell.gene_ids <- bitr(data$gene, fromType = "SYMBOL",
                      toType = "ENTREZID", OrgDb = org.Dm.eg.db)

data <- merge(data, cell.gene_ids, 
              by.x = "gene",        # Match based on gene column in all.degs
              by.y = "SYMBOL",      # Match based on SYMBOL column in cell.gene_ids
              all.x = TRUE)         # Keep all rows from all.degs, even if no match is found

data$Entrez <- data$ENTREZID

##If some genes fail to map, look them up in PANGEA and export a csv file
data2 <- read_csv(paste0(outputdir,"missingids.csv"))

merged_data <- data %>%
  left_join(data2, by = c("gene" = "Search Term"))

# Update the "entrez" column with the "ID" from data2
merged_data$Entrez <- ifelse(!is.na(merged_data$`Mapped to Entrez Gene`), merged_data$`Mapped to Entrez Gene`, merged_data$Entrez)

merged_data <- merged_data %>% 
  select(-`Mapped to Entrez Gene`, -`Gene Type`, -`Flybase ID`, -Symbol, - Messages, -Duplicates)

write.csv(merged_data, paste0(outputdir,"All_celltypes_DEGs.csv"))

data <- merged_data

##Perform GO term analysis =========
all.degs <- data
broadcell <- unique(data$celltype)
for (cell in broadcell){
  up.degs <- subset(all.degs, celltype == cell & significance == 'Up-regulated')
  
  ##Perform GO enrichment analysis for Cluster 2
  cell.ego <- enrichGO(gene = up.degs$Entrez, 
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
outputdir <- paste0(basedir,"DEG/soupX_SCT_data2/GO/GO_down/")
dir.create(outputdir,showWarnings = FALSE)
#dir.create(outputdir,showWarnings = FALSE)

for (cell in broadcell){
  down.degs <- subset(all.degs, celltype == cell & significance == 'Down-regulated')
  
  ##Perform GO enrichment analysis for Cluster 2
  cell.ego <- enrichGO(gene = down.degs$Entrez, 
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


# 5. GO analysis 2 ===================================================
## https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html

## Show GO term for both up- and down-regulated genes in the same plot for each celltype ===========================
outputdir <- paste0(basedir,"DEG/soupX_SCT_data2/GO/GO_updown/")
dir.create(outputdir,showWarnings = FALSE)

all.degs <- subset(data, significance != "Non-significant")

##Define dotplot theme
theme_dotplot <- theme(
  axis.text.y = element_text(size = 18),
  axis.text.x = element_text(size = 18),
  axis.title.x = element_text(size = 18),
  plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
)

for (cell in broadcell){
  message(paste("Processing cell type:", cell))
  
  cell_degs <- subset(all.degs, celltype == cell)
  bp_go_updown <- compareCluster(Entrez~significance, data=cell_degs, fun="enrichGO", OrgDb='org.Dm.eg.db',ont="BP")
  cc_go_updown <- compareCluster(Entrez~significance, data=cell_degs, fun="enrichGO", OrgDb='org.Dm.eg.db',ont="CC")
  mf_go_updown <- compareCluster(Entrez~significance, data=cell_degs, fun="enrichGO", OrgDb='org.Dm.eg.db',ont="MF")

  results_list <- list()
  
  if (!is.null(bp_go_updown) && nrow(bp_go_updown) > 0) {
    bp_go_updown@compareClusterResult$GO <- "BP"
    results_list[["BP"]] <- bp_go_updown@compareClusterResult
  }
  if (!is.null(cc_go_updown) && nrow(cc_go_updown) > 0) {
    cc_go_updown@compareClusterResult$GO <- "CC"
    results_list[["CC"]] <- cc_go_updown@compareClusterResult
  }
  if (!is.null(mf_go_updown) && nrow(mf_go_updown) > 0) {
    mf_go_updown@compareClusterResult$GO <- "MF"
    results_list[["MF"]] <- mf_go_updown@compareClusterResult
  }
  
  # Merge all non-NULL results
  if (length(results_list) > 0) {
    merged_go_updown <- do.call(rbind, results_list)
    # Save the merged results to a CSV file
    write.csv(merged_go_updown, file = paste0(outputdir, cell, "_GO_results.csv"), row.names = FALSE)
  } else {
    message(paste("No GO enrichment results for cell type:", cell))
  }

  if (!is.null(bp_go_updown) && nrow(bp_go_updown) > 0) {
    bp <- dotplot(bp_go_updown) + 
      ggtitle(str_wrap(paste("Biological process enrichment in", cell), width = 30)) +
      scale_y_discrete(guide = guide_axis(check.overlap = TRUE)) +
      scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
      theme_dotplot
    ggsave(bp, filename = paste0(outputdir, cell, "_BP_dotplot.png"), width = 14, height = 10)
  } else {
    message(paste("No BP enrichment found for", cell))
  }
  
  if (!is.null(cc_go_updown) && nrow(cc_go_updown) > 0) {
    cc <- dotplot(cc_go_updown) + 
      ggtitle(str_wrap(paste("Cellular component enrichment in", cell), width = 30)) +
      scale_y_discrete(guide = guide_axis(check.overlap = TRUE)) +
      scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
      theme_dotplot
    ggsave(cc, filename = paste0(outputdir, cell, "_CC_dotplot.png"), width = 14, height = 10)
  } else {
    message(paste("No CC enrichment found for", cell))
  }
  
  if (!is.null(mf_go_updown) && nrow(mf_go_updown) > 0) {
    mf <- dotplot(mf_go_updown) + 
      ggtitle(str_wrap(paste("Molecular function enrichment in", cell), width = 30)) +
      scale_y_discrete(guide = guide_axis(check.overlap = TRUE)) +
      scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
      theme_dotplot
    ggsave(mf, filename = paste0(outputdir, cell, "_MF_dotplot.png"), width = 14, height = 10)
  } else {
    message(paste("No MF enrichment found for", cell))
  }
  
}

kegg_updown <- compareCluster(Entrez~significance+celltype, data=cell.gene_ids, fun="enrichKEGG")







## GO analysis for genes in each category and generate figure for multiple groups ===========================

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






