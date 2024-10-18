##Load packages
library(Seurat)
library(patchwork)

#Load data (expression matrix was downloaded from GSE132465)
expr_matrix <- read.table("GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt", header = T)
rownames(expr_matrix) <- expr_matrix[,1]
expr_matrix <- expr_matrix[,-1]
GSE132465 <- CreateSeuratObject(expr_matrix, assay = "RNA", project = "GSE132465")

rm(expr_matrix)
# Load the Seurat object
seurat_obj <- GSE132465

# 1. Quality control - Filter low-quality cells and features
# Set minimum features (genes) and cells as per your data
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA > 500)

# 2. Normalize the data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# 3. Identify highly variable features (genes)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# 4. Scale the data
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

# 5. (Optional) Run PCA or further analysis
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Visualize PCA (optional)
ElbowPlot(seurat_obj)

# 2. Run t-SNE
# Set the number of dimensions based on PCA (10-20 is common)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:15)  # Change dims based on your elbow plot

# t-SNE Plot
DimPlot(seurat_obj, reduction = "tsne", group.by = "ident")  # 'ident' is by default cell type or cluster identities

# 3. Run UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)  # Adjust dims if needed

# UMAP Plot
DimPlot(seurat_obj, reduction = "umap", group.by = "ident")

# Optionally, save the t-SNE and UMAP results
saveRDS(seurat_obj, file = "GSE132465_tsne_umap.rds")

seurat_obj <- FindNeighbors(seurat_obj)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)

celltype.marker <- read.table("~/Cellmarker.txt", sep = "\t",header = T)
p <- VlnPlot(seurat_obj, features = celltype.marker$marker, stack = T, group.by = "seurat_clusters")+NoLegend()
ggsave(p, filename = "Feature_Marker_VlnPlot.png", width = 20, height = 8)



#Use cell markers for cell type identification. 
#T cells: 0, 1, 2, 16
#B cells: 3, 4, 10
#M1 macrophage: 5, 13
#M2 macrophage: 7
#Epithelial cells: 6, 8, 14
#Fibroblasts: 9, 11, 15
#Endothelial cells: 12

library(Seurat)

# Ensure that the active identities in the Seurat object correspond to clusters
# If clusters are not set as the active identity, set it
Idents(seurat_obj) <- "seurat_clusters"  # This assumes your clusters are stored in 'seurat_clusters'

# Create a mapping of clusters to cell types
cell_type_mapping <- c(
  "0" = "T cells", "1" = "T cells", "2" = "T cells", "16" = "T cells",
  "3" = "B cells", "4" = "B cells", "10" = "B cells",
  "5" = "M1 macrophage", "13" = "M1 macrophage",
  "7" = "M2 macrophage",
  "6" = "Epithelial cells", "8" = "Epithelial cells", "14" = "Epithelial cells",
  "9" = "Fibroblasts", "11" = "Fibroblasts", "15" = "Fibroblasts",
  "12" = "Endothelial cells"
)

# Add the cell type annotation to the metadata
# Ensure the clusters are converted to character to match the mapping
seurat_obj$cell_type <- plyr::mapvalues(x = as.character(Idents(seurat_obj)),
                                        from = names(cell_type_mapping),
                                        to = cell_type_mapping)

# Verify the assignment
table(seurat_obj$cell_type)

# t-SNE plot with cell type annotations
DimPlot(seurat_obj, reduction = "tsne", group.by = "cell_type", label = TRUE) + NoLegend()

# UMAP plot with cell type annotations
DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", label = TRUE) + NoLegend()
VlnPlot(seurat_obj, features = celltype.marker$marker, stack = T, group.by = "cell_type")+NoLegend()
marker <- c("CD8A", "CD3E", "CD3D", "VSIG4", "MRC1", "CD163", "PTGS2", "CD86", "COL1A1", "ACTA2", "EPCAM", "CEACAM5", "KRT19","VWF", "PECAM1", "ENG", "CD79A", "IGHA1", "MZB1")
VlnPlot(seurat_obj, features = marker, stack = T, group.by = "cell_type")+NoLegend()
Idents(seurat_obj) = seurat_obj$cell_type
FeaturePlot(seurat_obj, features = "DDR1", reduction = "tsne", label = T)
VlnPlot(seurat_obj, features = c("DDR1", "CXCL9", "CXCL10"), stack = T, group.by = "cell_type", pt.size = 1)+NoLegend()

#Coexpression
# Create a new metadata column with summed expression of DDR1, CXCL9, and CXCL10
Epithelial <- subset(seurat_obj, idents = "Epithelial cells")
FeaturePlot(Epithelial, features = c("DDR1", "CXCL9, CXCL10"), slot = "count", reduction = "tsne")

rm(expr_matrix)


# Remove prefixes and keep only "Normal" or "Tumor"
seurat_obj$tissue_type <- gsub(".*\\.N", "Normal", seurat_obj$orig.ident)  # Replace any prefix ending in .N with "Normal"
seurat_obj$tissue_type <- gsub(".*\\.T", "Tumor", seurat_obj$tissue_type)   # Replace any prefix ending in .T with "Tumor"

# Check the updated names
unique(seurat_obj$tissue_type)

Idents(Epithelial) = Epithelial$tissue_type
VlnPlot(Epithelial, "DDR1")+NoLegend()

#DDR1 GSEA
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Fetch DDR1 expression directly from the RNA assay (assuming default is "RNA")
DDR1_expr <- GetAssayData(Epithelial, assay = "RNA", layer = "data")["DDR1", ]

# Check if the expression data was retrieved correctly
head(DDR1_expr)

# Calculate the 75th percentile (top 25% cells)
quantile_DDR1 <- quantile(DDR1_expr, probs = 0.9)

# Classify cells into DDR1_H (high) and DDR1_L (low)
Epithelial$DDR1_group <- ifelse(DDR1_expr >= quantile_DDR1, "DDR1_H", "DDR1_L")

# Check the distribution of cells in each group
table(Epithelial$DDR1_group)

Idents(Epithelial) <- Epithelial$DDR1_group

# Identify differentially expressed genes (DEGs) between DDR1_H and DDR1_L groups
DEGs <- FindMarkers(Epithelial, ident.1 = "DDR1_H", ident.2 = "DDR1_L", logfc.threshold = 0.25)

# Convert gene symbols to Entrez IDs
DEGs$entrez <- mapIds(org.Hs.eg.db, keys = rownames(DEGs), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# Remove NAs
DEGs <- na.omit(DEGs)

# Rank genes based on log fold change
gene_list <- DEGs$avg_log2FC
names(gene_list) <- DEGs$entrez
gene_list <- sort(gene_list, decreasing = TRUE)

# Perform GSEA on KEGG pathways
gsea_results <- gseKEGG(geneList = gene_list, organism = "hsa", minGSSize = 10, pvalueCutoff = 0.05, verbose = FALSE)

# View top pathways
head(gsea_results@result)

# Select the top 10 enriched pathways
histone_pathways <- gsea_results@result[grep("istone", gsea_results@result$Description, ignore.case = TRUE), ]

# Create a bubble plot
ggplot(histone_pathways, aes(x = reorder(Description, NES), y = -log10(p.adjust))) +
  geom_point(aes(size = setSize, color = NES)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Top 10 Enriched KEGG Pathways in DDR1_H Group", x = "KEGG Pathway", y = "-log10(p.adjust)") +
  coord_flip()

# Get the MSigDB gene sets related to histone modification from "C2" category
msig_histone_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

# Filter for gene sets related to histone modifications (e.g., "HISTONE")
histone_sets <- msig_histone_sets[grep("HISTONE", msig_histone_sets$gs_name, ignore.case = TRUE), ]

# View the histone-related gene sets
unique(histone_sets$gs_name)

# Perform differential expression analysis for DDR1_H vs DDR1_L (or your specific groups)
DEGs <- FindMarkers(Epithelial, ident.1 = "DDR1_H", ident.2 = "DDR1_L", logfc.threshold = 0.25)

# Convert gene symbols to Entrez IDs
DEGs$entrez <- mapIds(org.Hs.eg.db, keys = rownames(DEGs), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# Remove NAs and rank genes by log fold change (logFC)
DEGs <- na.omit(DEGs)
gene_list <- DEGs$avg_log2FC
names(gene_list) <- DEGs$entrez
gene_list <- sort(gene_list, decreasing = TRUE)


# Perform GSEA on the histone-related gene sets
gsea_histone <- GSEA(geneList = gene_list, TERM2GENE = histone_sets[, c("gs_name", "entrez_gene")], pvalueCutoff = 0.1)

# View the GSEA results
head(gsea_histone@result)

ggplot(gsea_histone, aes(x = reorder(Description, NES), y = -log10(p.adjust))) +
  geom_point(aes(size = setSize, color = NES)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Histone Modification", 
       x = "KEGG Pathway", 
       y = "p.adj") +
  coord_flip()
