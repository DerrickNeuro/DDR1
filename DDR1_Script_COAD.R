# Load necessary libraries
library(limma)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)

#TCGA-COAD RNA-seq datasets were downloaded from UCSC Xena.. And the expression matrix was loaded in this environment named "COAD". 

# Extract DDR1 expression
DDR1_expr <- COAD["DDR1",]

# Create a group based on high/low DDR1 expression
cutoff <- median(DDR1_expr)
group <- ifelse(DDR1_expr > cutoff, "High", "Low")

# Make sure 'group' is a factor
group <- as.factor(group)

# Check if group aligns with the samples
length(group)  # This should be 499 to match the number of samples

# Differential expression analysis using limma
# Create a design matrix
design <- model.matrix(~ group)


# Fit the linear model
fit <- lmFit(COAD, design)
fit <- eBayes(fit)
deg_results <- topTable(fit, coef=2, number=Inf)

# Volcano plot of DEGs
deg_results$logP <- -log10(deg_results$P.Value)
deg_results$threshold <- as.factor(deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1)
volcanoplot(fit, coef=2, highlight=100, names=rownames(deg_results))

# GSVA analysis with KEGG pathways
# Load KEGG gene sets
kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
kegg_gene_sets <- split(x = kegg$gene_symbol, f = kegg$gs_name)

# Ensure that COAD has rownames (gene names)
# Assuming the gene names are in the first column or are already rownames in your data frame
if (is.null(rownames(COAD))) {
  # If gene names are in the first column, set them as rownames
  rownames(COAD) <- COAD[, 1]  # Adjust this if your data is structured differently
  COAD <- COAD[, -1]  # Remove the gene names column from the data frame
}

# Convert COAD to a numeric matrix (again ensuring it's properly formatted)
COAD_matrix <- as.matrix(COAD)
COAD_matrix <- apply(COAD_matrix, 2, as.numeric)

# Remove genes with constant expression (if needed)
gene_variances <- apply(COAD_matrix, 1, var)
COAD_matrix_filtered <- COAD_matrix[gene_variances > 0, ]

# Ensure rownames are properly set in COAD_matrix_filtered
rownames(COAD_matrix_filtered) <- rownames(COAD)[gene_variances > 0]

# Perform GSVA analysis
gsva_results <- gsva(COAD_matrix_filtered, kegg_gene_sets, method="ssgsea", kcdf="Gaussian", parallel.sz=1)

# Create design matrix for high/low DDR1 expression groups
gsva_design <- model.matrix(~ group)
gsva_fit <- lmFit(gsva_results, gsva_design)
gsva_fit <- eBayes(gsva_fit)
gsva_deg <- topTable(gsva_fit, coef=2, number=Inf)

# Filter significant GSVA results
sig_gsva <- gsva_deg[gsva_deg$adj.P.Val < 0.05, ]

# Export heatmap of significant GSVA results
heatmap_data <- gsva_results[rownames(sig_gsva), ]
pheatmap(heatmap_data, cluster_rows=TRUE, cluster_cols=TRUE)

# Optionally save the heatmap
pheatmap(heatmap_data, cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=data.frame(group),
         filename="GSVA_significant_heatmap.png")


# Extract the expression of cytokines across samples. 
DDR1_expr <- as.numeric(COAD["DDR1", ])
HMGB1_expr <- as.numeric(COAD["HMGB1", ])
IDO1_expr <- as.numeric(COAD["IDO1", ])
IFNG_expr <- as.numeric(COAD["IFNG", ])
IFNGR1_expr <- as.numeric(COAD["IFNGR1", ])
TGFB1_expr <- as.numeric(COAD["TGFB1", ])
CCL3_expr <- as.numeric(COAD["CCL3", ])
CXCL1_expr <- as.numeric(COAD["CXCL1", ])
CXCL2_expr <- as.numeric(COAD["CXCL2", ])
CXCL3_expr <- as.numeric(COAD["CXCL3", ])
CXCL5_expr <- as.numeric(COAD["CXCL5", ])
CXCL6_expr <- as.numeric(COAD["CXCL6", ])
CXCL7_expr <- as.numeric(COAD["CXCL7", ])
CXCL8_expr <- as.numeric(COAD["CXCL8", ])
CXCL9_expr <- as.numeric(COAD["CXCL9", ])
CXCL10_expr <- as.numeric(COAD["CXCL10", ])
CXCL11_expr <- as.numeric(COAD["CXCL11", ])
CXCL12_expr <- as.numeric(COAD["CXCL12", ])
CXCL13_expr <- as.numeric(COAD["CXCL13", ])
CXCR1_expr <- as.numeric(COAD["CXCR1", ])
CXCR2_expr <- as.numeric(COAD["CXCR2", ])
CXCR3_expr <- as.numeric(COAD["CXCR3", ])


# 加载必要的库
library(ggplot2)
library(ggpubr)

# Creating co-expression scatter plot. Writing a function for plotting
plot_correlation <- function(gene_x_expr, gene_y_expr, gene_x_name, gene_y_name) {
  df <- data.frame(x = gene_x_expr, y = gene_y_expr)
  
  plot <- ggscatter(df, x = "x", y = "y",
                    add = "reg.line", 
                    conf.int = TRUE,  
                    cor.coef = TRUE,  
                    cor.method = "pearson", 
                    xlab = paste0(gene_x_name, " Expression"), 
                    ylab = paste0(gene_y_name, " Expression"),
                    color = "skyblue") +  
    stat_cor(method = "pearson", label.x = min(gene_x_expr) * 0.8) + 
    geom_smooth(method = "lm", color = "grey", fill = "lightgrey", alpha = 0.5) +
    theme_minimal() +  
    theme(panel.grid = element_blank(), 
          axis.line = element_line(color = "black"), 
          axis.title = element_text(face = "bold"),
          axis.ticks = element_line(color = "black"))
  
  return(plot)
}

# Generate scatterplots
plot_list <- list()
gene_names <- c("HMGB1" , "IDO1", "IFNG", "IFNGR1", "TGFB1", "CCL3", "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", 
                "CXCL8", "CXCL9", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCR1", "CXCR2", "CXCR3")
gene_exprs <- list(HMGB1_expr, IDO1_expr, IFNG_expr, IFNGR1_expr, TGFB1_expr, CCL3_expr, CXCL1_expr, CXCL2_expr, CXCL3_expr, 
                   CXCL5_expr, CXCL6_expr, CXCL8_expr, CXCL9_expr, CXCL10_expr, CXCL11_expr, 
                   CXCL12_expr, CXCL13_expr, CXCR1_expr, CXCR2_expr, CXCR3_expr)

for (i in 1:length(gene_names)) {
  plot_list[[i]] <- plot_correlation(DDR1_expr, gene_exprs[[i]], "DDR1", gene_names[i])
}

# Combine plots
combined_plot <- ggarrange(plotlist = plot_list, ncol = 4, nrow = 6)
print(combined_plot)
ggsave( combined_plot, filename = "Correlation_DDR1_Chemokine.png", width = 15, height = 18)

