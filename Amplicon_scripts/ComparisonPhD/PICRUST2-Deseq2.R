library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggplot2)
#Import the table from QIIME2
pathway_table <- read.table("feature-table.tsv", 
                            header=TRUE, http://127.0.0.1:38141/graphics/plot_zoom_png?width=1468&height=687
                            sep="\t", 
                            skip=1, 
                            row.names=1, 
                            comment.char="")

pathway_table <- t(pathway_table)
metadata <- read.table("../../Coatemeta.txt",header = T,sep = "\t",row.names = 1)

# Convert periods in pathway_table row names to hyphens to match metadata
fixed_rownames <- gsub("\\.", "-", rownames(pathway_table))
rownames(pathway_table) <- fixed_rownames
fixed_rownames <- gsub("\\-1-", "-1.", rownames(pathway_table))
rownames(pathway_table) <- fixed_rownames
fixed_rownames <- gsub("\\Coatza-1.2018", "Coatza-1-2018", rownames(pathway_table))
rownames(pathway_table) <- fixed_rownames

#############Find common samples#############

# Find common samples
common_samples <- intersect(rownames(pathway_table), rownames(metadata))
cat("\n\nNumber of common samples found:", length(common_samples), "\n")

pathway_subset <- pathway_table[common_samples, ]
metadata_subset <- metadata[common_samples, , drop=FALSE]

# Print dimensions to verify
cat("\nDimensions of pathway_subset:", dim(pathway_subset)[1], "samples x", dim(pathway_subset)[2], "pathways\n")
cat("Dimensions of metadata_subset:", dim(metadata_subset)[1], "samples x", dim(metadata_subset)[2], "variables\n")

# Round to integers (DESeq2 requires count data)
pathway_subset <- round(pathway_subset)

# Create DESeq2 object - need to transpose back for DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = t(pathway_subset),  # Transpose back to get features as rows, samples as columns
  colData = metadata_subset,
  design = ~ Sediment  # Use your metadata column of interest
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results
res <- results(dds)

# Add pathway names
res$pathway <- rownames(res)

# Sort by adjusted p-value
res_sorted <- res[order(res$padj), ]

# Get significant pathways (adjust p-value threshold as needed)
sig_pathways <- subset(res_sorted, padj < 0.05)

# Write results to file
write.csv(as.data.frame(res_sorted), "deseq2_all_results.csv", row.names=FALSE)
write.csv(as.data.frame(sig_pathways), "deseq2_significant_pathways.csv", row.names=FALSE)

# Plot significant pathways
if(nrow(sig_pathways) > 0) {
  # Convert to data frame for plotting
  sig_df <- as.data.frame(sig_pathways)
  
  # If there are many significant pathways, limit to top 20
  if(nrow(sig_df) > 20) {
    sig_df <- sig_df[1:20, ]
  }
  
  # Create volcano plot
  ggplot(as.data.frame(res_sorted), aes(x=log2FoldChange, y=-log10(pvalue))) +
    geom_point(aes(color=padj<0.05), size=1) +
    scale_color_manual(values=c("grey", "red")) +
    labs(title="Volcano plot of differential abundance",
         x="Log2 Fold Change",
         y="-Log10 p-value") +
    theme_minimal() +
    geom_text(data=sig_df, 
              aes(label=pathway), 
              vjust=1.5, 
              size=3)
  
  ggsave("volcano_plot.png", width=10, height=8)
  
  # Create fold change plot for top significant pathways
  ggplot(sig_df, aes(x=reorder(pathway, log2FoldChange), y=log2FoldChange)) +
    geom_bar(stat="identity", aes(fill=log2FoldChange>0)) +
    scale_fill_manual(values=c("darkblue", "darkred"), 
                      labels=c("Decreased", "Increased"),
                      name="Change") +
    coord_flip() +
    labs(title="Significant Differentially Abundant Pathways",
         x="Pathway", 
         y="Log2 Fold Change") +
    theme_minimal()
  
  ggsave("fold_change_plot.png", width=10, height=8)
  
  # Create heatmap of significant pathways
  # Extract normalized counts
  normalized_counts <- counts(dds, normalized=TRUE)
  sig_norm_counts <- normalized_counts[rownames(sig_pathways), ]
  
  # Create metadata annotation
  anno <- as.data.frame(metadata_subset$Sediment)
  rownames(anno) <- rownames(metadata_subset)
  colnames(anno) <- "Sediment"
  
  # Plot heatmap
  PLOTH=pheatmap(log2(sig_norm_counts + 1),
           annotation_col = anno,
           scale = "row",
           show_rownames = TRUE,
           fontsize_row = 8,
           main = "Heatmap of Significant Pathways",
           filename = "significant_pathways_heatmap.png",
           width = 10,
           height = 20)
  
  # Print a summary of the results
  cat("Total pathways tested:", nrow(res_sorted), "\n")
  cat("Significant pathways (FDR < 0.05):", nrow(sig_pathways), "\n")
  
  # Print top 5 significant pathways
  cat("\nTop 5 significant pathways:\n")
  print(sig_pathways[1:min(5, nrow(sig_pathways)), c("pathway", "log2FoldChange", "padj")])
} else {
  cat("No significant pathways found at FDR < 0.05\n")
}
#############################################

# After running the DESeq2 analysis and getting sig_pathways

# Create a function to make a ggplot2 heatmap with better size control
create_heatmap <- function(norm_counts, metadata_info, num_pathways, comparison_name, reference_name) {
  # Get the top pathways based on adjusted p-value
  top_pathways <- head(rownames(sig_pathways), num_pathways)
  
  # Extract normalized counts for these pathways
  top_counts <- norm_counts[top_pathways, ]
  
  # Convert to long format for ggplot
  heatmap_data <- as.data.frame(log2(top_counts + 1))
  heatmap_data$pathway <- rownames(heatmap_data)
  
  # Melt to long format
  heatmap_long <- reshape2::melt(heatmap_data, id.vars = "pathway",
                                 variable.name = "sample", value.name = "abundance")
  
  # Add metadata information
  heatmap_long$sediment <- metadata_info$Sediment[match(heatmap_long$sample, rownames(metadata_info))]
  
  # Calculate the z-score for scaling (similar to pheatmap's "scale = row")
  heatmap_long <- heatmap_long %>%
    group_by(pathway) %>%
    mutate(z_score = (abundance - mean(abundance)) / sd(abundance)) %>%
    ungroup()
  
  # Create the plot
  p <- ggplot(heatmap_long, aes(x = sample, y = pathway, fill = z_score)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                         name = "Z-score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          panel.grid = element_blank()) +
    labs(title = paste0("Top ", num_pathways, " Differential Pathways"),
         subtitle = paste0(comparison_name, " vs ", reference_name),
         x = "Sample", y = "Pathway") +
    # Add sample group information at the top
    facet_grid(. ~ sediment, scales = "free_x", space = "free_x")
  
  # Return the plot
  return(p)
}

# First, determine the reference and comparison levels
sediment_levels <- levels(factor(metadata_subset$Sediment))
reference_level <- sediment_levels[1]
comparison_level <- sediment_levels[2]

# Get normalized counts from DESeq2
normalized_counts <- counts(dds, normalized = TRUE)

# Create and save heatmaps for different numbers of top pathways
# Use more reasonable sizes based on the number of pathways
top_numbers <- c(10, 20, 30, 50, 100, min(203, nrow(sig_pathways)))

for (num in top_numbers) {
  # Create heatmap
  p <- create_heatmap(normalized_counts, metadata_subset, num, comparison_level, reference_level)
  
  # Calculate reasonable dimensions - scale height with number of pathways
  # but cap at reasonable limits
  width <- 10  # Fixed width in inches
  height <- min(20, max(7, num/5 + 4))  # Scale height with pathway count, but cap between 7-20 inches
  
  # Save the plot with limitsize=FALSE to override size limits
  ggsave(paste0("heatmap_top_", num, "_pathways.png"), p, 
         width = width, height = height, limitsize = FALSE)
  
  # If you want to display the plot in R, uncomment the next line
  # print(p)
  
  cat("Created heatmap for top", num, "pathways with dimensions", width, "x", height, "inches\n")
}

# For the top 20 plot specifically (which you might want to customize further)
p20 <- create_heatmap(normalized_counts, metadata_subset, 20, comparison_level, reference_level)

# Customize p20 further
p20_custom <- p20 + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold")
  )

# Save the customized plot with reasonable dimensions
ggsave("heatmap_top_20_custom.png", p20_custom, width = 10, height = 8, limitsize = FALSE)

# If you want to display this plot
print(p20_custom)

####################################
# After you've run your DESeq2 analysis and have sig_pathways

# First, determine the reference and comparison levels
sediment_levels <- levels(factor(metadata_subset$Sediment))
reference_level <- sediment_levels[1]
comparison_level <- sediment_levels[2]
cat("Reference level:", reference_level, "\n")
cat("Comparison level:", comparison_level, "\n")
cat("Log2 fold change interpretation: positive values mean higher in", comparison_level, "compared to", reference_level, "\n")

# Create tables of top enriched pathways for each sediment type
# Get normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
sig_norm_counts <- normalized_counts[rownames(sig_pathways), ]

# Find top pathways enriched in each group
enriched_in_comparison <- sig_pathways[sig_pathways$log2FoldChange > 0, ]
enriched_in_reference <- sig_pathways[sig_pathways$log2FoldChange < 0, ]

# Sort by absolute log2 fold change
enriched_in_comparison <- enriched_in_comparison[order(-enriched_in_comparison$log2FoldChange), ]
enriched_in_reference <- enriched_in_reference[order(enriched_in_reference$log2FoldChange), ]

# Save top enriched pathways for each group
write.csv(as.data.frame(enriched_in_comparison), 
          paste0("pathways_enriched_in_", comparison_level, ".csv"), 
          row.names=FALSE)
write.csv(as.data.frame(enriched_in_reference), 
          paste0("pathways_enriched_in_", reference_level, ".csv"), 
          row.names=FALSE)

# Print summary
cat("\nPathways significantly enriched in", comparison_level, ":", nrow(enriched_in_comparison), "\n")
cat("Pathways significantly enriched in", reference_level, ":", nrow(enriched_in_reference), "\n")

# Print top 10 pathways for each group
cat("\nTop 10 pathways enriched in", comparison_level, ":\n")
print(head(enriched_in_comparison[, c("pathway", "log2FoldChange", "padj")], 10))

cat("\nTop 10 pathways enriched in", reference_level, ":\n")
print(head(enriched_in_reference[, c("pathway", "log2FoldChange", "padj")], 10))

# Create visualization of top 20 enriched pathways in each group
# For comparison group
top_comparison <- head(enriched_in_comparison, 20)
p_comparison <- ggplot(top_comparison, aes(x=reorder(pathway, log2FoldChange), y=log2FoldChange)) +
  geom_bar(stat="identity", fill="darkred") +
  coord_flip() +
  labs(title=paste("Top 20 Pathways Enriched in", comparison_level),
       subtitle=paste("Compared to", reference_level),
       x="Pathway", 
       y="Log2 Fold Change") +
  theme_minimal()

print(p_comparison)
ggsave(paste0("top_pathways_in_", comparison_level, ".png"), p_comparison, width=10, height=8)

# For reference group
top_reference <- head(enriched_in_reference, 20)
p_reference <- ggplot(top_reference, aes(x=reorder(pathway, -log2FoldChange), y=-log2FoldChange)) +
  geom_bar(stat="identity", fill="darkblue") +
  coord_flip() +
  labs(title=paste("Top 20 Pathways Enriched in", reference_level),
       subtitle=paste("Compared to", comparison_level),
       x="Pathway", 
       y="Log2 Fold Change (absolute value)") +
  theme_minimal()

print(p_reference)
ggsave(paste0("top_pathways_in_", reference_level, ".png"), p_reference, width=10, height=8)

# Create heatmaps of top pathways in each group
# For top pathways in comparison group
top_comp_pathways <- rownames(top_comparison)
comp_counts <- sig_norm_counts[top_comp_pathways, ]

anno <- data.frame(Sediment = metadata_subset$Sediment)
rownames(anno) <- colnames(comp_counts)

# Create and display heatmap for comparison group
pheatmap(log2(comp_counts + 1),
         annotation_col = anno,
         scale = "row",
         cluster_cols = FALSE,  # Don't cluster samples to keep groups together
         cluster_rows = TRUE,   # Cluster pathways by similarity
         show_rownames = TRUE,
         fontsize_row = 9,
         main = paste("Top 20 Pathways Enriched in", comparison_level),
         filename = paste0("heatmap_top_in_", comparison_level, ".png"),
         width = 10,
         height = 10)

# For top pathways in reference group
top_ref_pathways <- rownames(top_reference)
ref_counts <- sig_norm_counts[top_ref_pathways, ]

# Create and display heatmap for reference group
pheatmap(log2(ref_counts + 1),
         annotation_col = anno,
         scale = "row",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         fontsize_row = 9,
         main = paste("Top 20 Pathways Enriched in", reference_level),
         filename = paste0("heatmap_top_in_", reference_level, ".png"),
         width = 10,
         height = 10)

# For direct viewing in R, you need to explicitly print the pheatmap
# Create pheatmap object without saving file
plt_comp <- pheatmap(log2(comp_counts + 1),
                     annotation_col = anno,
                     scale = "row",
                     cluster_cols = FALSE,
                     cluster_rows = TRUE,
                     show_rownames = TRUE,
                     fontsize_row = 9,
                     main = paste("Top 20 Pathways Enriched in", comparison_level))

# Display it in R
print(plt_comp)

# Similarly for reference group
plt_ref <- pheatmap(log2(ref_counts + 1),
                    annotation_col = anno,
                    scale = "row", 
                    cluster_cols = FALSE,
                    cluster_rows = TRUE,
                    show_rownames = TRUE,
                    fontsize_row = 9,
                    main = paste("Top 20 Pathways Enriched in", reference_level))

# Display it in R
print(plt_ref)


#######ORDINATIONS##############


# First, load the required packages
# Install if needed: install.packages(c("vegan", "ggplot2", "umap", "Rtsne"))
library(vegan)
library(ggplot2)
library(umap)
library(Rtsne)

# Use the normalized counts from DESeq2
normalized_counts <- counts(dds, normalized=TRUE)

# Transpose so samples are rows (required for ordination)
count_matrix <- t(normalized_counts)

# Add sample metadata information (for plotting)
count_data <- as.data.frame(count_matrix)
count_data$Sediment <- metadata_subset$Sediment

# PCA (Principal Component Analysis)
pca_result <- prcomp(count_matrix, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x[,1:2])
pca_df$Sediment <- metadata_subset$Sediment

# Calculate variance explained
var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# Plot PCA
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Sediment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2) +
  labs(title = "PCA of MetaCyc Pathway Abundances",
       x = paste0("PC1 (", var_explained[1], "% variance)"),
       y = paste0("PC2 (", var_explained[2], "% variance)")) +
  theme_minimal()

print(pca_plot)
ggsave("pca_plot.png", pca_plot, width = 8, height = 6)

# PCoA (Principal Coordinates Analysis) with Bray-Curtis distance
# First calculate distance matrix
dist_matrix <- vegdist(count_matrix, method = "bray")
pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)

# Create data frame for plotting
pcoa_df <- as.data.frame(pcoa_result$points)
colnames(pcoa_df) <- c("Axis1", "Axis2")
pcoa_df$Sediment <- metadata_subset$Sediment

# Calculate variance explained
pcoa_var_explained <- round(100 * pcoa_result$eig / sum(pcoa_result$eig), 1)

# Plot PCoA
pcoa_plot <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Sediment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2) +
  labs(title = "PCoA of MetaCyc Pathway Abundances (Bray-Curtis Distance)",
       x = paste0("Axis 1 (", pcoa_var_explained[1], "% variance)"),
       y = paste0("Axis 2 (", pcoa_var_explained[2], "% variance)")) +
  theme_minimal()

print(pcoa_plot)
ggsave("pcoa_plot.png", pcoa_plot, width = 8, height = 6)

# NMDS (Non-metric Multidimensional Scaling)
nmds_result <- metaMDS(dist_matrix, k = 2)
nmds_df <- as.data.frame(nmds_result$points)
colnames(nmds_df) <- c("NMDS1", "NMDS2")
nmds_df$Sediment <- metadata_subset$Sediment

# Plot NMDS
nmds_plot <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Sediment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2) +
  labs(title = "NMDS of MetaCyc Pathway Abundances",
       subtitle = paste("Stress =", round(nmds_result$stress, 3))) +
  theme_minimal()

print(nmds_plot)
ggsave("nmds_plot.png", nmds_plot, width = 8, height = 6)

# UMAP (Uniform Manifold Approximation and Projection)
# Set random seed for reproducibility
set.seed(123)
umap_result <- umap(count_matrix)
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Sediment <- metadata_subset$Sediment

# Plot UMAP
umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Sediment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2) +
  labs(title = "UMAP of MetaCyc Pathway Abundances") +
  theme_minimal()

print(umap_plot)
ggsave("umap_plot.png", umap_plot, width = 8, height = 6)

# t-SNE (t-Distributed Stochastic Neighbor Embedding)
# Set random seed for reproducibility
set.seed(123)
tsne_result <- Rtsne(count_matrix, perplexity = min(30, nrow(count_matrix)-1), check_duplicates = FALSE)
tsne_df <- as.data.frame(tsne_result$Y)
colnames(tsne_df) <- c("tSNE1", "tSNE2")
tsne_df$Sediment <- metadata_subset$Sediment

# Plot t-SNE
tsne_plot <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Sediment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2) +
  labs(title = "t-SNE of MetaCyc Pathway Abundances") +
  theme_minimal()

print(tsne_plot)
ggsave("tsne_plot.png", tsne_plot, width = 8, height = 6)

umap_plot + 
  +     geom_text(aes(label = rownames(nmds_df)), 
                  +               nudge_x = 0.1, nudge_y = 0.15, 
                  +               size = 4)

# PERMANOVA (statistical test for group differences)
permanova_result <- adonis2(dist_matrix ~ Sediment, data = data.frame(Sediment = metadata_subset$Sediment))
print(permanova_result)

