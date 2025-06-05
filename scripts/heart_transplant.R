# Load Libraries ===============================================================
# This section loads all necessary R packages for the analysis.
# Ensure these packages are installed in your R environment.
# If using renv, these will be managed automatically upon `renv::restore()`.
library(GEOquery)         # For downloading gene expression data from GEO.
library(tidyverse)        # A collection of R packages for data manipulation and visualization (ggplot2, dplyr, etc.).
library(ggpubr)           # For creating publication-ready plots with statistical annotations.
library(pheatmap)         # For drawing pretty heatmaps.
library(limma)            # For differential expression analysis of microarray data.
library(ggfortify)        # Provides `autoplot` for ggplot2, useful for visualizing PCA results.
library(EnhancedVolcano)  # For creating customizable volcano plots.
library(hugene10sttranscriptcluster.db) # Annotation package for the HuGene-1_0-st-v1 array (GSE150059).
library(pd.hg.u133.plus.2) # Platform design annotation for Affymetrix Human Genome U133 Plus 2.0 Array (GSE87301).
library(hgu133plus2.db)    # Annotation package for the Affymetrix Human Genome U133 Plus 2.0 Array (GSE87301).
library(annotate)         # Provides tools for annotation data (e.g., getSYMBOL).
library(org.Hs.eg.db)     # Genome-wide annotation for Human (for gene symbol mapping in blood data).
library(cluster)          # For clustering algorithms like `clara` in PCA.
library(rstatix)          # Provides pipe-friendly framework for statistical tests (e.g., t-test).
library(cowplot)          # For arranging multiple ggplot2 plots into a grid.
library(combiroc)         # For calculating and visualizing ROC curves for single and combined markers.
library(reshape2)         # For reshaping data frames (e.g., melt).


# Define Custom Plotting Theme =================================================
# This custom theme is applied to ggplot2 plots for consistent styling,
# enhancing readability and publication readiness.
theme_custom_volcano <- theme_bw(base_size = 24) +
  theme(
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'plain', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'plain', vjust = 1),
    axis.text.x = element_text(angle = 0, size = 18, vjust = 1),
    axis.text.y = element_text(angle = 0, size = 18, vjust = 0.5),
    axis.title = element_text(size = 18),
    legend.position = "right",
    legend.key = element_blank(),
    legend.key.size = unit(0.5, 'cm'),
    legend.text = element_text(size = 14),
    title = element_text(size = 14),
    legend.title = element_blank()
  )

# Define Common Plotting Elements (to avoid repetition)
# This function applies common `ggpubr` and theme elements to plots.
apply_common_theme <- function(plot_obj, angle = 0, legend_pos = "right") {
  plot_obj +
    theme_pubr(base_size = 24, legend = legend_pos, x.text.angle = angle) +
    labs_pubr() +
    theme(axis.line = element_line(size = 1, colour = "black"))
}


# 1. Data Acquisition and Initial Processing ===================================

# Download gene expression datasets from GEO.
# Data will be saved in the 'data/' directory.
# `getGEO` returns a list of ExpressionSet objects; we select the first one if only one exists.
message("Downloading GSE150059 (Discovery Cohort)...")
gse150059 <- getGEO("GSE150059", destdir = "data/")
data_discovery <- gse150059[[1]]

message("Downloading GSE272655 (Validation Cohort 1)...")
gse272655 <- getGEO("GSE272655", destdir = "data/")
data_validation1 <- gse272655[[1]]

message("Downloading GSE87301 (Validation Cohort 2 - Blood)...")
gse87301 <- getGEO("GSE87301", destdir = "data/")
data_validation2_blood <- gse87301[[1]]


# Extract phenotype (clinical) data for the discovery cohort.
pheno_data_discovery <- pData(data_discovery)

# Standardize the column name for rejection status.
colnames(pheno_data_discovery)[which(colnames(pheno_data_discovery) == "mmdx (nr-normal, nr-minor, early injury, ptcmr, pabmr, abmr, tcmr, mixed):ch1")] <- "mmdx"

# Create a simplified 'group' variable for rejection status.
pheno_data_discovery <- pheno_data_discovery %>%
  mutate(group = if_else(str_detect(mmdx, "NR"), "No_rejection", "Rejection"))


# Plot sample counts by MMDx group for the discovery cohort.
message("Generating sample count plot...")
pheno_data_discovery %>%
  count(mmdx) %>%
  ggbarplot(x = "mmdx", y = "n", fill = "mmdx", color = "mmdx",
            xlab = "", ylab = "Sample count", label = TRUE, label.pos = "out") %>%
  apply_common_theme(angle = 30, legend_pos = "none") %>%
  ggsave("results/samples.png", dpi = 600, width = 8, height = 6)


# Extract feature (probe annotation) and expression data.
feature_data_discovery <- fData(data_discovery)
expression_data_discovery <- exprs(data_discovery)


# 2. Quality Control (QC) Analysis =============================================

# Boxplot of raw log2 expression values.
message("Generating expression boxplot...")
reshape2::melt(expression_data_discovery) %>%
  left_join(pheno_data_discovery %>% dplyr::select(Var2 = geo_accession, group), by = "Var2") %>%
  ggboxplot(x = "Var2", y = "value", fill = "group", xlab = "", ylab = "Log2 expression") %>%
  apply_common_theme(angle = 30, base_size = 7) %>%
  ggsave("results/boxplotggplot.png", dpi = 600, width = 10, height = 4)


# Sample correlation heatmap.
message("Generating correlation heatmap...")
pdf("results/CorHeatmap.pdf", width = 15, height = 15)
pheatmap(cor(expression_data_discovery),
         labels_row = pheno_data_discovery$group,
         labels_col = pheno_data_discovery$group,
         main = "Sample Correlation Heatmap") # Added a title for clarity
dev.off()


# Principal Component Analysis (PCA) for raw data.
message("Performing PCA on raw expression data...")
pca_raw <- prcomp(t(expression_data_discovery), scale. = FALSE)

percent_var_raw <- round(100 * pca_raw$sdev^2 / sum(pca_raw$sdev^2), 1)

data_pca_plot <- data.frame(PC1 = pca_raw$x[, 1], PC2 = pca_raw$x[, 2],
                            Group = pheno_data_discovery$group,
                            mmdx = pheno_data_discovery$mmdx)

ggplot(data_pca_plot, aes(PC1, PC2, color = Group)) +
  geom_point() +
  stat_ellipse() +
  ggtitle("PCA plot of the log-transformed expression data") +
  xlab(paste0("PC1, VarExp: ", percent_var_raw[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percent_var_raw[2], "%")) %>%
  apply_common_theme(legend_pos = "right") %>%
  ggsave("results/pca.png", dpi = 600, width = 8, height = 8)


# 3. Differential Expression (DE) Analysis =====================================

message("Starting differential expression analysis...")

# Average duplicated probes for the same gene symbol.
# This is a common step before DE analysis when multiple probes map to the same gene.
# Original data: probes as rows, samples as columns.
# `avereps` requires probe IDs as row names.
expression_data_discovery_averaged <- avereps(expression_data_discovery, ID = rownames(expression_data_discovery))

# Join with feature data to get gene symbols.
# Ensure 'ID' in `feature_data_discovery` corresponds to row names in `expression_data_discovery_averaged`.
expression_data_discovery_annotated <- expression_data_discovery_averaged %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(feature_data_discovery %>% dplyr::select(ID, `Gene Symbol`), by = "ID")

# Clean gene symbols: remove extra information after a space or slash.
# Example: "Gene1 // Description" becomes "Gene1".
expression_data_discovery_annotated$`Gene Symbol` <- gsub(" \\/.*", "", expression_data_discovery_annotated$`Gene Symbol`)

# Remove the original probe ID column and move 'Gene Symbol' to the first column.
expression_data_discovery_annotated <- expression_data_discovery_annotated %>%
  dplyr::select(-ID) %>%
  relocate(`Gene Symbol`)

# Handle duplicated gene symbols by averaging their expression.
# This is typically done if multiple probes still map to the same (cleaned) gene symbol.
expression_data_discovery_final <- expression_data_discovery_annotated %>%
  group_by(`Gene Symbol`) %>%
  summarize(across(starts_with("GSM"), mean), .groups = 'drop') %>%
  column_to_rownames("Gene Symbol")

# Remove first 28 rows (likely control probes or unannotated entries).
# This is a specific cleanup step based on the dataset's characteristics.
expression_data_discovery_final <- expression_data_discovery_final[-c(1:28), ]

# Prepare design matrix for limma.
# `group` variable from phenotype data is used to define experimental groups.
pheno_data_discovery$group <- factor(pheno_data_discovery$group, levels = c("No_rejection", "Rejection")) # Ensure correct order for contrast
design_matrix <- model.matrix(~ 0 + group, pheno_data_discovery)
colnames(design_matrix) <- levels(pheno_data_discovery$group)

# Fit linear model to expression data.
fit <- lmFit(expression_data_discovery_final, design_matrix)

# Define contrast matrix for differential expression: Rejection vs. No_rejection.
contrast_matrix <- makeContrasts(Rejection - No_rejection, levels = design_matrix)

# Apply contrasts to the fitted model.
fit2 <- contrasts.fit(fit, contrast_matrix)

# Apply eBayes to moderate gene-wise variances.
fit2 <- eBayes(fit2)

# Extract top differentially expressed genes.
# `adjust="fdr"` applies Benjamini-Hochberg for multiple testing correction.
# `sort.by="B"` sorts by B-statistic (log-odds of differential expression).
results_table_de <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)


# Generate Volcano Plot.
message("Generating volcano plot...")
EnhancedVolcano(results_table_de,
                x = "logFC",
                y = "adj.P.Val",
                lab = rownames(results_table_de),
                pCutoff = 0.05,
                FCcutoff = 1, # Fold change cutoff (log2FC)
                title = "Differential Expression: Rejection vs. No Rejection",
                caption = paste0('Total genes = ', nrow(results_table_de)),
                subtitle = "Filtered by |log2FC| > 1 and adj.P.Val < 0.05") %>%
  # Custom x-axis breaks for better visualization
  scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3), labels = c("-3", "-2", "-1", "0", "1", "2", "3")) %>%
  apply_common_theme(legend_pos = "none") %>%
  ggsave("results/volcano.png", dpi = 600, width = 8, height = 8)


# Filter and save significant genes.
sig_genes_df <- results_table_de %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  rownames_to_column("Genes")

write.csv(sig_genes_df, file = "results/significant_genes.csv", row.names = FALSE)

# Save expression data for significant genes.
expression_data_sig <- expression_data_discovery_final %>%
  filter(rownames(.) %in% sig_genes_df$Genes) %>%
  rownames_to_column("Genes")

write.csv(expression_data_sig, file = "results/expression.csv", row.names = FALSE)

# Save relevant clinical data.
# Selecting specific columns for easier sharing and analysis.
pheno_data_save <- pheno_data_discovery[, c("geo_accession", "mmdx", "group")] %>%
  rownames_to_column("geo_accesion_full") # Renamed to avoid confusion with internal GEO IDs

write.csv(pheno_data_save, file = "results/clinical_data.csv", row.names = FALSE)


# 4. Overrepresentation Analysis using Enrichr (Assumes pre-generated files) ===
# This section assumes that gene lists from DE analysis (e.g., significant_genes.csv)
# were uploaded to Enrichr (https://maayanlab.cloud/Enrichr/) and the results
# (e.g., "*_table.txt" files) have been downloaded into the 'results/' directory.
message("Performing overrepresentation analysis (requires pre-downloaded Enrichr results)...")

# List all Enrichr result files.
enrichr_files <- list.files(path = "results/", pattern = "_table.txt", recursive = TRUE, full.names = TRUE)

# Read and process each Enrichr result file.
enrich_results_list <- map(enrichr_files, read.delim)
names(enrich_results_list) <- gsub("_2.*", "", basename(enrichr_files)) # Clean up names for plotting

# Select and sort top 10 enriched terms.
enrich_sig_processed <- map(enrich_results_list, function(x) {
  temp <- arrange(x, Adjusted.P.value)
  temp <- head(temp, 10)
  temp <- arrange(temp, desc(Adjusted.P.value)) # Order for bar plot
  temp$Term <- factor(temp$Term, levels = unique(temp$Term)) # Ensure factor order for plotting
  return(temp)
})

# Pad shorter labels for consistent plot aesthetics across multiple plots.
# This helps align terms when `cowplot` is used.
max_length_term <- max(str_length(bind_rows(enrich_sig_processed)$Term)) # Find max length across all terms
enrich_sig_padded <- map(enrich_sig_processed, function(x) {
  mutate(x, Term = str_pad(Term, max_length_term, side = "right")) %>%
    arrange(desc(Adjusted.P.value)) %>%
    mutate(Term = factor(Term, levels = unique(Term))) # Re-factor after padding
})


# Generate bar plots for enriched terms.
enrichment_plots <- map2(enrich_sig_padded, names(enrich_sig_padded), function(x, y) {
  ggplot(x, aes(-log10(Adjusted.P.value), Term)) +
    geom_bar(aes(fill = -log10(Adjusted.P.value)), stat = "identity") +
    scale_fill_gradientn(colors = gplots::colorpanel(256, "#00468BFF", "#ED0000FF")) + # Custom color scale
    labs(x = "-log10(adjusted p-value)", y = "", fill = "-log10(adjusted p-value)") +
    ggtitle(y) +
    scale_y_discrete(position = "right") %>% # Place Y-axis labels on the right
    apply_common_theme(legend_pos = "right") # Apply common theme for consistency
})

# Combine all enrichment plots into a single figure.
message("Combining enrichment plots...")
figure_enrichment_combined <- plot_grid(plotlist = enrichment_plots, ncol = 2, align = "v")
ggsave("results/figure_enrichment_combined.png", figure_enrichment_combined, dpi = 500, width = 20, height = 10)


# 5. Validation Analysis (Internal and External) ===============================

message("Starting validation analysis...")

# Reload data saved from previous steps for consistency and modularity.
clinical_data_loaded <- read_csv("results/clinical_data.csv")
expression_loaded <- read_csv("results/expression.csv")
significant_genes_loaded <- read_csv("results/significant_genes.csv")


# PCA with clustering on significant genes (Internal Validation)
message("Performing PCA with clustering on significant genes...")
pca_res_deg <- prcomp(expression_loaded %>% column_to_rownames("Genes") %>% t(), scale. = TRUE)

# Perform clustering (e.g., clara for robust clustering).
# Number of clusters (3) is an assumption; may need justification or optimization in a full study.
set.seed(123) # for reproducibility of clustering
clust_pca <- clara(pca_res_deg$x, 3)

# Prepare data for PCA plotting with cluster assignments.
pca_df_deg <- pca_res_deg$x[, c(1, 2)] %>%
  as.data.frame() %>%
  rownames_to_column("geo_accession_full") %>%
  left_join(data.frame(geo_accession_full = names(clust_pca$clustering),
                       clusters = clust_pca$clustering), by = "geo_accession_full") %>%
  left_join(clinical_data_loaded %>% dplyr::select(geo_accession_full, group, mmdx), by = "geo_accession_full")


percent_var_deg <- round(100 * pca_res_deg$sdev^2 / sum(pca_res_deg$sdev^2), 1)

# Plot PCA coloured by rejection group and outlined by clusters.
ggplot(pca_df_deg, aes(PC1, PC2)) +
  geom_point(aes(color = group)) +
  stat_ellipse(aes(color = factor(clusters))) + # Ellipses for clusters
  xlab(paste0("PC1, VarExp: ", percent_var_deg[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percent_var_deg[2], "%")) %>%
  apply_common_theme() %>%
  ggsave("results/pca_deg.png", dpi = 600, width = 8, height = 8)

# Plot PCA for each cluster, colored by MMDx sub-group (for finer-grained analysis).
message("Generating individual PCA plots for clusters...")
for (i in 1:3) {
  ggplot(pca_df_deg %>% filter(clusters == i), aes(PC1, PC2)) +
    geom_point(aes(color = mmdx)) +
    xlab(paste0("PC1")) + # Simplified axis labels for sub-plots
    ylab(paste0("PC2")) %>%
    apply_common_theme() %>%
    ggsave(paste0("results/pca_deg_cluster_", i, ".png"), dpi = 600, width = 8, height = 8)
}


# Boxplots for top significant genes (Internal Validation)
# Selected top 10 genes based on preliminary analysis or specific criteria.
top_genes <- c("TYMS", "WARS", "AIM2", "CXCL9", "TRAT1", "HLA-DRB3", "TNFRSF9", "GZMH", "IL32", "AIF1")

expression_subset_topgenes <- expression_loaded %>%
  filter(Genes %in% top_genes) %>%
  column_to_rownames("Genes") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("geo_accesion_full")

boxplot_data_internal <- clinical_data_loaded %>%
  left_join(expression_subset_topgenes, by = "geo_accesion_full") %>%
  dplyr::select(-geo_accesion_full, -mmdx) %>% # Remove unneeded columns for melt
  reshape2::melt(id.vars = "group") # Melt for ggplot2 facetting

colnames(boxplot_data_internal)[colnames(boxplot_data_internal) == "variable"] <- "gene_name"

# Perform t-tests for each gene and adjust p-values.
stat_test_internal <- boxplot_data_internal %>%
  group_by(gene_name) %>%
  t_test(value ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Add significance annotation positions for plotting.
stat_test_internal_plot <- stat_test_internal %>%
  add_xy_position(x = "gene_name", dodge = 0.8) # Adjust dodge as needed


# Plot boxplots for top genes in the discovery cohort.
message("Generating boxplots for top genes (discovery cohort)...")
ggplot(boxplot_data_internal, aes(x = gene_name, y = value)) +
  geom_boxplot(aes(fill = group), outlier.size = 0.01) +
  stat_pvalue_manual(stat_test_internal_plot,
                     label = "p.adj.signif", hide.ns = TRUE,
                     tip.length = 0.03) +
  ylab("Expression") +
  xlab("") %>%
  apply_common_theme(angle = 30) %>%
  ggsave("results/rna_valid.png", dpi = 1000, width = 18, height = 12, units = "cm")


# Validation on External Dataset 1 (GSE272655 - Kidney biopsies)
message("Processing Validation Cohort 1 (GSE272655)...")

exp_valid1 <- exprs(data_validation1)
feature_data_valid1 <- fData(data_validation1)
pheno_data_valid1 <- pData(data_validation1)

# Average duplicated probes (if any) and map to gene symbols for validation data.
exp_data_valid1_averaged <- avereps(exp_valid1, ID = rownames(exp_valid1))

# Note: `annot_data` is from the discovery cohort. It's assumed that probe IDs
# are consistent enough between platforms for `left_join` to work for gene symbol mapping.
# For robust mapping, consider `annotate::mapIds` or similar specific to the platform.
exp_data_valid1_annotated <- exp_data_valid1_averaged %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(feature_data_discovery %>% dplyr::select(ID, `Gene Symbol`), by = "ID") # Using discovery feature data for mapping

# Clean gene symbols and process similar to discovery data.
exp_data_valid1_annotated$`Gene Symbol` <- gsub(" \\/.*", "", exp_data_valid1_annotated$`Gene Symbol`)
exp_data_valid1_final <- exp_data_valid1_annotated %>%
  dplyr::select(-ID) %>%
  relocate(`Gene Symbol`) %>%
  group_by(`Gene Symbol`) %>%
  summarize(across(starts_with("GSM"), mean), .groups = 'drop') %>%
  column_to_rownames("Gene Symbol")

exp_data_valid1_final <- exp_data_valid1_final[-c(1:28),] # Remove control probes based on observation from discovery data

# Prepare phenotype data for validation cohort 1.
pheno_data_valid1$group <- ifelse(pheno_data_valid1$`mmdx (abmr, tcmr, mixed, pabmr, nr):ch1` == "NR", "No_rejection", "Rejection")
pheno_data_valid1$group <- factor(pheno_data_valid1$group, levels = c("No_rejection", "Rejection"))

# Subset expression data for the top selected genes.
exp_data_valid1_subset_topgenes <- exp_data_valid1_final[which(rownames(exp_data_valid1_final) %in% top_genes), ]

# Prepare data for boxplots.
boxplot_data_valid1 <- exp_data_valid1_subset_topgenes %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("geo_accession") %>%
  left_join(pheno_data_valid1 %>% dplyr::select(geo_accession, group), by = "geo_accession") %>%
  dplyr::select(-geo_accession) %>% # Remove geo_accession before melt
  reshape2::melt(id.vars = "group")

colnames(boxplot_data_valid1)[colnames(boxplot_data_valid1) == "variable"] <- "gene_name"

# Perform t-tests and adjust p-values for validation cohort 1.
stat_test_valid1 <- boxplot_data_valid1 %>%
  group_by(gene_name) %>%
  t_test(value ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

stat_test_valid1_plot <- stat_test_valid1 %>%
  add_xy_position(x = "gene_name", dodge = 0.8)


# Plot boxplots for top genes in validation cohort 1.
message("Generating boxplots for top genes (validation cohort 1)...")
ggplot(boxplot_data_valid1, aes(x = gene_name, y = value)) +
  geom_boxplot(aes(fill = group), outlier.size = 0.01) +
  stat_pvalue_manual(stat_test_valid1_plot,
                     label = "p.adj.signif", hide.ns = TRUE,
                     tip.length = 0.03) +
  ylab("Expression") +
  xlab("") %>%
  apply_common_theme(angle = 30) %>%
  ggsave("results/rna_valid_valid.png", dpi = 1000, width = 18, height = 12, units = "cm")


# Validation on External Dataset 2 (GSE87301 - Blood)
message("Processing Validation Cohort 2 (GSE87301 - Blood)...")

exp_blood <- exprs(data_validation2_blood)
pheno_blood <- pData(data_validation2_blood)

# Map probes to gene symbols for the blood dataset (different platform).
# This uses platform-specific annotation packages (`hgu133plus2.db`, `org.Hs.eg.db`).
mapped_probes_blood <- mappedkeys(hgu133plus2SYMBOL)
gene_symbols_blood <- getSYMBOL(mapped_probes_blood, "hgu133plus2")
probe_to_symbol_df_blood <- data.frame(ProbeID = mapped_probes_blood, GeneSymbol = gene_symbols_blood) %>%
  # Remove NAs which occur when a probe cannot be mapped to a gene symbol
  filter(!is.na(GeneSymbol))

# Process expression data for blood.
exp_data_blood_processed <- exp_blood %>%
  as.data.frame() %>%
  rownames_to_column("ProbeID") %>%
  left_join(probe_to_symbol_df_blood, by = "ProbeID") %>%
  dplyr::select(-ProbeID) %>% # Remove probe ID as we have gene symbol
  filter(!is.na(GeneSymbol)) %>% # Remove rows where GeneSymbol is NA after join
  group_by(GeneSymbol) %>%
  dplyr::summarize(across(starts_with("GSM"), mean), .groups = 'drop') %>%
  column_to_rownames("GeneSymbol")

# Subset for top genes.
exp_data_blood_subset_topgenes <- exp_data_blood_processed %>%
  filter(rownames(.) %in% top_genes)

# Prepare phenotype data for blood cohort.
pheno_blood_processed <- pheno_blood %>%
  dplyr::select(geo_accession, group = `rejection status:ch1`) %>%
  mutate(group = if_else(group == "NR", "No_rejection", "Rejection")) %>%
  mutate(group = factor(group, levels = c("No_rejection", "Rejection")))


# Prepare data for boxplots for blood cohort.
boxplot_data_blood <- exp_data_blood_subset_topgenes %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("geo_accession") %>%
  left_join(pheno_blood_processed, by = "geo_accession") %>%
  dplyr::select(-geo_accession) %>% # Remove geo_accession before melt
  reshape2::melt(id.vars = "group")

colnames(boxplot_data_blood)[colnames(boxplot_data_blood) == "variable"] <- "gene_name"

# Perform t-tests and adjust p-values for blood cohort.
stat_test_blood <- boxplot_data_blood %>%
  group_by(gene_name) %>%
  t_test(value ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

stat_test_blood_plot <- stat_test_blood %>%
  add_xy_position(x = "gene_name", dodge = 0.8)

# Plot boxplots for top genes in blood validation cohort.
message("Generating boxplots for top genes (blood validation cohort)...")
ggplot(boxplot_data_blood, aes(x = gene_name, y = value)) +
  geom_boxplot(aes(fill = group), outlier.size = 0.01) +
  stat_pvalue_manual(stat_test_blood_plot,
                     label = "p.adj.signif", hide.ns = TRUE,
                     tip.length = 0.03) +
  ylab("Expression") +
  xlab("") +
  theme_classic() + # Using theme_classic for this plot as in original.
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) %>%
  ggsave("results/rna_valid_blood.png", dpi = 1000, width = 18, height = 10, units = "cm")


# 6. ROC Analysis (CombiROC) ===================================================

message("Starting ROC analysis using CombiROC...")

# The `combiroc` package expects `Patient.ID` and `Class` columns.
# Adjusting column names for the discovery cohort data.
roc_info_discovery <- expression_subset_topgenes %>%
  rownames_to_column("Patient.ID") %>%
  left_join(clinical_data_loaded %>% dplyr::select(geo_accesion_full, group),
            by = c("Patient.ID" = "geo_accesion_full")) %>%
  dplyr::rename(Class = group) %>%
  dplyr::select(Patient.ID, Class, everything()) # Reorder columns

# Fix the gene symbol for HLA-DRB3 if it was changed for R variable naming.
top_genes_roc <- c("TYMS","WARS","AIM2","CXCL9","TRAT1","HLA.DRB3","TNFRSF9","GZMH","IL32","AIF1")


## 6.1 Single Markers ROC Analysis (Discovery Cohort) ----

# Prepare data for combiroc.
data_long_roc_discovery <- combiroc_long(roc_info_discovery)

# Calculate single marker statistics to get optimal threshold.
sms_discovery <- single_markers_statistics(data_long_roc_discovery)
distr_discovery <- markers_distribution(data_long_roc_discovery,
                                        case_class = 'Rejection',
                                        signalthr_prediction = TRUE)
# Get the Youden's J threshold.
pr_discovery <- distr_discovery$Coord[distr_discovery$Coord$Youden == max(distr_discovery$Coord$Youden), "threshold"][1]

# Compute single marker ROCs.
tab_single_discovery <- combi(roc_info_discovery,
                              signalthr = pr_discovery,
                              combithr = 1, # Only single markers
                              case_class = 'Rejection',
                              max_length = 1)

# Generate ROC reports for selected single markers.
reports_single_discovery <- roc_reports(roc_info_discovery,
                                        markers_table = tab_single_discovery,
                                        case_class = 'Rejection',
                                        single_markers = top_genes_roc,
                                        selected_combinations = NULL)

# Extract ROC metrics and sort by AUC.
roc_stats_single_discovery <- reports_single_discovery$Metrics %>% arrange(desc(AUC))

write_csv(roc_stats_single_discovery %>% rownames_to_column("Marker"),
          "results/roc_single_discovery.csv")

# Function to generate distinct colors.
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Plot single marker ROC curves.
reports_single_discovery$Plot$data$name <- factor(reports_single_discovery$Plot$data$name,
                                                  levels = rownames(roc_stats_single_discovery))

plot_single_markers_discovery <- reports_single_discovery$Plot +
  coord_equal() + # Ensure 1:1 aspect ratio
  ggtitle("Single Markers ROC (Discovery Cohort)") +
  labs(color = "Markers") +
  scale_color_manual(labels = paste0(rownames(roc_stats_single_discovery), " (AUC = ",
                                     round(roc_stats_single_discovery$AUC * 100, 2), "%)"),
                     values = gg_color_hue(nrow(roc_stats_single_discovery))) +
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5))) %>%
  apply_common_theme(legend_pos = "right")

ggsave("results/single_markers_discovery.png", plot_single_markers_discovery,
       dpi = 1000, units = "in", height = 4, width = 8)


## 6.2 Dual Markers ROC Analysis (Discovery Cohort) ----

# Compute dual marker ROCs.
tab_dual_discovery <- combi(roc_info_discovery,
                            signalthr = pr_discovery,
                            combithr = 1, # Base threshold for combination
                            case_class = 'Rejection',
                            max_length = 2) # Combinations of up to two markers

# Generate ROC reports for dual markers.
# Selecting top N combinations for visualization.
top_n_combinations <- 10 # You might adjust this number
if (nrow(tab_dual_discovery) < top_n_combinations) {
  top_n_combinations <- nrow(tab_dual_discovery)
}

# Ensure `tab_dual_discovery` is sorted by AUC or other metric if needed for selection.
# `combi` function usually orders them, but explicit sort can be added.
ranked_combs_dual_discovery <- ranked_combs(tab_dual_discovery, min_SE = 40, min_SP = 70) # Optional filtering
selected_comb_indices <- as.numeric(gsub("Combination ", "", ranked_combs_dual_discovery$Combinations[1:top_n_combinations])) # Get indices of top N combinations

reports_dual_discovery <- roc_reports(roc_info_discovery,
                                      markers_table = tab_dual_discovery,
                                      case_class = 'Rejection',
                                      single_markers = NULL,
                                      selected_combinations = selected_comb_indices)

# Extract ROC metrics and combine with marker names for dual combinations.
roc_stats_dual_discovery <- reports_dual_discovery$Metrics %>%
  arrange(desc(AUC)) %>%
  rownames_to_column("Combination_ID") %>%
  left_join(tab_dual_discovery %>% rownames_to_column("Combination_ID") %>% dplyr::select(Combination_ID, Markers),
            by = "Combination_ID") %>%
  relocate(Combination_ID, Markers)

write_csv(roc_stats_dual_discovery, "results/roc_dual_discovery.csv")

# Plot dual marker ROC curves.
reports_dual_discovery$Plot$data$name <- factor(reports_dual_discovery$Plot$data$name,
                                                levels = roc_stats_dual_discovery$Combination_ID)

plot_dual_markers_discovery <- reports_dual_discovery$Plot +
  coord_equal() +
  ggtitle("Dual Markers ROC (Discovery Cohort)") +
  labs(color = "Combined Markers") +
  scale_color_manual(labels = paste0(roc_stats_dual_discovery$Markers, " (AUC = ",
                                     round(roc_stats_dual_discovery$AUC * 100, 2), "%)"),
                     values = gg_color_hue(top_n_combinations)) +
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5))) %>%
  apply_common_theme(legend_pos = "right", angle = 30)

ggsave("results/dual_markers_discovery.png", plot_dual_markers_discovery,
       dpi = 1000, units = "in", height = 4, width = 8)


## 6.3 Single Markers ROC Analysis (Validation Cohort 1) ----

message("Starting ROC analysis for Validation Cohort 1...")

# Prepare validation cohort 1 data for combiroc.
# Ensure `HLA.DRB3` renaming is consistent if needed.
# The original `top_genes` list may contain `HLA-DRB3` which needs to be `HLA.DRB3` for column names.
roc_info_valid1 <- exp_data_valid1_subset_topgenes %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Patient.ID") %>%
  left_join(pheno_data_valid1 %>% dplyr::select(geo_accession, group),
            by = c("Patient.ID" = "geo_accession")) %>%
  dplyr::rename(Class = group) %>%
  dplyr::select(Patient.ID, Class, everything())

# Ensure `HLA-DRB3` in `top_genes_roc` matches the column name in `roc_info_valid1` if it got converted.
# If `exp_data_valid1_subset_topgenes` already has `HLA.DRB3` as a column name, no change needed.
if ("HLA-DRB3" %in% colnames(roc_info_valid1)) {
  roc_info_valid1 <- roc_info_valid1 %>% dplyr::rename(`HLA.DRB3` = `HLA-DRB3`)
  top_genes_roc[which(top_genes_roc == "HLA-DRB3")] <- "HLA.DRB3" # Update list used for plotting
}

data_long_roc_valid1 <- combiroc_long(roc_info_valid1)
sms_valid1 <- single_markers_statistics(data_long_roc_valid1)
distr_valid1 <- markers_distribution(data_long_roc_valid1,
                                     case_class = 'Rejection',
                                     signalthr_prediction = TRUE)
pr_valid1 <- distr_valid1$Coord[distr_valid1$Coord$Youden == max(distr_valid1$Coord$Youden), "threshold"][1]

tab_single_valid1 <- combi(roc_info_valid1,
                           signalthr = pr_valid1,
                           combithr = 1,
                           case_class = 'Rejection',
                           max_length = 1)

reports_single_valid1 <- roc_reports(roc_info_valid1,
                                     markers_table = tab_single_valid1,
                                     case_class = 'Rejection',
                                     single_markers = top_genes_roc,
                                     selected_combinations = NULL)

roc_stats_single_valid1 <- reports_single_valid1$Metrics %>% arrange(desc(AUC))

write_csv(roc_stats_single_valid1 %>% rownames_to_column("Marker"),
          "results/roc_single_valid1.csv")

reports_single_valid1$Plot$data$name <- factor(reports_single_valid1$Plot$data$name,
                                               levels = rownames(roc_stats_single_valid1))

plot_single_markers_valid1 <- reports_single_valid1$Plot +
  coord_equal() +
  ggtitle("Single Markers ROC (Validation Cohort 1)") +
  labs(color = "Markers") +
  scale_color_manual(labels = paste0(rownames(roc_stats_single_valid1), " (AUC = ",
                                     round(roc_stats_single_valid1$AUC * 100, 2), "%)"),
                     values = gg_color_hue(nrow(roc_stats_single_valid1))) +
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5))) %>%
  apply_common_theme(legend_pos = "right")

ggsave("results/single_markers_valid1.png", plot_single_markers_valid1,
       dpi = 1000, units = "in", height = 4, width = 8)


## 6.4 Dual Markers ROC Analysis (Validation Cohort 1) ----

tab_dual_valid1 <- combi(roc_info_valid1,
                         signalthr = pr_valid1,
                         combithr = 1,
                         case_class = 'Rejection',
                         max_length = 2)

# Select top N combinations, similar to discovery cohort.
if (nrow(tab_dual_valid1) < top_n_combinations) {
  top_n_combinations_valid1 <- nrow(tab_dual_valid1)
} else {
  top_n_combinations_valid1 <- top_n_combinations
}

ranked_combs_dual_valid1 <- ranked_combs(tab_dual_valid1, min_SE = 40, min_SP = 70)
selected_comb_indices_valid1 <- as.numeric(gsub("Combination ", "", ranked_combs_dual_valid1$Combinations[1:top_n_combinations_valid1]))


reports_dual_valid1 <- roc_reports(roc_info_valid1,
                                   markers_table = tab_dual_valid1,
                                   case_class = 'Rejection',
                                   single_markers = NULL,
                                   selected_combinations = selected_comb_indices_valid1)

roc_stats_dual_valid1 <- reports_dual_valid1$Metrics %>%
  arrange(desc(AUC)) %>%
  rownames_to_column("Combination_ID") %>%
  left_join(tab_dual_valid1 %>% rownames_to_column("Combination_ID") %>% dplyr::select(Combination_ID, Markers),
            by = "Combination_ID") %>%
  relocate(Combination_ID, Markers)

write_csv(roc_stats_dual_valid1, "results/roc_dual_valid1.csv")

reports_dual_valid1$Plot$data$name <- factor(reports_dual_valid1$Plot$data$name,
                                             levels = roc_stats_dual_valid1$Combination_ID)

plot_dual_markers_valid1 <- reports_dual_valid1$Plot +
  coord_equal() +
  ggtitle("Dual Markers ROC (Validation Cohort 1)") +
  labs(color = "Combined Markers") +
  scale_color_manual(labels = paste0(roc_stats_dual_valid1$Markers, " (AUC = ",
                                     round(roc_stats_dual_valid1$AUC * 100, 2), "%)"),
                     values = gg_color_hue(top_n_combinations_valid1)) +
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5))) %>%
  apply_common_theme(legend_pos = "right", angle = 30)

ggsave("results/dual_markers_valid1.png", plot_dual_markers_valid1,
       dpi = 1000, units = "in", height = 4, width = 8)

message("Analysis complete. All results saved in the 'results/' directory.")