#!/usr/bin/Rscript

library(openxlsx)
library(scMINER)


# Do not create SparseExpressionSet object or add meta
pbmc.12k <- readscRNAseqData(file="../PBMC12k_input/", is.10x = T, CreateSparseEset = F, add.meta = F)


# set add.meta=T to run the quality control and store the info in the object
pbmc.12k.eset <- CreateSparseEset(data = pbmc.12k$raw.data, feature.data = pbmc.12k$feature.data, add.meta = T)


# Do not create SparseExpressionSet object or add meta
pbmc.12k.eset2 <- readscRNAseqData(file = "/Users/lding/Documents/scMINER/PBMC12k_input", is.10x = T, CreateSparseEset = T, add.meta = T)


cutoffs <- draw.scRNAseq.QC(SparseEset = pbmc.12k.eset,
                            project.name = "PBMC12k",
                            plot.dir = "/Users/lding/Documents/scMINER/PBMC12k_input/QC/",
                            group = "group",      # this indicate which meta data information will be use in x axis to group violin plots
                            output.cutoff = TRUE) # whether or not to output suggested cutoffs


# Perform the actual filtering
pbmc.12k.eset.filter <- preMICA.filtering(SparseEset = pbmc.12k.eset, cutoffs = cutoffs)


norm = 1e6 
exp.norm <- sweep(exprs(pbmc.12k.eset.filter), 2, norm/unname(Matrix::colSums(exprs(pbmc.12k.eset.filter))), '*')


# log transformation (required by MICA)
exp.log2 <- log(exp.norm + 1, base = 2)

# save as SparseEset
pbmc.12k.eset.log2 <- CreateSparseEset(data = exp.log2, meta.data = pData(pbmc.12k.eset.filter), feature.data = fData(pbmc.12k.eset.filter), add.meta = F)


# prepare MICA input
generateMICAinput(eset = pbmc.12k.eset.log2 , filepath = "/Users/lding/Documents/scMINER/PBMC12k_input/PBMC12k_MICA_input.h5ad")


# clean working environment
rm(exp.log2)
rm(exp.norm)

# load clustering results
pbmc.12k.eset.log2 <- readMICAoutput(eset = pbmc.12k.eset.log2, load_ClusterRes = TRUE, 
                       output_file = "/Users/lding/Documents/scMINER/PBMC12k_input/clustering_UMAP_euclidean_20_0.82.txt")

# X, Y, color_by specify coordinate and clustering label entries in the eset phenoData; pct is the size of the point
MICAplot(input_eset = pbmc.12k.eset.log2, X = "X", Y = "Y", color_by = "ClusterRes", pct = 0.5)


# Marker gene visualization
genes_of_interest <-c("CD3D", "CD27", "IL7R", "SELL", "CCR7", "IL32", "GZMA",
                      "GZMK", "DUSP2", "CD8A", "GZMH", "GZMB", "CD79A", "CD79B", "CD86", "CD14")

# Marker UMAP scatter plot
feature_highlighting(input_eset = pbmc.12k.eset.log2, target = genes_of_interest, feature = "geneSymbol", 
                     ylabel = "log2Exp", x = "X", y = "Y", pct.size = 0.3)

# Marker violoin plot
feature_vlnplot(input_eset = pbmc.12k.eset.log2, target = genes_of_interest, feature = "geneSymbol", 
                group_by = "ClusterRes", ylabel = "log2Exp", ncol = 4)

# Marker heatmap plot
feature_heatmap(input_eset = pbmc.12k.eset.log2, target = genes_of_interest, group_name = "ClusterRes", save_plot = FALSE, 
                width = 6, height = 6, name = "log2Exp")

# Read marker genes from an excel file
markers <- read.xlsx("/Users/lding/Git/scMINER/tech-docs/docs/tables/Immune_signatures.xlsx")

# Maker gene file has 3 required columns
head(markers)

# Draw a bubble plot
draw.marker.bbp(ref = markers, input_eset = pbmc.12k.eset.log2, width = 6, height = 4, feature = "geneSymbol",
                group_name = "ClusterRes", save_plot = FALSE)

