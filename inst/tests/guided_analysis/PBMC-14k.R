#!/usr/bin/Rscript

library(openxlsx)
library(anndata)
library(scMINER)


# Do not create SparseExpressionSet object or add meta
pbmc.14k <- readscRNAseqData(file="/Users/lding/Documents/scMINER/PBMC14k_input/PBMC_20k_MICA_input_filter_14k.txt", 
                             is.10x = F, CreateSparseEset = F, add.meta = F, sep='\t')

# set add.meta=T to run the quality control and store the info in the object
pbmc.14k.eset <- CreateSparseEset(data = t(as.matrix(pbmc.14k)), add.meta = T)

cutoffs <- draw.scRNAseq.QC(SparseEset = pbmc.14k.eset,
                            project.name = "PBMC14k",
                            plot.dir = "/Users/lding/Documents/scMINER/PBMC14k_input/QC/",
                            group = "group",      # this indicate which meta data information will be use in x axis to group violin plots
                            output.cutoff = TRUE) # whether or not to output suggested cutoffs

# Perform the actual filtering
pbmc.14k.eset.filter <- preMICA.filtering(SparseEset = pbmc.14k.eset, cutoffs = cutoffs)

norm = 1e6 
exp.norm <- sweep(exprs(pbmc.14k.eset.filter), 2, norm/unname(Matrix::colSums(exprs(pbmc.14k.eset.filter))), '*')

# log transformation (required by MICA)
exp.log2 <- log(exp.norm + 1, base = 2)

# save as SparseEset
pbmc.14k.eset.log2 <- CreateSparseEset(data = exp.log2, meta.data = pData(pbmc.14k.eset.filter), 
                                       feature.data = fData(pbmc.14k.eset.filter), add.meta = F)

# prepare MICA input
generateMICAinput(eset = pbmc.14k.eset.log2 , filepath = "/Users/lding/Documents/scMINER/PBMC14k_input/PBMC14k_MICA_input.h5ad")

# clean working environment
rm(exp.log2)
rm(exp.norm)

# load clustering results
pbmc.14k.eset.log2 <- readMICAoutput(eset = pbmc.14k.eset.log2, load_ClusterRes = TRUE, 
                      output_file = "/Users/lding/Git/scMINER/docs/docs/images/clustering_UMAP_euclidean_24_1.82212.txt")

# X, Y, color_by specify coordinate and clustering label entries in the eset phenoData; pct is the size of the point
MICAplot(input_eset = pbmc.14k.eset.log2, X = "X", Y = "Y", color_by = "ClusterRes", pct = 0.5)


# Marker gene visualization
genes_of_interest <-c("CD3D", "CD27", "IL7R", "SELL", "CCR7", "IL32", "GZMA",
                      "GZMK", "DUSP2", "CD8A", "GZMH", "GZMB", "CD79A", "CD79B", "CD86", "CD14")

# Marker UMAP scatter plot
feature_highlighting(input_eset = pbmc.14k.eset.log2, target = genes_of_interest, feature = "geneSymbol",
                     ylabel = "log2Exp", x = "X", y = "Y", pct.size = 0.3)

# Marker violoin plot
feature_vlnplot(input_eset = pbmc.14k.eset.log2, target = genes_of_interest, feature = "geneSymbol", 
                group_by = "ClusterRes", ylabel = "log2Exp", ncol = 4)

# Marker heatmap plot
feature_heatmap(input_eset = pbmc.14k.eset.log2, target = genes_of_interest, group_name = "ClusterRes", 
                save_plot = FALSE, width = 6, height = 6, name = "log2Exp")

# Read marker genes from an excel file
markers <- read.xlsx("/Users/lding/Git/scMINER/docs/docs/tables/Immune_signatures.xlsx")

# Maker gene file has 3 required columns
head(markers)

# Draw a bubble plot
draw.marker.bbp(ref = markers, input_eset = pbmc.14k.eset.log2, width = 6, height = 4, feature = "geneSymbol",
                group_name = "ClusterRes", save_plot = FALSE)


indx <- factor(x=c("Monocyte", "CD8NaiveCTL", "NK", "Bcell", "CD4TReg", "CD4TCM", "CD4NaiveT"),
               levels=c("Monocyte", "CD8NaiveCTL", "NK", "Bcell", "CD4TReg", "CD4TCM", "CD4NaiveT"))
pbmc.14k.eset.log2$celltype <- indx[pbmc.14k.eset.log2$ClusterRes]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             


generateSJARACNeInput(
  input_eset = pbmc.14k.eset.log2, funcType = "TF", ref = "hg",  # human
  wd.src = "/Users/lding/Documents/scMINER/PBMC14k_input/SJARACNe",  # output directory
  group_name = "celltype")


generateSJARACNeInput(
  input_eset = pbmc.14k.eset.log2, funcType = "SIG", ref = "hg",  # human
  wd.src = "/Users/lding/Documents/scMINER/PBMC14k_input/SJARACNe",  # output directory
  group_name = "celltype")


