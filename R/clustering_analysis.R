
#' Generate the standard input files for MICA from sparse eset object
#'
#' @description
#' This function is used to generate the standard input files for MICA (Mutual Information-based Clustering Analysis) from a sparse eset object. It
#' supports two file formats, "**.txt**" or "**.h5ad**". To generate a "**.h5ad**" file, the "**anndata**" package is required.
#'
#' @param input_eset The sparse eset object to generate MICA input from. It must be normalized and log-transformed.
#' @param output_file The output file, or MICA input file. Should be in either "**.txt**" or "**.h5ad**" format.
#' @param overwrite Logical, whether to overwrite the output_file if it already exists. Default: `FALSE`.
#' @param downSample_N A non-negative integer or `NULL`, number of cells to downsample to. Default: `NULL`.
#' @param seed Integer or `NULL`, the seed for sampling. Default: 1. Ignored if `downSample_N` = `NULL`.
#'
#' @return A .txt or .h5ad file that can be used as the MICA input
#' @export
#'
#' @examples
#' generateMICAinput(input_eset = log2cpm.ese, output_file = "./MICA/micaInput.txt")
generateMICAinput <- function(input_eset,
                              output_file,
                              overwrite = F,
                              downSample_N = NULL,
                              seed = 1)
{
  if(is.null(output_file)) {
    if (base::exists('scminer.par'))  {
      output_file = paste0(scminer.par$out.dir.MICA, "MICAinput.txt")
    } else {
      wd <- getwd()
      output_file = paste0(wd, "/MICAinput.txt")
    }
  }

  if (length(grep("\\/", output_file)) == 0) {
    if (base::exists('scminer.par')) {
      output_file = paste0(scminer.par$out.dir.MICA, output_file)
    } else {
      wd <- getwd()
      output_file = paste0(wd, output_file)
    }
  }

  if (is.null(overwrite)) {overwrite <- FALSE}

  if (file.exists(output_file)) {
    if (!overwrite) {
      stop(paste0('The output file exists: ', output_file, '.\nUse overwrite = TRUE to overwrite it.'))
    } else {
      file.remove(output_file)
    }
  }

  exp_mat <- Biobase::exprs(input_eset)
  base::suppressWarnings(exp_mat <- as.matrix(exp_mat))
  obs_dat <- Biobase::pData(input_eset)
  var_dat <- Biobase::fData(input_eset)
  cell_size <- ncol(exp_mat)

  if (is.null(downSample_N) == FALSE) {
    if (cell_size > downSample_N) {
      set.seed(seed = seed)
      s1 <- sample(1:cell_size, size = downSample_N, replace = FALSE)
      exp_mat <- exp_mat[, s1, drop = F]
      obs_dat <- obs_dat[s1, , drop = F]
      message(sprintf('original %s cells are downsampled to %s cells for MICA input', cell_size, length(s1)))
    }
  }

  if (length(grep(".txt$", output_file)) != 0) {
    mica_input <- as.data.frame(exp_mat)
    cat("Writing MICA input to:", output_file, "\n")
    write.table(mica_input, file = output_file, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE, append = FALSE)
  } else if (length(grep(".h5ad$", output_file)) != 0) {
    avail_anndata <- require(package = "anndata", quietly = TRUE)
    if (!avail_anndata) {
      stop('The package, anndata, is not installed. Please run install.packages("anndata") to install it and re-try!')
    }

    ad <- anndata::AnnData(X = t(exp_mat), obs = obs_dat, var = var_dat)
    cat("Writing MICA input to:", output_file, "\n")
    anndata::write_h5ad(ad, output_file)
  } else {
    stop("The output file should be ended with .txt or .h5ad.\n")
  }

  if (cell_size >= 5000) {
    cat('\nFor dataset with more than 5k cells, MICA GE mode is recommended.\n')
    recommend_cmd <- sprintf('With the default resolusion:\n\tmica ge -i %s -o project_space/MICA/mica_output_dir -res 1.822 -nw 8\n\n Or with a range of resolutions:\n\tmica ge -i %s -o project_space/MICA/mica_output_dir -minr 1.0 -maxr 2.0 -ss 0.1 -nw 8', output_file, output_file)
    cat("Suggested command line is:\n\n", recommend_cmd, "\n
    Where options represent:
        -res:   the number of communities (default: 1.822)
        -minr:  the minimum number of communities to sweep (default: 1.822)
        -maxr:  the maximum number of communities to sweep (default: 1.822)
        -ss:    the step size to sweep number of communities from the mimimum to the maximum (default: 1)
        -nw:    the number of workers to run in parallel (default: 1, suggested: 25)
        -annef: the ef value of hnsw
        -annm:  the M value of hnsw
    Use mica ge -h to find more options of MICA GE mode.")
  } else {
    cat('\nFor dataset with less than 5k cells, MICA MDS mode is recommended.\n')
    recommend_cmd <- sprintf('mica mds -i %s -o project_space/MICA/mica_output_dir -nck 4 5 6 7 8 9 10 -dd 20', output_file)
    cat("Suggested command line is:\n\n", recommend_cmd, "\n
    Where options represent:
        -pn:    specifies a project name for naming the output files;
        -nck:   an array of integers delimited by a single space, where each integer specifies a k to perform a k-mean clustering;
        -dd:    can be an integer or an array of integers delimited by a single space (default is 19), it specifies the number of dimensions used in k-mean clusterings.
    Use mica mds -h to see more options of MICA MDS mode.")
  }
  #return(recommend_cmd)
}


#' Add the MICA output (cluster labels and UMAP/tSNE coordinates) to sparse eset object
#'
#' @description
#' This function is used to add the clustering results by MICA into the sparse eset object. Two types of resluts would be added to the phenoData slot of sparse eset object:
#' - cluster ID of each cell
#' - dimension reduction coordinates of each cell, either `UMAP` ot `tSNE`.
#'
#' @param input_eset The sparse eset object to add the MICA output into
#' @param mica_output_file The .txt file generated by MICA. It includes 4 columns: "ID" (Cell ID), "X" (UMAP_1 or tSNE_1), "Y" (UMAP_2 or tSNE_2), "label" (ClusterID)
#' @param visual_method Character, method used for visualizing the clustering results: `umap` (the default) or `tsne`.
#'
#' @return A sparse eset object with clustering results added.
#' @export
#'
#' @examples
#' clustered.eset <- addMICAoutput(input_eset = input_eset, mica_output_file = "/path-to-mica-input/clustering_UMAP_euclidean_20_2.22554.txt", visual_method = "umap")
addMICAoutput <- function(input_eset,
                          mica_output_file,
                          visual_method = "umap")
{
  if (!file.exists(mica_output_file)) {
    stop("The MICA output file is not found: ", mica_output_file, ".")
  }

  micaResult <- read.table(mica_output_file, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, row.names = 1)
  if (setequal(colnames(input_eset), rownames(micaResult))) {
    micaResult <- micaResult[colnames(input_eset),]
  } else {
    stop('The eset and the mica_output_file do not match. Please check and re-try.')
  }

  visual_method <- base::toupper(visual_method)
  if (visual_method == "UMAP") {
    input_eset$UMAP_1 <- micaResult[,1]; input_eset$UMAP_2 <- micaResult[,2]
  } else if (visual_method == "TSNE" | visual_method == "T-SNE") {
    input_eset$tSNE_1 <- micaResult[,1]; input_eset$tSNE_2 <- micaResult[,2]
  } else {
    stop('The visual_method can not be recognized. Should be [umap | tsne].')
  }

  input_eset$clusterID <- as.factor(micaResult$label)
  return(input_eset)
}


#' Draw a scatter plot showing the coordinates and cluster id of each cell
#'
#' @description
#' This function is used to visualize the clustering results generated by MICA.
#'
#' @param input_eset The sparse eset object
#' @param color_by Factor, character or numeric, name of the column of MICA cluster labels. Default: "`clusterID`".
#' @param colors A character vector or `NULL`, colors of the MICA cluster labels. The length of this vector should be same as the number of groups in color_by column. If `NULL`, the ggplot default colors will be use. Default: `NULL`.
#' @param do.logTransform Logical, whether to do the log2(value + 1) transformation. Only valid when color_by is numeric Default: `TRUE`.
#' @param X,Y Character, name of the columns of x-axis and y-axis coordinates. Default: "`UMAP_1`", and "`UMAP_2`".
#' @param point.size Numeric, size of points. Default: 0.3.
#' @param point.alpha Numeric, transparency of points, ranging from 0 (more transparent) to 1 (less transparent). Default: 1.
#' @param name.plot_title Character or NULL, title of the plot. Default: `NULL`.
#' @param fontsize.plot_table Numeric, font size of the title. Default: 20.
#' @param show.cluster_label Logical, whether to show labels on the plot. Ignored when color_by is numeric. Default: `TRUE`.
#' @param fontsize.cluster_label Numeric, font size of the labels. Ignored when color_by is numeric. Default: 12.
#' @param legend.position Character, position of legend: "`right`", "`left`", "`top`", "`bottom`" or "`none`". Default: "`right`".
#' @param fontsize.legend_title Integer, font size of the legend title. Default: 10.
#' @param fontsize.legend_text Integer, font size of the legend text. Default: 8.
#' @param fontsize.axis_title Integer, font size of the axis title. Default: 10.
#' @param fontsize.axis_text Integer, font size of the axis text. Default: 8.
#'
#' @return A UMAP or T-SNE plot. It also print the plot to screen.
#' @export
#'
#' @examples
#' ## 1. color-coded by factor or character variable
#' p_umap <- MICAplot(input_eset = pbmc14k_log2cpm.eset, color_by = 'clusterID', point.size = 0.1)
#'
#' ## 2. color-coded by numeric variable
#' p_umap <- MICAplot(input_eset = pbmc14k_log2cpm.eset, color_by = 'nUMI', do.logTransform = TRUE)
MICAplot <- function(input_eset,
                     color_by = "clusterID", colors = NULL, do.logTransform = TRUE,
                     X = "UMAP_1", Y = "UMAP_2",
                     point.size = 0.3, point.alpha = 1,
                     name.plot_title = NULL, fontsize.plot_table = 20,
                     show.cluster_label = TRUE, fontsize.cluster_label = 12,
                     legend.position = "right", fontsize.legend_title = 10, fontsize.legend_text = 8,
                     fontsize.axis_title = 10, fontsize.axis_text = 8)
{
  input <- Biobase::pData(input_eset)
  if ((color_by %in% colnames(input)) == FALSE) {stop('The column', colore_by, 'was not found in the input eset.')}
  if ((X %in% colnames(input)) == FALSE) {stop('The column', X, 'was not found in the input eset.')}
  if ((Y %in% colnames(input)) == FALSE) {stop('The column', Y, 'was not found in the input eset.')}

  if (is.null(colors) == FALSE) {
    if (length(unique(input[, color_by])) != length(colors)) {
      stop('The length of [colors] must be same with the group numbers of "', color_by, '", which is ', length(unique(input[, color_by])), '.')
    }
  }

  if (is.numeric(input[,color_by]) == FALSE) {
    p <- ggplot(input, aes(x = input[, X], y = input[, Y], color = input[, color_by])) +
      geom_point(size = point.size, alpha = point.alpha)

    ## add cluster label
    if (show.cluster_label == TRUE) {
      loc <- stats::aggregate(input[, c(X,Y)], by = list(GroupID = input[, color_by]), mean)
      p <- p + geom_text(data = loc, aes(x = loc[, X], y = loc[, Y], label = loc[, "GroupID"]),
                         check_overlap = TRUE, na.rm = FALSE, size = fontsize.cluster_label, show.legend = F, inherit.aes = FALSE)
    }

    # add legend
    grps <- input[, color_by]
    lgnd <- paste0(names(table(grps))," (",table(grps),")")

    if (is.null(colors) == FALSE) {
      p <- p + scale_color_manual(values = colors, labels = lgnd)
    } else {
      p <- p + scale_color_discrete(labels = lgnd)
    }

    p <- p + guides(col = guide_legend(override.aes = list(size = 10), title = paste0(color_by, "\n(",dim(input)[1] ,")")), ncol = ceiling(length(unique(grps))/10))
  } else {
    if (do.logTransform == TRUE) {
      input[, color_by] <- log2(input[, color_by] + 1)
      message('The values in "', color_by, '" have been transformed by log2(value + 1). To turn transformation off, set do.logTransform = FALSE.')
    }

    p <- ggplot(input, aes(x = input[, X], y = input[, Y], color = input[, color_by])) +
      geom_point(size = point.size, alpha = point.alpha) +
      scale_colour_gradient(low = "lightgrey", high = "blue")
  }


  p <- p + theme_classic() +
    theme(
      legend.position = legend.position,
      legend.title = element_text(size = fontsize.legend_title, face = "bold"),
      legend.text = element_text(size = fontsize.legend_text),
      axis.title = element_text(size = fontsize.axis_title, face = "bold", hjust = 0.5, color = "black"),
      axis.text.x = element_text(size = fontsize.axis_text, hjust = 0.5, color = "black"),
      axis.text.y = element_text(size = fontsize.axis_text, hjust = 1, color = "black"))

  return(p)
}


