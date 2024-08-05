
#' Define the class: SparseExpressionSet
#'
#' @title SparseExpressionSet
#' @exportClass SparseExpressionSet
#' @import Biobase
methods::setClass(Class = "SparseExpressionSet",
                  contains = "ExpressionSet",
                  prototype = methods::prototype(methods::new("VersionedBiobase", versions = c(Biobase::classVersion("ExpressionSet"), SparseExpressionSet = "1.0.0" )))
)


#' Create a sparse expression set object from a data matrix
#'
#' @description
#' This function is used to create a pre-defined sparse expression set object from a data matrix of differnt classes: **`"dgCMatrix"`**, **`"dgTMatrix"`**, **`"dgeMatrix"`**,
#' **`"matrix"`**, **`"data.frame"`**. It allows the users to provide self-customized meta data for both cells (parameter **`cellData`**) and genes (parameter **`featureData`**).
#' It can also generate the meta data for both automatically, if **`addMetaData`** = `TRUE`. The automatically generated meta data includes:
#' - **"nUMI"**: number of total UMIs in each cell, only valid when the values in data matrix are raw UMI counts;
#' - **"nFeature"**: number of expressed features/genes in each cell;
#' - **"pctMito"**: percentage of UMIs of mitochondrial genes (defined by "^mt-|^MT-") in each cell;
#' - **"pctSpikeIn"**: percentage of UMIs of spike-in RNAs (defined by "^ERCC-|^Ercc-") in each cell;
#' - **"nCell"**: number of cells that each feature/gene was identified in.
#'
#' @param input_matrix A data matrix with Features/Genes as the rows and Cells as the columns. It should be one of: 'dgCMatrix', 'dgTMatrix', 'dgeMatrix', 'matrix', 'data.frame'.
#' @param do.sparseConversion Logical, whether to convert the **`input_matrix`** to a sparse matrix if it's not. Default: `TRUE`.
#' @param cellData A data frame containing meta data of cells or `NULL`. It's row.names should be consistent with the colnames of **`input_matrix`**. Default: `NULL`.
#' @param featureData A data frame containing meata data of features or `NULL`. It's row.names should be consistent with the row.names of **`input_matrix`**. Default: `NULL`.
#' @param annotation Character, a character describing the project properties. It's highly recommended to use the path to project space. Default: "".
#' @param projectID Character or `NULL`, the project name of the sparse eset object. Default: `NULL`.
#' @param addMetaData Logical, whether to calculate and add extra statistics (a.k.a. meta data) to cells and features. Default: `TRUE`.
#'
#' @return A sparse eset object with three slot: 1) gene by cell matrix; 2) data frame of cell information; 3) data frame of feature/gene information.
#' @export
#'
#' @examples
#' expression_raw.eset <- createSparseEset(input_matrix = sparseMatrix, projectID = "demoSample", addMetaData = T)
createSparseEset <- function(input_matrix,
                             do.sparseConversion = TRUE,
                             cellData = NULL,
                             featureData = NULL,
                             annotation = "",
                             projectID = NULL,
                             addMetaData = TRUE)
{
  ## data matrix
  if (class(input_matrix)[1] %in% c("dgCMatrix", "dgTMatrix", "dgeMatrix", "matrix", "data.frame")) {
    if (class(input_matrix)[1] %in% c("matrix", "data.frame") & do.sparseConversion == TRUE) {
      expression_data <- Matrix::Matrix(as.matrix(input_matrix), sparse = TRUE)
    } else {
      expression_data <- input_matrix
    }
  } else {
    stop("The format of input_matrix is not supported!\nscMINER supports: 'dgCMatrix', 'dgTMatrix', 'dgeMatrix', 'matrix', 'data.frame'. Use class() to check.")
  }

  cat("Creating sparse eset from the input_matrix ...\n")
  ## feature data
  if (is.null(featureData) == TRUE) {
    feature_data <- data.frame(row.names = rownames(expression_data), GeneSymbol = rownames(expression_data), stringsAsFactors = F)
  } else {
    if_involved <- row.names(expression_data) %in% row.names(featureData)
    if (all(if_involved) == TRUE) {
      feature_data <- featureData[row.names(expression_data), , drop = FALSE]
    } else {
      stop("Some features of the input_matrix were not found in the featureData: ", paste0(row.names(expression_data)[!if_involved], collapse = ", "), ".\n")
    }
  }

  ## phenotype data
  if (is.null(cellData) == TRUE) {
    cell_data <- data.frame(row.names = colnames(expression_data), CellID = colnames(expression_data), stringsAsFactors = F)
  } else {
    if_involved <- colnames(expression_data) %in% row.names(cellData)
    if (all(if_involved) == TRUE) {
      cell_data <- cellData[colnames(expression_data), , drop = FALSE]
    } else {
      stop("Some cells of the input_matrix were not found in the cellData: ", paste0(colnames(expression_data)[!if_involved], collapse = ", "), ".\n")
    }
  }

  if (is.null(annotation) == FALSE) {annotation <- as.character(annotation)}

  if (is.null(projectID) == FALSE) {cell_data$projectID <- projectID}

  ## meta data of each gene or cell
  if (addMetaData == TRUE) {
    cat("\tAdding meta data based on input_matrix ...\n")

    if_integers <- expression_data %% 1 == 0
    if (all(if_integers) == FALSE) {
      warning("The input_matrix might not be a raw count matrix, since non-integer values were found. The 'nUMI' column is only valid when the values in input_matrix are raw UMI counts.")
    }

    if ("nCell" %in% colnames(feature_data)) {warning("The 'nCell' column in featureData will be updated based on input_matrix To keep it, please rename it in featureData.")}
    feature_data$nCell <- Matrix::rowSums(expression_data != 0) # calculate number of cells each gene is identified in

    if ("nUMI" %in% colnames(cell_data)) {warning("The 'nUMI' column in cellData will be updated based on input_matrix. To keep it, please rename it in cellData.")}
    cell_data$nUMI <- Matrix::colSums(expression_data)

    if ("nFeature" %in% colnames(cell_data)) {warning("The 'nFeature' column in cellData will be updated based on input_matrix. To keep it, please rename it in cellData.")}
    cell_data$nFeature <- Matrix::colSums(expression_data != 0)

    if ("pctMito" %in% colnames(cell_data)) {warning("The 'pctMito' column in cellData will be updated based on input_matrix. To keep it, please rename it in cellData.")}
    mito_genes <- row.names(expression_data)[grepl(pattern = "^mt-|^MT-", x = row.names(expression_data))]
    cell_data$pctMito = round(Matrix::colSums(expression_data[row.names(expression_data) %in% mito_genes, ]) / Matrix::colSums(expression_data), 8)

    if ("pctSpikeIn" %in% colnames(cell_data)) {warning("The 'pctSpikeIn' column in cellData will be updated based on input_matrix. To keep it, please rename it in cellData.")}
    spikeIn_genes <- row.names(expression_data)[grepl(pattern = "^ERCC-|^Ercc-", x = row.names(expression_data))]
    cell_data$pctSpikeIn = round(Matrix::colSums(expression_data[row.names(expression_data) %in% spikeIn_genes, ]) / Matrix::colSums(expression_data), 8)

    if (("CellID" %in% colnames(cell_data)) == FALSE) {cell_data$CellID <- row.names(cell_data)}
  }

  eset <- new( "SparseExpressionSet",
               assayData = assayDataNew("environment", exprs = expression_data),
               phenoData = new("AnnotatedDataFrame", data = cell_data),
               featureData = new("AnnotatedDataFrame", data = feature_data),
               annotation = annotation)
  cat('Done! The sparse eset has been generated:', dim(eset)[1], 'genes,', dim(eset)[2], 'cells.\n')

  return(eset)
}


#' Combine multiple sparse expression set objects
#'
#' @description
#' This function is used to combine the sparse expression set objects. The combined eset object contains all cells and features of all input eset objects. If the eset objects are of
#' different features, `NA` values will be generated and automatically imputed by the minimum value of the combined gene expression matrix.
#'
#' @param eset_list A vector of sparse expression set objects to combine
#' @param projectID A character vector or `NULL`, set the project names of the eset objects to combine. Default: `NULL`.
#' @param addPrefix A character vector or `NULL`, add a **prefix** to the cell barcodes of each eset object to combine. It is highly recommended to use a prefix containing letters and/or numbers only, and not starting with numbers. Default: `NULL`.
#' @param addSurfix A character vector or `NULL`, add a **surfix** to the cell barcodes of each eset object to combine. It is highly recommended to use a surfix containing letters and/or numbers only, and not starting with numbers. Default: `NULL`.
#' @param addMetaData Logical, whether to update the meta data of cells and features after combination. Default: `TRUE`.
#' @param imputeNA Logical, whether to impute NA values in combined matrix. If `TRUE`, the min value of the matrix will be used to replace the NAs. If `FALSE`, NA values will retain. Default: `TRUE`.
#'
#' @return A sparse eset object with combined features and cells of multiple eset objects.
#' @export
#'
#' @examples
#' combined.eset <- combineSparseEset(c(sample_1.eset, sample_2.eset, sample_3.eset), projectID = c("sample1", "sample2", "sample3"), addPrefix = c("tag1", "tag2", "tag3"), addMetaData = TRUE)
combineSparseEset <- function(eset_list,
                              projectID = NULL,
                              addPrefix = NULL,
                              addSurfix = NULL,
                              addMetaData = TRUE,
                              imputeNA = TRUE)
{
  l_eset <- eset_list
  if (length(l_eset) <= 1) {stop("At least 2 esets are required for eset_list")}

  if (is.null(projectID) == FALSE) {
    if (length(projectID) != length(l_eset)) {stop("The length of projectID is different from eset_list")}
  }

  ## add cellPrefix and/or cellSurfix
  if (is.null(addPrefix) == FALSE) {
    if (length(addPrefix) == length(l_eset)) {
      for (i in 1:length(l_eset)) {
        colnames(l_eset[[i]]) <- paste0(addPrefix[i], "_", colnames(l_eset[[i]]))
      }
    } else {
      stop("The length of addPrefix is different from eset_list")
    }
  }

  if (is.null(addSurfix) == FALSE) {
    if (length(addSurfix) == length(l_eset)) {
      for (i in 1:length(l_eset)) {
        colnames(l_eset[[i]]) <- paste0(colnames(l_eset[[i]]), "_", addSurfix[i])
      }
    } else {
      stop("The length of addSurfix is different from eset_list")
    }
  }

  cat("Combining the input sparse eSets ...\n")
  ## merge matrix, barcodes, features
  if_duplicatedBarcodes <- lapply(l_eset, colnames) %>% unlist() %>% duplicated() %>% any()
  if (if_duplicatedBarcodes == TRUE) {
    stop("Duplicated barcode IDs were found among the esets. It's highly recommended to use addPrefix and/or addSurfix to avoid duplicated barcodes.")
  } else {
    for (i in 1:length(l_eset)) {
      suppressWarnings(exp_tmp <- as.matrix(exprs(l_eset[[i]])))
      pd_tmp <- pData(l_eset[[i]])

      if (is.null(projectID) == FALSE) {pd_tmp$projectID <- projectID[i]}

      if (i == 1) {
        exp_merged <- exp_tmp
        pd_merged <- pd_tmp
      } else {
        exp_merged <- base::merge(exp_merged, exp_tmp, by = "row.names", all = T)
        row.names(exp_merged) <- exp_merged[,1]
        exp_merged <- exp_merged[,-1]

        colnames_shared <- intersect(colnames(pd_merged), colnames(pd_tmp))
        if (length(colnames_shared) < length(colnames(pd_merged)) | length(colnames_shared) < length(colnames(pd_tmp))) {
          cat("Different colnames of phenotype data were found. Only the shared ones were kept.\n")
        }
        pd_merged <- pd_merged[, colnames_shared]
        pd_tmp <- pd_tmp[, colnames_shared]
        pd_merged <- rbind(pd_merged, pd_tmp)
      }
    }
  }

  fd_merged <- data.frame(row.names = row.names(exp_merged), geneSymbol = row.names(exp_merged), nCell = Matrix::rowSums(exp_merged != 0))

  if (any(is.na(exp_merged))) {
    if (imputeNA == TRUE) {
      min_v <- min(exp_merged, na.rm = TRUE)
      exp_merged[is.na(exp_merged)] <- min_v
      cat('NA values were found in the merged matrix and have been replaced by the minimum value: ', min_v, '.\n')
    } else {
      cat('NA values were found and retained. To impute NAs, please set imputeNA = TRUE.\n')
    }
  } else {
    cat('No NA values was found in the combined gene expression matrix!\n')
  }
  exp_merged <- Matrix::Matrix(as.matrix(exp_merged), sparse = TRUE)

  ## meta data of each gene or cell
  if (addMetaData == TRUE) {
    cat("Adding meta data based on merged data matrix ...\n")

    if_integers <- exp_merged %% 1 == 0
    if (all(if_integers) == FALSE) {
      warning("The data matrix might not be a raw count matrix, since non-integer values were found. The 'nUMI' column is only valid when the values in data matrix are raw UMI counts.")
    }

    fd_merged$nCell <- Matrix::rowSums(exp_merged != 0) # calculate number of cells each gene is identified in

    pd_merged$nUMI <- Matrix::colSums(exp_merged)
    pd_merged$nFeature <- Matrix::colSums(exp_merged != 0)

    mito_genes <- row.names(exp_merged)[grepl(pattern = "^mt-|^MT-", x = row.names(exp_merged))]
    pd_merged$pctMito = round(Matrix::colSums(exp_merged[mito_genes, ]) / Matrix::colSums(exp_merged), 8)

    spikeIn_genes <- row.names(exp_merged)[grepl(pattern = "^ERCC-|^Ercc-", x = row.names(exp_merged))]
    pd_merged$pctSpikeIn = round(Matrix::colSums(exp_merged[spikeIn_genes, ]) / Matrix::colSums(exp_merged), 8)

    pd_merged$CellID <- row.names(pd_merged)
  }

  eset <- new( "SparseExpressionSet",
               assayData = assayDataNew("environment", exprs = exp_merged),
               phenoData = new("AnnotatedDataFrame", data = pd_merged),
               featureData= new("AnnotatedDataFrame", data = fd_merged))
  cat('Done! The combined sparse eset has been generated:', dim(eset)[1], 'genes,', dim(eset)[2], 'cells.\n')

  return(eset)
}


#' Update the slots and/or meta data of the sparse eset object
#'
#' @description
#' This function is used to update the three slots ('**assayData**', '**phenoData**', '**featureData**') and/or '**meta data**' of sparse eset object.
#'
#' @param input_eset The sparse eset object to update
#' @param dataMatrix A data matrix with Features/Genes as the rows and Cells as the columns. It's row.names and colnames must be consistent with the **`input_eset`**. Default: `NULL`.
#' @param cellData A data frame containing meta data of cells or `NULL`. It's row.names should be consistent with the colnames of **`input_eset`**. Default: `NULL`.
#' @param featureData A data frame containing meata data of features or `NULL`. It's row.names should be consistent with the row.names of **`input_eset`**. Default: `NULL`.
#' @param addMetaData Logical, whether to update the meta data of features and cells based on the expression matrix. Default: `FALSE`.
#'
#' @return A sparse eset object with updated information
#' @import Biobase
#' @export
#'
#' @examples
#' updated.eset <- updateSparseEset(input_eset = input.eset, cellData = data.frame(pData(input.eset), cellType = "B_cells"), addMetaData = TRUE)
updateSparseEset <- function(input_eset,
                             dataMatrix = NULL,
                             cellData = NULL,
                             featureData = NULL,
                             addMetaData = FALSE)
{
  if (is.null(dataMatrix) == FALSE) {
    if (setequal(row.names(dataMatrix), row.names(input_eset))) {
      dataMatrix <- dataMatrix[row.names(input_eset), , drop = FALSE]
    } else {
      stop("The rownames (Features) of dataMatrix and input_eset do NOT match.")
    }

    if (setequal(colnames(dataMatrix), colnames(input_eset))) {
      dataMatrix <- dataMatrix[,colnames(input_eset)]
    } else {
      stop("The colnames (Cells) of dataMatrix and input_eset do NOT match.")
    }

    input_eset <- new( "SparseExpressionSet",
                 assayData = assayDataNew("environment", exprs = dataMatrix),
                 phenoData = new("AnnotatedDataFrame", data = Biobase::pData(input_eset)),
                 featureData= new("AnnotatedDataFrame", data = Biobase::fData(input_eset)))
    cat("The data matrix of input_eset has been updated!\n")
  }

  if (is.null(cellData) == TRUE) {
    cell_data <- data.frame(row.names = colnames(input_eset), CellID = colnames(input_eset), stringsAsFactors = F)
  } else {
    if_involved <- colnames(input_eset) %in% row.names(cellData)
    if (all(if_involved) == TRUE) {
      cell_data <- cellData[colnames(input_eset), , drop = FALSE]
    } else {
      stop("Some cells of the dataMatrix were not found in the cellData: ", paste0(colnames(input_eset)[if_involved], collapse = ", "), ".\n")
    }
  }

  if (is.null(cellData) == FALSE) {
    if_involved <- colnames(input_eset) %in% row.names(cellData)
    if (all(if_involved) == TRUE) {
      pData(input_eset) <- cellData[colnames(input_eset), , drop = FALSE]
      cat("The cell information of input_eset has been updated!\n")
    } else {
      stop("Some cells of the input_eset were not found in the cellData: ", paste0(colnames(input_eset)[!if_involved], collapse = ", "), ". Please check and re-try.")
    }
  }

  if (is.null(featureData) == FALSE) {
    if_involved <- row.names(input_eset) %in% row.names(featureData)
    if (all(if_involved) == TRUE) {
      fData(input_eset) <- featureData[row.names(input_eset), , drop = FALSE]
      cat("The feature information of input_eset has been updated!\n")
    } else {
      stop("Some features of the input_eset were not found in the featureData: ", paste0(row.names(input_eset)[!if_involved], collapse = ", "), ". Please check and re-try.")
    }
  }

  if (addMetaData == TRUE) {
    cat("Updating meta data based on data matrix ...\n")

    if_integers <- exprs(input_eset) %% 1 == 0
    if (all(if_integers) == FALSE) {
      warning("The data matrix might not be a raw count matrix, since non-integer values were found. The 'nUMI' column is only valid when the values in data matrix are raw UMI counts.")
    }

    Biobase::fData(input_eset)$nCell <- Matrix::rowSums(exprs(input_eset) != 0) # calculate number of cells each gene is identified in

    Biobase::pData(input_eset)$nUMI <- Matrix::colSums(exprs(input_eset))
    Biobase::pData(input_eset)$nFeature <- Matrix::colSums(exprs(input_eset) != 0)

    mito_genes <- row.names(exprs(input_eset))[grepl(pattern = "^mt-|^MT-", x = row.names(exprs(input_eset)))]
    Biobase::pData(input_eset)$pctMito = round(Matrix::colSums(exprs(input_eset)[mito_genes, ]) / Matrix::colSums(exprs(input_eset)), 8)

    spikeIn_genes <- row.names(exprs(input_eset))[grepl(pattern = "^ERCC-|^Ercc-", x = row.names(exprs(input_eset)))]
    Biobase::pData(input_eset)$pctSpikeIn = round(Matrix::colSums(exprs(input_eset)[spikeIn_genes, ]) / Matrix::colSums(exprs(input_eset)), 8)

    if (("CellID" %in% colnames(Biobase::pData(input_eset))) == FALSE) {Biobase::pData(input_eset)$CellID <- colnames(input_eset)}
  }

  return(input_eset)
}


#' Filter the cells and/or features of sparse eset object using automatic or self-customized cutoffs
#'
#' @description
#' This function is used to remove the cells and features of low quality. It provides two modes to define the cutoffs:
#' - **"auto"**: in this mode, scMINER will estimate the cutoffs based on Median ± 3*MAD (maximum absolute deviation). This mode works well for the matrix of raw UMI counts or TPM (Transcripts Per Million) values.
#' - **"manual"**: in this mode, the users can manually specify the cutoffs, both low and high, of all 5 metrics: **nUMI**, **nFeature**, **pctMito**, **pctSpikeIn** for cells, and **nCell** for genes. No cells or
#' features would be removed under the default cutoffs of each metrics.
#'
#' @param input_eset The sparse eset object to be filtered
#' @param filter_mode Character, mode to apply the filtration cutoffs: **"auto"** (the default) or **"manual"**.
#' @param filter_type Character, objective type to be filtered: **"both"** (the default), **"cell"** or **"feature"** .
#' @param gene.nCell_min Numeric, the minimum number of cells that the qualified genes are identified in. Default: 1.
#' @param gene.nCell_max Numeric, the maximum number of cells that the qualified genes are identified in. Default: Inf.
#' @param cell.nUMI_min Numeric, the minimum number of total UMI counts per cell that the qualified cells carry. Default: 1.
#' @param cell.nUMI_max Numeric, the maximum number of total UMI counts per cell that the qualified cells carry. Default: Inf.
#' @param cell.nFeature_min Numeric, the minimum number of non-zero Features per cell that the qualified cells carry. Default: 1.
#' @param cell.nFeature_max Numeric, the maximum number of non-zero Features per cell that the qualified cells carry. Default: Inf.
#' @param cell.pctMito_min Numeric, the minimum percentage of UMI counts of mitochondrial genes that the qualified cells carry. Default: 0.
#' @param cell.pctMito_max Numeric, the maximum percentage of UMI counts of mitochondrial genes that the qualified cells carry. Default: 1.
#' @param cell.pctSpikeIn_min Numeric, the minimum percentage of UMI counts of spike-in that the qualified cells carry. Default: 0.
#' @param cell.pctSpikeIn_max Numeric, the maximum percentage of UMI counts of spike-in that the qualified cells carry. Default: 1.
#'
#' @return A filtered sparse eset object. It also prints the summary of filtration to the screen.
#' @export
#'
#' @examples
#' filtered.eset <- filterSparseEset(raw.eset) ## filter the input eset using the cutoffs calculated by scMINER.
#' filtered.eset <- filterSparseEset(raw.eset, gene.nCell_min = 10, cell.nUMI_min = 500, cell.nFeature_min = 100, cell.nFeature_max = 5000, cell.pctMito_max = 0.15)
filterSparseEset <- function(input_eset,
                             filter_mode = "auto",
                             filter_type = "both",
                             gene.nCell_min = NULL, gene.nCell_max = NULL,
                             cell.nUMI_min = NULL, cell.nUMI_max = NULL,
                             cell.nFeature_min = NULL, cell.nFeature_max = NULL,
                             cell.pctMito_min = NULL, cell.pctMito_max = NULL,
                             cell.pctSpikeIn_min = NULL, cell.pctSpikeIn_max = NULL)
                             #gene.nCell_min = 1, gene.nCell_max = Inf,
                             #cell.nUMI_min = 1, cell.nUMI_max = Inf,
                             #cell.nFeature_min = 1, cell.nFeature_max = Inf,
                             #cell.pctMito_min = 0, cell.pctMito_max = 1,
                             #cell.pctSpikeIn_min = 0, cell.pctSpikeIn_max = 1)
{
  ## check parameters
  if ((filter_mode %in% c("auto", "manual")) == FALSE) {stop('The filter_mode was not recognized. Please specify one of them: [auto | manual].')}
  if ((filter_type %in% c("both", "cell", "feature")) == FALSE) {stop('The filter_type was not recognized. Please specify one of them: [both | cell | feature].')}

  ## check the availability of 5 QC metrics
  cat("Checking the availability of the 5 metrics ('nCell', 'nUMI', 'nFeature', 'pctMito', 'pctSpikeIn') used for filtration ...\n")
  if (("nCell" %in% colnames(Biobase::fData(input_eset))) & all(c("nUMI", "nFeature", "pctMito", "pctSpikeIn") %in% colnames(Biobase::pData(input_eset)))) {
    cat("Checking passed! All 5 metrics are available.\n")
  } else {
    cat("Part of the 5 metrics are not available. Updating the input eset automatically...\n")
    input_eset <- updateSparseEset(input_eset, addMetaData = T)
    cat("Checking passed! All 5 metrics are available after updating.\n")
  }

  ## prepare the cutoffs
  if (filter_mode == "auto") {
    if (is.null(gene.nCell_min) == TRUE) {gene.nCell_min = max(floor(0.005 * dim(input_eset)[2]), 1)}
    if (is.null(gene.nCell_max) == TRUE) {gene.nCell_max = Inf}
    if (is.null(cell.nUMI_min) == TRUE) {cell.nUMI_min = max(floor(exp(median(log(input_eset$nUMI)) - 3 * mad(log(input_eset$nUMI)))), 100)}
    if (is.null(cell.nUMI_max) == TRUE) {cell.nUMI_max = ceiling(exp(median(log(input_eset$nUMI)) + 3 * mad(log(input_eset$nUMI))))}
    if (is.null(cell.nFeature_min) == TRUE) {cell.nFeature_min = max(floor(exp(median(log(input_eset$nFeature)) - 3 * mad(log(input_eset$nFeature)))), 50)}
    if (is.null(cell.nFeature_max) == TRUE) {cell.nFeature_max = Inf}
    if (is.null(cell.pctMito_min) == TRUE) {cell.pctMito_min = 0}
    if (is.null(cell.pctMito_max) == TRUE) {cell.pctMito_max = round(median(input_eset$pctMito) + 3 * mad(input_eset$pctMito), 4)}
    if (is.null(cell.pctSpikeIn_min) == TRUE) {cell.pctSpikeIn_min = 0}
    if (is.null(cell.pctSpikeIn_max) == TRUE) {cell.pctSpikeIn_max = round(median(input_eset$pctSpikeIn) + 3 * mad(input_eset$pctSpikeIn), 4)}
  } else if (filter_mode == "manual") {
    if (is.null(gene.nCell_min) == TRUE) {gene.nCell_min = 1}
    if (is.null(gene.nCell_max) == TRUE) {gene.nCell_max = Inf}
    if (is.null(cell.nUMI_min) == TRUE) {cell.nUMI_min = 1}
    if (is.null(cell.nUMI_max) == TRUE) {cell.nUMI_max = Inf}
    if (is.null(cell.nFeature_min) == TRUE) {cell.nFeature_min = 1}
    if (is.null(cell.nFeature_max) == TRUE) {cell.nFeature_max = Inf}
    if (is.null(cell.pctMito_min) == TRUE) {cell.pctMito_min = 0}
    if (is.null(cell.pctMito_max) == TRUE) {cell.pctMito_max = 1}
    if (is.null(cell.pctSpikeIn_min) == TRUE) {cell.pctSpikeIn_min = 0}
    if (is.null(cell.pctSpikeIn_max) == TRUE) {cell.pctSpikeIn_max = 1}
  }

  ## filtering: genes
  fd_pre <- Biobase::fData(input_eset)
  if (filter_type %in% c("both", "feature")) {
    gene_qualified <- row.names(input_eset)[fd_pre$nCell >= gene.nCell_min & fd_pre$nCell <= gene.nCell_max]
  } else {
    gene_qualified <- row.names(input_eset)
  }
  msg_gene <- sprintf("\tGene filtration statistics:\n\t\tMetrics\t\tnCell\n\t\tCutoff_Low\t%.0f\n\t\tCutoff_High\t%f\n\t\tGene_total\t%d\n\t\tGene_passed\t%d(%.2f%%)\n\t\tGene_failed\t%d(%.2f%%)\n",
                      gene.nCell_min, gene.nCell_max, nrow(input_eset), length(gene_qualified), (length(gene_qualified)/nrow(input_eset)*100), (nrow(input_eset)-length(gene_qualified)), ((nrow(input_eset)-length(gene_qualified))/nrow(input_eset)*100))

  ## filtering: cells
  if (filter_type %in% c("both", "cell")) {
    cell_qualified.nUMI <- which((input_eset$nUMI >= cell.nUMI_min) & (input_eset$nUMI <= cell.nUMI_max))
    cell_qualified.nFeature <- which((input_eset$nFeature >= cell.nFeature_min) & (input_eset$nFeature <= cell.nFeature_max)); cell_qualified.combined <- intersect(cell_qualified.nUMI, cell_qualified.nFeature)
    cell_qualified.pctMito <- which((input_eset$pctMito >= cell.pctMito_min) & (input_eset$pctMito <= cell.pctMito_max)); cell_qualified.combined <- intersect(cell_qualified.combined, cell_qualified.pctMito)
    cell_qualified.pctSpikeIn <- which((input_eset$pctSpikeIn >= cell.pctSpikeIn_min) & (input_eset$pctSpikeIn <= cell.pctSpikeIn_max)); cell_qualified.combined <- intersect(cell_qualified.combined, cell_qualified.pctSpikeIn)
    cell_qualified <- colnames(input_eset)[cell_qualified.combined]
  } else {
    cell_qualified.nUMI <- which((input_eset$nUMI >= cell.nUMI_min) & (input_eset$nUMI <= cell.nUMI_max))
    cell_qualified.nFeature <- which((input_eset$nFeature >= cell.nFeature_min) & (input_eset$nFeature <= cell.nFeature_max)); cell_qualified.combined <- intersect(cell_qualified.nUMI, cell_qualified.nFeature)
    cell_qualified.pctMito <- which((input_eset$pctMito >= cell.pctMito_min) & (input_eset$pctMito <= cell.pctMito_max)); cell_qualified.combined <- intersect(cell_qualified.combined, cell_qualified.pctMito)
    cell_qualified.pctSpikeIn <- which((input_eset$pctSpikeIn >= cell.pctSpikeIn_min) & (input_eset$pctSpikeIn <= cell.pctSpikeIn_max)); cell_qualified.combined <- intersect(cell_qualified.combined, cell_qualified.pctSpikeIn)
    cell_qualified <- colnames(input_eset)
  }

  msg_cell <- sprintf("\tCell filtration statistics:\n\t\tMetrics\t\tnUMI\t\tnFeature\tpctMito\t\tpctSpikeIn\tCombined\n\t\tCutoff_Low\t%.0f\t\t%.0f\t\t%.0f\t\t%.0f\t\t%s\n\t\tCutoff_High\t%.0f\t\t%.0f\t\t%.4f\t\t%.4f\t\t%s\n\t\tCell_total\t%d\t\t%d\t\t%d\t\t%d\t\t%d\n\t\tCell_passed\t%d(%.2f%%)\t%d(%.2f%%)\t%d(%.2f%%)\t%d(%.2f%%)\t%d(%.2f%%)\n\t\tCell_failed\t%d(%.2f%%)\t%d(%.2f%%)\t%d(%.2f%%)\t%d(%.2f%%)\t%d(%.2f%%)\n",
                  cell.nUMI_min, cell.nFeature_min, cell.pctMito_min, cell.pctSpikeIn_min, "NA",
                  cell.nUMI_max, cell.nFeature_max, cell.pctMito_max, cell.pctSpikeIn_max, "NA",
                  ncol(input_eset), ncol(input_eset), ncol(input_eset), ncol(input_eset), ncol(input_eset),
                  length(cell_qualified.nUMI), (length(cell_qualified.nUMI)/ncol(input_eset)*100),
                  length(cell_qualified.nFeature), (length(cell_qualified.nFeature)/ncol(input_eset)*100),
                  length(cell_qualified.pctMito), (length(cell_qualified.pctMito)/ncol(input_eset)*100),
                  length(cell_qualified.pctSpikeIn), (length(cell_qualified.pctSpikeIn)/ncol(input_eset)*100),
                  length(cell_qualified.combined), (length(cell_qualified.combined)/ncol(input_eset)*100),
                  (ncol(input_eset)-length(cell_qualified.nUMI)), ((ncol(input_eset)-length(cell_qualified.nUMI))/ncol(input_eset)*100),
                  (ncol(input_eset)-length(cell_qualified.nFeature)), ((ncol(input_eset)-length(cell_qualified.nFeature))/ncol(input_eset)*100),
                  (ncol(input_eset)-length(cell_qualified.pctMito)), ((ncol(input_eset)-length(cell_qualified.pctMito))/ncol(input_eset)*100),
                  (ncol(input_eset)-length(cell_qualified.pctSpikeIn)), ((ncol(input_eset)-length(cell_qualified.pctSpikeIn))/ncol(input_eset)*100),
                  (ncol(input_eset)-length(cell_qualified.combined)), ((ncol(input_eset)-length(cell_qualified.combined))/ncol(input_eset)*100))

  eset_filtered <- input_eset[gene_qualified, cell_qualified]
  cat("Filtration is done!\n")

  message("Filtration Summary:")
  message("\t", nrow(eset_filtered), "/", nrow(input_eset), " genes passed!")
  message("\t", ncol(eset_filtered), "/", ncol(input_eset), " cells passed!")
  if (filter_type == "both") {
    message("\nFor more details:\n", msg_gene, "\n", msg_cell)
  } else if (filter_type == "cell") {
    message("\nFor more details:\n", msg_cell)
  } else if (filter_type == "feature") {
    message("\nFor more details:\n", msg_gene)
  }

  return(eset_filtered)
}


#' Normalize and log-transform the sparse eset object
#'
#' @description
#' This function is used to normalize and log-transform the sparse eset object. The default method is "log21p".
#'
#' @param input_eset The sparse eset object for normalization
#' @param scale_factor Numeric, the library size to normalize to. Default: 1000000.
#' @param do.logTransform Logical, whether to do log-transformation. Default: `TRUE`.
#' @param log_base Numeric, the base of log-transformation. Usually **2**, **exp(1)** or **10**. Default: 2.
#' @param log_pseudoCount Numeric, the pseudo count to add to avoid "`-Inf`" in log-transformation. Default: 1.
#'
#' @return A sparse eset object that has been normalized and log-transformed
#' @export
#'
#' @examples normalized.eset <- normalizeSparseEset(input_eset = filtered.eset, scale_factor = 1000000, do.logTransform = TRUE)
normalizeSparseEset <- function(input_eset,
                                scale_factor = 1000000,
                                do.logTransform = TRUE,
                                log_base = 2,
                                log_pseudoCount = 1
                                )
{
  exp_mat <- Biobase::exprs(input_eset)

  suppressWarnings(exp_mat.normalized <- base::sweep(exp_mat, 2, scale_factor/unname(Matrix::colSums(exp_mat)), '*'))

  if (do.logTransform == TRUE) {
    base::suppressWarnings(exp_mat.normalized <- log(exp_mat.normalized + log_pseudoCount, base = log_base))
    cat("Done! The data matrix of eset has been normalized and log-transformed!\n")
  } else {
    cat("Done! The data matrix of eset has been normalized but NOT log-transformed!\n")
  }

  eset <- new( "SparseExpressionSet",
               assayData = assayDataNew("environment", exprs = exp_mat.normalized),
               phenoData = new("AnnotatedDataFrame", data = Biobase::pData(input_eset)),
               featureData= new("AnnotatedDataFrame", data = Biobase::fData(input_eset)))
  cat('The returned eset contains:', dim(eset)[1], 'genes,', dim(eset)[2], 'cells.\n')

  return(eset)
}

#' Generate a quality control report from sparse eset object
#'
#' @description
#' This function is used to generate a html quality control report from a sparse eset object. Compared with the summary table return by filterSparseEset(), the output report
#' contains more comprehensive and detailed QC results and can be used to estimate the cutoffs to filter the eset object. It also contains some plots for presentation purpose.
#'
#' @param input_eset The sparse eset object for quality control analysis
#' @param output_html_file The path of the output .html file
#' @param overwrite Logical, whether to overwrite the output .html file if it already exists. Default: `FALSE`.
#' @param group_by Character or `NULL`, Name of column in pData(eset) used for grouping. Default: `NULL`.
#'
#' @return A html-formatted quality control report of sparse eset object
#' @export
#'
#' @examples
#' drawSparseEsetQC(input_eset, output_html_file = "./QC/esetQCreport.html", overwrite = FALSE, group = "projectID")
drawSparseEsetQC <- function(input_eset,
                             output_html_file,
                             overwrite = FALSE,
                             group_by = NULL)
{
  ## check the output file
  cat("Checkinig the output html file ...\n")
  if (grepl("html$", output_html_file, ignore.case = TRUE) == FALSE) {stop('The output_html_file should be in .html format.')}

  if (file.exists(output_html_file) == TRUE) {
    if (overwrite == TRUE) {
      cat("\tThe output_html_file", output_html_file, "already exists, and will be overwritten sincel overwrite is set TRUE.\n")
      file.remove(output_html_file)
    } else {
      stop("The output_html_file", output_html_file, "already exists. Please set overwrite to TRUE to overwrite it, or use another name.\n")
    }
  } else {
    output_dir <- dirname(output_html_file)
    if (dir.exists(output_dir) == FALSE) {
      dir.create(output_dir, recursive = TRUE)
      cat("\tThe directory of output html file does not exists. Now it has been created!\n")
    }
  }

  ## check the input eset
  if (is.null(group_by) == FALSE) {
    cat("Checkinig the group_by information in eset ...\n")
    if (group_by %in% colnames(Biobase::pData(input_eset))) {
      grps <- unique(Biobase::pData(input_eset)[, group_by])
      cat("\t", length(grps), "groups are found in the input eset!\n")
    } else {
      stop("The group_by was not found in the input eset. Please check and re-try.")
    }
  }

  cat("Checking the availability of the 5 metrics ('nCell', 'nUMI', 'nFeature', 'pctMito', 'pctSpikeIn') used for quanlity control ...\n")
  if (("nCell" %in% colnames(Biobase::fData(input_eset))) & all(c("nUMI", "nFeature", "pctMito", "pctSpikeIn") %in% colnames(Biobase::pData(input_eset)))) {
    cat("\tChecking passed! All 5 metrics are available.\n")
  } else {
    cat("\tPart of the 5 metrics are not available. Updating the eset automatically ...\n")
    input_eset <- updateSparseEset(input_eset, addMetaData = T)
    cat("\tChecking passed! All 5 metrics are available after updating.\n")
  }

  # generate the html QC report
  output_rmd_file <- gsub("html$", "Rmd", output_html_file)
  file.copy(from = system.file("rmd", "SparseEset_QC.Rmd", package = "scMINER"), to = output_rmd_file)
  rmarkdown::render(output_rmd_file, clean = TRUE, quiet = FALSE)
  file.remove(output_rmd_file)
}


