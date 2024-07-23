
#' Combine P values using Fisher's method or Stouffer's method
#'
#' @description
#' This function is used to calculate the combined P value and Z score from multiple P values. it can also be used to convert P value to Z score.
#'
#' @param pvals A vector of numeric, P values to be combined or converted to Z scores
#' @param method Character, method used to combine P values: "`Stouffer`" (the default) or "`Fisher`".
#' @param signed Logical, whether the input P values are signed or not. Usually they are signed by folder change. Default: `TRUE`.
#' @param twosided Logical, whether the input P values are two-sided. If FALSE, the input P values will be treated as one-tailed. Default: `TRUE`.
#'
#' @return A vector containing "Z-statistics" and "P.Value".
#' @export
#'
#' @examples
#' ## 1. combine P values from a vector
#' combinePvalVector(c(0.1,1e-3,1e-5), method = 'Stouffer')
#'
#' ## 2. combine P values from a table
#' df_Pcombined <- sapply(df$Pval*sign(df$log2FC), function(x) {combinePvalVector(x, twosided = TRUE)[2]})
#' df_Zscore <- sapply(df$Pval*sign(df$log2FC), function(x) {combinePvalVector(x, twosided = TRUE)[1]})
combinePvalVector <- function(pvals,
                              method = 'Stouffer',
                              signed = TRUE,
                              twosided = TRUE) {

  #remove NA pvalues
  pvals <- pvals[!is.na(pvals) & !is.null(pvals)]
  pvals[which(abs(pvals) <= 0)] <- .Machine$double.xmin
  if (sum(is.na(pvals)) >= 1) {
    stat <- NA
    pval <- NA
  } else {
    if (twosided & (sum(pvals > 1 | pvals < -1) >= 1)) {stop('pvalues must between 0 and 1!\n')}
    if (!twosided & (sum(pvals > 0.5 | pvals < -0.5) >= 1)) {stop('One-sided pvalues must between 0 and 0.5!\n')}
    if (!signed) {pvals <- abs(pvals)}

    signs <- sign(pvals)
    signs[signs == 0] <- 1

    if (grepl('Fisher', method, ignore.case = TRUE)) {
      if (twosided & signed) {
        neg.pvals <- pos.pvals <- abs(pvals) / 2
        pos.pvals[signs < 0] <- 1 - pos.pvals[signs < 0]
        neg.pvals[signs > 0] <- 1 - neg.pvals[signs > 0]
      } else {
        neg.pvals <- pos.pvals <- abs(pvals)
      }

      pvals <-
        c(1, -1) * c(
          pchisq(
            -2 * sum(log(as.numeric(pos.pvals))),
            df = 2 * base::length(pvals),
            lower.tail = FALSE
          ) / 2,
          pchisq(
            -2 * sum(log(as.numeric(neg.pvals))),
            df = 2 * base::length(pvals),
            lower.tail = FALSE
          ) / 2
        )
      pval <- base::min(abs(pvals))[1] #if two pvals are equal, pick up the first one
      stat <- sign(pvals[abs(pvals) == pval])[1] * qnorm(pval, lower.tail = F)[1]
      pval <- 2 * pval
    }
    else if (grepl('Stou', method, ignore.case = TRUE)) {
      if (twosided) {
        zs <- signs * qnorm(abs(pvals) / 2, lower.tail = FALSE)
        stat <- sum(zs) / sqrt(base::length(zs))
        pval <- 2 * pnorm(abs(stat), lower.tail = FALSE)
      }
      else{
        zs <- signs * qnorm(abs(pvals), lower.tail = FALSE)
        stat <- sum(zs) / sqrt(base::length(zs))
        pval <- pnorm(abs(stat), lower.tail = FALSE)
      }
    }
    else{
      stop('Only \"Fisher\" or \"Stouffer\" method is supported!!!\n')
    }
  }
  return(c(`Z-statistics` = stat, `P.Value` = pval))
}


#' Perform differential analysis between two groups
#'
#' @description
#' This function is used to perform the differential analysis between two groups. It supports there methods: "`limma`", "`wilcoxon`", and "`t.test`". This is an internal function for `getDE()` and `getDA()`.
#'
#' @param input_eset The expression set object that filtered, normalized and log-transformed
#' @param group_by Character, name of the column for grouping, usually the column of cell types or clusters. Default: "`clusterID`".
#' @param g1 A vector of character defining the fore-ground group or `NULL`. Default: `NULL`.
#' @param g0 A vector of character defining the back-ground group or `NULL`. Default: `NULL`.
#' @param use_method Character, method used for differential analysis: "`limma`" (the default), "`wilcoxon`", and "`t.test`".
#'
#' @return A data frame. Rows are genes/drivers, and columns are 11 statistics of differential analysis.
#' @export
#' @noRd
#'
#' @examples
#' ## to call this function
#' res <- compare2groups(input_eset = input_eset, group_by = group_by, g1 = g1_tmp, g0 = g0_tmp, use_method = use_method)
compare2groups <- function(input_eset,
                           group_by = "clusterID",
                           g1 = NULL, g0 = NULL,
                           use_method = "limma")
{
  ## check parameters
  if ((group_by %in% colnames(Biobase::pData(input_eset))) == FALSE) {stop('The group_by column,', group_by, ', was not found in the input set. Please check and re-try.')}

  pd <- Biobase::pData(input_eset)
  g1.cells <- row.names(pd[pd[,group_by] %in% g1,])
  if (length(g1.cells) > 0) {cat("\t", length(g1.cells), "cells were found for g1.\n")} else {stop('No cell was found for g1. Please check and re-try.')}
  g0.cells <- row.names(pd[pd[,group_by] %in% g0,])
  if (length(g0.cells) > 0) {cat("\t", length(g0.cells), "cells were found for g0.\n")} else {stop('No cell was found for g0. Please check and re-try.')}

  ## prepare matrix for limma
  use_samples <- c(g1.cells, g0.cells)
  exp_mat <- Biobase::exprs(input_eset[, use_samples])

  ## calculate the average and percentage
  exp_mat.g1 <- as.matrix(exp_mat[, g1.cells, drop = FALSE])
  exp_mat.g0 <- as.matrix(exp_mat[, g0.cells, drop = FALSE])
  g1_tag <- paste0(g1, collapse = ",")
  g0_tag <- paste0(g0, collapse = ",")
  g1_avg <- base::rowMeans(exp_mat.g1)
  g0_avg <- base::rowMeans(exp_mat.g0)
  g1_pct <- base::apply(exp_mat.g1, 1, function(x) {sum(x > 0) / length(x)})
  g0_pct <- base::apply(exp_mat.g0, 1, function(x) {sum(x > 0) / length(x)})
  df_stats <- data.frame(feature = row.names(exp_mat), g1_tag = as.character(g1_tag), g0_tag = g0_tag, g1_avg = g1_avg, g0_avg = g0_avg, g1_pct = g1_pct, g0_pct = g0_pct, log2FC = g1_avg - g0_avg)

  ## calculate P value using different methods
  if (use_method == "limma") {
    # prepare data
    group_id <- factor(c(rep("g1", length(g1.cells)), rep("g0", length(g0.cells))))
    design <- model.matrix(~ 0 + group_id)
    colnames(design) <- levels(group_id); rownames(design) <- colnames(exp_mat)

    # calculate P
    fit <- limma::lmFit(exp_mat, design) # fit linear model
    contrasts <- limma::makeContrasts(g1 - g0, levels = design) # create contrast matrix
    fit2 <- limma::contrasts.fit(fit, contrasts = contrasts) # fit contrasts
    fit2 <- limma::eBayes(fit2, trend = TRUE) # apply empirical Bayes moderation to improve the estimates
    topTable <- limma::topTable(fit2, adjust.method = "fdr", number = Inf, coef = 1)
    df_pval <- data.frame(feature = row.names(topTable), Pval = topTable$P.Value, FDR = topTable$adj.P.Val)
  } else if (use_method == "wilcoxon") {
    group_id <- factor(c(rep("g1", length(g1.cells)), rep("g0", length(g0.cells)))) # prepare group vector
    p_values <- numeric(nrow(exp_mat)) # Initialize a vector to store p-values
    for (i in 1:nrow(exp_mat)) {
      gene_expr <- exp_mat[i, ]
      test_result <- stats::wilcox.test(gene_expr ~ group_id)
      p_values[i] <- test_result$p.value
    }
    df_pval <- data.frame(feature = row.names(exp_mat), Pval = p_values, FDR = stats::p.adjust(p_values, method = "fdr"))
  } else if (use_method == "t.test") {
    group_id <- factor(c(rep("g1", length(g1.cells)), rep("g0", length(g0.cells)))) # prepare group vector
    p_values <- numeric(nrow(exp_mat)) # Initialize a vector to store p-values
    for (i in 1:nrow(exp_mat)) {
      gene_expr <- exp_mat[i, ]
      gene_expr.g1 <- gene_expr[group_id == "g1"]
      gene_expr.g0 <- gene_expr[group_id == "g0"]
      if (stats::var(gene_expr.g1) == 0 | stats::var(gene_expr.g1) == 0) {
        gene_expr[group_id == "g1"] <- gene_expr[group_id == "g1"] + stats::rnorm(length(gene_expr[group_id == "g1"]), mean = 0, sd = 1e-8)
        gene_expr[group_id == "g0"] <- gene_expr[group_id == "g0"] + stats::rnorm(length(gene_expr[group_id == "g0"]), mean = 0, sd = 1e-8)
      }
      test_result <- stats::t.test(gene_expr ~ group_id,)
      p_values[i] <- test_result$p.value
    }
    df_pval <- data.frame(feature = row.names(exp_mat), Pval = p_values, FDR = stats::p.adjust(p_values, method = "fdr"))
  } else {
    stop('The use_method was not recognized. Should be [limma | wilcoxon].')
  }

  ## merge two data frames
  df_res <- merge(df_stats, df_pval, by = "feature")

  df_res$Pval[df_res$Pval <= 0] <- .Machine$double.xmin # change 0 to 2.225074e-308 to avoid the failure in z score transformation
  df_res$Zscore <- sapply(df_res$Pval*sign(df_res$log2FC), function(x) {combinePvalVector(x, twosided = TRUE)[1]})
  df_res <- df_res[order(df_res$log2FC, decreasing = TRUE), ]

  return(df_res)
}


#' Perform differential expression analysis on expression set object
#'
#' @description
#' This function is used to perform the differential expression analysis on sparse eset object. It supports there methods: "`limma`", "`wilcoxon`", and "`t.test`".
#'
#'
#' @param input_eset The expression set object that filtered, normalized and log-transformed
#' @param group_by Character, name of the column for grouping, usually the column of cell types or clusters. Default: "`clusterID`".
#' @param g1 A vector of character defining the fore-ground group or `NULL`. Default: `NULL`.
#' @param g0 A vector of character defining the back-ground group or `NULL`. Default: `NULL`.
#' @param use_method Character, method used for differential analysis: "`limma`" (the default), "`wilcoxon`", and "`t.test`".
#'
#' @return A data frame. Rows are genes/drivers, and columns are 11 statistics of differential analysis.
#' @export
#'
#' @examples
#' ## 1. To perform differential expression analysis in a 1-vs-rest manner for all groups in "clusterID" column
#' de_res <- getDE(input_eset = clustered.eset, group_by = "clusterID", use_method = "limma")
#'
#' ## 2. To perform differential expression analysis in a 1-vs-rest manner for one specific group in "clusterID" column
#' de_res <- getDE(input_eset = clustered.eset, group_by = "clusterID", g1 = c("1"), use_method = "limma")
#'
#' ## 3. To perform differential expression analysis in a rest-vs-1 manner for one specific group in "clusterID" column
#' de_res <- getDE(input_eset = clustered.eset, group_by = "clusterID", g0 = c("1"), use_method = "limma")
#'
#' ## 4. To perform differential expression analysis in a 1-vs-1 manner for groups in "clusterID" column
#' de_res <- getDE(input_eset = clustered.eset, group_by = "clusterID", g1 = c("1"), g0 = c("3"), use_method = "limma")
getDE <- function(input_eset,
                  group_by = "clusterID",
                  g1 = NULL, g0 = NULL,
                  use_method = "limmma")
{
  ## check parameters
  if ((group_by %in% colnames(Biobase::pData(input_eset))) == FALSE) {stop('The group_by column [', group_by, '] was not found in the input set. Please check and re-try.')}
  grps <- sort(unique(Biobase::pData(input_eset)[, group_by]))
  cat(length(grps), "groups were found in group_by column [", group_by, "].\n")

  if (is.null(g1) == TRUE) {
    if (is.null(g0) == TRUE) {
      cat("Since no group was specified, the differential analysis will be conducted among all groups in the group_by column [", group_by, "] in the 1-vs-rest manner.\n")
      de_res <- data.frame()
      for (i in 1:length(grps)) {
        cat("\t", i, "/", length(grps), ": group 1 (", grps[i], ") vs the rest...\n")
        g1_tmp <- grps[i]; g0_tmp <- grps[! grps %in% g1_tmp]
        res <- compare2groups(input_eset = input_eset, group_by = group_by, g1 = g1_tmp, g0 = g0_tmp, use_method = use_method)

        if (i == 1) {de_res <- res} else {de_res <- rbind(de_res, res)}
      }
      return(de_res)
    } else {
      if (all(g0 %in% grps)) {
        cat("Since g0 was specified but g1 was not, all the other cells except those of g0 will be defined as g1.\n")
        g0_tmp <- g0; g1_tmp <- grps[! grps %in% g0_tmp]
        cat("\t", "1 / 1 : group 1 (", paste0(g1_tmp, collapse = ", "), ") vs group 0 (", paste0(g0_tmp, collapse = ", "), ") ...\n")
        res <- compare2groups(input_eset = input_eset, group_by = group_by, g1 = g1_tmp, g0 = g0_tmp, use_method = use_method)
        return(res)
      } else {
        g0.bad <- g0[! g0 %in% grps]
        stop(length(g0.bad), '/', length(g0), 'were not found in g0: ', paste0(g0.bad, collapse = ", "), ".")
      }
    }
  } else {
    if (is.null(g0) == TRUE) {
      if (all(g1 %in% grps)) {
        cat("Since g1 was specified but g0 was not, all the other cells except those of g1 will be defined as g0.\n")
        g1_tmp <- g1; g0_tmp <- grps[! grps %in% g1_tmp]
        cat("\t", "1 / 1 : group 1 (", paste0(g1_tmp, collapse = ", "), ") vs group 0 (", paste0(g0_tmp, collapse = ", "), ") ...\n")
        res <- compare2groups(input_eset = input_eset, group_by = group_by, g1 = g1_tmp, g0 = g0_tmp, use_method = use_method)
        return(res)
      } else {
        g1.bad <- g1[! g1 %in% grps]
        stop(length(g1.bad), '/', length(g1), 'were not found in g1: ', paste0(g1.bad, collapse = ", "), ".")
      }
    } else {
      if (all(c(g1, g0) %in% grps)) {
        g1_tmp <- g1; g0_tmp <- g0
        cat("\t", "1 / 1 : group 1 (", paste0(g1_tmp, collapse = ", "), ") vs group 0 (", paste0(g0_tmp, collapse = ", "), ") ...\n")
        res <- compare2groups(input_eset = input_eset, group_by = group_by, g1 = g1_tmp, g0 = g0_tmp, use_method = use_method)
        return(res)
      } else {
        g1.bad <- g1[! g1 %in% grps]; g0.bad <- g0[! g0 %in% grps]
        if (length(g1.bad) > 0) {
          if (length(g0.bad) > 0) {
            stop(length(g1.bad), '/', length(g1), 'were not found in g1: ', paste0(g1.bad, collapse = ", "), ".\n", length(g0.bad), '/', length(g0), 'were not found in g0: ', paste0(g0.bad, collapse = ", "), ".")
          } else {
            stop(length(g1.bad), '/', length(g1), 'were not found in g1: ', paste0(g1.bad, collapse = ", "), ".")
          }
        } else {
          stop(length(g0.bad), '/', length(g0), 'were not found in g0: ', paste0(g0.bad, collapse = ", "), ".")
        }
      }
    }
  }
}

#' Perform differential activity analysis on expression set
#'
#' @param input_eset The expression set object that filtered, normalized and log-transformed
#' @param group_by Character, name of the column for grouping, usually the column of cell types or clusters. Default: "`clusterID`".
#' @param g1 A vector of character defining the fore-ground group or `NULL`. Default: `NULL`.
#' @param g0 A vector of character defining the back-ground group or `NULL`. Default: `NULL`.
#' @param use_method Character, method used for differential analysis: "`limma`", "`wilcoxon`", and "`t.test`" (the default).
#'
#' @return A data frame. Rows are genes/drivers, and columns are 11 statistics of differential analysis.
#' @export
#'
#' @examples
#' ## 1. To perform differential activity analysis in a 1-vs-rest manner for all groups in "clusterID" column
#' da_res <- getDA(input_eset = activity_clustered.eset, group_by = "clusterID", use_method = "t.test")
#'
#' ## 2. To perform differential activity analysis in a 1-vs-rest manner for one specific group in "clusterID" column
#' da_res <- getDA(input_eset = activity_clustered.eset, group_by = "clusterID", g1 = c("1"), use_method = "t.test")
#'
#' ## 3. To perform differential activity analysis in a rest-vs-1 manner for one specific group in "clusterID" column
#' da_res <- getDA(input_eset = activity_clustered.eset, group_by = "clusterID", g0 = c("1"), use_method = "t.test")
#'
#' ## 4. To perform differential activity analysis in a 1-vs-1 manner for groups in "clusterID" column
#' da_res <- getDA(input_eset = activity_clustered.eset, group_by = "clusterID", g1 = c("1"), g0 = c("3"), use_method = "t.test")
getDA <- function(input_eset,
                  group_by = "clusterID",
                  g1 = NULL, g0 = NULL,
                  use_method = "t.test")
{
  ## check parameters
  if ((group_by %in% colnames(Biobase::pData(input_eset))) == FALSE) {stop('The group_by column [', group_by, '] was not found in the input set. Please check and re-try.')}
  grps <- sort(unique(Biobase::pData(input_eset)[, group_by]))
  cat(length(grps), "groups were found in group_by column [", group_by, "].\n")

  if (is.null(g1) == TRUE) {
    if (is.null(g0) == TRUE) {
      cat("Since no group was specified, the differential analysis will be conducted among all groups in the group_by column [", group_by, "] in the 1-vs-rest manner.\n")
      de_res <- data.frame()
      for (i in 1:length(grps)) {
        cat("\t", i, "/", length(grps), ": group 1 (", grps[i], ") vs the rest...\n")
        g1_tmp <- grps[i]; g0_tmp <- grps[! grps %in% g1_tmp]
        res <- compare2groups(input_eset = input_eset, group_by = group_by, g1 = g1_tmp, g0 = g0_tmp, use_method = use_method)

        if (i == 1) {de_res <- res} else {de_res <- rbind(de_res, res)}
      }
      return(de_res)
    } else {
      if (all(g0 %in% grps)) {
        cat("Since g0 was specified but g1 was not, all the other cells except those of g0 will be defined as g1.\n")
        g0_tmp <- g0; g1_tmp <- grps[! grps %in% g0_tmp]
        cat("\t", "1 / 1 : group 1 (", paste0(g1_tmp, collapse = ", "), ") vs group 0 (", paste0(g0_tmp, collapse = ", "), ") ...\n")
        res <- compare2groups(input_eset = input_eset, group_by = group_by, g1 = g1_tmp, g0 = g0_tmp, use_method = use_method)
        return(res)
      } else {
        g0.bad <- g0[! g0 %in% grps]
        stop(length(g0.bad), '/', length(g0), 'were not found in g0: ', paste0(g0.bad, collapse = ", "), ".")
      }
    }
  } else {
    if (is.null(g0) == TRUE) {
      if (all(g1 %in% grps)) {
        cat("Since g1 was specified but g0 was not, all the other cells except those of g1 will be defined as g0.\n")
        g1_tmp <- g1; g0_tmp <- grps[! grps %in% g1_tmp]
        cat("\t", "1 / 1 : group 1 (", paste0(g1_tmp, collapse = ", "), ") vs group 0 (", paste0(g0_tmp, collapse = ", "), ") ...\n")
        res <- compare2groups(input_eset = input_eset, group_by = group_by, g1 = g1_tmp, g0 = g0_tmp, use_method = use_method)
        return(res)
      } else {
        g1.bad <- g1[! g1 %in% grps]
        stop(length(g1.bad), '/', length(g1), 'were not found in g1: ', paste0(g1.bad, collapse = ", "), ".")
      }
    } else {
      if (all(c(g1, g0) %in% grps)) {
        g1_tmp <- g1; g0_tmp <- g0
        cat("\t", "1 / 1 : group 1 (", paste0(g1_tmp, collapse = ", "), ") vs group 0 (", paste0(g0_tmp, collapse = ", "), ") ...\n")
        res <- compare2groups(input_eset = input_eset, group_by = group_by, g1 = g1_tmp, g0 = g0_tmp, use_method = use_method)
        return(res)
      } else {
        g1.bad <- g1[! g1 %in% grps]; g0.bad <- g0[! g0 %in% grps]
        if (length(g1.bad) > 0) {
          if (length(g0.bad) > 0) {
            stop(length(g1.bad), '/', length(g1), 'were not found in g1: ', paste0(g1.bad, collapse = ", "), ".\n", length(g0.bad), '/', length(g0), 'were not found in g0: ', paste0(g0.bad, collapse = ", "), ".")
          } else {
            stop(length(g1.bad), '/', length(g1), 'were not found in g1: ', paste0(g1.bad, collapse = ", "), ".")
          }
        } else {
          stop(length(g0.bad), '/', length(g0), 'were not found in g0: ', paste0(g0.bad, collapse = ", "), ".")
        }
      }
    }
  }
}


#' Pick the top genes/drivers for each group from differential analysis results
#'
#' @description
#' This function is used to pick the top genes/drivers for each group from differential analysis results based on either `fold change`, `p value` or `FDR`.
#'
#'
#' @param input_table The table generated by `getDE()` or `getDA()`, containing 11 statistics of differential analysis.
#' @param number An Integer, number of top genes/driver to pick from each group. Default: 10.
#' @param group_by Character, name of the column for grouping, usually the column of cell types or clusters. Default: "`g1_tag`".
#' @param sort_by Character, name of the column for sorting. Default: "`log2FC`".
#' @param sort_decreasing Logical, whether to sort the column specified by `sort_by` in decreasing order. If `FALSE`, the column will be sorted in increasing order. Default: `TRUE`.
#'
#' @return A data frame with top genes/driver of each group
#' @export
#'
#' @examples
#' top_drivers <- getTopFeature(da_res, number = 10, group_by, "g1_tag")
getTopFeatures <- function(input_table,
                           number = 10,
                           group_by = "g1_tag",
                           sort_by = "log2FC",
                           sort_decreasing = TRUE)
{
  ## check parameters
  if ((group_by %in% colnames(input_table)) == FALSE) {stop('The group_by column [', group_by, '] was not found in the input table. Please check and re-try.')}
  if ((sort_by %in% colnames(input_table)) == FALSE) {stop('The sort_by column [', sort_by, '] was not found in the input table. Please check and re-try.')}

  input_table.sorted <- input_table[order(input_table[, sort_by], decreasing = sort_decreasing), ]
  grps <- unique(input_table[, group_by])
  table_res <- data.frame()
  for (i in 1:length(grps)) {
    table.sel <- input_table.sorted[input_table.sorted[, group_by] %in% c(grps[i]), , drop = FALSE]
    if (nrow(table.sel) >= number) {
      if (i == 1) {table_res <- table.sel[c(1:number),]} else {table_res <- rbind(table_res, table.sel[number,])}
    } else {
      if (i == 1) {table_res <- table.sel} else {table_res <- rbind(table_res, table.sel)}
    }
  }
  return(table_res)
}

