
#' Violin plot showing the expression or activity of selected features by self-defined groups
#'
#' @description
#' This function is used to draw a violin plot of selected features among self-defined groups from a sparse eset object.
#'
#' @param input_eset The expression set object that filtered, normalized and log-transformed
#' @param features A vector of genes or drivers (row.names of the input eset) to plot
#' @param group_by Character, name of the column for grouping, usually the column of cell types or clusters. Default: "`clusterID`".
#' @param ncol Integer, number of columns when multiple plots are displayed. Default: 3.
#' @param colors A vector of colors for filling the violins. The length should be same as the number of groups. Default: `NULL` (ggplot default colors).
#' @param legend.position Character, position of legend: "`right`" (the default), "`left`", "`top`", "`bottom`" or "`none`".
#' @param fontsize.legend_title Integer, font size of the legend title. Default: 10.
#' @param fontsize.legend_text Integer, font size of the legend text. Default: 8.
#' @param fontsize.strip Integer, font size of the plot strip. Default: 10.
#' @param fontsize.axis_title Integer, font size of the axis title. Default: 10.
#' @param fontsize.axis_text Integer, font size of the axis text. Default: 8.
#' @param xlabel.angle Numeric, the angle of the a-axis title. When it's set not 0, the x-axis text will automatically right-justified. Default: 0.
#' @param ylabel.text Character, the title of y-axis. Default: "`Expression (log2CPM)`"
#' @param stat_method Character or `NULL`. method of the stat point to show: "`median`" (the default), "`mean`". If `NULL`, the stat point won't show up.
#' @param add_boxplot Logical, whether to add box plot. Default: `FALSE`.
#' @param boxplot.width Numeric, width of the box plot relative to the body of violin plot, ranging from 0 to 1. Default: 0.3. Ignored if `add_boxplot` = `FALSE`.
#' @param boxplot.fill Character, color used to fill the box plots. Default: "`white`". Ignored if `add_boxplot` = `FALSE`.
#' @param boxplot.alpha Numerical, transparency of box plots, ranging from 0 (more transparent) to 1 (less transparent). Default: 0.8. Ignored if `add_boxplot` = `FALSE`.
#' @param add_jitter Logical, whether to add jittered points. Default: `FALSE`.
#' @param jitter.height Numeric, amount of vertical jitter. Default: 0.
#' @param jitter.width Numeric, amount of horizontal jitter. Default: 0.3.
#' @param jitter.size Numeric, size of the jittered points. Default: 0.1.
#'
#' @return A ggplot object with one or multiple violin plots
#' @import ggplot2
#' @export
#'
#' @examples
#' ## 1. violin plots grouped by clusters (say the column name is 'clusterID')
#' p_vln <- feature_vlnplot(input_eset = clustered.eset, features = c("CD14", "CD19", "CD8A"), group_by = "clusterID")
#'
#' ## 2. violin plots grouped by cell types (say the column name is 'cellType')
#' p_vln <- feature_vlnplot(input_eset = clustered.eset, features = c("CD14", "CD19", "CD8A"), group_by = "cellType")
#'
#' ## 3. customize the colors to fill the violin plots
#' p_vln <- feature_vlnplot(input_eset = clustered.eset, features = c("CD14", "CD19", "CD8A"), group_by = "clusterID", colors = c("blue", "red", "green"))
#'
#' ## 4. add jittered points
#' p_vln <- feature_vlnplot(input_eset = clustered.eset, features = c("CD14", "CD19", "CD8A"), group_by = "clusterID", add_jitter = TRUE, jitter.width = 0.5, jitter.size = 0.5)
#'
#' ## 5. using activity data
#' p_vln <- feature_vlnplot(input_eset = activity_clustered.eset, features = c("CD14_SIG", "CD19_SIG", "CD8A_SIG"), group_by = "clusterID", ylabel_text = "Activity")
feature_vlnplot <- function(input_eset,
                            features = NULL,
                            group_by = "clusterID",
                            ncol = 3, colors = NULL,
                            legend.position = "right", fontsize.legend_title = 10, fontsize.legend_text = 8,
                            fontsize.strip = 10, fontsize.axis_title = 10, fontsize.axis_text = 8, xlabel.angle = 0, ylabel.text = "Expression (log2CPM)",
                            stat_method = "median",
                            add_boxplot = FALSE, boxplot.width = 0.3, boxplot.fill = "white", boxplot.alpha = 0.8,
                            add_jitter = FALSE, jitter.height = 0, jitter.width = 0.3, jitter.size = 0.1)
{
  ## check parameters
  if (is.null(features) == TRUE) {
    stop("The features can not be none. Please specify the features/genes and re-try.")
  } else {
    features_in <- features[features %in% row.names(input_eset)]
    if (length(features_in) == 0) {stop('None of the specified features was found in the input eset. Please check and re-try!')}

    features_out <- features[! features %in% row.names(input_eset)]
    if (length(features_out) > 0) {warning(length(features_out), '/', length(features), ' specified features were not found in the input eset: ', paste0(features_out, collapse = ", "), '.')}
  }

  if ((group_by %in% colnames(Biobase::pData(input_eset))) == FALSE) {stop('The group_by column,', group_by, ', was not found in the input set. Please check and re-try.')}

  ## prepare the master table for visualization
  exp_mat <- as.matrix(Biobase::exprs(input_eset))
  master <- data.frame(groups = as.factor(Biobase::pData(input_eset)[, group_by]), t(exp_mat[features_in, , drop = FALSE]))
  master_melt <- reshape2::melt(master, id.vars = c("groups"))

  ## visualize
  p <- ggplot(master_melt, aes(x = groups, y = value, fill = groups)) + geom_violin(trim = TRUE, scale = "width", na.rm = TRUE) + theme_classic() + guides(fill=guide_legend(title = group_by))

  # add stat dot, supporting median or mean
  if (is.null(stat_method) == FALSE) {
    if (stat_method == "median") {p <- p + stat_summary(fun = median, geom = "point", size = 1.2, color = "black", position = position_dodge(width = 1))}
    else if (stat_method == "mean") {p <- p + stat_summary(fun = mean, geom = "point", size = 1.2, color = "black", position = position_dodge(width = 1))}
    else {stop('The stat_method can not be recognized. Should be [median | mean]. Set it to NULL to skip.')}
  }

  # add box plot
  if (add_boxplot == TRUE) {p <- p + geom_boxplot(fill = boxplot.fill, alpha = boxplot.alpha, width = boxplot.width, outlier.shape = NA, show.legend = FALSE, na.rm = TRUE)}

  # add jitter point
  if (add_jitter == TRUE) {p <- p + geom_jitter(height = jitter.height, width = jitter.width, size = jitter.size)}

  hjust_x = base::ifelse (xlabel.angle == 0, 0.5, 1)
  p <- p + facet_wrap(~variable, scales = "free", ncol = ncol) + labs(x = group_by, y = ylabel.text) +
    theme(
      legend.position = legend.position,
      legend.title = element_text(size = fontsize.legend_title, face = "bold"),
      legend.text = element_text(size = fontsize.legend_text),
      strip.text = element_text(size = fontsize.strip, face = "bold"),
      axis.title = element_text(size = fontsize.axis_title, face = "bold", hjust = 0.5, color = "black"),
      axis.text.x = element_text(size = fontsize.axis_text, hjust = hjust_x, color = "black", angle = xlabel.angle),
      axis.text.y = element_text(size = fontsize.axis_text, hjust = 1, color = "black"))

  if (is.null(colors) == FALSE) {p <- p + scale_fill_manual(values = colors)}

  return(p)
}


#' Box plot showing the expression or activity of selected features by self-defined groups
#'
#' @description
#' This function is used to draw a box plot of selected features among self-defined groups from a sparse eset object.
#'
#' @param input_eset The expression set object that filtered, normalized and log-transformed
#' @param features A vector of genes or drivers (row.names of the input eset) to plot
#' @param group_by Character, name of the column for grouping, usually the column of cell types or clusters. Default: "`clusterID`".
#' @param ncol Integer, number of columns when multiple plots are displayed. Default: 3.
#' @param colors A vector of colors for filling the violins. The length should be same as the number of groups. Default: `NULL` (ggplot default colors).
#' @param legend.position Character, position of legend: "`right`" (the default), "`left`", "`top`", "`bottom`" or "`none`".
#' @param fontsize.legend_title Integer, font size of the legend title. Default: 10.
#' @param fontsize.legend_text Integer, font size of the legend text. Default: 8.
#' @param fontsize.strip Integer, font size of the plot strip. Default: 10.
#' @param fontsize.axis_title Integer, font size of the axis label and text. Default: 10.
#' @param fontsize.axis_text Integer, font size of the axis label and text. Default: 8.
#' @param xlabel.angle Numeric, the angle of the a-axis title. When it's set not 0, the x-axis text will automatically right-justified. Default: 0.
#' @param ylabel.text Character, the title of y-axis. Default: "`Expression (log2CPM)`"
#' @param stat_method Character or NULL. method of the stat point to show: "`median`" (the default), "`mean`". If `NULL`, the stat point won't show up.
#' @param add_jitter Logical, whether to add jittered points. Default: `FALSE`.
#' @param jitter.height Numeric, amount of vertical jitter. Default: 0.
#' @param jitter.width Numeric, amount of horizontal jitter. Default: 0.3.
#' @param jitter.size Numeric, size of the jittered points. Default: 0.1.
#'
#' @return A ggplot object with one or multiple box plots
#' @export
#'
#' @examples
#' ## 1. violin plots grouped by clusters (say the column name is 'clusterID')
#' p_box <- feature_boxplot(input_eset = clustered.eset, features = c("CD14", "CD19", "CD8A"), group_by = "clusterID")
#'
#' ## 2. violin plots grouped by cell types (say the column name is 'cellType')
#' p_box <- feature_boxplot(input_eset = clustered.eset, features = c("CD14", "CD19", "CD8A"), group_by = "cellType")
#'
#' ## 3. customize the colors to fill the violin plots
#' p_box <- feature_boxplot(input_eset = clustered.eset, features = c("CD14", "CD19", "CD8A"), group_by = "clusterID", colors = c("blue", "red", "green"))
#'
#' ## 4. add jittered points
#' p_box <- feature_boxplot(input_eset = clustered.eset, features = c("CD14", "CD19", "CD8A"), group_by = "clusterID", add_jitter = TRUE, jitter.width = 0.5, jitter.size = 0.5)
#'
#' ## 5. using activity data
#' p_box <- feature_boxplot(input_eset = activity_clustered.eset, features = c("CD14_SIG", "CD19_SIG", "CD8A_SIG"), group_by = "clusterID", ylabel_text = "Activity")
feature_boxplot <- function(input_eset,
                            features = NULL,
                            group_by = "clusterID",
                            ncol = 3, colors = NULL,
                            legend.position = "right", fontsize.legend_title = 10, fontsize.legend_text = 8,
                            fontsize.strip = 10, fontsize.axis_title = 10, fontsize.axis_text = 8, xlabel.angle = 0, ylabel.text = "Expression (log2CPM)",
                            stat_method = "median",
                            add_jitter = FALSE, jitter.height = 0, jitter.width = 0.3, jitter.size = 0.1)
{
  ## check parameters
  if (is.null(features) == TRUE) {
    stop("The features can not be none. Please specify the features/genes and re-try.")
  } else {
    features_in <- features[features %in% row.names(input_eset)]
    if (length(features_in) == 0) {stop('None of the specified features was found in the input eset. Please check and re-try!')}

    features_out <- features[! features %in% row.names(input_eset)]
    if (length(features_out) > 0) {warning(length(features_out), '/', length(features), ' specified features were not found in the input eset: ', paste0(features_out, collapse = ", "), '.')}
  }

  if ((group_by %in% colnames(Biobase::pData(input_eset))) == FALSE) {stop('The group_by column,', group_by, ', was not found in the input set. Please check and re-try.')}

  ## prepare the master table for visualization
  exp_mat <- as.matrix(Biobase::exprs(input_eset))
  master <- data.frame(groups = as.factor(Biobase::pData(input_eset)[, group_by]), t(exp_mat[features_in, , drop = FALSE]))
  master_melt <- reshape2::melt(master, id.vars = c("groups"))

  ## visualize
  p <- ggplot(master_melt, aes(x = groups, y = value, fill = groups)) +
    geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
    theme_classic() + guides(fill=guide_legend(title = group_by))

  # add stat dot, supporting median or mean
  if (is.null(stat_method) == FALSE) {
    if (stat_method == "median") {p <- p + stat_summary(fun = median, geom = "point", size = 1.2, color = "black", position = position_dodge(width = 1))}
    else if (stat_method == "mean") {p <- p + stat_summary(fun = mean, geom = "point", size = 1.2, color = "black", position = position_dodge(width = 1))}
    else {stop('The stat_method can not be recognized. Should be [median | mean]. Set it to NULL to skip.')}
  }

  # add jitter point
  if (add_jitter == TRUE) {p <- p + geom_jitter(height = jitter.height, width = jitter.width, size = jitter.size)}

  hjust_x = base::ifelse (xlabel.angle == 0, 0.5, 1)
  p <- p + facet_wrap(~variable, scales = "free", ncol = ncol) + labs(x = group_by, y = ylabel.text) +
    theme(
      legend.position = legend.position,
      legend.title = element_text(size = fontsize.legend_title, face = "bold"),
      legend.text = element_text(size = fontsize.legend_text),
      strip.text = element_text(size = fontsize.strip, face = "bold"),
      axis.title = element_text(size = fontsize.axis_title, face = "bold", hjust = 0.5, color = "black"),
      axis.text.x = element_text(size = fontsize.axis_text, hjust = hjust_x, color = "black", angle = xlabel.angle),
      axis.text.y = element_text(size = fontsize.axis_text, hjust = 1, color = "black"))

  if (is.null(colors) == FALSE) {p <- p+ scale_fill_manual(values = colors)}

  return(p)
}


#' Scatter plot showing the expression or activity of selected features on UMAP or t-SNE coordinates
#'
#' @description
#' This function is used to draw a scatter plot of selected features on UMAP or t-SNE coordinates from a sparse eset object.
#'
#' @param input_eset The expression set object that filtered, normalized and log-transformed
#' @param features A vector of genes or drivers (row.names of the input eset) to plot
#' @param location_x Character, name of the column of x-axis coordinates. Default: "`UMAP_1`".
#' @param location_y Character, name of the column of y-axis coordinates. Default: "`UMAP_2`".
#' @param colors A vector of two colors indicating the low and high values respectively. Default: c("`lightgrey`", "`red`").
#' @param ncol Integer, number of columns when multiple plots are displayed. Default: 3.
#' @param point.size Numeric, size of the scatter points. Default: 0.5.
#' @param legend.position Character, position of legend: "`right`" (the default), "`left`", "`top`", "`bottom`" or "`none`".
#' @param legend.key_height Numeric, height of the legend key in unit of "cm". Default: 3.
#' @param legend.key_width Numeric, width of the legend key in unit of "cm". Default: 1.
#' @param fontsize.legend_title Numeric, font size of the legend title. Default: 10.
#' @param fontsize.legend_text Numeric, font size of the legend text. Default: 8.
#' @param fontsize.strip Numeric, font size of the plot strip. Default: 10.
#' @param fontsize.axis_title Numeric, font size of the axis label and text. Default: 10.
#' @param fontsize.axis_text Numeric, font size of the axis label and text. Default: 8.
#'
#' @return Print a plot to screen and return a gtable containing a list of plots, can be visualized by plot(g), and saved by ggsave(file = "output.pdf", g)
#' @export
#'
#' @examples
#' ## 1. scatter plots with UMAP projections
#' p_scatter <- feature_scatterplot(input_eset = clustered.eset, features = c("CD14", "CD19", "CD8A"), location_x = "UMAP_1", location_y = "UMAP_2")
#'
#' ## 2. scatter plots with t-SNE projections
#' p_scatter <- feature_scatterplot(input_eset = clustered.eset, features = c("CD14", "CD19", "CD8A"), location_x = "tSNE_1", location_y = "tSNE_2")
#'
#' ## 3. change the point size and font size
#' p_scatter <- feature_scatterplot(input_eset = clustered.eset, features = c("CD14", "CD19", "CD8A"), location_x = "UMAP_1", location_y = "UMAP_2", point.size = 1, fontsize.strip = 12, fontsize.axis = 10)
feature_scatterplot <- function(input_eset,
                                features = NULL,
                                location_x = "UMAP_1", location_y = "UMAP_2",
                                colors = NULL,
                                ncol = 3, point.size = 0.5,
                                legend.position = "right",
                                legend.key_height = 3, legend.key_width = 1,
                                fontsize.legend_title = 10, fontsize.legend_text = 8,
                                fontsize.strip = 10, fontsize.axis_title = 10, fontsize.axis_text = 8)
{
  ## check parameters
  if (is.null(features) == TRUE) {
    stop("The features can not be none. Please specify the features/genes and re-try.")
  } else {
    features_in <- features[features %in% row.names(input_eset)]
    if (length(features_in) == 0) {stop('None of the specified features was found in the input eset. Please check and re-try!')}

    features_out <- features[! features %in% row.names(input_eset)]
    if (length(features_out) > 0) {warning(length(features_out), '/', length(features), ' specified features were not found in the input eset: ', paste0(features_out, collapse = ", "), '.')}
  }

  if ((location_x %in% colnames(Biobase::pData(input_eset))) == FALSE) {stop('The location_x column,', location_x, ', was not found in the input set. Please check and re-try.')}
  if ((location_y %in% colnames(Biobase::pData(input_eset))) == FALSE) {stop('The location_y column,', location_y, ', was not found in the input set. Please check and re-try.')}

  ## prepare the master table for visualization
  exp_mat <- as.matrix(Biobase::exprs(input_eset))
  master <- data.frame(Biobase::pData(input_eset)[, c(location_x, location_y)], t(exp_mat[features_in, , drop = FALSE]))
  master_melt <- reshape2::melt(master, id.vars = c(location_x, location_y))

  ps <- list()
  if (is.null(colors) == TRUE) {
    color.low <- "lightgrey"; color.high <- "red"
  } else {
      if (length(colors) == 2) {color.low <- colors[1]; color.high <- colors[2]} else {stop('The length of "colors" must be of 2.')}
  }

  for (i in 1:length(features_in)) {
    master_melt.sel <- master_melt[master_melt$variable == features_in[i],]

    p <- ggplot(master_melt.sel, aes(x = master_melt.sel[, location_x], y = master_melt.sel[, location_y])) +
      geom_point(aes(colour = value), size = point.size, pch = 16, stroke = NA) +
      scale_colour_gradient(low = color.low, high = color.high) + theme_classic() +
      facet_wrap(~variable, scales = "fixed") + labs(x = location_x, y = location_y) +
      theme(
        legend.position = legend.position,
        legend.key.height = unit(legend.key_height, 'cm'), legend.key.width = unit(legend.key_width, 'cm'),
        legend.title = element_text(size = fontsize.legend_title, face = "bold"),
        legend.text = element_text(size = fontsize.legend_text),
        strip.text = element_text(size = fontsize.strip, face = "bold"),
        axis.title = element_text(size = fontsize.axis_title, face = "bold", hjust = 0.5, color = "black"),
        axis.text.x = element_text(size = fontsize.axis_text, hjust = 0.5, color = "black"),
        axis.text.y = element_text(size = fontsize.axis_text, hjust = 1, color = "black"))

    ps[[i]] <- p
  }

  g <- gridExtra::grid.arrange(grobs = ps, ncol = ncol)
}


#' Bubble blot showing the expression or activity of selected features by self-defined groups
#'
#' @description
#' This function is used to draw a bubble plot of selected features among self-defined groups from a sparse eset object.
#'
#' @param input_eset The expression set object that filtered, normalized and log-transformed
#' @param features A vector of genes or drivers (row.names of the input eset) to plot
#' @param group_by Character, name of the column for grouping, usually the column of cell types or clusters. Default: "`clusterID`".
#' @param colors A vector of two colors indicating the low and high values respectively. Default: c("`lightgrey`", "`red`").
#' @param legend.position Character, position of legend: "`right`" (the default), "`left`", "`top`", "`bottom`" or "`none`".
#' @param fontsize.legend_title Integer, font size of the legend title. Default: 10.
#' @param fontsize.legend_text Integer, font size of the legend text. Default: 8.
#' @param fontsize.axis_title Integer, font size of the axis label and text. Default: 10.
#' @param fontsize.axis_text Integer, font size of the axis label and text. Default: 8.
#' @param xlabel.angle Numeric, the angle of the a-axis title. When it's set not 0, the x-axis text will automatically right-justified. Default: 0.
#'
#' @return A ggplot object of bubble plot
#' @export
#'
#' @examples
#' features_of_interest <- c("CD3D","CD27","IL7R","SELL","CCR7","IL32","GZMA","GZMK","DUSP2","CD8A","GZMH","GZMB","CD79A","CD79B","CD86","CD14")
#' ## 1. the most commonly used command
#' p_bubble <- feature_bubbleplot(input_eset = clustered.eset, features = features_of_interest, group_by = "clusterID")
#'
#' ## 2. customize the colors
#' p_bubble <- feature_bubbleplot(input_eset = clustered.eset, features = features_of_interest, group_by = "clusterID", colors = c("lightgrey", "red"))
feature_bubbleplot <- function(input_eset,
                               features = NULL,
                               group_by = "clusterID",
                               colors = NULL,
                               legend.position = "right", fontsize.legend_title = 10, fontsize.legend_text = 8,
                               fontsize.axis_title = 10, fontsize.axis_text = 8, xlabel.angle = 0)
{
  ## check parameters
  if (is.null(features) == TRUE) {
    stop("The features can not be none. Please specify the features/genes and re-try.")
  } else {
    features_in <- features[features %in% row.names(input_eset)]
    if (length(features_in) == 0) {stop('None of the specified features was found in the input eset. Please check and re-try!')}

    features_out <- features[! features %in% row.names(input_eset)]
    if (length(features_out) > 0) {warning(length(features_out), '/', length(features), ' specified features were not found in the input eset: ', paste0(features_out, collapse = ", "), '.')}
  }

  if ((group_by %in% colnames(Biobase::pData(input_eset))) == FALSE) {stop('The group_by column,', group_by, ', was not found in the input set. Please check and re-try.')}

  ## prepare the master table for visualization
  exp <- t(as.matrix(Biobase::exprs(input_eset)[features_in, , drop = FALSE]))
  pd <- Biobase::pData(input_eset)[, group_by, drop = FALSE]

  # calculate the mean and percentage
  groups <- sort(unique(pd[, group_by]))
  markers <- colnames(exp)
  df_mean <- matrix(NA, nrow = length(groups), ncol = length(markers), dimnames = list(groups, markers))
  df_pct <- matrix(NA, nrow = length(groups), ncol = length(markers), dimnames = list(groups, markers))
  for (i in 1:length(groups)) {
    pd.sel <- pd[pd[, group_by] %in% c(groups[i]), , drop = FALSE]
    exp.sel <- exp[row.names(pd.sel), , drop = FALSE]
    df_mean[groups[i],] <- colMeans(exp.sel, na.rm = TRUE)
    df_pct[groups[i],] <- apply(exp.sel, 2, function(x) {sum(x > 0) / length(x) * 100})
  }

  df_mean.long <- reshape2::melt(df_mean); colnames(df_mean.long) <- c(group_by, "Features", "Mean_value")
  df_pct.long <- reshape2::melt(df_pct); colnames(df_pct.long) <- c(group_by, "Features", "Percentage_value")
  master <- dplyr::left_join(df_mean.long, df_pct.long, by = c(group_by, "Features"))

  # visualize
  hjust_x = base::ifelse (xlabel.angle == 0, 0.5, 1)
  if (is.null(colors) == TRUE) {color.low <- "lightgrey"; color.high <- "red"}
  else {if (length(colors) == 2) {color.low <- colors[1]; color.high <- colors[2]} else {stop('The length of "colors" must be of 2.')}}
  p <- ggplot(master, aes(x = Features, y = as.factor(master[, group_by]))) + geom_point(aes(color = Mean_value, size = Percentage_value), pch = 16, stroke = NA) +
    scale_colour_gradient(low = color.low, high = color.high) + theme_classic() + labs(x = "", y = group_by) +
    theme(
      legend.position = legend.position,
      legend.title = element_text(size = fontsize.legend_title, face = "bold"),
      legend.text = element_text(size = fontsize.legend_text),
      axis.title = element_text(size = fontsize.axis_title, face = "bold", hjust = 0.5, color = "black"),
      axis.text.x = element_text(size = fontsize.axis_text, hjust = hjust_x, color = "black", angle = xlabel.angle),
      axis.text.y = element_text(size = fontsize.axis_text, hjust = 1, color = "black"))

  return(p)
}

#' Heatmap showing the expression or activity of selected features by self-defined groups
#'
#' @description
#' This function is used to draw a heatmap of selected features among self-defined groups from a sparse eset object. By default, the groups are sorted by size, from largest to smallest. Within each group, the cells are sorted alphabetically.
#'
#' @param input_eset The expression set object that filtered, normalized and log-transformed
#' @param features A vector of genes or drivers (row.names of the input eset) to plot
#' @param group_by Character, name of the column for grouping, usually the column of cell types or clusters. Default: "`clusterID`".
#' @param scale_method Character, method for data scaling: "`none`" (the default), "`column`", "`row`".
#' @param annotation_columns Character, name(s) of the column(s) to add for cell annotation. Default: `NULL`.
#' @param use_gaps.column Logical, whether to put a gap between cell groups. Default: `FALSE`.
#' @param use_gaps.row Logical, whether to put a gap between features. Default: `FALSE`.
#' @param show_rownames Logical, whether to show the rownames. Default: `TRUE`.
#' @param fontsize.row Numeric, font size of the rownames. Defualt: 10.
#' @param cluster_rows Logical, whether to cluster the rows. If `TRUE`, the rows will be clustered. If `FALSE`, the rows are displays following the order in '`features`'. Default: `FALSE`.
#'
#' @return Print the heatmap to screen
#' @import pheatmap
#' @export
#'
#' @examples
#' features_of_interest <- c("CD3D","CD27","IL7R","SELL","CCR7","IL32","GZMA","GZMK","DUSP2","CD8A","GZMH","GZMB","CD79A","CD79B","CD86","CD14")
#' ## 1. the most commonly used command
#' feature_heatmap(input_eset = clustered.eset, features = features_of_interest, group_by = "clusterID")
#'
#' ## 2. add one more column ('true_label') for cell annotation
#' feature_heatmap(input_eset = clustered.eset, features = features_of_interest, group_by = "clusterID", annotation_columns = c("true_label"))
#'
#' ## 3. scale the data by row
#' feature_heatmap(input_eset = clustered.eset, features = features_of_interest, group_by = "clusterID", scale_method = "row")
#'
#' ## 4. cluster the rows
#' feature_heatmap(input_eset = clustered.eset, features = features_of_interest, group_by = "clusterID", cluster_rows = TRUE)
#'
#' ## 5. add gaps
#' feature_heatmap(input_eset = clustered.eset, features = features_of_interest, group_by = "clusterID", use_gaps.column = TRUE, use_gaps.row = TRUE)
feature_heatmap <- function(input_eset,
                            features = NULL,
                            group_by = "clusterID",
                            scale_method = "none", annotation_columns = NULL, use_gaps.column = FALSE,
                            cluster_rows = FALSE, show_rownames = TRUE, fontsize.row = 10, use_gaps.row = FALSE)
{
  ## check parameters
  if (is.null(features) == TRUE) {
    stop("The features can not be none. Please specify the features/genes and re-try.")
  } else {
    features_in <- features[features %in% row.names(input_eset)]
    if (length(features_in) == 0) {stop('None of the specified features was found in the input eset. Please check and re-try!')}

    features_out <- features[! features %in% row.names(input_eset)]
    if (length(features_out) > 0) {warning(length(features_out), '/', length(features), ' specified features were not found in the input eset: ', paste0(features_out, collapse = ", "), '.')}
  }

  if (is.null(annotation_columns) == TRUE) {
    annotation_columns <- group_by
  } else {
    annotation_columns <- unique(c(group_by, annotation_columns))
  }
  annotation_columns.not_found <- annotation_columns[!annotation_columns %in% colnames(Biobase::pData(input_eset))]
  if (length(annotation_columns.not_found) > 0) {stop('These columns used for groupping or annotation were not found: ', paste0(annotation_columns.not_found, collapse = ", "), '. Please check and re-try.')}

  ## prepare the master table for visualization
  exp_mat <- as.matrix(Biobase::exprs(input_eset))
  master <- exp_mat[features_in, , drop = FALSE]

  # sort the groups from largest to smallest
  groupID_sorted <- names(sort(table(Biobase::pData(input_eset)[, group_by]), decreasing = TRUE))
  groupSize_sorted <- as.numeric(sort(table(Biobase::pData(input_eset)[, group_by]), decreasing = TRUE))

  # sort cells from the largest group to the smallest one, in each group, sort cells alphabetically by cell ID
  cellID_sorted <- c()
  for (i in 1:length(groupID_sorted)) {
    cells.tmp <- colnames(input_eset)[Biobase::pData(input_eset)[, group_by] == groupID_sorted[i]]
    cells.tmp <- sort(cells.tmp, decreasing = TRUE)
    if (i == 1) {cellID_sorted <- cells.tmp} else {cellID_sorted <- c(cellID_sorted, cells.tmp)}
  }
  master_sorted <- master[, cellID_sorted]

  ## prepare the cell annotation matrix
  cellAnnotation <- Biobase::pData(input_eset)[, unique(c(group_by, annotation_columns)), drop = FALSE]

  ## make the heatmap
  if (use_gaps.column == TRUE) {use_gaps.column <- base::cumsum(groupSize_sorted)}
  if (use_gaps.row == TRUE) {use_gaps.row <- 1:length(features_in)}
  p <- pheatmap::pheatmap(mat = as.matrix(master_sorted),
                          color = rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(256)),
                          border_color = NA, scale = scale_method,
                          cluster_rows = cluster_rows, treeheight_row = 0, show_rownames = show_rownames, fontsize_row = fontsize.row,
                          gaps_col = use_gaps.column, gaps_row = use_gaps.row,
                          cluster_cols = FALSE, show_colnames = FALSE, annotation_col = cellAnnotation)

  return(p)
}


#' Bar plot showing the cell composition of self-defined groups
#'
#' @description
#' This function is used to draw a bar plot showing the cell composition of self-defined groups.
#'
#' @param input_eset The expression set object that filtered, normalized and log-transformed
#' @param group_by Character, name of the column for grouping, usually the column of cell types or clusters. Default: "`clusterID`".
#' @param color_by Character, name of the column for color-coding, usually the column of cell types or clusters. Default: "`cell_type`".
#' @param colors A vector of colors for filling the violins. The length should be same as the number of groups. Default: `NULL` (ggplot default colors).
#' @param legend.position Character, position of legend: "`right`" (the default), "`left`", "`top`", "`bottom`" or "`none`".
#' @param xlabel.angle Numeric, the angle of the a-axis title. When it's set not 0, the x-axis text will automatically right-justified. Default: 0.
#' @param fontsize.legend_title Integer, font size of the legend title. Default: 10.
#' @param fontsize.legend_text Integer, font size of the legend text. Default: 8.
#' @param fontsize.axis_title Integer, font size of the axis label and text. Default: 10.
#' @param fontsize.axis_text Integer, font size of the axis label and text. Default: 8.
#'
#' @return A ggplot object that can be visualized by "`p`" or `ggsave(file = "output.pdf", p)`
#' @export
#'
#' @examples
#' ## 1. bar plot grouped by clusters ("clusterID") and colored by true labels ("true_label)
#' p_bar <- draw_barplot(input_eset = clustered.eset, group_by = "clusterID", color_by = "true_label")
#'
#' ## 2. customize the colors
#' p_bar <- draw_barplot(input_eset = clustered.eset, group_by = "clusterID", color_by = "true_label", colors = c("green", "red", "blue", "grey", "orange", "purple", "yellow"))
draw_barplot <- function(input_eset,
                         group_by = "clusterID",
                         color_by = "cell_type",
                         colors = NULL,
                         legend.position = "right",
                         xlabel.angle = 0,
                         fontsize.legend_title = 12, fontsize.legend_text = 10,
                         fontsize.axis_title = 12, fontsize.axis_text = 10)
{
  ## check parameters
  if ((group_by %in% colnames(Biobase::pData(input_eset))) == FALSE) {stop('The group_by column,', group_by, ', was not found in the input set. Please check and re-try.')}
  if ((color_by %in% colnames(Biobase::pData(input_eset))) == FALSE) {stop('The color_by column,', color_by, ', was not found in the input set. Please check and re-try.')}

  pd <- Biobase::pData(input_eset)
  pd[,group_by] <- as.factor(pd[,group_by])
  pd[,color_by] <- as.factor(pd[,color_by])

  hjust_x = base::ifelse (xlabel.angle == 0, 0.5, 1)
  p <- ggplot(data = pd, aes(pd[,group_by])) + theme_classic() + geom_bar(aes(fill = pd[,color_by]), position = "fill") +
    labs(x = group_by, y = "Percent Proportion") + guides(fill = guide_legend(title = color_by)) +
    theme(
      legend.position = legend.position,
      legend.title = element_text(size = fontsize.legend_title, face = "bold"),
      legend.text = element_text(size = fontsize.legend_text),
      axis.title = element_text(size = fontsize.axis_title, face = "bold", hjust = 0.5, color = "black"),
      axis.text.x = element_text(size = fontsize.axis_text, hjust = hjust_x, color = "black", angle = xlabel.angle),
      axis.text.y = element_text(size = fontsize.axis_text, hjust = 1, color = "black"))

  if (is.null(colors) == FALSE) {p <- p + scale_fill_manual(values = colors)}

  return(p)
}


#' Bubble plot showing the signature scores by self-defined groups
#'
#' @description
#' This function is used to draw a bubble plot of signature scores among self-defined groups.
#'
#' @param input_eset The expression set object that filtered, normalized and log-transformed
#' @param signature_table A matrix or data frame containing three columns: signature_name, signature_feature, weight. Default: `NULL`.
#' @param group_by Character, name of the column for grouping, usually the column of cell types or clusters. Default: "`clusterID`".
#' @param colors A vector of two colors indicating the low and high values respectively. Default: c("`lightgrey`", "`red`").
#' @param legend.position Character, position of legend: "`right`" (the default), "`left`", "`top`", "`bottom`" or "`none`".
#' @param fontsize.legend_title Integer, font size of the legend title. Default: 10.
#' @param fontsize.legend_text Integer, font size of the legend text. Default: 8.
#' @param fontsize.axis_title Integer, font size of the axis label and text. Default: 10.
#' @param fontsize.axis_text Integer, font size of the axis label and text. Default: 8.
#' @param xlabel.angle Numeric, the angle of the a-axis title. When it's set not 0, the x-axis text will automatically right-justified. Default: 0.
#'
#' @return A ggplot object of bubble plot
#' @export
#'
#' @examples
#' marker_file <- system.file('PBMC14KDS_DemoDataSet/DATA/', 'Immune_signatures.xlsx', package = "scMINER")
#' signature_table <- openxlsx::read.xlsx(marker_file)
#' head(signature_table)
#' ## 1. the most commonly used command
#' p_bubbleplot <- draw_bubbleplot(input_eset = clustered.eset, signature_table = signature_table, group_by = "clusterID")
#'
#' ## 2. customize the colors
#' p_bubbleplot <- draw_bubbleplot(input_eset = clustered.eset, signature_table = signature_table, group_by = "clusterID", colors = c("lightgrey", "red"))
draw_bubbleplot <- function(input_eset,
                            signature_table = NULL,
                            group_by = "clusterID",
                            colors = NULL,
                            legend.position = "right", fontsize.legend_title = 10, fontsize.legend_text = 8,
                            fontsize.axis_title = 10, fontsize.axis_text = 8, xlabel.angle = 0)
{
  ## check parameters
  if (is.null(signature_table) == TRUE) {
    stop('The "signature_table" is not found. Please check and re-try.')
  } else {
    if (all(c("signature_name","signature_feature","weight") %in% colnames(signature_table))) {
      signature_table <- signature_table[, c("signature_name","signature_feature","weight")]
      signature_table <- signature_table[!duplicated(signature_table),]
      features_in <- unique(signature_table$signature_feature[signature_table$signature_feature %in% row.names(input_eset)])
      if (length(features_in) == 0) {stop('None of the features specified in "signatures" was found in the input eset. Please check and re-try!')}

      signature_table <- signature_table[signature_table$signature_feature %in% features_in,]
      cat(length(features_in), "features of", length(unique(signature_table$signature_name)), "signatures were found in the input eset and will be used in calculation.\n")
    } else {
      stop('The matrix specified by "signatures" must contain these three columns: signature_name, signature_feature, and weight. Please check and re-try.')
    }
  }

  if ((group_by %in% colnames(Biobase::pData(input_eset))) == FALSE) {stop('The group_by column,', group_by, ', was not found in the input set. Please check and re-try.')}

  ## prepare the master table for visualization
  exp <- base::apply(Biobase::exprs(input_eset), 2, z_normalization)
  exp <- as.matrix(exp[features_in, , drop = FALSE])
  exp <- exp[, colSums(is.na(exp)) == 0] # when all gene expressions of one cell is exactly same, NaN happens.
  pd <- Biobase::pData(input_eset)[, group_by, drop = FALSE]
  pd <- pd[colnames(exp), , drop = FALSE]

  # get the raw signature scores
  signatures <- unique(signature_table$signature_name)
  signature_score <- matrix(NA, nrow = ncol(exp), ncol = length(signatures), dimnames = list(colnames(exp), signatures))
  for (i in 1:length(signatures)) {
    signature_table.sel <- signature_table[signature_table$signature_name %in% c(signatures[i]), , drop = FALSE]
    n <- nrow(signature_table.sel)

    if (n > 1) {
      mat <- t(exp[signature_table.sel$signature_feature,]) %*% as.vector(signature_table.sel$weight)
      signature_score[,signatures[i]] <- mat[,1]/n
    } else if (n == 1) {
      signature_score[,signatures[i]] <- exp[signature_table.sel$signature_feature,] * signature_table.sel$weight
    }
  }

  # get the scaled signature scores
  signature_score.scaled <- base::apply(signature_score, 2, z_normalization)

  groups <- sort(unique(pd[, group_by]))
  df_mean <- matrix(NA, nrow = length(groups), ncol = ncol(signature_score.scaled), dimnames = list(groups, colnames(signature_score.scaled)))
  df_pct <- matrix(NA, nrow = length(groups), ncol = ncol(signature_score.scaled), dimnames = list(groups, colnames(signature_score.scaled)))
  for (i in 1:length(groups)) {
    pd.sel <- pd[pd[, group_by] %in% c(groups[i]), , drop = FALSE]
    signature_score.scaled.sel <- signature_score.scaled[row.names(pd.sel), , drop = FALSE]
    df_mean[groups[i],] <- colMeans(signature_score.scaled.sel, na.rm = TRUE)
    df_pct[groups[i],] <- apply(signature_score.scaled.sel, 2, function(x) {sum(x > 0) / length(x) * 100})
  }

  df_mean.long <- reshape2::melt(df_mean); colnames(df_mean.long) <- c(group_by, "Signatures", "Mean_value")
  df_pct.long <- reshape2::melt(df_pct); colnames(df_pct.long) <- c(group_by, "Signatures", "Percentage_value")
  master <- dplyr::left_join(df_mean.long, df_pct.long, by = c(group_by, "Signatures"))

  # visualize
  hjust_x = base::ifelse (xlabel.angle == 0, 0.5, 1)
  if (is.null(colors) == TRUE) {color.low <- "lightgrey"; color.high <- "red"} else {if (length(colors) == 2) {color.low <- colors[1]; color.high <- colors[2]} else {stop('The length of "colors" must be of 2.')}}
  p <- ggplot(master, aes(x = as.factor(master[, group_by]), y = Signatures)) + geom_point(aes(color = Mean_value, size = Percentage_value), pch = 16) +
    scale_colour_gradient(low = color.low, high = color.high) + theme_classic() + labs(x = group_by, y = "Signatures") +
    theme(
      legend.position = legend.position,
      legend.title = element_text(size = fontsize.legend_title, face = "bold"),
      legend.text = element_text(size = fontsize.legend_text),
      axis.title = element_text(size = fontsize.axis_title, face = "bold", hjust = 0.5, color = "black"),
      axis.text.x = element_text(size = fontsize.axis_text, hjust = hjust_x, color = "black", angle = xlabel.angle),
      axis.text.y = element_text(size = fontsize.axis_text, hjust = 1, color = "black"))

  return(p)
}


#' Prepare the standard input files for scMINER Portal
#'
#' @description
#' This function is used to generated the standard input files that can be directly uploaded to scMINER Portal for visualization and beyond.
#'
#' @param input_expression.eset A eset object of expression data which has been filerted, normalized and log-transformed. Default: `NULL`.
#' @param input_expression.seuratObj A Seurat object with reduction results (umap and/or tsne). Default: `NULL`.
#' @param group_by Character, name of the column in phenoData of eset object or meta.data of Seurat Obj that defines the groups for colorcoding in scMINER Portal. Default: `NULL`.
#' @param input_activity.eset A eset object of activity data calculated by scMINER. Default: `NULL`.
#' @param input_network.dir Character, the path to the SJARACNe directory. Default: `NULL`.
#' @param input_network.table A table of network information, it contains at least three column: "`CellGroup`" (name of groups), "`NetworkType`" (type of network, TF or SIG) and "`NetworkFile`" (path to network files). Default: `NULL`.
#' @param output_dir Character, the path to the output directory. Default: `NULL`.
#'
#' @return This function generated 1-3 standard files which can be uploaded to scMINER Portal directly.
#' @export
#'
#' @examples
#' ## 1. the most commonly used command
#' generatePortalInputs(input_expression.eset = expression_clustered.eset, group_by = "cellType", input_activity.eset = activity_clustered.eset, input_network.dir = "./SJARACNe", output_dir = "./scMINERportal")
#'
#' ## 2. prepare expression data from Seurat object ("pbmc14.obj")
#' generatePortalInputs(input_expression.seuratObj = pbmc14.obj, output_dir = "./path-to-output_dir")
#'
#' ## 3. prepare network data from a table
#' network.table <- data.frame(CellGroup = c("CD4__CD25_T_Reg", "CD4__CD25_T_Reg", "CD19__B", "CD19__B"),NetworkType = c("SIG", "TF", "SIG", "TF"),
#'                             NetworkFile = c("./sjaracne/CD4__CD25_T_Reg/SIG/b100_pce-3/sjaracne_workflow-1474c41b-067b-4f86-ab99-09f73dadb16g/consensus_network_ncol_.txt",
#'                                             "./sjaracne/CD4__CD25_T_Reg/TF/b100_pce-3/sjaracne_workflow-a93cd6db-7253-4ffb-ae4e-633b9dedf11d/consensus_network_ncol_.txt",
#'                                             "./sjaracne/CD19__B/SIG/sjaracne_workflow-da0c3c72-7afb-44fa-973b-e4d767e20b6f/consensus_network_ncol_.txt",
#'                                             "./sjaracne/CD19__B/TF/sjaracne_workflow-0426ea12-10bf-428c-b199-d5bd1a7aab5f/consensus_network_ncol_.txt"))
#' generatePortalInputs(input_expression.eset = expression_clustered.eset, group_by = "cellType", input_network.table = network.table, output_dir = "./path-to-output_dir")
generatePortalInputs <- function(input_expression.eset = NULL,
                                 input_expression.seuratObj = NULL,
                                 group_by = NULL,
                                 input_activity.eset = NULL,
                                 input_network.dir = NULL,
                                 input_network.table = NULL,
                                 output_dir = NULL)
{
  ## check output directory
  if (is.null(output_dir)) {
    stop('Please specify output_dir and re-try.')
  } else {
    if (dir.exists(output_dir) == FALSE) {
      dir.create(output_dir, recursive = TRUE)
      cat("The output directory has been created: ", output_dir, ".\n")
    }
  }

  ## prepare expression data
  if (is.null(input_expression.eset) == TRUE) {
    if (is.null(input_expression.seuratObj) == TRUE) {
      cat("The expression data preparation for scMINER Portal has been skipped, since no inputs was provided.\n")
    } else {
      cat("The expression data for scMINER Portal will be prepared from the input Seurat object...\n")

      # umap or tsne
      cat("\tChecking the reduction results ...")
      if (is.null(input_expression.seuratObj@reductions$umap) == TRUE) {
        if (is.null(input_expression.seuratObj@reductions$tsne) == TRUE) {
          stop('No reduction results was found in the input Seurat object. Both "input_expression.seuratObj@reductions$umap" and "input_expression.seuratObj@reductions$tsne" are NULL. Please check and re-try.')
        } else {
          seurat_reductions <- data.frame(CellID = row.names(input_expression.seuratObj@reductions$tsne@cell.embeddings), input_expression.seuratObj@reductions$tsne@cell.embeddings)
        }
      } else {
        if (is.null(input_expression.seuratObj@reductions$tsne) == TRUE) {
          seurat_reductions <- data.frame(CellID = row.names(input_expression.seuratObj@reductions$umap@cell.embeddings), input_expression.seuratObj@reductions$umap@cell.embeddings)
        } else {
          umap <- data.frame(CellID = row.names(input_expression.seuratObj@reductions$umap@cell.embeddings), input_expression.seuratObj@reductions$umap@cell.embeddings)
          tsne <- data.frame(CellID = row.names(input_expression.seuratObj@reductions$tsne@cell.embeddings), input_expression.seuratObj@reductions$tsne@cell.embeddings)
          seurat_reductions <- merge(umap, tsne, by = "CellID")
        }
      }
      cat("\t\tDone! The reduction results are ready.\n")

      # meta data
      cat("\tChecking the meta.data ...")
      seurat_metadata <- data.frame(CellID = row.names(input_expression.seuratObj@meta.data), input_expression.seuratObj@meta.data, CellGroup = Seurat::Idents(input_expression.seuratObj))
      if (is.null(group_by) == FALSE) {
        if (group_by %in% colnames(seurat_metadata)) {
          seurat_metadata[, "CellGroup"] <- seurat_metadata[, group_by]
        } else {
          stop("The group_by was not found in meta.data. Please check and re-try.")
        }
      }
      cat("\t\tDone! The meta.data are ready.\n")

      # expression matrix
      cat("\tChecking the normalized expression matrix...\n")
      if (all(input_expression.seuratObj[['RNA']]@counts@x == input_expression.seuratObj[['RNA']]@data@x)) {
        cat("\tThe data matrix in input_expression.seuratObj@assays$RNA@data is not normalized. scMINER will normalize and log-transforme it...")
        expression.eset <- createSparseEset(input_matrix = input_expression.seuratObj[['RNA']]@counts, cellData = merge(seurat_metadata, seurat_reductions, by = "CellID"), addMetaData = FALSE)
        expression.eset <- normalizeSparseEset(input_eset = expression.eset)
        cat("\t\tDone! The data matrix has been normalized.\n")
        saveRDS(expression.eset, file = paste0(output_dir, "/expression.eset"))
        cat("The expression data for scMINER Portal has been generated:", paste0(output_dir, "/expression.eset"),".\n")
      }
    }
  } else {
    if (is.null(input_expression.seuratObj)) {
      cat("The expression data for scMINER Portal will be prepared from the input eset object ...\n")

      ## check the expression matrix
      if (all(Biobase::exprs(input_expression.eset) %% 1 == 0)) {
        cat("All values in the expression matrix of the input eset are integers. Please confirm the input eset is normalized and log-transformed.\n")
      }

      ## check pheno data
      pd <- Biobase::pData(input_expression.eset)
      if (("CellID" %in% colnames(pd)) == FALSE) {pd$CellID <- row.names(pd)}

      cat("\tChecking the reduction results ...")
      reductions_col <- grepl("umap|tsne", colnames(pd), ignore.case = TRUE)
      if (length(reductions_col) > 0) {
        cat("Done! These columns will be used as the reductions results:", paste0(reductions_col, collapse = ", "), ".\n")
      } else {
        stop('No reduction results was found in the input eset object. Failed to retrieve the columns of reduction results by "umap|tsne". Please check and retry.')
      }

      cat("\tChecking the group info ...")
      if (is.null(group_by) == FALSE) {
        if (group_by %in% colnames(pd)) {
          pd[, "CellGroup"] <- pd[, group_by]
        } else {
          stop("The group_by was not found in input eset object. Please check and re-try.")
        }
        cat("\t\tDone! The CellGroup column had been set by:", group_by, ".\n")
      } else {
        if ("clusterID" %in% colnames(pd)) {
          pd[, "CellGroup"] <- pd[, "clusterID"]
          cat("\t\tDone! The 'group_by' is not specified and the default value is used: 'clusterID'. To change it, please specify 'group_by' argument.\n")
        } else {
          stop('Please specify "group_by" to define the group info.')
        }
      }

      expression.eset <- updateSparseEset(input_expression.eset, cellData = pd, addMetaData = FALSE)
      saveRDS(expression.eset, file = paste0(output_dir, "/expression.rds"))
      cat("The expression data for scMINER Portal has been generated:", paste0(output_dir, "/expression.rds"),".\n")
    } else {
      stop('Please use ONLY ONE of "input_expression.eset" or "input_expression.seuratObj" to specify the input expression data.')
    }
  }

  ## prepare activity data
  if (is.null(input_activity.eset) == FALSE) {
    cat("The activity data for scMINER Portal will be prepared from the input eset object ...\n")

    fd <- Biobase::fData(input_activity.eset)
    if (("GeneSymbol" %in% colnames(fd)) == FALSE) {fd$GeneSymbol <- gsub("[_TF|_SIG]", "", row.names(fd))}

    pd <- Biobase::pData(input_activity.eset)
    if (("CellID" %in% colnames(pd)) == FALSE) {pd$CellID <- row.names(pd)}

    cat("\tChecking the reduction results ...")
    reductions_col <- grepl("umap|tsne", colnames(pd), ignore.case = TRUE)
    if (length(reductions_col) > 0) {
      cat("Done! These columns will be used as the reductions results:", paste0(reductions_col, collapse = ", "), ".\n")
    } else {
      stop('No reduction results was found in the input eset object. Failed to retrieve the columns of reduction results by "umap|tsne". Please check and retry.')
    }

    cat("\tChecking the group info ...")
    if (is.null(group_by) == FALSE) {
      if (group_by %in% colnames(pd)) {
        pd[, "CellGroup"] <- pd[, group_by]
      } else {
        stop("The group_by was not found in input eset object. Please check and re-try.")
      }
      cat("\t\tDone! The CellGroup column had been set by:", group_by, ".\n")
    } else {
      if ("clusterID" %in% colnames(pd)) {
        pd[, "CellGroup"] <- pd[, "clusterID"]
        cat("\t\tDone! The CellGroup column had been set by: clusterID.\n")
      } else {
        stop('Please specify "group_by" to define the group info.')
      }
    }

    activity.eset <- updateSparseEset(input_activity.eset, cellData = pd, featureData = fd, addMetaData = FALSE)
    saveRDS(activity.eset, file = paste0(output_dir, "/activity.rds"))
    cat("The activity data for scMINER Portal has been generated:", paste0(output_dir, "/activity.rds"),".\n")
  } else {
    cat("The activity data preparation for scMINER Portal has been skipped, since no inputs was provided.\n")
  }

  ## prepare network data
  if (is.null(input_network.dir) == FALSE) {
    if (is.null(input_network.table) == FALSE) {
      stop('Please use ONLY ONE of "input_network.dir" or "input_network.table" to specify the input network data.')
    } else {
      cat("The network data for scMINER Portal will be prepared from the specified directory ...\n")
      if (dir.exists(input_network.dir) == FALSE) {stop('The input directory for network files was not found: ', input_network.dir, '.\n')}

      grps <- list.dirs(input_network.dir, recursive = FALSE, full.names = FALSE)
      cat("\t", length(grps), "groups were found in the input network directory.\n")

      network_merged <- data.frame()
      for (i in 1:length(grps)) {
        dir.tmp <- paste0(input_network.dir, "/", grps[i])
        network.tf <- list.files(path = paste0(dir.tmp, "/TF"), pattern = "consensus_network_ncol_.txt", recursive = TRUE, full.names = TRUE)
        if (length(network.tf) == 1) {
          net.tf <- read.delim(file = network.tf, stringsAsFactors = FALSE)
          net.tf$CellGroup <- grps[i]; net.tf$NetworkType <- "TF";
        } else if (length(network.tf) < 1) {
          cat("\t\tNo TF network file was found for group:", grps[i], ".\n")
        } else {
          stop('Multiple TF network files were found for group:', grps[i], '.')
        }

        network.sig <- list.files(path = paste0(dir.tmp, "/SIG"), pattern = "consensus_network_ncol_.txt", recursive = TRUE, full.names = TRUE)
        if (length(network.sig) == 1) {
          net.sig <- read.delim(file = network.sig, stringsAsFactors = FALSE)
          net.sig$CellGroup <- grps[i]; net.sig$NetworkType <- "SIG";
        } else if (length(network.sig) < 1) {
          cat("\t\tNo SIG network file was found for group:", grps[i], ".\n")
        } else {
          stop('Multiple SIG network files were found for group:', grps[i], '.')
        }

        network.grp <- rbind(net.tf, net.sig)

        if (i == 1) {network_merged <- network.grp} else {network_merged <- rbind(network_merged, network.grp)}
      }

      write.table(network_merged, file = paste0(output_dir, "/networks.txt"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
      cat("The network data for scMINER Portal has been generated:", paste0(output_dir, "/networks.txt"),".\n")
    }
  } else {
    if (is.null(input_network.table) == FALSE) {
      cat("The network data for scMINER Portal will be prepared from the specified table ...\n")

      col.missed <- c("CellGroup", "NetworkType", "NetworkFile")[! c("CellGroup", "NetworkType", "NetworkFile") %in% colnames(input_network.table)]
      if (length(col.missed) > 0) {stop('Part of the required columns were not found in the input network table: ', paste0(col.missed, ", "), '.')}

      type.unrecognized <- unique(input_network.table$NetworkType[! input_network.table$NetworkType %in% c("TF", "SIG")])
      if (length(type.unrecognized) > 0) {stop('Part of the network type was not recognized: ', paste0(type.unrecognized, ", "), '. Should be [TF | SIG].')}

      network.duplicated <- input_network.table$NetworkFile[duplicated(input_network.table$NetworkFile)]
      if (length(network.duplicated) > 0) {stop('These network files are duplicated: ', paste0(network.duplicated, ", "), '.')}

      grps <- unique(input_network.table$CellGroup)
      cat("\t", length(grps), "groups were found in the input network table.\n")
      network_merged <- data.frame()
      for (i in 1:nrow(input_network.table)) {
        file.tmp <- input_network.table$NetworkFile[i]
        if (file.exists(file.tmp)) {
          net.tmp <- read.delim(file = file.tmp, stringsAsFactors = FALSE)
          net.tmp$CellGroup <- input_network.table$CellGroup[i]; net.tmp$NetworkType <- input_network.table$NetworkType[i]
        } else {
          cat("\t\tThis network file is not found:", file.tmp, ".\n")
        }

        if (i == 1) {network_merged <- net.tmp} else {network_merged <- rbind(network_merged, net.tmp)}
      }

      write.table(network_merged, file = paste0(output_dir, "/networks.txt"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
      cat("The network data for scMINER Portal has been generated:", paste0(output_dir, "/networks.txt"),".\n")
    } else {
      cat("The network data preparation for scMINER Portal has been skipped, since no inputs was provided.\n")
    }
  }
}
