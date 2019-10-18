#'@import Biobase ggplot2 kableExtra knitr openxlsx NetBID2 ComplexHeatmap
#'@importFrom reshape2 melt
#'@importFrom rmarkdown render pandoc_available html_document
#'@importFrom utils read.delim write.table
#'@importFrom Matrix rowSums colSums readMM t
#'@importFrom dplyr filter select starts_with arrange desc left_join
#'@importFrom plyr ddply
#'@importFrom methods as new
#'@importFrom stats IQR aggregate as.dendrogram as.dist cutree density dist fisher.test gaussian glm mad t.test hclust kmeans ks.test lm median model.matrix na.omit order.dendrogram p.adjust pchisq pnorm prcomp pt qnorm quantile sd splinefun complete.cases
#'@importFrom utils read.delim write.table read.table read.delim2
#'@importFrom NetBID2 get.SJAracne.network
#'@importFrom grDevices colorRampPalette dev.off png
#'@importFrom RColorBrewer brewer.pal
#'@importFrom scales hue_pal
#'
#############
library(scales)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)## for colors coding
library(grDevices)

library(Matrix)
library(stats) ##  t.test
library(methods)
library(openxlsx)
library(dplyr) #for easy filtering and apply function
library(plyr)

library(Biobase) ## basic functions for bioconductor

library(rmarkdown)
library(kableExtra) #for Rmarkdown
library(knitr) #for Rmarkdown

library(NetBID2)
