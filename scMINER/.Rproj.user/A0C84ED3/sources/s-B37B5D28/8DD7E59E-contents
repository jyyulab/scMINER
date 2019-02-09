\name{gene_vlnplot}
\alias{gene_vlnplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gene visualization on violinplots}
\description{
Visualize gene expression/activity in different groups using vlnplots.
}
\usage{
gene_vlnplot(eset, group_tag = "celltype", target, ylabel = "Expression", ncol = 3, drawquantiles = FALSE, title.size = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{A ExpressionSet,which stores group information in pData.}
  \item{group_tag}{A string, Column name in pData(eset) that indicates group info.}
  \item{target}{One or a set of gene name(s) to be visualized.}
  \item{ylabel}{A string, indicating name of value on y axis.}
  \item{ncol}{number of plot per row.}
  \item{drawquantiles}{logical, whether to draw quantile info on violinplot or not.}
  \item{title.size}{title text size}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
ggplot object
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
YuLab
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##
gene_vlnplot (eset=eset.demo, group_tag = "label", target=c("GZMK"))
gene_vlnplot (eset=eset.demo, group_tag ="celltype",target=c("GZMK","GZMA","CD8A"),draw_quantiles=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ vlnplots }% use one of  RShowDoc("KEYWORDS")
\keyword{ scMINER }% __ONLY ONE__ keyword per line
