\name{gene_heatmap}
\alias{gene_heatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Heatmap visualization for genes}
\description{
Run Heatmap to visualize gene expression/activity profile
}
\usage{
gene_heatmap(eset, target, group_tag = "label", name = "log2_expression",
              save_plot = TRUE, width = 4, height = 8, plot_name = "GeneHeatmap.png")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{ExpressionSet with group information}
  \item{target}{A vector of gene name(s) to visualize.}
  \item{group_tag}{A string, column name in pData(eset) that indicates group info.}
  \item{name}{Name of color scale, default as "log2 Expression".}
  \item{save_plot}{logical, TRUE/FALSE.}
  \item{width}{width of plot to save in "inch".}
  \item{height}{height of plot to save in "inch".}
  \item{plot_name}{plot title, plot device as ".png"}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a ggplot object
}
\references{
%% ~put references to the literature/web site here ~
}
\author{YuLab}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
gene_vlnplot(eset=eset.demo, group_tag="label",target=c("GZMK,CD8A,TCF7"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ heatmap }% use one of  RShowDoc("KEYWORDS")
\keyword{ scMINER }% __ONLY ONE__ keyword per line
