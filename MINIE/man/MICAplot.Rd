\name{MICAplot}
\alias{MICAplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot visualization with sample info
}
\description{
ggplot2 function to visualize phenotype data in expressionset
}
\usage{
MICAplot(input_eset = eset, label = metaName, visualize = "UMAP",
          title.size = 5, title.name = "", pct = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_eset}{Input expressionSet with visualization coordinates stored in pData}
  \item{label}{String, Name of Metadata using as color label}
  \item{visualize}{String, used as x,y label,default="UMAP"}
  \item{title.size}{Size of plot title}
  \item{title.name}{Plot title name}
  \item{pct}{Point size}
}
\details{

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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
