\name{pre.MICA}
\alias{pre.MICA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
scRNA-seq data pre-processing
}
\description{
Conduct quality control,normalizaiton and log2 transformation of scRNA-seq data
}
\usage{
pre.MICA (data.input= NULL, gene_filter=TRUE,cell_filter=TRUE,
          cell_percentage=0.005,ERCC_filter=TRUE,Mito_filter=TRUE,
          UMI_filter=c("both","low","high"),plotting=TRUE,
          plot.dir=".",norm=10e6,sampleID="PBMC_12k",log2Transform=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.input}{Input raw UMI count scRNA-seq expression matrix}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
