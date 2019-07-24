\name{readscRNAseqData}
\alias{readscRNAseqData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read raw UMI count scRNA-seq data file
}
\description{
Read data from 10x genomics or other common file type via function read.delim2()
}
\usage{
readscRNAseqData(file,is.10x=TRUE,sep="\t",stringsAsFactors=FALSE,quote="",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{filename/filepath}
  \item{is.10x}{logical, if TRUE, then 'file' need to be the path to 10x genomics outputs}
  \item{...}{other attributes that could be pass to read.delim2() function}
}
\details{

}
\value{
a data frame
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Yu Lab
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
