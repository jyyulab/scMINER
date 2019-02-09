\name{readMICAoutput}
\alias{readMICAoutput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{readMICAoutput}
%%  ~~function to do ... ~~

\description{
Read MICA I/O and create an Expression set:
Read MICA input expression matrix and result file that MICA outputs,
including gene expression matrix, t-SNE coordinates and clustering labels(default as TRUE).
}
\usage{
readMICAoutput(input_file,output_file,load_clust_label=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_file}{Cell by gene matrix in .txt.}
  \item{output_file}{Clustering result .txt file including tSNE coordinates and clustering label.}
  \item{load_clust_label}{logical, define whether to store clustering label under eset$label.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
\item{eset}{Expression Set with tSNE and clustering label information stored in pData.}

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

\examples{
## Only takes txt files
eset<-readMICAoutput(input_file="demo_input.txt",
                     load_clust_label=TRUE,
                     output_file="demo.ggplot.txt")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MICA}% use one of  RShowDoc("KEYWORDS")
\keyword{scMINER}% __ONLY ONE__ keyword per line
