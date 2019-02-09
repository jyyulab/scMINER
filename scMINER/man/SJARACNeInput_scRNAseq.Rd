\name{SJARACNeInput_scRNAseq}
\alias{SJARACNeInput_scRNAseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SJAracneInput_scRNAseq}
\description{Inner function to generate SJArance input file from ExpressionSet data for one subgroup.}
\usage{
SJARACNeInput_scRNAseq(eset.sel=eset.demo,tf.ref=ref,wd.src="./",grp.tag="Cluster1")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset.sel}{ExpressionSet to generate SJaracne input}
  \item{tf.ref}{A vector of reference transcription factors}
  \item{wd.src}{path to store SjAracne input.}
  \item{group.tag}{name of group for identification}
}
}
\details{
Non-expressed genes in subgroups are filtered.
tf.ref should be coordinate with featureNames(eset.sel).
}
\value{
  \item{Expression file}{Expression matrix ending with .exp extension.}
  \item{TF list}{A list of Transcription factors selected as hub genes for network generation in .txt}

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
\code{\link{help(generateSJARACNeInput)}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SJARACNe }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
