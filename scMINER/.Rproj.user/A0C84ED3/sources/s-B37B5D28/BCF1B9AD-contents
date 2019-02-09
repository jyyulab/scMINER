\name{generateSJARACNeInput}
\alias{generateSJARACNeInput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bulk generating of SJARACNe Input across all clusters}
\description{
Generate separate SJARACNe Input files for all clusters/groups.
}
\usage{
generateSJARACNeInput(eset = eset.demo, tf.ref = ref, wd.src = "./", group_tag = "celltype")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{ExpressionSet with group/cluster info.}
  \item{tf.ref}{A vector of reference TF.}
  \item{wd.src}{path to store SJAracne input files.}
  \item{group_tag}{String, name of group to separate expression matrix.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
tf.ref should be cordinated with FeatureNames of eset}
\value{
SJARACNe input files for each subgroups
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
  \code{\link{SJARACNeInput_scRNAseq}}}
\examples{
#transcription factor list from knowledge base

tf.ref<-c("TCF7","LEF1","FOXO1","FOXP1")
generateSJARACNeInput(eset = eset.demo,
                      tf.ref = ref,
                      wd.src = "./",
                      group_tag = "celltype")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SJARACNe }% use one of  RShowDoc("KEYWORDS")

