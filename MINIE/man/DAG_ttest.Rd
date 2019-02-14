\name{DAG_ttest}
\alias{DAG_ttest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inner function of FindDAG
}
\description{
inner function of FindDAG, compare mean activity of two unbalanced group
}
\usage{
DAG_ttest(d, group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{A vector of gene expression}
  \item{group}{Vector, Group information}
}
\details{
Inner function of FindDAG
}
\value{
  \item{rs.t}{A vector including id, t.statistics, p value, z value, 0.95 conficence interval, mean activity of two groups}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
YuLab
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
 \code{\link{help(FindDAG)}}
 \code{\link{help(t.test)}}
}
\examples{
group<-as.factor(c(1,1,2,2))
d<-data.frame(id="gene1",c(1,3,4,5),stringAsFactors=FALSE)
DAG_ttest(d,group)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ DAG }% use one of  RShowDoc("KEYWORDS")

