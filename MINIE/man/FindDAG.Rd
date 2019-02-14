\name{FindDAG}
\alias{FindDAG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find differential activity genes from activity matrix
}
\description{
Conduct two_sided t.test on all genes in specific group VS Others to find differential activity genes, a table with essential statistics will be outputted.
}
\usage{
#Try find DAG for all groups
FindDAG(eset, group_tag="group")

#Try find DAG for only group 1
FindDAG(eset, group_tag="group",group_case="group1")

}

\arguments{
  \item{eset}{ExpressionSet that stores group information in pData}
  \item{group_tag}{a character string, column name in pData(eset) that indicates group info.}
  \item{group_case}{a character string, particular group name stored under group_tag}
}
\details{
For scRNA-seq data, due to its large number of cells, p value is often refered as "0" due to the calculation capacity in R. This directly affects the z score calculation. Thus, if pvalue = 0, z scores are mannually set to +/- 40.}

\value{
Output will be a table containing following information:
  \item{id}{GeneID(gene symbol or isoform ID)}
  \item{FuncType}{whether it is TF or not}
  \item{degree}{Number of targets of this gene in each group specific network}
  \item{t}{t statistics outputted from two-sided t.test}
  \item{pval}{p value outputted from two-sided t.test}
  \item{z}{z score}
  \item{MeanAct}{Mean activity value of this gene in specific group}
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
\code{\link{help}}
}
\examples{
#eset should be a Expressionset Obj

FindDAG(eset, group_tag="group")

#Try find DAG for only group 1
FindDAG(eset, group_tag="group",group_case="group1")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FindDAG }% use one of  RShowDoc("KEYWORDS")
