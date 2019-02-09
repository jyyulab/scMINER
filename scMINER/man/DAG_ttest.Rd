\name{DAG_ttest}
\alias{DAG_ttest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
inner function of FindDAG
}
\usage{
DAG_ttest(d, group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{group}{
%%     ~~Describe \code{group} here~~
}
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (d, group)
{
    d.tmp <- unlist(d[, -1])
    d.tmp <- data.frame(acs = d.tmp, group = group, stringsAsFactors = FALSE)
    d.sel <- d.tmp[complete.cases(d.tmp), ]
    d.sel$group <- as.factor(d.sel$group)
    n.cases <- nlevels(d.sel$group)
    if (n.cases == 2) {
        res <- t.test(dplyr::filter(d.tmp, group == "Aim")$acs,
            dplyr::filter(d.tmp, group == "Ctrl")$acs)
        pval.t <- res$p.value
        t.stat <- res$statistic
        rs.t <- c(id = d$acs.id, t.stat, pval = pval.t, z = sign(pval.t) *
            abs(qnorm(pval.t/2)), res$parameter, CI.low = res$conf.int[1],
            CI.high = res$conf.int[2], MeanAct.Aim = unname(res$estimate[1]),
            MeanAct.Ctrl = unname(res$estimate[2]))
    }
    else {
        rs.t <- c(id = d$acs.id, t.stat = NA, pval = NA, z = NA,
            df = NA, CI.low = NA, CI.high = NA, MeanAct.Aim = mean(filter(d.tmp,
                group == "Aim")$acs), MeanAct.Ctrl = mean(filter(d.tmp,
                group == "Ctrl")$acs))
    }
    return(rs.t)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
