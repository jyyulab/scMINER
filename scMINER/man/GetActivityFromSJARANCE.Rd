\name{GetActivityFromSJARANCe}
\alias{GetActivityFromSJARANCe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate gene activity from SJARACNe
}
\description{
Aggregate network information from SJARACNe and calculate activity score for each hub genes.
}
\usage{
GetActivityFromSJARANCE(SJARACNe_output_path = NA,
                        SJARACNe_input_eset = NA,
                        group_tag = NA,
                        activity.method = "unweighted",
                        activity.norm = TRUE,
                        save_network_file = FALSE,
                        save_path = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SJARACNe_output_path}{Path to SJARACNe output folder(s).}
  \item{SJARACNe_input_eset}{Expressionset that you generate input from.}
  \item{group_tag}{Group names that defines expression matrix separation.}
  \item{activity.method}{c("weighted,unweighted), default to "unweighted".}
  \item{activity.norm}{logical,default to TRUE.}
  \item{save_network_file}{logical,default to FALSE}
  \item{save_path}{Path to save network file}
}
\details{
Group information stored in group_tag should be coordinates with the one used for SJARACNe input generation.
activity method is default to unweighted mean. It calculates mean of all normalized expression from targets.
}
\value{
An expressionset with activity values
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
acs.eset<-GetActivityFromSJARANCE(SJARACNe_output_path = ,
                        SJARACNe_input_eset = eset.demo,
                        group_tag = NA,
                        activity.method = "unweighted",
                        activity.norm = TRUE,
                        save_network_file = FALSE,
                        save_path = NA)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Activity}% use one of  RShowDoc("KEYWORDS")

