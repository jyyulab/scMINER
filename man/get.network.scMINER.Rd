% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functionsFromNetBID2.R
\name{get.network.scMINER}
\alias{get.network.scMINER}
\title{Read SJARACNe Network Result and Return it as List Object(adapted from NetBID2)}
\usage{
get.network.scMINER(network_file = NULL)
}
\arguments{
\item{network_file}{character, the path for storing network file. For the output of SJAracne, the name of the network file will be "consensus_network_ncol_.txt" under the output directory.}
}
\value{
Return a list containing three elements, \code{network_dat}, \code{target_list} and \code{igraph_obj}.
}
\description{
\code{get.network.scMINER} reads SJARACNe network construction result and returns a list object
with network data frame, driver-to-target list and igraph object wrapped inside.
}
\details{
In the demo, "consensus_network_ncol_.txt" file will be read and convert into a list object.
This list contains three elements, \code{network_data}, \code{target_list} and \code{igraph_obj}.
\code{network_dat} is a data.frame, contains all the information of the network SJARACNe constructed.
\code{target_list} is a driver-to-target list object. Please check details in \code{get_net2target_list}.
\code{igraph_obj} is an igraph object used to save this directed and weighted network.
Each edge of the network has two attributes, \code{weight} and \code{sign}.
\code{weight} is the "MI (mutual information)" value and \code{sign} is the sign of the spearman
correlation coefficient (1, positive regulation; -1, negative regulation).
}
\examples{
\dontrun{
tf.network  <- get.network.scMINER(network_file="consensus_network_ncol_.txt")
}
}
