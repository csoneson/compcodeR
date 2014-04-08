\name{check_compData}
\alias{check_compData}
\title{Check the validity of a \code{compData} object}
\usage{
check_compData(object)
}
\arguments{
  \item{object}{A \code{compData} object}
}
\description{
Check the validity of a \code{compData} object. An object
that passes the check can be used as the input for the
differential expression analysis methods interfaced by
\code{compcodeR}.
}
\examples{
mydata <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
                                samples.per.cond = 5, n.diffexp = 100)
check_compData(mydata)
}
\author{
Charlotte Soneson
}

