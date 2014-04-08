\name{checkDataObject}
\alias{checkDataObject}
\title{Check a list or a \code{compData} object for compatibility with the differential expression functions interfaced by \code{compcodeR}}
\usage{
checkDataObject(data.obj)
}
\arguments{
  \item{data.obj}{A list containing data and condition
  information, or a \code{compData} object.}
}
\description{
Check if a list or a \code{compData} object contains the
necessary slots for applying the differential expression
functions interfaced by the \code{compcodeR} package. This
function is provided for backward compatibility, see also
\code{\link{check_compData}} and
\code{\link{check_compData_results}}.
}
\examples{
mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
                                    samples.per.cond = 5, n.diffexp = 100)
checkDataObject(mydata.obj)
}
\author{
Charlotte Soneson
}

