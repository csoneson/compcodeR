% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\name{check_compData_results}
\alias{check_compData_results}
\title{Check the validity of a \code{compData} result object}
\usage{
check_compData_results(object)
}
\arguments{
\item{object}{A \code{compData} object}
}
\description{
Check the validity of a \code{compData} object containing differential expression results. An object that passes the check can be used as the input for the method comparison functions in \code{compcodeR}.
}
\examples{
tmpdir <- normalizePath(tempdir(), winslash = "/")
mydata <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
                                samples.per.cond = 5, n.diffexp = 100,
                                output.file = file.path(tmpdir, "mydata.rds"))
## Check an object without differential expression results
check_compData_results(mydata)

runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), 
           result.extent = "voom.limma", 
           Rmdfunction = "voom.limma.createRmd", 
           output.directory = tmpdir, norm.method = "TMM")
resdata <- readRDS(file.path(tmpdir, "mydata_voom.limma.rds"))
## Check an object containing differential expression results
check_compData_results(resdata)
}
\author{
Charlotte Soneson
}
