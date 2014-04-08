\name{SAMseq.createRmd}
\alias{SAMseq.createRmd}
\title{Generate a \code{.Rmd} file containing code to perform differential expression analysis with SAMseq}
\usage{
SAMseq.createRmd(data.path, result.path, codefile)
}
\arguments{
  \item{data.path}{The path to a .rds file containing the
  \code{compData} object that will be used for the
  differential expression analysis.}

  \item{result.path}{The path to the file where the result
  object will be saved.}

  \item{codefile}{The path to the file where the code will
  be written.}
}
\value{
The function generates a \code{.Rmd} file containing the
code for performing the differential expression analysis.
This file can be executed using e.g. the \code{knitr}
package.
}
\description{
A function to generate code that can be run to perform
differential expression analysis of RNAseq data (comparing
two conditions) using SAMseq. The code is written to a
\code{.Rmd} file. This function is generally not called by
the user, the main interface for performing differential
expression analysis is the \code{\link{runDiffExp}}
function.
}
\details{
For more information about the methods and the
interpretation of the parameters, see the \code{SAMseq}
package and the corresponding publications.
}
\examples{
try(
if (require(samr)) {
tmpdir <- normalizePath(tempdir(), winslash = "/")
mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
                                    samples.per.cond = 5, n.diffexp = 100,
                                    output.file = file.path(tmpdir, "mydata.rds"))
runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "SAMseq",
           Rmdfunction = "SAMseq.createRmd",
           output.directory = tmpdir)
})
}
\author{
Charlotte Soneson
}
\references{
Li J and Tibshirani R (2011): Finding consistent patterns:
a nonparametric approach for identifying differential
expression in RNA-Seq data. Statistical Methods in Medical
Research
}

