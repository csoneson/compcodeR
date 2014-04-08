\name{HTSDiff.createRmd}
\alias{HTSDiff.createRmd}
\title{Generate a \code{.Rmd} file containing code to perform differential expression analysis with HTSDiff}
\usage{
HTSDiff.createRmd(data.path, result.path, codefile, norm.method,
  DEclusters = 4, epsilon = 0.8)
}
\arguments{
  \item{data.path}{The path to a .rds file containing the
  \code{compData} object that will be used for the
  differential expression analysis.}

  \item{result.path}{The path to the file where the result
  object will be saved.}

  \item{codefile}{The path to the file where the code will
  be written.}

  \item{norm.method}{The between-sample normalization
  method used to compensate for varying library sizes and
  composition in the differential expression analysis.
  Possible values are \code{"TMM"}, \code{"TC"},
  \code{"UQ"}, \code{"Median"} and \code{"DESeq"}.}

  \item{epsilon}{Cutoff on the log-fold change to consider
  a gene differentially expressed. Default is 0.8.}

  \item{DEclusters}{The number of clusters used to
  represent differentially expressed genes. Default is 4.}
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
two conditions) with the \code{HTSDiff} package. The code
is written to a \code{.Rmd} file. This function is
generally not called by the user, the main interface for
performing differential expression analysis is the
\code{\link{runDiffExp}} function.
}
\details{
For more information about the methods and the
interpretation of the parameters, see the \code{HTSDiff}
package and the corresponding publications.
}
\examples{
try(
if (require(HTSDiff)) {
tmpdir <- normalizePath(tempdir(), winslash = "/")
mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
                                    samples.per.cond = 5, n.diffexp = 100,
                                    output.file = file.path(tmpdir, "mydata.rds"))
runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "HTSDiff",
           Rmdfunction = "HTSDiff.createRmd",
           output.directory = tmpdir, norm.method = "TMM")
})
}
\author{
Charlotte Soneson
}
\references{
Balzergue S, Rigaill G, Brunaud V, Blondet E, Rau A, Rogier
O, Caius J, Maugis-Rabusseau C, Soubigou-Taconnat L,
Aubourg S, Lurin C, Delannoy E and Martin-Magniette M-L.
HTSDiff: a sensitive analysis of differential gene
expression by RNA-seq. Submitted.
}

