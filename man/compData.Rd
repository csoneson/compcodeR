\name{compData}
\alias{compData}
\title{Create a \code{compData} object}
\usage{
compData(count.matrix, sample.annotations, info.parameters,
  variable.annotations = data.frame(), filtering = "no info",
  analysis.date = "", package.version = "", method.names = list(),
  code = "", result.table = data.frame())
}
\arguments{
  \item{count.matrix}{A count matrix, with genes as rows
  and observations as columns.}

  \item{sample.annotations}{A data frame, containing at
  least one column named 'condition', encoding the grouping
  of the observations into two groups. The row names should
  be the same as the column names of the
  \code{count.matrix}.}

  \item{info.parameters}{A list containing information
  regarding simulation parameters etc. The only mandatory
  entries are \code{dataset} and \code{uID}, but it may
  contain entries such as the ones listed below (see
  \code{generateSyntheticData} for more detailed
  information about each of these entries). \itemize{ \item
  \code{dataset}: an informative name or identifier of the
  data set (e.g., summarizing the simulation settings).
  \item \code{samples.per.cond} \item \code{n.diffexp}
  \item \code{repl.id} \item \code{seqdepth} \item
  \code{minfact} \item \code{maxfact} \item
  \code{fraction.upregulated} \item
  \code{between.group.diffdisp} \item
  \code{filter.threshold.total} \item
  \code{filter.threshold.mediancpm} \item
  \code{fraction.non.overdispersed} \item
  \code{random.outlier.high.prob} \item
  \code{random.outlier.low.prob} \item
  \code{single.outlier.high.prob} \item
  \code{single.outlier.low.prob} \item \code{effect.size}
  \item \code{uID}: a unique ID for the data set. In
  contrast to \code{dataset}, the \code{uID} is unique e.g.
  for each instance of replicated data sets generated with
  the same simulation settings. }}

  \item{variable.annotations}{A data frame with variable
  annotations (with number of rows equal to the number of
  rows in \code{count.matrix}, that is, the number of
  variables in the data set). Not mandatory, but may
  contain columns such as the ones listed below. If
  present, the row names should be the same as the row
  names of the \code{count.matrix}. \itemize{ \item
  \code{truedispersions.S1}: the true dispersion for each
  gene in condition S1. \item \code{truedispersions.S2}:
  the true dispersion for each gene in condition S2. \item
  \code{truemeans.S1}: the true mean value for each gene in
  condition S1. \item \code{truemeans.S2}: the true mean
  value for each gene in condition S2. \item
  \code{n.random.outliers.up.S1}: the number of 'random'
  outliers with extremely high counts for each gene in
  condition S1. \item \code{n.random.outliers.up.S2}: the
  number of 'random' outliers with extremely high counts
  for each gene in condition S2. \item
  \code{n.random.outliers.down.S1}: the number of 'random'
  outliers with extremely low counts for each gene in
  condition S1. \item \code{n.random.outliers.down.S2}: the
  number of 'random' outliers with extremely low counts for
  each gene in condition S2. \item
  \code{n.single.outliers.up.S1}: the number of 'single'
  outliers with extremely high counts for each gene in
  condition S1. \item \code{n.single.outliers.up.S2}: the
  number of 'single' outliers with extremely high counts
  for each gene in condition S2. \item
  \code{n.single.outliers.down.S1}: the number of 'single'
  outliers with extremely low counts for each gene in
  condition S1. \item \code{n.single.outliers.down.S2}: the
  number of 'single' outliers with extremely low counts for
  each gene in condition S2. \item \code{M.value}: the
  M-value (observed log2 fold change between condition S1
  and condition S2) for each gene. \item \code{A.value}:
  the A-value (observed average expression level across
  condition S1 and condition S2) for each gene. \item
  \code{truelog2foldchanges}: the true (simulated) log2
  fold changes between condition S1 and condition S2. \item
  \code{upregulation}: a binary vector indicating which
  genes are simulated to be upregulated in condition S2
  compared to condition S1. \item \code{downregulation}: a
  binary vector indicating which genes are simulated to be
  downregulated in condition S2 compared to condition S1.
  \item \code{differential.expression}: a binary vector
  indicating which genes are simulated to be differentially
  expressed in condition S2 compared to condition S1. }}

  \item{filtering}{A character string containing
  information about the filtering that has been applied to
  the data set.}

  \item{analysis.date}{If a differential expression
  analysis has been performed, a character string detailing
  when it was performed.}

  \item{package.version}{If a differential expression
  analysis has been performed, a character string giving
  the version of the differential expression packages that
  were applied.}

  \item{method.names}{If a differential expression analysis
  has been performed, a list with entries \code{full.name}
  and \code{short.name}, giving the full name of the
  differential expression method (may including version
  number and parameter settings) and a short name or
  abbreviation.}

  \item{code}{If a differential expression analysis has
  been performed, a character string containing the code
  that was run to perform the analysis. The code should be
  in R markdown format, and can be written to an HTML file
  using the \code{\link{generateCodeHTMLs}} function.}

  \item{result.table}{If a differential expression analysis
  has been performed, a data frame containing the results
  of the analysis. The number of rows should be equal to
  the number of rows in \code{count.matrix} and if present,
  the row names should be identical. The only mandatory
  column is \code{score}, which gives a score for each
  gene, where a higher score suggests a "more highly
  differentially expressed" gene. Different comparison
  functions use different columns of this table, if
  available. The list below gives the columns that are used
  by the interfaced methods. \itemize{ \item \code{pvalue}
  nominal p-values \item \code{adjpvalue} p-values adjusted
  for multiple comparisons \item \code{logFC} estimated
  log-fold changes between the two conditions \item
  \code{score} the score that will be used to rank the
  genes in order of significance. Note that high scores
  always signify differential expression, that is, a strong
  association with the predictor. For example, for methods
  returning a nominal p-value the score can be defined as 1
  - pvalue. \item \code{FDR} false discovery rate estimates
  \item \code{posterior.DE} posterior probabilities of
  differential expression \item \code{prob.DE} conditional
  probabilities of differential expression \item
  \code{lfdr} local false discovery rates \item
  \code{statistic} test statistics from the differential
  expression analysis \item \code{dispersion.S1} dispersion
  estimates in condition S1 \item \code{dispersion.S2}
  dispersion estimates in condition S2 }}
}
\value{
A \code{compData} object.
}
\description{
The \code{compData} class is used to store information
about the experiment, such as the count matrix, sample and
variable annotations, information regarding the generation
of the data and results from applying a differential
expression analysis to the data. This constructor function
creates a \code{compData} object.
}
\examples{
count.matrix <- round(matrix(1000*runif(4000), 1000))
sample.annotations <- data.frame(condition = c(1, 1, 2, 2))
info.parameters <- list(dataset = "mydata", uID = "123456")
cpd <- compData(count.matrix, sample.annotations, info.parameters)
}
\author{
Charlotte Soneson
}

