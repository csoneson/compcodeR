#' List available *.createRmd functions
#'
#' Print a list of all \code{*.createRmd} functions that are available in the search path. These functions can be used together with the \code{\link{runDiffExp}} function to perform differential expression analysis. Consult the help pages for the respective functions for more information.
#' @export
#' @author Charlotte Soneson
#' @examples
#' listcreateRmd()
#'
#' @importFrom stats cor hclust as.dist runif rexp median loess predict na.omit rnbinom rpois rnorm sd lm pf pt qf qt
#' @importFrom grDevices heat.colors
#' @importFrom graphics par lines legend title axis
#'
listcreateRmd <- function() {
  s <- unlist(lapply(search(), ls, all.names = TRUE))
  print(unname(s[grep("\\.createRmd$", s)]))
}

#' Generate a .Rmd file containing code to perform differential expression analysis with EBSeq
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the \code{EBSeq} package. The code is written to a .Rmd file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the meaning of the parameters, see the \code{EBSeq} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. Possible values are \code{"median"} and \code{"quantile"}.
#' @references
#' Leng N, Dawson JA, Thomson JA, Ruotti V, Rissman AI, Smits BMG, Haag JD, Gould MN, Stewart RM and Kendziorski C (2013): EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics
#' @author Charlotte Soneson
#' @return The function generates a .Rmd file containing the differential expression code. This file can be executed using e.g. the \code{knitr} package.
#' @export
#' @examples
#' try(
#' if (require(EBSeq)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
#'                                     samples.per.cond = 5, n.diffexp = 100,
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "EBSeq",
#'            Rmdfunction = "EBSeq.createRmd",
#'            output.directory = tmpdir, norm.method = "median")
#' }
#' )
#' @importFrom utils packageVersion
EBSeq.createRmd <- function(data.path, result.path, codefile,
                            norm.method) {
  codefile <- file(codefile, open = 'w')
  writeLines("### EBSeq", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = TRUE, error = TRUE, warning = TRUE}",
               "require(EBSeq)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')"),
             codefile)
  if (norm.method == "median") {
    writeLines("sizes <- MedianNorm(count.matrix(cdata))", codefile)
  } else if (norm.method == "quantile") {
    writeLines("sizes <- QuantileNorm(count.matrix(cdata), .75)", codefile)
  }
  writeLines(c("EBSeq.test <- EBTest(Data = count.matrix(cdata), Conditions = factor(sample.annotations(cdata)$condition), sizeFactors = sizes, maxround = 10)",
               "EBSeq.ppmat <- GetPPMat(EBSeq.test)",
               "EBSeq.posteriors.DE <- EBSeq.ppmat[, 'PPDE']",
               "EBSeq.lfdr <- 1 - EBSeq.ppmat[, 'PPDE']",
               "EBSeq.FDR <- rep(NA, length(EBSeq.lfdr))",
               "for (l in 1:length(EBSeq.FDR)) {",
               "EBSeq.FDR[l] <- mean(EBSeq.lfdr[which(EBSeq.lfdr <= EBSeq.lfdr[l])])",
               "}"), codefile)
  writeLines(c("EBSeq.score <- 1 - EBSeq.FDR",
               "result.table <- data.frame('FDR' = EBSeq.FDR, 'lfdr' = EBSeq.lfdr, 'score' = EBSeq.score, 'posterior.DE' = EBSeq.posteriors.DE)",
               "result.table <- result.table[match(rownames(count.matrix(cdata)), rownames(result.table)), ]",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('EBSeq,', packageVersion('EBSeq'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'EBSeq', 'full.name' = '", paste('EBSeq.', utils::packageVersion('EBSeq'), '.', norm.method, sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with the edgeR exact test
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the exact test functionality from the edgeR package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{edgeR} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}.
#' @param trend.method The method used to estimate the trend in the mean-dispersion relationship. Possible values are \code{"none"}, \code{"movingave"} and \code{"loess"}
#' @param disp.type The type of dispersion estimate used. Possible values are \code{"common"}, \code{"trended"} and \code{"tagwise"}.
#' @references
#' Robinson MD, McCarthy DJ and Smyth GK (2010): edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#' @export
#' @author Charlotte Soneson
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
#'                                     samples.per.cond = 5, n.diffexp = 100,
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "edgeR.exact",
#'            Rmdfunction = "edgeR.exact.createRmd",
#'            output.directory = tmpdir, norm.method = "TMM",
#'            trend.method = "movingave", disp.type = "tagwise")
#' @importFrom utils packageVersion
edgeR.exact.createRmd <- function(data.path, result.path, codefile,
                                  norm.method, trend.method, disp.type) {
  ## Write the code for applying edgeR exact test to an Rmd file
  codefile <- file(codefile, open = "w")
  writeLines("### edgeR", codefile)
  writeLines(paste("Data file: ", data.path, sep = ""), codefile)

  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(edgeR)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               "edgeR.dgelist <- edgeR::DGEList(counts = count.matrix(cdata), group = factor(sample.annotations(cdata)$condition))",
               paste("edgeR.dgelist <- edgeR::calcNormFactors(edgeR.dgelist, method = '", norm.method, "')", sep = ''),
               "edgeR.dgelist <- edgeR::estimateCommonDisp(edgeR.dgelist)"), codefile)
  if (disp.type == 'tagwise') {
    writeLines(c(paste("edgeR.dgelist <- edgeR::estimateTagwiseDisp(edgeR.dgelist, trend = '", trend.method, "')", sep = ''),
                 "dispersion.S1 <- edgeR.dgelist$tagwise.dispersion",
                 "dispersion.S2 <- edgeR.dgelist$tagwise.dispersion"), codefile)
  } else {
    if (disp.type == 'common') {
      writeLines(c("dispersion.S1 <- rep(edgeR.dgelist$common.dispersion, nrow(count.matrix(cdata)))",
                   "dispersion.S2 <- rep(edgeR.dgelist$common.dispersion, nrow(count.matrix(cdata)))"), codefile)
    } else if (disp.type == 'trended') {
      writeLines(c("edgeR.dgelist <- edgeR::estimateTrendedDisp(edgeR.dgelist)",
                   "dispersion.S1 <- edgeR.dgelist$trended.dispersion",
                   "dispersion.S2 <- edgeR.dgelist$trended.dispersion"), codefile)
    }
  }
  writeLines(c("edgeR.exacttest <- edgeR::exactTest(edgeR.dgelist)",
               "edgeR.pvalues <- edgeR.exacttest$table$PValue",
               "edgeR.adjpvalues <- p.adjust(edgeR.pvalues, method = 'BH')",
               "edgeR.logFC <- edgeR.exacttest$table$logFC",
               "edgeR.score <- 1 - edgeR.pvalues",
               "result.table <- data.frame('pvalue' = edgeR.pvalues, 'adjpvalue' = edgeR.adjpvalues, 'logFC' = edgeR.logFC, 'score' = edgeR.score, 'dispersion.S1' = dispersion.S1, 'dispersion.S2' = dispersion.S2)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('edgeR,', packageVersion('edgeR'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'edgeR', 'full.name' = '", paste('edgeR.', utils::packageVersion('edgeR'), '.exact.', norm.method, '.', trend.method, '.', disp.type, sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with the edgeR GLM approach
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the GLM functionality from the edgeR package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{edgeR} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}.
#' @param disp.type The type of dispersion estimate used. Possible values are \code{"common"}, \code{"trended"} and \code{"tagwise"}.
#' @param disp.method The method used to estimate the dispersion. Possible values are \code{"CoxReid"}, \code{"Pearson"} and \code{"deviance"}.
#' @param trended Logical parameter indicating whether or not a trended dispersion estimate should be used.
#' @references
#' Robinson MD, McCarthy DJ and Smyth GK (2010): edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#' @export
#' @author Charlotte Soneson
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
#'                                     samples.per.cond = 5, n.diffexp = 100,
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "edgeR.GLM",
#'            Rmdfunction = "edgeR.GLM.createRmd",
#'            output.directory = tmpdir, norm.method = "TMM",
#'            disp.type = "tagwise", disp.method = "CoxReid",
#'            trended = TRUE)
#' @importFrom utils packageVersion
edgeR.GLM.createRmd <- function(data.path, result.path, codefile,
                                norm.method, disp.type, disp.method, trended) {
  codefile <- file(codefile, open = "w")
  writeLines("### edgeR", codefile)
  writeLines(paste("Data file: ", data.path, sep = ""), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(edgeR)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               "edgeR.dgelist <- edgeR::DGEList(counts = count.matrix(cdata), group = factor(sample.annotations(cdata)$condition))",
               paste("edgeR.dgelist <- edgeR::calcNormFactors(edgeR.dgelist, method = '", norm.method, "')", sep = ''),
               paste("edgeR.dgelist <- edgeR::estimateGLMCommonDisp(edgeR.dgelist, design = model.matrix(~factor(sample.annotations(cdata)$condition)), method = '", disp.method, "')", sep = '')), codefile)
  if ((disp.type == 'tagwise' & trended == TRUE) | disp.type == 'trended') {
    writeLines("edgeR.dgelist <- edgeR::estimateGLMTrendedDisp(edgeR.dgelist, design = model.matrix(~factor(sample.annotations(cdata)$condition)), method = 'auto')", codefile)
  }
  if (disp.type == 'tagwise') {
    writeLines(paste("edgeR.dgelist <- edgeR::estimateGLMTagwiseDisp(edgeR.dgelist, design = model.matrix(~factor(sample.annotations(cdata)$condition)), trend = ", trended, ")", sep = ''), codefile)
  }
  if (disp.type == 'common') {
    writeLines(c("dispersion.S1 <- rep(edgeR.dgelist$common.dispersion, nrow(count.matrix(cdata)))",
                 "dispersion.S2 <- rep(edgeR.dgelist$common.dispersion, nrow(count.matrix(cdata)))"), codefile)
  } else {
    if (disp.type == 'trended') {
      writeLines(c("dispersion.S1 <- edgeR.dgelist$trended.dispersion",
                   "dispersion.S2 <- edgeR.dgelist$trended.dispersion"), codefile)
    } else {
      writeLines(c("dispersion.S1 <- edgeR.dgelist$tagwise.dispersion",
                   "dispersion.S2 <- edgeR.dgelist$tagwise.dispersion"), codefile)
    }
  }
  writeLines(c("edgeR.glmfit <- edgeR::glmFit(edgeR.dgelist, design = model.matrix(~factor(sample.annotations(cdata)$condition)))",
               "edgeR.glmlrt <- edgeR::glmLRT(edgeR.glmfit, coef = 2)",
               "edgeR.pvalues <- edgeR.glmlrt$table$PValue",
               "edgeR.adjpvalues <- p.adjust(edgeR.pvalues, method = 'BH')",
               "edgeR.logFC <- edgeR.glmlrt$table$logFC",
               "edgeR.score <- 1 - edgeR.pvalues",
               "result.table <- data.frame('pvalue' = edgeR.pvalues, 'adjpvalue' = edgeR.adjpvalues, 'logFC' = edgeR.logFC, 'score' = edgeR.score, 'dispersion.S1' = dispersion.S1, 'dispersion.S2' = dispersion.S2)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('edgeR,', packageVersion('edgeR'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'edgeR', 'full.name' = '", paste('edgeR.', utils::packageVersion('edgeR'), '.GLM.', norm.method, '.', as.vector(ifelse(trended == TRUE, 'trend', 'notrend')), '.', disp.method, '.', disp.type, sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}


#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with DESeq2
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the DESeq2 package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{DESeq2} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param fit.type The fitting method used to get the dispersion-mean relationship. Possible values are \code{"parametric"}, \code{"local"} and \code{"mean"}.
#' @param test The test to use. Possible values are \code{"Wald"} and \code{"LRT"}.
#' @param beta.prior Whether or not to put a zero-mean normal prior on the non-intercept coefficients. Default is \code{TRUE}.
#' @param independent.filtering Whether or not to perform independent filtering of the data. With independent filtering=TRUE, the adjusted p-values for genes not passing the filter threshold are set to NA.
#' @param cooks.cutoff The cutoff value for the Cook's distance to consider a value to be an outlier. Set to Inf or FALSE to disable outlier detection. For genes with detected outliers, the p-value and adjusted p-value will be set to NA.
#' @param impute.outliers Whether or not the outliers should be replaced by a trimmed mean and the analysis rerun.
#' @param nas.as.ones Whether or not adjusted p values that are returned as \code{NA} by \code{DESeq2} should be set to \code{1}. This option is useful for comparisons with other methods. For more details, see section "I want to benchmark DESeq2 comparing to other DE tools" from the \code{DESeq2} vignette (available by running \code{vignette("DESeq2", package = "DESeq2")}). Default to \code{FALSE}.
#'
#' @export
#' @author Charlotte Soneson
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references
#' Anders S and Huber W (2010): Differential expression analysis for sequence count data. Genome Biology 11:R106
#' @examples
#' try(
#' if (require(DESeq2)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
#'                                     samples.per.cond = 5, n.diffexp = 100,
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "DESeq2",
#'            Rmdfunction = "DESeq2.createRmd",
#'            output.directory = tmpdir, fit.type = "parametric",
#'            test = "Wald")
#' })
#' @importFrom utils packageVersion
DESeq2.createRmd <- function(data.path, result.path, codefile,
                             fit.type, test, beta.prior = TRUE,
                             independent.filtering = TRUE, cooks.cutoff = TRUE,
                             impute.outliers = TRUE,
                             nas.as.ones = FALSE) {
  codefile <- file(codefile, open = 'w')
  writeLines("### DESeq2", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = TRUE, error = TRUE, warning = TRUE}",
               "require(DESeq2)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               "DESeq2.ds <- DESeq2::DESeqDataSetFromMatrix(countData = count.matrix(cdata), colData = data.frame(condition = factor(sample.annotations(cdata)$condition)), design = ~ condition)",
               paste("DESeq2.ds <- DESeq2::DESeq(DESeq2.ds, fitType = '", fit.type, "', test = '", test, "', betaPrior = ", beta.prior, ")", sep = "")), codefile)
  if (impute.outliers == TRUE) {
    writeLines(c("DESeq2.ds.clean <- DESeq2::replaceOutliersWithTrimmedMean(DESeq2.ds)",
                 paste("DESeq2.ds.clean <- DESeq2::DESeq(DESeq2.ds.clean, fitType = '", fit.type, "', test = '", test, "', betaPrior = ", beta.prior, ")", sep = ""),
                 "DESeq2.ds <- DESeq2.ds.clean"), codefile)
  }
  writeLines(paste("DESeq2.results <- DESeq2::results(DESeq2.ds, independentFiltering = ", independent.filtering, ", cooksCutoff = ", cooks.cutoff, ")", sep = ""),
             codefile)
  if (nas.as.ones) {
    # see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#i-want-to-benchmark-deseq2-comparing-to-other-de-tools.
    message("As `nas.as.ones=TRUE`, all NAs in adjusted p values are replaced by 1 to allow for benchmarking with other methods. For more details, see section 'I want to benchmark DESeq2 comparing to other DE tools' from the `DESeq2` vignette (available by running `vignette('DESeq2', package = 'DESeq2')`)")
    writeLines("DESeq2.results$padj <- ifelse(is.na(DESeq2.results$padj), 1, DESeq2.results$padj)", codefile)
  } else {
    message("As `nas.as.ones=FALSE`, there might be some NAs in the adjusted p values computed by DESeq2. This might bias the comparison of the results with other methods. For more details, see section 'I want to benchmark DESeq2 comparing to other DE tools' from the `DESeq2` vignette (available by running `vignette('DESeq2', package = 'DESeq2')`)")
  }
  writeLines(c("DESeq2.pvalues <- DESeq2.results$pvalue",
               "DESeq2.adjpvalues <- DESeq2.results$padj",
               "DESeq2.logFC <- DESeq2.results$log2FoldChange",
               "DESeq2.score <- 1 - DESeq2.pvalues",
               "result.table <- data.frame('pvalue' = DESeq2.pvalues, 'adjpvalue' = DESeq2.adjpvalues, 'logFC' = DESeq2.logFC, 'score' = DESeq2.score)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('DESeq2,', packageVersion('DESeq2'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'DESeq2', 'full.name' = '", paste('DESeq2.', utils::packageVersion('DESeq2'), '.', fit.type, '.', test, '.', ifelse(beta.prior == TRUE, 'bp', 'nobp'), '.', ifelse(independent.filtering == TRUE, 'indf', 'noindf'), paste(".cook_", cooks.cutoff, sep = ""), ifelse(impute.outliers, ".imp", ".noimp"), sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with voom+limma
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) by applying the voom transformation (from the limma package) followed by differential expression analysis with limma. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{limma} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. The normalization factors are calculated using the \code{calcNormFactors} of the \code{edgeR} package. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}
#' @export
#' @author Charlotte Soneson
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references
#' Smyth GK (2005): Limma: linear models for microarray data. In: 'Bioinformatics and Computational Biology Solutions using R and Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds), Springer, New York, pages 397-420
#'
#'  Law CW, Chen Y, Shi W and Smyth GK (2014): voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology 15, R29
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
#'                                     samples.per.cond = 5, n.diffexp = 100,
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "voom.limma",
#'            Rmdfunction = "voom.limma.createRmd",
#'            output.directory = tmpdir, norm.method = "TMM")
#' @importFrom utils packageVersion
voom.limma.createRmd <- function(data.path, result.path, codefile, norm.method) {
  codefile <- file(codefile, open = 'w')
  writeLines("### voom + limma", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(limma)",
               "require(edgeR)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               paste("nf <- edgeR::calcNormFactors(count.matrix(cdata), method = '", norm.method, "')", sep = ''),
               "voom.data <- limma::voom(count.matrix(cdata), design = model.matrix(~factor(sample.annotations(cdata)$condition)), lib.size = colSums(count.matrix(cdata)) * nf)",
               "voom.data$genes <- rownames(count.matrix(cdata))",
               "voom.fitlimma <- limma::lmFit(voom.data, design = model.matrix(~factor(sample.annotations(cdata)$condition)))",
               "voom.fitbayes <- limma::eBayes(voom.fitlimma)",
               "voom.pvalues <- voom.fitbayes$p.value[, 2]",
               "voom.adjpvalues <- p.adjust(voom.pvalues, method = 'BH')",
               "voom.logFC <- voom.fitbayes$coefficients[, 2]",
               "voom.score <- 1 - voom.pvalues",
               "result.table <- data.frame('pvalue' = voom.pvalues, 'adjpvalue' = voom.adjpvalues, 'logFC' = voom.logFC, 'score' = voom.score)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('limma,', packageVersion('limma'), ';', 'edgeR,', packageVersion('edgeR'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'voom', 'full.name' = '",
                     paste('voom.', utils::packageVersion('limma'), '.limma.', norm.method, sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with baySeq
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the \code{baySeq} package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{baySeq} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. Possible values are \code{"quantile"}, \code{"total"} and \code{"edgeR"}.
#' @param equaldisp Logical parameter indicating whether or not equal dispersion should be assumed across all conditions.
#' @param sample.size The size of the sample used to estimate the priors (default 5000).
#' @param estimation The approach used to estimate the priors. Possible values are \code{"QL"} (default), \code{"ML"} and \code{"edgeR"}.
#' @param pET The method used to re-estimate the priors. Possible values are \code{"BIC"} (default), \code{"none"} and \code{"iteratively"}.
#' @noRd
#' @author Charlotte Soneson
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references
#' Hardcastle TJ (2012): baySeq: Empirical Bayesian analysis of patterns of differential expression in count data. R package
#'
#' Hardcastle TJ and Kelly KA (2010): baySeq: Empirical Bayesian methods for identifying differential expression in sequence count data. BMC Bioinformatics 11:422
#' @examples
#' try(
#' if (require(baySeq)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
#'                                     samples.per.cond = 5, n.diffexp = 100,
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' ## Note! In the interest of speed, we set sample.size=10 in this example.
#' ## In a real analysis, much larger sample sizes are recommended (the default is 5000).
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "baySeq",
#'            Rmdfunction = "baySeq.createRmd",
#'            output.directory = tmpdir, norm.method = "edgeR",
#'            equaldisp = TRUE, sample.size = 10)
#' })
#' @importFrom utils packageVersion
baySeq.createRmd <- function(data.path, result.path, codefile,
                             norm.method, equaldisp, sample.size = 5000,
                             estimation = "QL", pET = "BIC") {
  codefile <- file(codefile, open = 'w')
  writeLines("### baySeq", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = TRUE, error = TRUE, warning = TRUE}",
               "require(baySeq)",
               "require(edgeR)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               "baySeq.cd <- new('countData', data = count.matrix(cdata), replicates = sample.annotations(cdata)$condition, groups = list(NDE = rep(1, length(sample.annotations(cdata)$condition)), DE = sample.annotations(cdata)$condition))",
               paste("libsizes(baySeq.cd) <- baySeq::getLibsizes(baySeq.cd, estimationType = '", norm.method, "')", sep = '')), codefile)
  writeLines(c(paste("baySeq.cd <- baySeq::getPriors.NB(baySeq.cd, samplesize =", sample.size, ", equalDispersions = ", equaldisp, ", estimation = '", estimation, "', cl = NULL)", sep = ''),
               paste("baySeq.cd <- baySeq::getLikelihoods(baySeq.cd, prs = c(0.5, 0.5), pET = '", pET, "', cl = NULL)", sep = '')), codefile)
  writeLines(c("baySeq.posteriors.DE <- exp(baySeq.cd@posteriors)[, 'DE']",
               "baySeq.FDR <- baySeq::topCounts(baySeq.cd, group = 'DE', FDR = 1)$FDR.DE[match(rownames(count.matrix(cdata)), rownames(baySeq::topCounts(baySeq.cd, group = 'DE', FDR = 1)))]",
               "baySeq.score <- 1 - baySeq.FDR",
               "result.table <- data.frame('FDR' = baySeq.FDR, 'score' = baySeq.score, 'posterior.DE' = baySeq.posteriors.DE)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('baySeq,', packageVersion('baySeq'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'baySeq', 'full.name' = '", paste('baySeq.', utils::packageVersion('baySeq'), '.', norm.method, '.', ifelse(equaldisp == TRUE, 'equaldisp', 'nonequaldisp'), '.samplesize', sample.size, '.', estimation, '.', pET, sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}


#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with voom+t-test
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) by applying the voom transformation (from the \code{limma} package) followed by differential expression analysis with a t-test. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{limma} and \code{edgeR} packages and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. The normalization factors are calculated using the \code{calcNormFactors} function from the \code{edgeR} package. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}.
#' @export
#' @author Charlotte Soneson
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references
#' Smyth GK (2005): Limma: linear models for microarray data. In: 'Bioinformatics and Computational Biology Solutions using R and Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds), Springer, New York, pages 397-420
#'
#' Law CW, Chen Y, Shi W and Smyth GK (2014): voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology 15, R29
#' @examples
#' try(
#' if (require(genefilter)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
#'                                     samples.per.cond = 5, n.diffexp = 100,
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "voom.ttest",
#'            Rmdfunction = "voom.ttest.createRmd",
#'            output.directory = tmpdir, norm.method = "TMM")
#' })
#' @importFrom utils packageVersion
voom.ttest.createRmd <- function(data.path, result.path, codefile, norm.method) {
  codefile <- file(codefile, open = 'w')
  writeLines("### voom (t-test)", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(limma)",
               "require(edgeR)",
               "require(genefilter)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               paste("nf <- edgeR::calcNormFactors(count.matrix(cdata), method = '", norm.method, "')", sep = ''),
               "voom.data <- limma::voom(count.matrix(cdata), design = model.matrix(~factor(sample.annotations(cdata)$condition)), lib.size = colSums(count.matrix(cdata)) * nf)",
               "voom.data$genes <- rownames(count.matrix(cdata))",
               "voom.ttest <- genefilter::rowttests(voom.data$E, factor(sample.annotations(cdata)$condition))",
               "voom.pvalues <- voom.ttest$p.value",
               "voom.adjpvalues <- p.adjust(voom.pvalues, method = 'BH')",
               "voom.logFC <- voom.ttest$dm",
               "voom.score <- 1 - voom.pvalues",
               "result.table <- data.frame('pvalue' = voom.pvalues, 'adjpvalue' = voom.adjpvalues, 'logFC' = voom.logFC, 'score' = voom.score)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('limma,', packageVersion('limma'), ';', 'edgeR,', packageVersion('edgeR'), ';', 'genefilter,', packageVersion('genefilter'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'voom', 'full.name' = '",
                     paste('voom.', utils::packageVersion('limma'), '.ttest.', norm.method, sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with limma after log-transforming the counts per million (cpm)
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using limma, after preprocessing the counts by computing the counts per million (cpm) and applying a logarithmic transformation. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{edgeR} and \code{limma} packages and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. The normalization factors are calculated using the \code{calcNormFactors} function from the \code{edgeR} package. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}
#' @export
#' @author Charlotte Soneson
#' @references
#' Smyth GK (2005): Limma: linear models for microarray data. In: 'Bioinformatics and Computational Biology Solutions using R and Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds), Springer, New York, pages 397-420
#'
#' Robinson MD, McCarthy DJ and Smyth GK (2010): edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#'
#' Robinson MD and Oshlack A  (2010): A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology 11:R25
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
#'                                     samples.per.cond = 5, n.diffexp = 100,
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "logcpm.limma",
#'            Rmdfunction = "logcpm.limma.createRmd",
#'            output.directory = tmpdir, norm.method = "TMM")
#' @importFrom utils packageVersion
logcpm.limma.createRmd <- function(data.path, result.path, codefile, norm.method) {
  codefile <- file(codefile, open = 'w')
  writeLines("### log(cpm) + limma", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(limma)",
               "require(edgeR)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               paste("nf <- edgeR::calcNormFactors(count.matrix(cdata), method = '", norm.method, "')", sep = ''),
               "logcpm.data <- log2(sweep((count.matrix(cdata) + 0.5), 2, colSums(count.matrix(cdata)) * nf, '/')*1e6)",
               "rownames(logcpm.data) <- rownames(count.matrix(cdata))",
               "logcpm.fitlimma <- limma::lmFit(logcpm.data, design = model.matrix(~factor(sample.annotations(cdata)$condition)))",
               "logcpm.fitbayes <- limma::eBayes(logcpm.fitlimma)",
               "logcpm.pvalues <- logcpm.fitbayes$p.value[, 2]",
               "logcpm.adjpvalues <- p.adjust(logcpm.pvalues, method = 'BH')",
               "logcpm.logFC <- logcpm.fitbayes$coefficients[, 2]",
               "logcpm.score <- 1 - logcpm.pvalues",
               "result.table <- data.frame('pvalue' = logcpm.pvalues, 'adjpvalue' = logcpm.adjpvalues, 'logFC' = logcpm.logFC, 'score' = logcpm.score)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('limma,', packageVersion('limma'), ';', 'edgeR,', packageVersion('edgeR'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'logcpm.limma', 'full.name' = '",
                     paste('logcpm.limma.', utils::packageVersion('limma'), '.', norm.method, sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with limma after square root-transforming the counts per million (cpm)
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using limma, after preprocessing the counts by computing the counts per million (cpm) and applying a square-root transformation. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{edgeR} and \code{limma} packages and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. The normalization factors are calculated using the \code{calcNormFactors} function from the \code{edgeR} package. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}.
#' @export
#' @author Charlotte Soneson
#' @references
#' Smyth GK (2005): Limma: linear models for microarray data. In: 'Bioinformatics and Computational Biology Solutions using R and Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds), Springer, New York, pages 397-420
#'
#' Robinson MD, McCarthy DJ and Smyth GK (2010): edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#'
#' Robinson MD and Oshlack A (2010): A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology 11:R25
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
#'                                     samples.per.cond = 5, n.diffexp = 100,
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "sqrtcpm.limma",
#'            Rmdfunction = "sqrtcpm.limma.createRmd",
#'            output.directory = tmpdir, norm.method = "TMM")
#' @importFrom utils packageVersion
sqrtcpm.limma.createRmd <- function(data.path, result.path, codefile, norm.method) {
  codefile <- file(codefile, open = 'w')
  writeLines("### sqrt(cpm) + limma", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(limma)",
               "require(edgeR)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               paste("nf <- edgeR::calcNormFactors(count.matrix(cdata), method = '", norm.method, "')", sep = ''),
               "sqrtcpm.data <- sqrt(sweep((count.matrix(cdata) + 0.5), 2, colSums(count.matrix(cdata)) * nf, '/')*1e6)",
               "rownames(sqrtcpm.data) <- rownames(count.matrix(cdata))",
               "sqrtcpm.fitlimma <- limma::lmFit(sqrtcpm.data, design = model.matrix(~factor(sample.annotations(cdata)$condition)))",
               "sqrtcpm.fitbayes <- limma::eBayes(sqrtcpm.fitlimma)",
               "sqrtcpm.pvalues <- sqrtcpm.fitbayes$p.value[, 2]",
               "sqrtcpm.adjpvalues <- p.adjust(sqrtcpm.pvalues, method = 'BH')",
               "sqrtcpm.logFC <- sqrtcpm.fitbayes$coefficients[, 2]",
               "sqrtcpm.score <- 1 - sqrtcpm.pvalues",
               "result.table <- data.frame('pvalue' = sqrtcpm.pvalues, 'adjpvalue' = sqrtcpm.adjpvalues, 'logFC' = sqrtcpm.logFC, 'score' = sqrtcpm.score)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('limma,', packageVersion('limma'), ';', 'edgeR,', packageVersion('edgeR'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'sqrtcpm.limma', 'full.name' = '",
                     paste('sqrtcpm.limma.', utils::packageVersion('limma'), '.', norm.method, sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with a t-test
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using a t-test, applied to the normalized counts. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{edgeR} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. The normalization factors are calculated using the \code{calcNormFactors} function from the \code{edgeR} package. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}
#' @export
#' @author Charlotte Soneson
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references
#' Robinson MD, McCarthy DJ and Smyth GK (2010): edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#'
#' Robinson MD and Oshlack A (2010): A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology 11:R25
#' @examples
#' try(
#' if (require(genefilter)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
#'                                     samples.per.cond = 5, n.diffexp = 100,
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "ttest",
#'            Rmdfunction = "ttest.createRmd",
#'            output.directory = tmpdir, norm.method = "TMM")
#' })
#' @importFrom utils packageVersion
ttest.createRmd <- function(data.path, result.path, codefile, norm.method) {
  codefile <- file(codefile, open = 'w')
  writeLines("### t-test", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(genefilter)",
               "require(edgeR)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               paste("nf <- edgeR::calcNormFactors(count.matrix(cdata), method = '", norm.method, "') * colSums(count.matrix(cdata))/exp(mean(log(colSums(count.matrix(cdata)))))", sep = ''),
               "count.matrix <- round(sweep(count.matrix(cdata), 2, nf, '/'))",
               "ttest.result <- genefilter::rowttests(count.matrix, factor(sample.annotations(cdata)$condition))",
               "ttest.pvalues <- ttest.result$p.value",
               "ttest.adjpvalues <- p.adjust(ttest.pvalues, method = 'BH')",
               "ttest.logFC <- ttest.result$dm",
               "ttest.score <- 1 - ttest.pvalues",
               "result.table <- data.frame('pvalue' = ttest.pvalues, 'adjpvalue' = ttest.adjpvalues, 'logFC' = ttest.logFC, 'score' = ttest.score)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('edgeR,', packageVersion('edgeR'), ';', 'genefilter,', packageVersion('genefilter'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'ttest', 'full.name' = 'ttest.", utils::packageVersion('genefilter'), '.', norm.method, "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with NOISeq
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using \code{NOISeq}. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{NOISeq} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. The normalization factors are calculated using the \code{calcNormFactors} function from the \code{edgeR} package. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}.
#' @export
#' @author Charlotte Soneson
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references
#' Robinson MD, McCarthy DJ and Smyth GK (2010): edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#'
#' Robinson MD and Oshlack A (2010): A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology 11:R25
#'
#' Tarazona S, Furio-Tari P, Ferrer A and Conesa A (2012): NOISeq: Exploratory analysis and differential expression for RNA-seq data. R package
#'
#' Tarazona S, Garcia-Alcalde F, Dopazo J, Ferrer A and Conesa A (2011): Differential expression in RNA-seq: a matter of depth. Genome Res 21(12), 2213-2223
#' @examples
#' try(
#' if (require(NOISeq)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
#'                                     samples.per.cond = 5, n.diffexp = 100,
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "NOISeq",
#'            Rmdfunction = "NOISeq.prenorm.createRmd",
#'            output.directory = tmpdir, norm.method = "TMM")
#' })
#' @importFrom utils packageVersion
NOISeq.prenorm.createRmd <- function(data.path, result.path, codefile, norm.method) {
  codefile <- file(codefile, open = 'w')
  cdata <- readRDS(data.path)
  writeLines("### NOISeq", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(edgeR)",
               "require(NOISeq)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(cdata)) {
    cdata <- convertListTocompData(cdata)
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               paste("nf <- edgeR::calcNormFactors(count.matrix(cdata), method = '", norm.method, "') * colSums(count.matrix(cdata))/exp(mean(log(colSums(count.matrix(cdata)))))", sep = ''),
               "count.matrix <- round(sweep(count.matrix(cdata), 2, nf, '/'))",
               "inp.data <- NOISeq::readData(count.matrix, factors = sample.annotations(cdata))",
               paste("NOISeq.test <- NOISeq::noiseqbio(inp.data, k = 0.5, norm = 'n', nclust = 15, factor = 'condition', conditions = c('", unique(as.character(sample.annotations(cdata)$condition))[1], "', '", unique(as.character(sample.annotations(cdata)$condition))[2], "'), filter = 0)@results", sep = ""),
               ##               paste("NOISeq.test = noiseqbio(inp.data, k = 0.5, norm = 'n', nclust = 15, factor = 'condition', conditions = c(unique(as.character(sample.annotations(cdata)$condition))[1], unique(as.character(sample.annotations(cdata)$condition))[2]), filter = 0)@results", sep = ""),
               "NOISeq.prob <- NOISeq.test[[1]]$prob",
               "NOISeq.statistic <- NOISeq.test[[1]]$theta",
               "NOISeq.score <- NOISeq.prob",
               "result.table <- data.frame('probabilities' = NOISeq.prob, 'statistic' = NOISeq.statistic, 'score' = NOISeq.score, 'FDR' = 1 - NOISeq.prob)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('NOISeq,', packageVersion('NOISeq'), ';', 'edgeR,', packageVersion('edgeR'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'NOISeq', 'full.name' = '",
                     paste('NOISeq.', utils::packageVersion('NOISeq'), '.', norm.method, sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with NBPSeq
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using \code{NBPSeq}. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{NBPSeq} and \code{edgeR} packages and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. The normalization factors are calculated using the \code{calcNormFactors} function from the \code{edgeR} package. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}.
#' @param disp.method The method to use to estimate the dispersion values. Possible values are \code{"NBP"} and \code{"NB2"}.
#' @export
#' @author Charlotte Soneson
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references
#' Robinson MD, McCarthy DJ and Smyth GK (2010): edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#'
#' Robinson MD and Oshlack A (2010): A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology 11:R25
#'
#' Di Y, Schafer DW, Cumbie JS, and Chang JH (2011): The NBP Negative Binomial Model for Assessing Differential Gene Expression from RNA-Seq. Statistical Applications in Genetics and Molecular Biology 10(1), 1-28
#' @examples
#' try(
#' if (require(NBPSeq)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
#'                                     samples.per.cond = 5, n.diffexp = 100,
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "NBPSeq",
#'            Rmdfunction = "NBPSeq.createRmd",
#'            output.directory = tmpdir, norm.method = "TMM", disp.method = "NBP")
#' })
#' @importFrom utils packageVersion
NBPSeq.createRmd <- function(data.path, result.path, codefile, norm.method, disp.method) {
  codefile <- file(codefile, open = 'w')
  writeLines("### NBPSeq", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(edgeR)",
               "require(NBPSeq)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               paste("nf <- edgeR::calcNormFactors(count.matrix(cdata), method = '", norm.method, "')", sep = ''),
               paste("NBPSeq.test <- NBPSeq::nbp.test(count.matrix(cdata), grp.ids = sample.annotations(cdata)$condition, grp1 = levels(factor(sample.annotations(cdata)$condition))[1], grp2 = levels(factor(sample.annotations(cdata)$condition))[2], norm.factors = nf, model.disp = '", disp.method, "', print.level = 5)", sep = ''),
               "NBPSeq.pvalues <- NBPSeq.test$p.values",
               "NBPSeq.adjpvalues <- NBPSeq.test$q.values",
               "NBPSeq.pvalues[which(is.na(NBPSeq.pvalues))] <- 1",
               "NBPSeq.adjpvalues[which(is.na(NBPSeq.adjpvalues))] <- 1",
               "NBPSeq.logFC <- NBPSeq.test$log.fc",
               "NBPSeq.score <- 1 - NBPSeq.pvalues",
               "result.table <- data.frame('pvalue' = NBPSeq.pvalues, 'adjpvalue' = NBPSeq.adjpvalues, 'logFC' = NBPSeq.logFC, 'score' = NBPSeq.score)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('NBPSeq,', packageVersion('NBPSeq'), ';', 'edgeR,', packageVersion('edgeR'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'NBPSeq', 'full.name' = '",
                     paste('NBPSeq.', utils::packageVersion('NBPSeq'), '.', norm.method, '.', disp.method, sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with DSS
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the DSS package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{DSS} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. Possible values are \code{"quantile"}, \code{"total"} and \code{"median"}.
#' @param disp.trend A logical parameter indicating whether or not to include a trend in the dispersion estimation.
#' @export
#' @author Charlotte Soneson
#' @references
#' Wu H, Wang C and Wu Z (2013): A new shrinkage estimator for dispersion improves differential expression detection in RNA-seq data. Biostatistics 14(2), 232-243
#' @examples
#' try(
#' if (require(DSS)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
#'                                     samples.per.cond = 5, n.diffexp = 100,
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "DSS",
#'            Rmdfunction = "DSS.createRmd",
#'            output.directory = tmpdir, norm.method = "quantile",
#'            disp.trend = TRUE)
#' })
#' @importFrom utils packageVersion
DSS.createRmd <- function(data.path, result.path, codefile, norm.method, disp.trend) {
  codefile <- file(codefile, open = 'w')
  writeLines("### DSS", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(DSS)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               "designs <- sample.annotations(cdata)$condition",
               "names(designs) <- rownames(sample.annotations(cdata))",
               "DSS.seqData <- DSS::newSeqCountSet(count.matrix(cdata), designs)",
               paste("DSS.seqData <- DSS::estNormFactors(DSS.seqData, method = '", norm.method, "')", sep = ''),
               paste("DSS.seqData <- DSS::estDispersion(DSS.seqData, trend = ", disp.trend, ")", sep = ''),
               "DSS.test <- DSS::waldTest(DSS.seqData, levels(factor(sample.annotations(cdata)$condition))[1], levels(factor(sample.annotations(cdata)$condition))[2])",
               "DSS.test <- DSS.test[match(rownames(count.matrix(cdata)), rownames(DSS.test)), ]",
               "DSS.pvalues <- DSS.test[, 'pval']",
               "DSS.adjpvalues <- DSS.test[, 'fdr']",
               "DSS.logFC <- DSS.test[, 'lfc']",
               "DSS.statistic <- DSS.test[, 'stats']",
               "DSS.score <- 1 - DSS.pvalues",
               "result.table <- data.frame('pvalue' = DSS.pvalues, 'adjpvalue' = DSS.adjpvalues, 'logFC' = DSS.logFC, 'statistic' = DSS.statistic, 'score' = DSS.score)",
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('DSS,', packageVersion('DSS'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'DSS', 'full.name' = '",
                     paste('DSS.', utils::packageVersion('DSS'), '.', norm.method, '.', ifelse(disp.trend, 'trend', 'notrend'), sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with TCC
#'
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the TCC package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#'
#' For more information about the methods and the interpretation of the parameters, see the \code{TCC} package and the corresponding publications.
#'
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. Possible values are \code{"tmm"}, and \code{"deseq"}.
#' @param test.method The method used in TCC to find differentially expressed genes. Possible values are \code{"edger"}, \code{"deseq"} and \code{"bayseq"}.
#' @param iteration The number of iterations used to find the normalization factors. Default value is 3.
#' @param normFDR The FDR cutoff for calling differentially expressed genes in the computation of the normalization factors. Default value is 0.1.
#' @param floorPDEG The minimum value to be eliminated as potential differentially expressed genes before performing step 3 in the TCC algorithm. Default value is 0.05.
#' @export
#' @author Charlotte Soneson
#' @references
#' Kadota K, Nishiyama T, and Shimizu K. A normalization strategy for comparing tag count data. Algorithms Mol Biol. 7:5, 2012.
#'
#' Sun J, Nishiyama T, Shimizu K, and Kadota K. TCC: an R package for comparing tag count data with robust normalization strategies. BMC Bioinformatics 14:219, 2013.
#' @examples
#' try(
#' if (require(TCC)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
#'                                     samples.per.cond = 5, n.diffexp = 100,
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "TCC",
#'            Rmdfunction = "TCC.createRmd",
#'            output.directory = tmpdir, norm.method = "tmm",
#'            test.method = "edger")
#' })
#' @importFrom utils packageVersion
TCC.createRmd <- function(data.path, result.path, codefile, norm.method,
                          test.method, iteration = 3,
                          normFDR = 0.1, floorPDEG = 0.05) {
  codefile <- file(codefile, open = 'w')
  writeLines("### TCC", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(TCC)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')",
               "TCC.data <- new('TCC', count.matrix(cdata), sample.annotations(cdata)$condition)",
               paste("TCC.data = TCC::calcNormFactors(TCC.data, norm.method = '", norm.method, "', test.method = '", test.method, "', iteration = ", iteration, ", FDR = ", normFDR, ", floorPDEG = ", floorPDEG, ")", sep = ""),
               paste("TCC.results <- TCC::estimateDE(TCC.data, test.method = '", test.method, "', FDR = 1)", sep = "")), codefile)
  if (test.method %in% c("edger", "deseq")) {
    writeLines(c("TCC.pvalues <- TCC.results$stat$p.value",
                 "TCC.adjpvalues <- TCC.results$stat$q.value",
                 "TCC.score <- 1 - TCC.pvalues",
                 "result.table <- data.frame('pvalue' = TCC.pvalues, 'adjpvalue' = TCC.adjpvalues, 'score' = TCC.score)"), codefile)
  } else if (test.method == "bayseq") {
    writeLines(c("TCC.posteriors.DE <- 1 - TCC.results$stat$p.value",
                 "TCC.FDR <- TCC.results$stat$q.value",
                 "TCC.score <- 1 - TCC.FDR",
                 "result.table <- data.frame('posterior.DE' = TCC.posteriors.DE, 'FDR' = TCC.FDR, 'score' = TCC.score)"), codefile)
  }
  writeLines(c("rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('TCC,', packageVersion('TCC'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'TCC', 'full.name' = '",
                     paste('TCC.', utils::packageVersion('TCC'), '.', norm.method, '.', test.method, '.iter', iteration, '.normFDR', normFDR, '.floorPDEG', floorPDEG, sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}
