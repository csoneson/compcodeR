#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with DESeq2 with custom model matrix
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
#' @param extraDesignFactors A vector containing the extra factors to be passed to the design matrix of \code{DESeq2}. All the factors need to be a \code{sample.annotations} from the \code{\link{compData}} object. It should not contain the "condition" factor column, that will be added automatically.
#' @param divByLengths If TRUE, the counts are divided by the sequence lengths. If FALSE, the normalizing method explained in the details section is used. Default to FALSE.
#' 
#' @details 
#' The length matrix are used as a normalization factor and applied to the \code{DESeq2}
#' model in the way explained in \code{\link[DESeq2]{normalizationFactors}}
#' (see details and examples of this function).
#' The provided matrix will be multiplied by the default normalization factor 
#' obtained through the \code{\link[DESeq2]{estimateSizeFactors}} function.
#' 
#' @export 
#' @author Charlotte Soneson, Paul Bastide, Mélina Gallopin
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references 
#' Anders S and Huber W (2010): Differential expression analysis for sequence count data. Genome Biology 11:R106
#' @examples
#' try(
#' if (require(DESeq2)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
#'                                     samples.per.cond = 5, n.diffexp = 100, 
#'                                     id.species = 1:10,
#'                                     lengths.relmeans = rpois(1000, 1000),
#'                                     lengths.dispersions = rgamma(1000, 1, 1),
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' ## Add annotations
#' sample.annotations(mydata.obj)$test_factor <- factor(rep(1:2, each = 5))
#' sample.annotations(mydata.obj)$test_reg <- 1:10
#' saveRDS(mydata.obj, file.path(tmpdir, "mydata.rds"))
#' ## Diff Exp
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "DESeq2", 
#'            Rmdfunction = "DESeq2.length.createRmd", 
#'            output.directory = tmpdir, fit.type = "parametric",
#'            test = "Wald",
#'            extraDesignFactors = c("test_factor", "test_reg"))
#' })
DESeq2.length.createRmd <- function(data.path, result.path, codefile, 
                                    fit.type, test, beta.prior = TRUE, 
                                    independent.filtering = TRUE, cooks.cutoff = TRUE, 
                                    impute.outliers = TRUE,
                                    extraDesignFactors = NULL,
                                    divByLengths = FALSE) {
  codefile <- file(codefile, open = 'w')
  writeLines("### DESeq2.length", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = TRUE, error = TRUE, warning = TRUE}", 
               "require(DESeq2)", 
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')"),
             codefile)
  if (divByLengths) {
    writeLines("count_matrix <- round(count.matrix(cdata) / length.matrix(cdata))", codefile)
  } else {
    writeLines("count_matrix <- count.matrix(cdata)", codefile)
  }
  if (is.null(extraDesignFactors)) {
    writeLines(c(
      paste("DESeq2.length.ds <- DESeq2::DESeqDataSetFromMatrix(countData = count_matrix,",
            "colData = data.frame(condition = factor(sample.annotations(cdata)$condition)),",
            "design = ~ condition)")),
      codefile)
  } else {
    writeLines(c(
      paste0("DESeq2.length.ds <- DESeq2::DESeqDataSetFromMatrix(countData = count_matrix,",
             " colData = cbind(sample.annotations(cdata)[, c('", paste(extraDesignFactors, collapse = "', '"), "'), drop = FALSE], data.frame(condition = factor(sample.annotations(cdata)$condition))),",
             " design = as.formula(paste(' ~ ', paste(c('", paste(extraDesignFactors, collapse = "', '"), "'), collapse= '+'), '+ condition')))")),
      codefile)
  }
  if (!divByLengths) {
    writeLines(c(
      "## Size Factors",
      "DESeq2.length.ds <-  estimateSizeFactors(DESeq2.length.ds)",
      "size_fac <- sizeFactors(DESeq2.length.ds)",
      "mat_size_fac <- matrix(size_fac, ncol = length(size_fac), nrow = nrow(count.matrix(cdata)), byrow = T)",
      "## Extra factors",
      "extraNormFactor <- length.matrix(cdata)",
      "normFactors <- (mat_size_fac * extraNormFactor) / exp(rowMeans(log(mat_size_fac * extraNormFactor)))",
      "normalizationFactors(DESeq2.length.ds) <- as.matrix(normFactors)"),
      codefile)
  }
  writeLines(paste("DESeq2.length.ds <- DESeq2::DESeq(DESeq2.length.ds, fitType = '", fit.type, "', test = '", test, "', betaPrior = ", beta.prior, ")", sep = ""),
             codefile)
  if (impute.outliers == TRUE) {
    writeLines(c("DESeq2.length.ds.clean <- DESeq2::replaceOutliersWithTrimmedMean(DESeq2.length.ds)",
                 paste("DESeq2.length.ds.clean <- DESeq2::DESeq(DESeq2.length.ds.clean, fitType = '", fit.type, "', test = '", test, "', betaPrior = ", beta.prior, ")", sep = ""), 
                 "DESeq2.length.ds <- DESeq2.length.ds.clean"), codefile)
  }
  writeLines(c(paste("DESeq2.length.results <- DESeq2::results(DESeq2.length.ds, independentFiltering = ", independent.filtering, ", cooksCutoff = ", cooks.cutoff, ")", sep = ""),
               "DESeq2.length.pvalues <- DESeq2.length.results$pvalue", 
               "DESeq2.length.adjpvalues <- DESeq2.length.results$padj", 
               "DESeq2.length.logFC <- DESeq2.length.results$log2FoldChange", 
               "DESeq2.length.score <- 1 - DESeq2.length.pvalues", 
               "result.table <- data.frame('pvalue' = DESeq2.length.pvalues, 'adjpvalue' = DESeq2.length.adjpvalues, 'logFC' = DESeq2.length.logFC, 'score' = DESeq2.length.score)", 
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table", 
               "package.version(cdata) <- paste('DESeq2,', packageVersion('DESeq2'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'DESeq2.length', 'full.name' = '", paste('DESeq2.length.', utils::packageVersion('DESeq2'), '.', fit.type, '.', test, '.', ifelse(beta.prior == TRUE, 'bp', 'nobp'), '.', ifelse(independent.filtering == TRUE, 'indf', 'noindf'), paste(".cook_", cooks.cutoff, sep = ""), ifelse(impute.outliers, ".imp", ".noimp"), ifelse(divByLengths, ".divByLengths", ""), ifelse(!is.null(extraDesignFactors), paste0(".", paste(extraDesignFactors, collapse = ".")), ""), sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)  
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with sqrtTPM+limma
#' 
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) by applying the sqrt(TPM) transformation followed by differential expression analysis with limma. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#' 
#' For more information about the methods and the interpretation of the parameters, see the \code{limma} package and the corresponding publications.
#' 
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. The normalization factors are calculated using the \code{calcNormFactors} of the \code{edgeR} package. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}
#' @param extraDesignFactors A vector containing the extra factors to be passed to the design matrix of \code{limma}. All the factors need to be a \code{sample.annotations} from the \code{\link{compData}} object. It should not contain the "condition" factor column, that will be added automatically.
#' @param lengthNormalization one of "none" (no correction), "TPM", "RPKM" (default) or "gwRPKM". See details.
#' @param dataTransformation one of "log2", "asin(sqrt)" or "sqrt". Data transformation to apply to the normalized data.
#' @param trend should an intensity-trend be allowed for the prior variance? Default to \code{FALSE}.
#' @param blockFactor Name of the factor specifying a blocking variable, to be passed to \code{\link[limma]{duplicateCorrelation}}. All the factors need to be a \code{sample.annotations} from the \code{\link{compData}} object. Default to null (no block structure).
#' 
#' @details 
#' The \code{length.matrix} field of the \code{compData} object 
#' is used to normalize the counts, by computing the square root of the TPM.
#' 
#' @export 
#' @author Charlotte Soneson, Paul Bastide, Mélina Gallopin
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references 
#' Smyth GK (2005): Limma: linear models for microarray data. In: 'Bioinformatics and Computational Biology Solutions using R and Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds), Springer, New York, pages 397-420
#'
#' Musser, JM, Wagner, GP. (2015): Character trees from transcriptome data: Origin and individuation of morphological characters and the so‐called “species signal”. J. Exp. Zool. (Mol. Dev. Evol.) 324B: 588– 604. 
#' 
#' @examples
#' try(
#' if (require(limma)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
#'                                     samples.per.cond = 5, n.diffexp = 100, 
#'                                     id.species = factor(1:10),
#'                                     lengths.relmeans = rpois(1000, 1000),
#'                                     lengths.dispersions = rgamma(1000, 1, 1),
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "length.limma", 
#'            Rmdfunction = "lengthNorm.limma.createRmd", 
#'            output.directory = tmpdir, norm.method = "TMM")
#' })
#' 
lengthNorm.limma.createRmd <- function(data.path, result.path, codefile, norm.method, 
                                       extraDesignFactors = NULL,
                                       lengthNormalization = "RPKM",
                                       dataTransformation = "log2",
                                       trend = FALSE,
                                       blockFactor = NULL) {
  codefile <- file(codefile, open = 'w')
  writeLines("###  limma + length", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(limma)", 
               "require(edgeR)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')"),
             codefile)
  
  writeLines(c("", "# Design"),codefile)
  if (is.null(extraDesignFactors)) {
    writeLines(c(
      "design_formula <- as.formula(~ condition)",
      "design_data <- sample.annotations(cdata)[, 'condition', drop = FALSE]"),
      codefile)
  } else {
    writeLines(c(
      paste0("design_formula <- as.formula(paste(' ~ ', paste(c('", paste(extraDesignFactors, collapse = "', '"), "'), collapse= '+'), '+ condition'))"),
      paste0("design_data <- sample.annotations(cdata)[, c('", paste(extraDesignFactors, collapse = "', '"), "', 'condition'), drop = FALSE]")),
      codefile)
  }
  writeLines(c(
    "design_data <- data.frame(apply(design_data, 2, as.factor))",
    "design <- model.matrix(design_formula, design_data)"),
    codefile)
  
  writeNormalization(norm.method, lengthNormalization, dataTransformation, codefile)
  
  if (!is.null(blockFactor)) {
    if (length(blockFactor) > 1) stop("Only one factor can be taken for block definition.")
    writeLines(c("# Fitting Block correlations",
                 paste0("block <- sample.annotations(cdata)[['", paste(blockFactor), "']]"),
                 "corfit <- duplicateCorrelation(data.trans, design = design, block = block, ndups = 1)"),
               codefile)
    writeLines(c("", "# Fit"), codefile)
    writeLines(c("length.fitlimma <- limma::lmFit(data.trans, design = design, correlation = corfit$consensus, block = block)"),
               codefile)
  } else {
    writeLines(c("", "# Fit"), codefile)
    writeLines(c("length.fitlimma <- limma::lmFit(data.trans, design = design)"),
               codefile)
  }
  if (!trend) {
    writeLines("length.fitbayes <- limma::eBayes(length.fitlimma)", codefile)
  } else {
    writeLines("length.fitbayes <- limma::eBayes(length.fitlimma, trend = TRUE)", codefile)
  }
  writeLines(c("length.pvalues <- length.fitbayes$p.value[, ncol(length.fitbayes$p.value)]", 
               "length.adjpvalues <- p.adjust(length.pvalues, method = 'BH')", 
               "length.logFC <- length.fitbayes$coefficients[, ncol(length.fitbayes$coefficients)]", 
               "length.score <- 1 - length.pvalues", 
               "result.table <- data.frame('pvalue' = length.pvalues, 'adjpvalue' = length.adjpvalues, 'logFC' = length.logFC, 'score' = length.score)",
               "rownames(result.table) <- rownames(count.matrix(cdata))", 
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('limma,', packageVersion('limma'), ';', 'edgeR,', packageVersion('edgeR'))", 
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'sqrtTPM', 'full.name' = '", 
                     paste('length.', utils::packageVersion('limma'), '.limma.', norm.method,
                           ".lengthNorm.", lengthNormalization, '.',
                           "dataTrans.", dataTransformation,
                           ifelse(trend, '.with_trend', ".no_trend"),
                           ifelse(!is.null(extraDesignFactors), paste0(".", paste(extraDesignFactors, collapse = ".")), ""),
                           ifelse(!is.null(blockFactor), paste0(".", paste(blockFactor, collapse = ".")), ""),
                           sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}