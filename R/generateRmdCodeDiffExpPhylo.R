#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with DESeq2 with custom model matrix
#' 
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the DESeq2 package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#' 
#' For more information about the methods and the interpretation of the parameters, see the \code{DESeq2} package and the corresponding publications. 
#' 
#' @inheritParams DESeq2.createRmd
#' @param data.path The path to a .rds file containing the \code{phyloCompData} object that will be used for the differential expression analysis.
#' @param extra.design.covariates A vector containing the names of extra control variables to be passed to the design matrix of \code{DESeq2}. All the covariates need to be a column of the \code{sample.annotations} data frame from the \code{\link{phyloCompData}} object, with a matching column name. The covariates can be a numeric vector, or a factor. Note that "condition" factor column is always included, and should not be added here. See Details.
#' 
#' @details 
#' The lengths matrix is used as a normalization factor and applied to the \code{DESeq2}
#' model in the way explained in \code{\link[DESeq2]{normalizationFactors}}
#' (see examples of this function).
#' The provided matrix will be multiplied by the default normalization factor 
#' obtained through the \code{\link[DESeq2]{estimateSizeFactors}} function.
#' 
#' The \code{design} model used in the \code{\link[DESeq2]{DESeqDataSetFromMatrix}}
#' uses the "condition" column of the \code{sample.annotations} data frame from the \code{\link{phyloCompData}} object
#' as well as all the covariates named in \code{extra.design.covariates}.
#' For example, if \code{extra.design.covariates = c("var1", "var2")}, then
#' \code{sample.annotations} must have two columns named "var1" and "var2", and the design formula
#' in the \code{\link[DESeq2]{DESeqDataSetFromMatrix}} function will be:
#' \code{~ condition + var1 + var2}.
#' 
#' @export 
#' @author Charlotte Soneson, Paul Bastide, Mélina Gallopin
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references 
#' Anders S and Huber W (2010): Differential expression analysis for sequence count data. Genome Biology 11:R106
#' 
#' Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8.
#' @examples
#' try(
#' if (require(DESeq2)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' ## Simulate data
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
#'                                     samples.per.cond = 5, n.diffexp = 100, 
#'                                     id.species = 1:10,
#'                                     lengths.relmeans = rpois(1000, 1000),
#'                                     lengths.dispersions = rgamma(1000, 1, 1),
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' ## Add covariates
#' ## Model fitted is count.matrix ~ condition + test_factor + test_reg
#' sample.annotations(mydata.obj)$test_factor <- factor(rep(1:2, each = 5))
#' sample.annotations(mydata.obj)$test_reg <- rnorm(10, 0, 1)
#' saveRDS(mydata.obj, file.path(tmpdir, "mydata.rds"))
#' ## Diff Exp
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "DESeq2", 
#'            Rmdfunction = "DESeq2.length.createRmd", 
#'            output.directory = tmpdir, fit.type = "parametric",
#'            test = "Wald",
#'            extra.design.covariates = c("test_factor", "test_reg"))
#' })
DESeq2.length.createRmd <- function(data.path, result.path, codefile, 
                                    fit.type, test, beta.prior = TRUE, 
                                    independent.filtering = TRUE, cooks.cutoff = TRUE, 
                                    impute.outliers = TRUE,
                                    extra.design.covariates = NULL,
                                    nas.as.ones = FALSE) {
  codefile <- file(codefile, open = 'w')
  writeLines("### DESeq2.length", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = TRUE, error = TRUE, warning = TRUE}", 
               "require(DESeq2)", 
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTophyloCompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_phyloCompData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid phyloCompData object.')"),
             codefile)
  writeLines("count_matrix <- count.matrix(cdata)", codefile)
  if (is.null(extra.design.covariates)) {
    writeLines(c(
      paste("DESeq2.length.ds <- DESeq2::DESeqDataSetFromMatrix(countData = count_matrix,",
            "colData = data.frame(condition = factor(sample.annotations(cdata)$condition)),",
            "design = ~ condition)")),
      codefile)
  } else {
    writeLines(c(
      paste0("DESeq2.length.ds <- DESeq2::DESeqDataSetFromMatrix(countData = count_matrix,",
             " colData = cbind(sample.annotations(cdata)[, c('", paste(extra.design.covariates, collapse = "', '"), "'), drop = FALSE], data.frame(condition = factor(sample.annotations(cdata)$condition))),",
             " design = as.formula(paste(' ~ ', paste(c('", paste(extra.design.covariates, collapse = "', '"), "'), collapse= '+'), '+ condition')))")),
      codefile)
  }
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
  writeLines(paste("DESeq2.length.ds <- DESeq2::DESeq(DESeq2.length.ds, fitType = '", fit.type, "', test = '", test, "', betaPrior = ", beta.prior, ")", sep = ""),
             codefile)
  if (impute.outliers == TRUE) {
    writeLines(c("DESeq2.length.ds.clean <- DESeq2::replaceOutliersWithTrimmedMean(DESeq2.length.ds)",
                 paste("DESeq2.length.ds.clean <- DESeq2::DESeq(DESeq2.length.ds.clean, fitType = '", fit.type, "', test = '", test, "', betaPrior = ", beta.prior, ")", sep = ""), 
                 "DESeq2.length.ds <- DESeq2.length.ds.clean"), codefile)
  }
  writeLines(paste("DESeq2.length.results <- DESeq2::results(DESeq2.length.ds, independentFiltering = ", independent.filtering, ", cooksCutoff = ", cooks.cutoff, ")", sep = ""),
             codefile)  
  if (nas.as.ones) {
    # see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#i-want-to-benchmark-deseq2-comparing-to-other-de-tools.
    message("As `nas.as.ones=TRUE`, all NAs in adjusted p values are replaced by 1 to allow for benchmarking with other methods. For more details, see section 'I want to benchmark DESeq2 comparing to other DE tools' from the `DESeq2` vignette (available by running `vignette('DESeq2', package = 'DESeq2')`)")
    writeLines("DESeq2.length.results$padj <- ifelse(is.na(DESeq2.length.results$padj), 1, DESeq2.length.results$padj)", codefile)  
  } else {
    message("As `nas.as.ones=FALSE`, there might be some NAs in the adjusted p values computed by DESeq2. This might bias the comparison of the results with other methods. For more details, see section 'I want to benchmark DESeq2 comparing to other DE tools' from the `DESeq2` vignette (available by running `vignette('DESeq2', package = 'DESeq2')`)")
  }
  writeLines(c("DESeq2.length.pvalues <- DESeq2.length.results$pvalue", 
               "DESeq2.length.adjpvalues <- DESeq2.length.results$padj", 
               "DESeq2.length.logFC <- DESeq2.length.results$log2FoldChange", 
               "DESeq2.length.score <- 1 - DESeq2.length.pvalues", 
               "result.table <- data.frame('pvalue' = DESeq2.length.pvalues, 'adjpvalue' = DESeq2.length.adjpvalues, 'logFC' = DESeq2.length.logFC, 'score' = DESeq2.length.score)", 
               "rownames(result.table) <- rownames(count.matrix(cdata))",
               "result.table(cdata) <- result.table", 
               "package.version(cdata) <- paste('DESeq2,', packageVersion('DESeq2'))",
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'DESeq2.length', 'full.name' = '", paste('DESeq2.length.', utils::packageVersion('DESeq2'), '.', fit.type, '.', test, '.', ifelse(beta.prior == TRUE, 'bp', 'nobp'), '.', ifelse(independent.filtering == TRUE, 'indf', 'noindf'), paste(".cook_", cooks.cutoff, sep = ""), ifelse(impute.outliers, ".imp", ".noimp"), ifelse(!is.null(extra.design.covariates), paste0(".", paste(extra.design.covariates, collapse = ".")), ""), sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid phyloCompData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)  
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with length normalized counts + limma
#' 
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) by applying a length normalizing transformation followed by differential expression analysis with limma. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#' 
#' For more information about the methods and the interpretation of the parameters, see the \code{limma} package and the corresponding publications.
#' 
#' @inheritParams DESeq2.length.createRmd
#' @inheritParams voom.limma.createRmd 
#' @param extra.design.covariates A vector containing the names of extra control variables to be passed to the design matrix of \code{limma}. All the covariates need to be a column of the \code{sample.annotations} data frame from the \code{\link{phyloCompData}} object, with a matching column name. The covariates can be a numeric vector, or a factor. Note that "condition" factor column is always included, and should not be added here. See Details.
#' @param length.normalization one of "none" (no length correction), "TPM", or "RPKM" (default). See details.
#' @param data.transformation one of "log2", "asin(sqrt)" or "sqrt". Data transformation to apply to the normalized data.
#' @param trend should an intensity-trend be allowed for the prior variance? Default to \code{FALSE}.
#' @param block.factor Name of the factor specifying a blocking variable, to be passed to \code{\link[limma]{duplicateCorrelation}} function of the \code{limma} package. All the factors need to be a \code{sample.annotations} from the \code{\link{phyloCompData}} object. Default to null (no block structure).
#' 
#' @details 
#' The \code{length.matrix} field of the \code{phyloCompData} object 
#' is used to normalize the counts, using one of the following formulas:
#' * \code{length.normalization="none"} : \eqn{CPM_{gi} = \frac{N_{gi} + 0.5}{NF_i \times \sum_{g} N_{gi} + 1} \times 10^6}
#' * \code{length.normalization="TPM"} : \eqn{TPM_{gi} = \frac{(N_{gi} + 0.5) / L_{gi}}{NF_i \times \sum_{g} N_{gi}/L_{gi} + 1} \times 10^6}
#' * \code{length.normalization="RPKM"} : \eqn{RPKM_{gi} = \frac{(N_{gi} + 0.5) / L_{gi}}{NF_i \times \sum_{g} N_{gi} + 1} \times 10^9}
#' 
#' where \eqn{N_{gi}} is the count for gene g and sample i,
#' where \eqn{L_{gi}} is the length of gene g in sample i,
#' and \eqn{NF_i} is the normalization for sample i,
#' normalized using \code{calcNormFactors} of the \code{edgeR} package.
#' 
#' The function specified by the \code{data.transformation} is then applied
#' to the normalized count matrix.
#' 
#' The "\eqn{+0.5}" and "\eqn{+1}" are taken from Law et al 2014,
#' and dropped from the normalization 
#' when the transformation is something else than \code{log2}.
#' 
#' The "\eqn{\times 10^6}" and "\eqn{\times 10^9}" factors are omitted when
#' the \code{asin(sqrt)} transformation is taken, as \eqn{asin} can only
#' be applied to real numbers smaller than 1.
#' 
#' The \code{design} model used in the \code{\link[limma]{lmFit}}
#' uses the "condition" column of the \code{sample.annotations} data frame from the \code{\link{phyloCompData}} object
#' as well as all the covariates named in \code{extra.design.covariates}.
#' For example, if \code{extra.design.covariates = c("var1", "var2")}, then
#' \code{sample.annotations} must have two columns named "var1" and "var2", and the design formula
#' in the \code{\link[limma]{lmFit}} function will be:
#' \code{~ condition + var1 + var2}.
#' 
#' @md
#' 
#' @export 
#' @author Charlotte Soneson, Paul Bastide, Mélina Gallopin
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references 
#' Smyth GK (2005): Limma: linear models for microarray data. In: 'Bioinformatics and Computational Biology Solutions using R and Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds), Springer, New York, pages 397-420
#' 
#' Smyth, G. K., Michaud, J., and Scott, H. (2005). The use of within-array replicate spots for assessing differential expression in microarray experiments. Bioinformatics 21(9), 2067-2075.
#' 
#' Law, C.W., Chen, Y., Shi, W. et al. (2014) voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol 15, R29.
#'
#' Musser, JM, Wagner, GP. (2015): Character trees from transcriptome data: Origin and individuation of morphological characters and the so‐called “species signal”. J. Exp. Zool. (Mol. Dev. Evol.) 324B: 588– 604. 
#' 
#' @examples
#' try(
#' if (require(limma)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' ## Simulate data
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
#'                                     samples.per.cond = 5, n.diffexp = 100, 
#'                                     id.species = factor(1:10),
#'                                     lengths.relmeans = rpois(1000, 1000),
#'                                     lengths.dispersions = rgamma(1000, 1, 1),
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' ## Add covariates
#' ## Model fitted is count.matrix ~ condition + test_factor + test_reg
#' sample.annotations(mydata.obj)$test_factor <- factor(rep(1:2, each = 5))
#' sample.annotations(mydata.obj)$test_reg <- rnorm(10, 0, 1)
#' saveRDS(mydata.obj, file.path(tmpdir, "mydata.rds"))
#' ## Diff Exp
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "length.limma", 
#'            Rmdfunction = "lengthNorm.limma.createRmd", 
#'            output.directory = tmpdir, norm.method = "TMM")
#' })
#' 
lengthNorm.limma.createRmd <- function(data.path, result.path, codefile, norm.method, 
                                       extra.design.covariates = NULL,
                                       length.normalization = "RPKM",
                                       data.transformation = "log2",
                                       trend = FALSE,
                                       block.factor = NULL) {
  if (!is.null(block.factor)) {
    if (!requireNamespace("statmod", quietly = TRUE)) stop("Package `statmod` is required for correlation modeling with `block.factor`.")
  }
  codefile <- file(codefile, open = 'w')
  writeLines("###  limma + length", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(limma)", 
               "require(edgeR)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTophyloCompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_phyloCompData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid phyloCompData object.')"),
             codefile)
  
  writeLines(c("", "# Design"),codefile)
  if (is.null(extra.design.covariates)) {
    writeLines(c(
      "design_formula <- as.formula(~ condition)",
      "design_data <- sample.annotations(cdata)[, 'condition', drop = FALSE]"),
      codefile)
  } else {
    writeLines(c(
      paste0("design_formula <- as.formula(paste(' ~ ', paste(c('", paste(extra.design.covariates, collapse = "', '"), "'), collapse= '+'), '+ condition'))"),
      paste0("design_data <- sample.annotations(cdata)[, c('", paste(extra.design.covariates, collapse = "', '"), "', 'condition'), drop = FALSE]")),
      codefile)
  }
  writeLines(c(
    "design_data$condition <- factor(design_data$condition)",
    "design <- model.matrix(design_formula, design_data)"),
    codefile)
  
  writeNormalization(norm.method, length.normalization, data.transformation, codefile)
  
  if (!is.null(block.factor)) {
    if (length(block.factor) > 1) stop("Only one factor can be taken for block definition.")
    writeLines(c("# Fitting Block correlations",
                 paste0("block <- sample.annotations(cdata)[['", paste(block.factor), "']]"),
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
                           ".lengthNorm.", length.normalization, '.',
                           "dataTrans.", data.transformation,
                           ifelse(trend, '.with_trend', ".no_trend"),
                           ifelse(!is.null(extra.design.covariates), paste0(".", paste(extra.design.covariates, collapse = ".")), ""),
                           ifelse(!is.null(block.factor), paste0(".", paste(block.factor, collapse = ".")), ""),
                           sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid phyloCompData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with length normalized counts + limma
#' 
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) by applying a length normalizing transformation followed by differential expression analysis with limma. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#' 
#' For more information about the methods and the interpretation of the parameters, see the \code{limma} package and the corresponding publications.
#' 
#' @inheritParams lengthNorm.limma.createRmd
#' @param n.sv The number of surogate variables to estimate (see \code{\link[sva]{sva}}). Default to \code{"auto"}: will be estimated with \code{\link[sva]{num.sv}}.
#' 
#' @details 
#' 
#' See the \code{details} section of \code{\link{lengthNorm.limma.createRmd}} for details
#' on the normalization and the extra design covariates.
#' 
#' @md
#' 
#' @export 
#' @author Charlotte Soneson, Paul Bastide, Mélina Gallopin
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references 
#' Smyth GK (2005): Limma: linear models for microarray data. In: 'Bioinformatics and Computational Biology Solutions using R and Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds), Springer, New York, pages 397-420
#' 
#' Smyth, G. K., Michaud, J., and Scott, H. (2005). The use of within-array replicate spots for assessing differential expression in microarray experiments. Bioinformatics 21(9), 2067-2075.
#' 
#' Law, C.W., Chen, Y., Shi, W. et al. (2014) voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol 15, R29.
#'
#' Musser, JM, Wagner, GP. (2015): Character trees from transcriptome data: Origin and individuation of morphological characters and the so‐called “species signal”. J. Exp. Zool. (Mol. Dev. Evol.) 324B: 588– 604. 
#' 
#' Leek JT, Johnson WE, Parker HS, Jaffe AE, and Storey JD. (2012) The sva package for removing batch effects and other unwanted variation in high-throughput experiments. Bioinformatics DOI:10.1093/bioinformatics/bts034
#' 
#' @examples
#' try(
#' if (require(limma) && require(sva)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' ## Simulate data
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
#'                                     samples.per.cond = 5, n.diffexp = 100, 
#'                                     id.species = factor(1:10),
#'                                     lengths.relmeans = rpois(1000, 1000),
#'                                     lengths.dispersions = rgamma(1000, 1, 1),
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' ## Add covariates
#' ## Model fitted is count.matrix ~ condition + test_factor + test_reg
#' sample.annotations(mydata.obj)$test_factor <- factor(rep(1:2, each = 5))
#' sample.annotations(mydata.obj)$test_reg <- rnorm(10, 0, 1)
#' saveRDS(mydata.obj, file.path(tmpdir, "mydata.rds"))
#' ## Diff Exp
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "lengthNorm.sva.limma", 
#'            Rmdfunction = "lengthNorm.sva.limma.createRmd", 
#'            output.directory = tmpdir, norm.method = "TMM")
#' })
#' 
lengthNorm.sva.limma.createRmd <- function(data.path, result.path, codefile, norm.method, 
                                           extra.design.covariates = NULL,
                                           length.normalization = "RPKM",
                                           data.transformation = "log2",
                                           trend = FALSE,
                                           # block.factor = NULL,
                                           n.sv = "auto") {
  # if (!is.null(block.factor)) {
  #   if (!requireNamespace("statmod", quietly = TRUE)) stop("Package `statmod` is required for correlation modeling with `block.factor`.")
  # }
  if (!requireNamespace("sva", quietly = TRUE)) stop("Package `sva` is required for the SVA analysis.")
  codefile <- file(codefile, open = 'w')
  writeLines("###  limma + length + SVA", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = FALSE, error = TRUE, warning = TRUE}",
               "require(limma)", 
               "require(edgeR)",
               "require(sva)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTophyloCompData(cdata)", codefile)
  }
  writeLines(c("is.valid <- check_phyloCompData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid phyloCompData object.')"),
             codefile)
  
  writeLines(c("", "# Design"),codefile)
  if (is.null(extra.design.covariates)) {
    writeLines(c(
      "design_formula <- as.formula(~ condition)",
      "design_formula_0 <- as.formula(~ 1)",
      "design_data <- sample.annotations(cdata)[, 'condition', drop = FALSE]"),
      codefile)
  } else {
    writeLines(c(
      paste0("design_formula <- as.formula(paste(' ~ ', paste(c('", paste(extra.design.covariates, collapse = "', '"), "'), collapse= '+'), '+ condition'))"),
      paste0("design_formula_0 <- as.formula(paste(' ~ ', paste(c('", paste(extra.design.covariates, collapse = "', '"), "'), collapse= '+')))"),
      paste0("design_data <- sample.annotations(cdata)[, c('", paste(extra.design.covariates, collapse = "', '"), "', 'condition'), drop = FALSE]")),
      codefile)
  }
  writeLines(c(
    "design_data$condition <- factor(design_data$condition)",
    "design <- model.matrix(design_formula, design_data)",
    "design_0 <- model.matrix(design_formula_0, design_data)"),
    codefile)
  
  writeNormalization(norm.method, length.normalization, data.transformation, codefile)
  
  writeLines(c("", "# SVA"), codefile)
  if (n.sv == "auto") {
    writeLines("n.sv <- sva::num.sv(data.trans, design, method = 'leek')", codefile)
  } else {
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    if (!is.wholenumber(n.sv)) stop("'n.sv' must be either 'auto' or a whole number.")
    writeLines(paste0("n.sv <- ", paste0(n.sv)), codefile)
  }
  n.sv <- ifelse(is.null(n.sv), "NULL", paste0(n.sv))
  writeLines(c("svobj <- sva::sva(data.trans, design, design_0, n.sv =  n.sv, method = 'irw')",
               "design <- cbind(design, svobj$sv)"),
             codefile)
  
  # if (!is.null(block.factor)) {
  #   if (length(block.factor) > 1) stop("Only one factor can be taken for block definition.")
  #   writeLines(c("# Fitting Block correlations",
  #                paste0("block <- sample.annotations(cdata)[['", paste(block.factor), "']]"),
  #                "corfit <- duplicateCorrelation(data.trans, design = design, block = block, ndups = 1)"),
  #              codefile)
  #   writeLines(c("", "# Fit"), codefile)
  #   writeLines(c("length.fitlimma <- limma::lmFit(data.trans, design = design, correlation = corfit$consensus, block = block)"),
  #              codefile)
  # } else {
  #   writeLines(c("", "# Fit"), codefile)
  writeLines(c("length.fitlimma <- limma::lmFit(data.trans, design = design)"),
             codefile)
  # }
  if (!trend) {
    writeLines("length.fitbayes <- limma::eBayes(length.fitlimma)", codefile)
  } else {
    writeLines("length.fitbayes <- limma::eBayes(length.fitlimma, trend = TRUE)", codefile)
  }
  writeLines(c("length.pvalues <- length.fitbayes$p.value[, 'condition2']", 
               "length.adjpvalues <- p.adjust(length.pvalues, method = 'BH')", 
               "length.logFC <- length.fitbayes$coefficients[, 'condition2']", 
               "length.score <- 1 - length.pvalues", 
               "length.n.sv <- svobj$n.sv", 
               "result.table <- data.frame('pvalue' = length.pvalues, 'adjpvalue' = length.adjpvalues, 'logFC' = length.logFC, 'score' = length.score, 'n.sv' = length.n.sv)",
               "rownames(result.table) <- rownames(count.matrix(cdata))", 
               "result.table(cdata) <- result.table",
               "package.version(cdata) <- paste('limma,', packageVersion('limma'), ';', 'edgeR,', packageVersion('edgeR'))", 
               "analysis.date(cdata) <- date()",
               paste("method.names(cdata) <- list('short.name' = 'sqrtTPM', 'full.name' = '", 
                     paste('length.',
                           utils::packageVersion('sva'), '.sva.',
                           utils::packageVersion('limma'), '.limma.', norm.method,
                           ".lengthNorm.", length.normalization, '.',
                           "dataTrans.", data.transformation,
                           ifelse(trend, '.with_trend', ".no_trend"),
                           ifelse(!is.null(extra.design.covariates), paste0(".", paste(extra.design.covariates, collapse = ".")), ""),
                           ".n.sv.", n.sv,
                           # ifelse(!is.null(block.factor), paste0(".", paste(block.factor, collapse = ".")), ""),
                           sep = ''), "')", sep = ''),
               "is.valid <- check_compData_results(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid phyloCompData result object.')",
               paste("saveRDS(cdata, '", result.path, "')", sep = "")), codefile)
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}