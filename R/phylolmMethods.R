#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with \code{\link[phylolm]{phylolm}}.
#' 
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the phylolm package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#' 
#' For more information about the methods and the interpretation of the parameters, see the \code{\link[phylolm]{phylolm}} package and the corresponding publications. 
#' 
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. The normalization factors are calculated using the \code{calcNormFactors} of the \code{edgeR} package. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}
#' @param model The model for trait evolution on the tree. Default to "BM".
#' @param measurement_error A logical value indicating whether there is measurement error. Default to TRUE.
#' @param extraDesignFactors A vector containing the extra factors to be passed to the design matrix of \code{DESeq2}. All the factors need to be a \code{sample.annotations} from the \code{\link{compData}} object. It should not contain the "condition" factor column, that will be added automatically.
#' @param lengthNormalization one of "none" (no correction), "TPM", "RPKM" (default) or "gwRPKM". See details.
#' @param ... Further arguments to be passed to function \code{\link[phylolm]{phylolm}}.
#' 
#' @details 
#' The \code{length.matrix} field of the \code{compData}
#' object are used in the \code{voom} method to normalize the counts. 
#' \describe{
#' \item{\code{none}:}{No length normalization.}
#' \item{\code{TPM}:}{The raw counts are divided by the length of their associated genes before normalization by \code{voom}.}
#' \item{\code{RPKM}:}{The log2 length is substracted to the log2 CPM computed by \code{voom} for each gene and sample.}
#' \item{\code{gwRPKM}:}{The log2 CPM computed by \code{voom} is regressed against the log2 lengths for each genes, and then the contribution of the length is substracted from the score.}
#' \item{\code{asPredictor}:}{No length normalization, but the log2 lengths are included, for each gene, to the (phylogenetic) linear model as a predictor.}
#' }
#' 
#' @export 
#' 
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references 
#' Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @examples
#' try(
#' if (require(ape) && require(phylolm)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' set.seed(20200317)
#' tree <- rphylo(10, 0.1, 0)
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
#'                                     samples.per.cond = 5, n.diffexp = 100, 
#'                                     tree = tree,
#'                                     id.species = 1:10,
#'                                     lengths.relmeans = rpois(1000, 1000),
#'                                     lengths.dispersions = rgamma(1000, 1, 1),
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' ## Add annotations
#' sample.annotations(mydata.obj)$test_factor <- factor(rep(1:2, 5))
#' saveRDS(mydata.obj, file.path(tmpdir, "mydata.rds"))
#' ## Diff Exp
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "DESeq2", 
#'            Rmdfunction = "voom.phylolm.createRmd", 
#'            output.directory = tmpdir,
#'            norm.method = "TMM",
#'            extraDesignFactors = c("test_factor", "test_reg"),
#'            lengthNormalization = "RPKM")
#' })
voom.phylolm.createRmd <- function(data.path, result.path, codefile, 
                                   norm.method,
                                   model = "BM", measurement_error = TRUE,
                                   extraDesignFactors = NULL,
                                   lengthNormalization = "RPKM",
                                   ...) {
  codefile <- file(codefile, open = 'w')
  writeLines("### phylolm", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = TRUE, error = TRUE, warning = TRUE}", 
               "require(phylolm)", 
               "require(limma)", 
               "require(edgeR)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')"),
             codefile)
  writeLines("count_data <- count.matrix(cdata)", codefile)
  ## Design for normalization
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
  writeLines(c("", "# Normalisation using voom"),codefile)
  lengthNormalization <- match.arg(lengthNormalization, c("gwRPKM", "RPKM", "TPM", "none", "asPredictor"))
  if (lengthNormalization == "none" || lengthNormalization == "asPredictor") {
    writeLines(c(paste("nf <- edgeR::calcNormFactors(count.matrix(cdata), method = '", norm.method, "')", sep = ''),
                 "voom.data <- limma::voom(count.matrix(cdata), design = design, lib.size = colSums(count.matrix(cdata)) * nf)"), 
               codefile)
  } else if (lengthNormalization == "TPM") {
    writeLines(c(paste("nf <- edgeR::calcNormFactors(count.matrix(cdata) / length.matrix(cdata), method = '", norm.method, "')", sep = ''),
                 "voom.data <- limma::voom(count.matrix(cdata) / length.matrix(cdata), design = design, lib.size = colSums(count.matrix(cdata) / length.matrix(cdata)) * nf)"), 
               codefile)
  } else if (lengthNormalization == "RPKM") {
    writeLines(c(paste("nf <- edgeR::calcNormFactors(count.matrix(cdata), method = '", norm.method, "')", sep = ''),
                 "voom.data <- limma::voom(count.matrix(cdata), design = design, lib.size = colSums(count.matrix(cdata)) * nf * t(length.matrix(cdata)))"), 
               codefile)
  } else if (lengthNormalization == "gwRPKM") {
    writeLines(c(paste("nf <- edgeR::calcNormFactors(count.matrix(cdata), method = '", norm.method, "')", sep = ''),
                 "lib.size <- colSums(count.matrix(cdata)) * nf",
                 "y <- t(log2(t(count.matrix(cdata)+0.5)/(lib.size+1)*1e6))",
                 "length_factor <- sapply(1:nrow(y), function(i) lm(y[i, ] ~ log2(length.matrix(cdata))[i, ])$coefficients[2])",
                 "length_factor <- sweep(length.matrix(cdata), 1, length_factor, '^')",
                 "voom.data <- limma::voom(count.matrix(cdata), design = design, lib.size = colSums(count.matrix(cdata)) * nf * t(length_factor))"),
               codefile)
  }
  writeLines(c(
    "voom.data$genes <- rownames(count_data)",
    "voom.data <- voom.data$E"), 
    codefile)
  ## Functions to apply phylolm
  writeLines(c("", "# Wrapper functions"),codefile)
  writeLines(c(
    "extract_results_phylolm <- function(phylo_lm_obj) {",
    "  res <- as.data.frame(summary(phylo_lm_obj)$coefficients)",
    "  result.table <- data.frame('pvalue' = res['condition2', 'p.value'],",
    "                             'logFC' = res['condition2', 'Estimate'],",
    "                             'score' = 1 - res['condition2', 'p.value'])",
    "  return(result.table)",
    "}"
  ),
  codefile)
  if (lengthNormalization != "asPredictor") {
    writeLines(c(
      "phylolm_analysis <- function(dat, model, measurement_error, ...) {",
      "  data_reg <- cbind(data.frame(expr = dat), design_data)",
      "  levels(data_reg$condition) <- c(1, 2)",
      "  res <- try(extract_results_phylolm(phylolm(paste('expr', paste(as.character(design_formula), collapse = '')),",
      "                                             data = data_reg,",
      "                                             phy = tree,",
      "                                             model = model,",
      "                                             measurement_error = measurement_error, ",
      "                                             ...)))",
      "  if (inherits(res, 'try-error')) {",
      "    if (model == 'BM' && measurement_error) {",
      "      res <- try(extract_results_phylolm(phylolm(paste('expr', paste(as.character(design_formula), collapse = '')),",
      "                                                 data = data_reg,",
      "                                                 phy = tree,",
      "                                                 model = 'lambda',",
      "                                                 measurement_error = FALSE, ",
      "                                                 ...)))",
      "    }",
      "  }",
      "  if (inherits(res, 'try-error')) {",
      "    res <- data.frame('pvalue' = 1.0,",
      "                      'logFC' = 0.0,",
      "                      'score' = 0.0)",
      "    warning(paste0(gene, 'produced an error.'))",
      "  }",
      "  return(res)",
      "}"
    ),
    codefile)
  } else {
    writeLines(c(
      "phylolm_analysis <- function(gene, model, measurement_error, ...) {",
      "  data_reg <- cbind(data.frame(expr = voom.data[gene, ], lengths = log2(length.matrix(cdata)[gene, ])), design_data)",
      "  levels(data_reg$condition) <- c(1, 2)",
      "  res <- try(extract_results_phylolm(phylolm(paste('expr', paste(as.character(design_formula), collapse = ''), '+lengths'),",
      "                                             data = data_reg,",
      "                                             phy = tree,",
      "                                             model = model,",
      "                                             measurement_error = measurement_error, ",
      "                                             ...)))",
      "  if (inherits(res, 'try-error')) {",
      "    if (model == 'BM' && measurement_error) {",
      "      res <- try(extract_results_phylolm(phylolm(paste('expr', paste(as.character(design_formula), collapse = '')),",
      "                                                 data = data_reg,",
      "                                                 phy = tree,",
      "                                                 model = 'lambda',",
      "                                                 measurement_error = FALSE, ",
      "                                                 ...)))",
      "    }",
      "  }",
      "  if (inherits(res, 'try-error')) {",
      "    res <- data.frame('pvalue' = 1.0,",
      "                      'logFC' = 0.0,",
      "                      'score' = 0.0)",
      "    warning(paste0(gene, 'produced an error.'))",
      "  }",
      "  return(res)",
      "}"
    ),
    codefile)
  }
  ## Apply analysis
  writeLines(c("", "# Analysis"),codefile)
  extra_args <- eval(substitute(alist(...)))
  extra_args <- sapply(extra_args, function(x) paste(" = ", x))
  extra_args <- paste(names(extra_args), extra_args, collapse = ", ")
  writeLines(c("tree <- getTree(cdata)"),codefile)
  if (lengthNormalization != "asPredictor") {
  writeLines(c(
    paste0("voom.phylolm.results_list <- apply(voom.data, 1, phylolm_analysis, model = '", model, "', measurement_error = ", measurement_error, ", ", extra_args, ")"),
    "result.table <- do.call(rbind, voom.phylolm.results_list)"),
    codefile)
  } else {
    writeLines(c(
      paste0("result.table <- t(sapply(1:nrow(voom.data), phylolm_analysis, model = '", model, "', measurement_error = ", measurement_error, ", ", extra_args, "))"),
      "result.table <- as.data.frame(result.table)"),
      codefile)
  }
  writeLines(c(
    "result.table$adjpvalue <- p.adjust(result.table$pvalue, 'BH')"),
    codefile)
  writeLines(c("", "# Save the results"),codefile)
  writeLines(c(
    "rownames(result.table) <- rownames(count.matrix(cdata))",
    "result.table(cdata) <- result.table", 
    "package.version(cdata) <- paste('phylolm,', packageVersion('phylolm'))",
    "package.version(cdata) <- paste('limma,', packageVersion('limma'))",
    "analysis.date(cdata) <- date()",
    paste("method.names(cdata) <- list('short.name' = 'voom.phylolm', 'full.name' = '", paste('voom.phylolm', packageVersion('limma'), packageVersion('phylolm'), '.', norm.method, '.', model, '.', ifelse(!is.null(measurement_error), 'me', 'nome'), '.', "lengthNorm.", lengthNormalization, ifelse(!is.null(extraDesignFactors), paste0(".", paste(extraDesignFactors, collapse = ".")), ""), sep = ''), "')", sep = ''),
    "is.valid <- check_compData_results(cdata)",
    "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
    paste("saveRDS(cdata, '", result.path, "')", sep = "")),
    codefile)  
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

getTree <- function(cdata) {
  if (is.null(cdata@info.parameters$tree)) {
    message("There were no tree in the data object. Using a star tree of unit height in phylolm.")
    ntaxa <- cdata@info.parameters$samples.per.cond * 2
    tree <- ape::stree(ntaxa, "star")
    tree$edge.length <- rep(1, nrow(tree$edge))
    tree$tip.label <- rownames(cdata@sample.annotations)
    return(tree)
  } else {
    tree <- cdata@info.parameters$tree
    if (!ape::is.ultrametric(tree)) stop("The tree must be ultrametric.")
    return(phytools::force.ultrametric(tree))
  }
}

