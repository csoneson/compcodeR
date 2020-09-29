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
#' object is used in the \code{voom} method to normalize the counts. 
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
    "voom.data$genes <- rownames(count.matrix(cdata))",
    "voom.data <- voom.data$E"), 
    codefile)
  ## Functions to apply phylolm
  writeLines(c("", "# Wrapper functions"),codefile)
  ff <- deparse(extract_results_phylolm)
  ff[1] <- paste0("extract_results_phylolm <- ", ff[1])
  writeLines(ff, codefile)
  if (lengthNormalization != "asPredictor") {
    ff <- deparse(phylolm_analysis)
    ff[1] <- paste0("phylolm_analysis <- ", ff[1])
    writeLines(ff, codefile)
  } else {
    ff <- deparse(phylolm_analysis_lengthAsPredictor)
    ff[1] <- paste0("phylolm_analysis_lengthAsPredictor <- ", ff[1])
    writeLines(ff, codefile)
  }
  ## Apply analysis
  writeLines(c("", "# Analysis"),codefile)
  extra_args <- eval(substitute(alist(...)))
  extra_args <- sapply(extra_args, function(x) paste(" = ", x))
  extra_args <- paste(names(extra_args), extra_args, collapse = ", ")
  writeLines(c("tree <- getTree(cdata)"),codefile)
  if (lengthNormalization != "asPredictor") {
    writeLines(c(
      paste0("voom.phylolm.results_list <- apply(voom.data, 1, phylolm_analysis, design_data = design_data, design_formula = design_formula, tree = tree, model = '", model, "', measurement_error = ", measurement_error, ", ", extra_args, ")"),
      "result.table <- do.call(rbind, voom.phylolm.results_list)"),
      codefile)
  } else {
    writeLines(c(
      paste0("result.table <- t(sapply(1:nrow(voom.data), phylolm_analysis_lengthAsPredictor, all_dat = voom.data, all_len = length.matrix(cdata), design_data = design_data, design_formula = design_formula, tree = tree, model = '", model, "', measurement_error = ", measurement_error, ", ", extra_args, "))"),
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
#' @param ... Further arguments to be passed to function \code{\link[phylolm]{phylolm}}.
#' 
#' @details 
#' The \code{length.matrix} field of the \code{compData} object 
#' is used to normalize the counts, by computing the square root of the TPM.
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
#'            Rmdfunction = "sqrtTPM.phylolm.createRmd", 
#'            output.directory = tmpdir,
#'            norm.method = "TMM",
#'            extraDesignFactors = c("test_factor", "test_reg"))
#' })
sqrtTPM.phylolm.createRmd <- function(data.path, result.path, codefile, 
                                      norm.method,
                                      model = "BM", measurement_error = TRUE,
                                      extraDesignFactors = NULL,
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
  writeLines(c("", "# Normalisation using sqrt TPM"),codefile)
  writeLines(c("countsLengths <- count.matrix(cdata) / length.matrix(cdata)",
               paste("nf <- edgeR::calcNormFactors(countsLengths, method = '", norm.method, "')", sep = ''),
               "lib.size <- colSums(countsLengths) * nf",
               "sqrtTPM.data <- t(sqrt(t(countsLengths) / lib.size))",
               "rownames(sqrtTPM.data) <- rownames(count.matrix(cdata))"
  ), 
  codefile)
  
  ## Functions to apply phylolm
  writeLines(c("", "# Wrapper functions"),codefile)
  ff <- deparse(extract_results_phylolm)
  ff[1] <- paste0("extract_results_phylolm <- ", ff[1])
  writeLines(ff, codefile)
  ff <- deparse(phylolm_analysis)
  ff[1] <- paste0("phylolm_analysis <- ", ff[1])
  writeLines(ff, codefile)
  ## Apply analysis
  writeLines(c("", "# Analysis"),codefile)
  extra_args <- eval(substitute(alist(...)))
  extra_args <- sapply(extra_args, function(x) paste(" = ", x))
  extra_args <- paste(names(extra_args), extra_args, collapse = ", ")
  writeLines(c("tree <- getTree(cdata)"),codefile)
  writeLines(c(
    paste0("sqrtTPM.phylolm.results_list <- apply(sqrtTPM.data, 1, phylolm_analysis, design_data = design_data, design_formula = design_formula, tree = tree, model = '", model, "', measurement_error = ", measurement_error, ", ", extra_args, ")"),
    "result.table <- do.call(rbind, sqrtTPM.phylolm.results_list)"),
    codefile)
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
    paste("method.names(cdata) <- list('short.name' = 'sqrtTPM.phylolm', 'full.name' = '", paste('sqrtTPM.phylolm', packageVersion('limma'), packageVersion('phylolm'), '.', norm.method, '.', model, '.', ifelse(!is.null(measurement_error), 'me', 'nome'), ifelse(!is.null(extraDesignFactors), paste0(".", paste(extraDesignFactors, collapse = ".")), ""), sep = ''), "')", sep = ''),
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
    return(tree)
  }
}

#' @title Extract phylom results
#'
#' @description
#' Extrat results from a phylolm object.
#' The coefficient of interest must be named "condition".
#' 
#' @param phylo_lm_obj a phylolm object.
#' 
#' @return A list, with:
#' \describe{
#' \item{pvalue}{the p value of the differential expression.}
#' \item{logFC}{the log fold change of the differential expression.}
#' \item{score}{1 - pvalue.}
#' }
#' 
#' @keywords internal
#' 
extract_results_phylolm <- function(phylo_lm_obj) {
  res <- as.data.frame(summary(phylo_lm_obj)$coefficients)
  result.table <- data.frame('pvalue' = res['condition2', 'p.value'],
                             'logFC' = res['condition2', 'Estimate'],
                             'score' = 1 - res['condition2', 'p.value'])
  return(result.table)
}

#' @title Perform the phylolm analysis
#'
#' @description
#' Perform the phylolm analysis for a given gene.
#' 
#' @param dat the data associated with a gene
#' @param design_data design matrix
#' @param design_formula design formula
#' @param tree phylogenetic tree
#' @param model the model to be used in phylolm
#' @param measurement_error boolean
#' 
#' @return A list, with:
#' \describe{
#' \item{pvalue}{the p value of the differential expression.}
#' \item{logFC}{the log fold change of the differential expression.}
#' \item{score}{1 - pvalue.}
#' }
#' 
#' @keywords internal
#' 
phylolm_analysis <- function(dat, design_data, design_formula, tree, model, measurement_error, ...) {
  data_reg <- cbind(data.frame(expr = dat), design_data)
  levels(data_reg$condition) <- c(1, 2)
  res <- try(extract_results_phylolm(phylolm::phylolm(paste('expr', paste(as.character(design_formula), collapse = '')),
                                                      data = data_reg,
                                                      phy = tree,
                                                      model = model,
                                                      measurement_error = measurement_error, 
                                                      ...)))
  if (inherits(res, 'try-error')) {
    if (model == 'BM' && measurement_error) {
      res <- try(extract_results_phylolm(phylolm::phylolm(paste('expr', paste(as.character(design_formula), collapse = '')),
                                                          data = data_reg,
                                                          phy = tree,
                                                          model = 'lambda',
                                                          measurement_error = FALSE, 
                                                          ...)))
    }
  }
  if (inherits(res, 'try-error')) {
    res <- data.frame('pvalue' = 1.0,
                      'logFC' = 0.0,
                      'score' = 0.0)
    warning(paste0('A gene produced an error.'))
  }
  return(res)
}

#' @title Perform the phylolm analysis
#'
#' @description
#' Perform the phylolm analysis for a given gene, with legnths as a predictor.
#' 
#' @param gene a gene
#' @param all_dat the matrix of all expression data
#' @param all_len the matrix of all length data
#' @param design_data design matrix
#' @param design_formula design formula
#' @param tree phylogenetic tree
#' @param model the model to be used in phylolm
#' @param measurement_error boolean
#' 
#' @return A list, with:
#' \describe{
#' \item{pvalue}{the p value of the differential expression.}
#' \item{logFC}{the log fold change of the differential expression.}
#' \item{score}{1 - pvalue.}
#' }
#' 
#' @keywords internal
#' 
phylolm_analysis_lengthAsPredictor <- function(gene, all_dat, all_len, design_data, design_formula, tree, model, measurement_error, ...) {
  data_reg <- cbind(data.frame(expr = all_dat[gene, ], lengths = log2(all_len[gene, ])), design_data)
  levels(data_reg$condition) <- c(1, 2)
  res <- try(extract_results_phylolm(phylolm::phylolm(paste('expr', paste(as.character(design_formula), collapse = ''), '+lengths'),
                                                      data = data_reg,
                                                      phy = tree,
                                                      model = model,
                                                      measurement_error = measurement_error, 
                                                      ...)))
  if (inherits(res, 'try-error')) {
    if (model == 'BM' && measurement_error) {
      res <- try(extract_results_phylolm(phylolm::phylolm(paste('expr', paste(as.character(design_formula), collapse = '')),
                                                          data = data_reg,
                                                          phy = tree,
                                                          model = 'lambda',
                                                          measurement_error = FALSE, 
                                                          ...)))
    }
  }
  if (inherits(res, 'try-error')) {
    res <- data.frame('pvalue' = 1.0,
                      'logFC' = 0.0,
                      'score' = 0.0)
    warning(paste0(gene, 'produced an error.'))
  }
  return(res)
}

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
#' @param dataTransformation one of "log2", "log2+1" or "sqrt". Data transformation to apply to the normalized data.
#' @param ... Further arguments to be passed to function \code{\link[phylolm]{phylolm}}.
#' 
#' @details 
#' The \code{length.matrix} field of the \code{compData}
#' object is used to normalize the counts. 
#' \describe{
#' \item{\code{none}:}{No length normalization.}
#' \item{\code{TPM}:}{The raw counts are divided by the length of their associated genes before normalization by \code{voom}.}
#' \item{\code{RPKM}:}{The log2 length is substracted to the log2 CPM computed by \code{voom} for each gene and sample.}
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
#'            Rmdfunction = "phylolm.createRmd", 
#'            output.directory = tmpdir,
#'            norm.method = "TMM",
#'            extraDesignFactors = c("test_factor", "test_reg"),
#'            lengthNormalization = "RPKM")
#' })
phylolm.createRmd <- function(data.path, result.path, codefile, 
                              norm.method,
                              model = "BM", measurement_error = TRUE,
                              extraDesignFactors = NULL,
                              lengthNormalization = "RPKM",
                              dataTransformation = "log2",
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
  writeNormalization(norm.method, lengthNormalization, dataTransformation, codefile)
  ## Functions to apply phylolm
  writeLines(c("", "# Wrapper functions"),codefile)
  ff <- deparse(extract_results_phylolm)
  ff[1] <- paste0("extract_results_phylolm <- ", ff[1])
  writeLines(ff, codefile)
  if (lengthNormalization != "asPredictor") {
    ff <- deparse(phylolm_analysis)
    ff[1] <- paste0("phylolm_analysis <- ", ff[1])
    writeLines(ff, codefile)
  } else {
    ff <- deparse(phylolm_analysis_lengthAsPredictor)
    ff[1] <- paste0("phylolm_analysis_lengthAsPredictor <- ", ff[1])
    writeLines(ff, codefile)
  }
  ## Apply analysis
  writeLines(c("", "# Analysis"),codefile)
  extra_args <- eval(substitute(alist(...)))
  extra_args <- sapply(extra_args, function(x) paste(" = ", x))
  extra_args <- paste(names(extra_args), extra_args, collapse = ", ")
  writeLines(c("tree <- getTree(cdata)"),codefile)
  if (lengthNormalization != "asPredictor") {
    writeLines(c(
      paste0("phylolm.results_list <- apply(data.trans, 1, phylolm_analysis, design_data = design_data, design_formula = design_formula, tree = tree, model = '", model, "', measurement_error = ", measurement_error, ", ", extra_args, ")"),
      "result.table <- do.call(rbind, phylolm.results_list)"),
      codefile)
  } else {
    writeLines(c(
      paste0("result.table <- t(sapply(1:nrow(data.trans), phylolm_analysis_lengthAsPredictor, all_dat = data.trans, all_len = length.matrix(cdata), design_data = design_data, design_formula = design_formula, tree = tree, model = '", model, "', measurement_error = ", measurement_error, ", ", extra_args, "))"),
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
    paste("method.names(cdata) <- list('short.name' = 'phylolm', 'full.name' = '",
          paste('phylolm', packageVersion('limma'), packageVersion('phylolm'), '.', norm.method, '.', model, '.',
                ifelse(!is.null(measurement_error), 'me', 'nome'), '.',
                "lengthNorm.", lengthNormalization, '.',
                "dataTrans.", dataTransformation,
                ifelse(!is.null(extraDesignFactors), paste0(".", paste(extraDesignFactors, collapse = ".")), ""),
                sep = ''),
          "')", sep = ''),
    "is.valid <- check_compData_results(cdata)",
    "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
    paste("saveRDS(cdata, '", result.path, "')", sep = "")),
    codefile)  
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}

#' Generate a \code{.Rmd} file containing code to normalize data.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. The normalization factors are calculated using the \code{calcNormFactors} of the \code{edgeR} package. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}
#' @param lengthNormalization one of "none" (no correction), "TPM", "RPKM" (default) or "gwRPKM". See details.
#' @param dataTransformation one of "log2", "log2+1" or "sqrt." Data transformation to apply to the normalized data.
#' @param codefile 
#' 
#' @details 
#' The \code{length.matrix} field of the \code{compData}
#' object is used to normalize the counts. 
#' \describe{
#' \item{\code{none}:}{No length normalization.}
#' \item{\code{TPM}:}{The raw counts are divided by the length of their associated genes before normalization by \code{voom}.}
#' \item{\code{RPKM}:}{The log2 length is substracted to the log2 CPM computed by \code{voom} for each gene and sample.}
#' }
#' 
#' @keywords internal
#' 
writeNormalization <- function(norm.method, lengthNormalization, dataTransformation, codefile) {
  writeLines(c("", "# Normalisation"),codefile)
  lengthNormalization <- match.arg(lengthNormalization, c("RPKM", "TPM", "none", "gwRPKM"))
  dataTransformation <- match.arg(dataTransformation, c("log2", "log2+1", "sqrt"))
  if (lengthNormalization == "none" || lengthNormalization == "asPredictor") {
    writeLines(c(paste("nf <- edgeR::calcNormFactors(count.matrix(cdata), method = '", norm.method, "')", sep = ''),
                 "lib.size <- colSums(count.matrix(cdata)) * nf"),
               codefile)
    if (dataTransformation == "log2") {
      writeLines("data.norm <- sweep(count.matrix(cdata) + 0.5, 2, lib.size + 1, '/')", codefile)
    } else {
      writeLines("data.norm <- sweep(count.matrix(cdata), 2, lib.size, '/')", codefile)
    }
    writeLines("data.norm <- data.norm * 1e6", codefile)
  } else if (lengthNormalization == "TPM") {
    writeLines(c(paste("nf <- edgeR::calcNormFactors(count.matrix(cdata) / length.matrix(cdata), method = '", norm.method, "')", sep = ''),
                 "lib.size <- colSums(count.matrix(cdata) / length.matrix(cdata)) * nf"),
               codefile)
    if (dataTransformation == "log2") {
      writeLines("data.norm <- sweep((count.matrix(cdata) + 0.5) / length.matrix(cdata), 2, lib.size + 1, '/')", codefile)
    } else {
      writeLines("data.norm <- sweep((count.matrix(cdata)) / length.matrix(cdata), 2, lib.size, '/')", codefile)
    }
    writeLines("data.norm <- data.norm * 1e6", codefile)
  } else if (lengthNormalization == "RPKM") {
    writeLines(c(paste("nf <- edgeR::calcNormFactors(count.matrix(cdata), method = '", norm.method, "')", sep = ''),
                 "lib.size <- colSums(count.matrix(cdata)) * nf"),
               codefile)
    if (dataTransformation == "log2") {
      writeLines("data.norm <- sweep((count.matrix(cdata) + 0.5) / length.matrix(cdata), 2, lib.size + 1, '/')", codefile)
    } else {
      writeLines("data.norm <- sweep((count.matrix(cdata)) / length.matrix(cdata), 2, lib.size, '/')", codefile)
    }
    writeLines("data.norm <- data.norm * 1e9", codefile)
  } else if (lengthNormalization == "gwRPKM") {
    if (dataTransformation != "log2") stop("gwRPKM normalisation is only available for the log2 transformation.")
    writeLines(c(paste("nf <- edgeR::calcNormFactors(count.matrix(cdata), method = '", norm.method, "')", sep = ''),
                 "lib.size <- colSums(count.matrix(cdata)) * nf",
                 "data.norm.none <- sweep(count.matrix(cdata) + 0.5, 2, lib.size + 1, '/')",
                 "data.norm.none <- data.norm.none * 1e9",
                 "get_gwRPKM <- function(i) {",
                 "    fitLength <- lm(log2(data.norm.none[i, ]) ~ log2(length.matrix(cdata)[i, ]))",
                 "    if (is.na(fitLength$coefficients[2])) return(log2(data.norm.none[i, ])) ## All lengths are the same",
                 "    return(log2(data.norm.none[i, ]) - fitLength$coefficients[2] * log2(length.matrix(cdata)[i, ]))",
                 "}",
                 "data.trans <- t(sapply(1:nrow(data.norm.none), get_gwRPKM))"),
               codefile)
  }
  if (lengthNormalization != "gwRPKM") {
    writeLines(c("", "# Transformation"),codefile)
    if (dataTransformation == "log2") {
      writeLines("data.trans <- log2(data.norm)", codefile)
    } else if (dataTransformation == "log2+1") {
      writeLines("data.trans <- log2(data.norm + 1)", codefile)
    } else if (dataTransformation == "sqrt") {
      writeLines("data.trans <- sqrt(data.norm)", codefile)
    }
  }
  writeLines("rownames(data.trans) <- rownames(count.matrix(cdata))", codefile)
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with \code{\link[phylolm]{phylolm}}.
#' 
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the phylolm package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#' 
#' For more information about the methods and the interpretation of the parameters, see the \code{\link[phylolm]{phylolm}} package and the corresponding publications. 
#' 
#' @param data.path The path to a .rds file containing the \code{compData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param fit.type The fitting method used to get the dispersion-mean relationship. Possible values are \code{"parametric"} and \code{"local"}.
#' @param model The model for trait evolution on the tree. Default to "BM".
#' @param measurement_error A logical value indicating whether there is measurement error. Default to TRUE.
#' @param extraDesignFactors A vector containing the extra factors to be passed to the design matrix of \code{DESeq2}. All the factors need to be a \code{sample.annotations} from the \code{\link{compData}} object. It should not contain the "condition" factor column, that will be added automatically.
#' @param ... Further arguments to be passed to function \code{\link[phylolm]{phylolm}}.
#' 
#' @details 
#' The \code{length.matrix} field of the \code{compData}
#' object is used to normalize the counts. 
#' \describe{
#' \item{\code{none}:}{No length normalization.}
#' \item{\code{TPM}:}{The raw counts are divided by the length of their associated genes before normalization by \code{voom}.}
#' \item{\code{RPKM}:}{The log2 length is substracted to the log2 CPM computed by \code{voom} for each gene and sample.}
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
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "vst.phylolm", 
#'            Rmdfunction = "vst.phylolm.createRmd", 
#'            output.directory = tmpdir,
#'            fit.type = "parametric",
#'            extraDesignFactors = c("test_factor", "test_reg"))
#' })
vst.phylolm.createRmd <- function(data.path, result.path, codefile, 
                                  fit.type,
                                  model = "BM", measurement_error = TRUE,
                                  extraDesignFactors = NULL,
                                  ...) {
  codefile <- file(codefile, open = 'w')
  writeLines("### phylolm", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = TRUE, error = TRUE, warning = TRUE}", 
               "require(phylolm)", 
               "require(limma)", 
               "require(DESeq2)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTocompData(cdata)", codefile)
  }
  
  writeLines(c("is.valid <- check_compData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid compData object.')"),
             codefile)
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
  
  ## Normalization with DESeq2
  writeLines(c("", "## Normalization with DESeq2"), codefile)
  writeLines(c(
    paste("DESeq2.length.ds <- DESeq2::DESeqDataSetFromMatrix(countData = count.matrix(cdata),",
          "colData = design_data,",
          "design = design_formula)")),
    codefile)
  writeLines(c(
    "# Size Factors",
    "DESeq2.length.ds <-  estimateSizeFactors(DESeq2.length.ds)",
    "size_fac <- sizeFactors(DESeq2.length.ds)",
    "mat_size_fac <- matrix(size_fac, ncol = length(size_fac), nrow = nrow(count.matrix(cdata)), byrow = T)",
    "# Extra factors",
    "extraNormFactor <- length.matrix(cdata)",
    "normFactors <- (mat_size_fac * extraNormFactor) / exp(rowMeans(log(mat_size_fac * extraNormFactor)))",
    "normalizationFactors(DESeq2.length.ds) <- as.matrix(normFactors)"),
    codefile)
  writeLines(c("# VST",
               paste("DESeq2.length.ds <- DESeq2::estimateDispersions(DESeq2.length.ds, fitType = '", fit.type, "')", sep = ""),
               "vst.data <- DESeq2::getVarianceStabilizedData(DESeq2.length.ds)"),
             codefile)
  
  ## Functions to apply phylolm
  writeLines(c("", "# Wrapper functions"),codefile)
  ff <- deparse(extract_results_phylolm)
  ff[1] <- paste0("extract_results_phylolm <- ", ff[1])
  writeLines(ff, codefile)
  ff <- deparse(phylolm_analysis)
  ff[1] <- paste0("phylolm_analysis <- ", ff[1])
  writeLines(ff, codefile)
  ## Apply analysis
  writeLines(c("", "# Analysis"),codefile)
  extra_args <- eval(substitute(alist(...)))
  extra_args <- sapply(extra_args, function(x) paste(" = ", x))
  extra_args <- paste(names(extra_args), extra_args, collapse = ", ")
  writeLines(c("tree <- getTree(cdata)"),codefile)
  writeLines(c(
    paste0("phylolm.results_list <- apply(vst.data, 1, phylolm_analysis, design_data = design_data, design_formula = design_formula, tree = tree, model = '", model, "', measurement_error = ", measurement_error, ", ", extra_args, ")"),
    "result.table <- do.call(rbind, phylolm.results_list)"),
    codefile)
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
    paste("method.names(cdata) <- list('short.name' = 'phylolm', 'full.name' = '",
          paste('vst.phylolm', packageVersion('DESeq2'), packageVersion('limma'), packageVersion('phylolm'), '.', fit.type, '.', model, '.',
                ifelse(!is.null(measurement_error), 'me', 'nome'), '.',
                ifelse(!is.null(extraDesignFactors), paste0(".", paste(extraDesignFactors, collapse = ".")), ""),
                sep = ''),
          "')", sep = ''),
    "is.valid <- check_compData_results(cdata)",
    "if (!(is.valid == TRUE)) stop('Not a valid compData result object.')",
    paste("saveRDS(cdata, '", result.path, "')", sep = "")),
    codefile)  
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}