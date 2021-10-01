#' Check the validity of a \code{compData} object
#' 
#' Check the validity of a \code{compData} object. An object that passes the check can be used as the input for the differential expression analysis methods interfaced by \code{compcodeR}.
#' 
#' @export 
#' @author Charlotte Soneson
#' @param object A \code{compData} object
#' @examples
#' mydata <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
#'                                 samples.per.cond = 5, n.diffexp = 100)
#' check_compData(mydata)
check_compData <- function(object) {
  if (!isS4(object)) {
    return("This is not an S4 object.")
  } else {
    if (length(count.matrix(object)) == 0) {
      return("Object must contain a non-empty count matrix.")
    }
    
    if (length(sample.annotations(object)) == 0) {
      return("Object must contain a non-empty sample annotation data frame.")
    }
    
    if (length(sample.annotations(object)$condition) == 0) {
      return("The sample.annotations must contain a column named condition.")
    }
    
    if (length(info.parameters(object)) == 0) {
      return("Object must contain a non-empty list called info.parameters.")
    }
    
    if (length(info.parameters(object)$dataset) == 0) {
      return("info.parameters list must contain an entry named 'dataset'.")
    }
    
    if (length(info.parameters(object)$uID) == 0) {
      return("info.parameters list must contain an entry named 'uID'.")
    }
    
    if (length(count.matrix(object)) != 0 && 
          length(variable.annotations(object)) != 0) {
      if (nrow(count.matrix(object)) != 
            nrow(variable.annotations(object))) {
        return("count.matrix and variable.annotations do not contain the same number of rows.")
      } else {
        if (length(rownames(count.matrix(object))) != 0 && 
            length(rownames(variable.annotations(object))) != 0 && 
            any(rownames(count.matrix(object)) != 
                rownames(variable.annotations(object)))) {
          return("The rownames of count.matrix and variable.annotations are not the same.")
        }
      } 
    } 
    
    if (length(count.matrix(object)) != 0 && 
          length(sample.annotations(object)) != 0) {
      if (ncol(count.matrix(object)) != nrow(sample.annotations(object))) {
        return("The number of columns of count.matrix is different from the number of rows of sample.annotations.")
      } else {
        if (length(colnames(count.matrix(object))) != 0 && 
            length(rownames(sample.annotations(object))) != 0 && 
            any(colnames(count.matrix(object)) != 
                rownames(sample.annotations(object)))) {
          return("The colnames of count.matrix are different from the rownames of sample.annotations.")
        }
      } 
    }
    
    if (length(result.table(object)) != 0 && length(count.matrix(object)) != 0) {
      if (nrow(result.table(object)) != nrow(count.matrix(object))) {
        return("result.table must have the same number of rows as count.matrix.")
      }
    
      if (!is.null(rownames(result.table(object))) && !is.null(rownames(count.matrix(object)))
          && !all(rownames(result.table(object)) == rownames(count.matrix(object)))) {
        return("If present, rownames of result.table must agree with those of count.matrix.")
      }
    }
  }
  return(TRUE)
}

#' Check the validity of a \code{compData} result object
#' 
#' Check the validity of a \code{compData} object containing differential expression results. An object that passes the check can be used as the input for the method comparison functions in \code{compcodeR}.
#' 
#' @export 
#' @author Charlotte Soneson
#' @param object A \code{compData} object
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
#'                                 samples.per.cond = 5, n.diffexp = 100,
#'                                 output.file = file.path(tmpdir, "mydata.rds"))
#' ## Check an object without differential expression results
#' check_compData_results(mydata)
#' 
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), 
#'            result.extent = "voom.limma", 
#'            Rmdfunction = "voom.limma.createRmd", 
#'            output.directory = tmpdir, norm.method = "TMM")
#' resdata <- readRDS(file.path(tmpdir, "mydata_voom.limma.rds"))
#' ## Check an object containing differential expression results
#' check_compData_results(resdata)
check_compData_results = function(object) {
  tmp <- check_phyloCompData(object)
  if (is.character(tmp)) {
    return(tmp)
  } else {
    if(check_phyloCompData(object)) {
      if (length(method.names(object)) == 0) {
        return("Object must contain a list named 'method.names' identifying the differential expression method used.")
      }
      
      if (length(result.table(object)) == 0) {
        return("Object must contain a data frame named 'result.table'.")
      }
      
      if (length(result.table(object)$score) == 0) {
        return("result.table must contain a column named 'score'.")
      }
    }
  }
  return(TRUE)
}

#' Class compData
#'
#' The \code{compData} class is used to store information about the experiment, such as the count matrix, sample and variable annotations, information regarding the generation of the data and results from applying a differential expression analysis to the data.
#' 
#' @section Slots:
#' \describe{
#' \item{\code{count.matrix}:}{The read count matrix, with genes as rows and samples as columns.  Class \code{matrix}}
#' \item{\code{sample.annotations}:}{A data frame containing sample annotation information for all samples in the data set. Must contain at least a column named \code{condition}, encoding the division of the samples into two classes. The row names should be the same as the column names of \code{count.matrix}. Class \code{data.frame}}
#' \item{\code{info.parameters}:}{A list of parameters detailing the simulation process used to generate the data. Must contain at least two entries, named \code{dataset} (an informative name for the data set/simulation setting) and \code{uID} (a unique ID for the specific data set instance). Class \code{list}}
#' \item{\code{filtering}:}{A character string detailing the filtering process that has been applied to the data. Class \code{character}}
#' \item{\code{variable.annotations}:}{Contains information regarding the variables, such as the differential expression status, the true mean, dispersion and effect sizes. If present, the row names should be the same as those of \code{count.matrix}. Class \code{data.frame}}
#' \item{\code{analysis.date}:}{(If a differential expression analysis has been performed and the results are included in the \code{compData} object). Gives the date when the differential expression analysis was performed. Class \code{character}}
#' \item{\code{package.version}:}{(If a differential expression analysis has been performed and the results are included in the \code{compData} object). Gives the version numbers of the package(s) used for the differential expression analysis. Class \code{character}}
#' \item{\code{method.names}:}{(If a differential expression analysis has been performed and the results are included in the \code{compData} object). A list, containing the name of the method used for the differential expression analysis. The list should have two entries: \code{full.name} and \code{short.name}, where the \code{full.name} is the full (potentially long) name identifying the method, and \code{short.name} may be an abbreviation. Class \code{list}}
#' \item{\code{code}:}{(If a differential expression analysis has been performed and the results are included in the \code{compData} object). A character string containing the code that was used to run the differential expression analysis. The code should be in R markdown format. Class \code{character}}
#' \item{\code{result.table}:}{(If a differential expression analysis has been performed and the results are included in the \code{compData} object). Contains the results of the differential expression analysis, in the form of a data frame with one row per gene. Must contain at least one column named \code{score}, where a higher value corresponds to 'more strongly differentially expressed genes'. Class \code{data.frame}}
#' }
#' 
#' @section Methods:
#' \describe{
#' \item{count.matrix}{\code{signature(x="compData")}}
#' \item{count.matrix<-}{\code{signature(x="compData",value="matrix")}:
#' Get or set the count matrix in a \code{compData} object. \code{value} should be a numeric matrix.}
#' \item{sample.annotations}{\code{signature(x="compData")}}
#' \item{sample.annotations<-}{\code{signature(x="compData",value="data.frame")}:
#' Get or set the sample annotations data frame in a \code{compData} object. \code{value} should be a data frame with at least a column named 'condition'.}
#' \item{info.parameters}{\code{signature(x="compData")}}
#' \item{info.parameters<-}{\code{signature(x="compData",value="list")}:
#' Get or set the list with info parameters in a \code{compData} object. \code{value} should be a list with at least elements named 'dataset' and 'uID'.}
#' \item{filtering}{\code{signature(x="compData")}}
#' \item{filtering<-}{\code{signature(x="compData",value="character")}:
#' Get or set the information about the filtering in a \code{compData} object. \code{value} should be a character string describing the filtering that has been performed.}
#' \item{variable.annotations}{\code{signature(x="compData")}}
#' \item{variable.annotations<-}{\code{signature(x="compData",value="data.frame")}:
#' Get or set the variable annotations data frame in a \code{compData} object. \code{value} should be a data frame.}
#' \item{analysis.date}{\code{signature(x="compData")}}
#' \item{analysis.date<-}{\code{signature(x="compData",value="character")}:
#' Get or set the analysis date in a \code{compData} object. \code{value} should be a character string describing when the differential expression analysis of the data was performed.}
#' \item{package.version}{\code{signature(x="compData")}}
#' \item{package.version<-}{\code{signature(x="compData",value="character")}:
#' Get or set the information about the package version in a \code{compData} object. \code{value} should be a character string detailing which packages and versions were used to perform the differential expression analysis of the data.}
#' \item{method.names}{\code{signature(x="compData")}}
#' \item{method.names<-}{\code{signature(x="compData",value="list")}:
#' Get or set the method names in a \code{compData} object. \code{value} should be a list with slots \code{full.name} and \code{short.name}, giving the full name and an abbreviation for the method that was used to perform the analysis of the data.}
#' \item{code}{\code{signature(x="compData")}}
#' \item{code<-}{\code{signature(x="compData",value="character")}:
#' Get or set the code slot in a \code{compData} object. \code{value} should be a character string in R markdown format, giving the code that was run to obtain the results from the differential expression analysis.}
#' \item{result.table}{\code{signature(x="compData")}}
#' \item{result.table<-}{\code{signature(x="compData",value="data.frame")}:
#' Get or set the result table in a \code{compData} object. \code{value} should be a data frame with one row per gene, and at least a column named 'score'.}
#' }
#' @section Construction:
#' An object of the class \code{compData} can be constructed using the \code{\link{compData}} function. 
#' @name compData-class
#' @rdname compData-class
#' @exportClass compData
#' @author Charlotte Soneson
setClass(
  Class = "compData", 
  
  representation = representation(count.matrix = "matrix",
                                  sample.annotations = "data.frame",
                                  info.parameters = "list",
                                  filtering = "character",
                                  variable.annotations = "data.frame",
                                  analysis.date = "character",
                                  package.version = "character",
                                  method.names = "list",
                                  code = "character",
                                  result.table = "data.frame"),
  
  validity = check_compData)

setGeneric("count.matrix", function(x) standardGeneric("count.matrix"))
setMethod("count.matrix", "compData", function(x) x@count.matrix)

setGeneric("sample.annotations", function(x, ...) standardGeneric("sample.annotations"))
setMethod("sample.annotations", "compData", function(x) x@sample.annotations)

setGeneric("info.parameters", function(x, ...) standardGeneric("info.parameters"))
setMethod("info.parameters", "compData", function(x) x@info.parameters)

setGeneric("filtering", function(x, ...) standardGeneric("filtering"))
setMethod("filtering", "compData", function(x) x@filtering)

setGeneric("variable.annotations", function(x, ...) standardGeneric("variable.annotations"))
setMethod("variable.annotations", "compData", function(x) x@variable.annotations)

setGeneric("analysis.date", function(x, ...) standardGeneric("analysis.date"))
setMethod("analysis.date", "compData", function(x) x@analysis.date)

setGeneric("package.version", function(x, ...) standardGeneric("package.version"))
setMethod("package.version", "compData", function(x) x@package.version)

setGeneric("method.names", function(x, ...) standardGeneric("method.names"))
setMethod("method.names", "compData", function(x) x@method.names)

setGeneric("code", function(x, ...) standardGeneric("code"))
setMethod("code", "compData", function(x) x@code)

setGeneric("result.table", function(x, ...) standardGeneric("result.table"))
setMethod("result.table", "compData", function(x) x@result.table)

setGeneric("count.matrix<-", function(x, value) standardGeneric("count.matrix<-"))
setReplaceMethod("count.matrix", "compData",
                 function(x, value) {x@count.matrix <- value; check_compData(x); x})

setGeneric("sample.annotations<-", function(x, value) standardGeneric("sample.annotations<-"))
setReplaceMethod("sample.annotations", "compData",
                 function(x, value) {x@sample.annotations <- value; check_compData(x); x})

setGeneric("info.parameters<-", function(x, value) standardGeneric("info.parameters<-"))
setReplaceMethod("info.parameters", "compData",
                 function(x, value) {x@info.parameters <- value; check_compData(x); x})

setGeneric("filtering<-", function(x, value) standardGeneric("filtering<-"))
setReplaceMethod("filtering", "compData",
                 function(x, value) {x@filtering <- value; check_compData(x); x})

setGeneric("variable.annotations<-", 
           function(x, value) standardGeneric("variable.annotations<-"))
setReplaceMethod("variable.annotations", "compData",
                 function(x, value) {x@variable.annotations <- value; check_compData(x); x})

setGeneric("analysis.date<-", function(x, value) standardGeneric("analysis.date<-"))
setReplaceMethod("analysis.date", "compData",
                 function(x, value) {x@analysis.date <- value; check_compData_results(x); x})

setGeneric("package.version<-", function(x, value) standardGeneric("package.version<-"))
setReplaceMethod("package.version", "compData",
                 function(x, value) {x@package.version <- value; check_compData_results(x); x})

setGeneric("method.names<-", function(x, value) standardGeneric("method.names<-"))
setReplaceMethod("method.names", "compData",
                 function(x, value) {x@method.names <- value; check_compData_results(x); x})

setGeneric("code<-", function(x, value) standardGeneric("code<-"))
setReplaceMethod("code", "compData",
                 function(x, value) {x@code <- value; check_compData_results(x); x})

setGeneric("result.table<-", function(x, value) standardGeneric("result.table<-"))
setReplaceMethod("result.table", "compData",
                 function(x, value) {x@result.table <- value; check_compData_results(x); x})

#' Create a \code{compData} object
#' 
#' The \code{compData} class is used to store information about the experiment, such as the count matrix, sample and variable annotations, information regarding the generation of the data and results from applying a differential expression analysis to the data. This constructor function creates a \code{compData} object.
#' 
#' @param count.matrix A count matrix, with genes as rows and observations as columns.
#' @param sample.annotations A data frame, containing at least one column named 'condition', encoding the grouping of the observations into two groups. The row names should be the same as the column names of the \code{count.matrix}.
#' @param info.parameters A list containing information regarding simulation parameters etc. The only mandatory entries are \code{dataset} and \code{uID}, but it may contain entries such as the ones listed below (see \code{generateSyntheticData} for more detailed information about each of these entries).
#' \itemize{
#' \item \code{dataset}: an informative name or identifier of the data set (e.g., summarizing the simulation settings).
#' \item \code{samples.per.cond}
#' \item \code{n.diffexp} 
#' \item \code{repl.id}
#' \item \code{seqdepth}
#' \item \code{minfact}
#' \item \code{maxfact}
#' \item \code{fraction.upregulated}
#' \item \code{between.group.diffdisp}
#' \item \code{filter.threshold.total}
#' \item \code{filter.threshold.mediancpm}
#' \item \code{fraction.non.overdispersed}
#' \item \code{random.outlier.high.prob}
#' \item \code{random.outlier.low.prob}
#' \item \code{single.outlier.high.prob}
#' \item \code{single.outlier.low.prob}
#' \item \code{effect.size}
#' \item \code{uID}: a unique ID for the data set. In contrast to \code{dataset}, the \code{uID} is unique e.g. for each instance of replicated data sets generated with the same simulation settings. 
#' }
#' @param variable.annotations A data frame with variable annotations (with number of rows equal to the number of rows in \code{count.matrix}, that is, the number of variables in the data set). Not mandatory, but may contain columns such as the ones listed below. If present, the row names should be the same as the row names of the \code{count.matrix}.
#' \itemize{
#' \item \code{truedispersions.S1}: the true dispersion for each gene in condition S1.
#' \item \code{truedispersions.S2}: the true dispersion for each gene in condition S2.
#' \item \code{truemeans.S1}: the true mean value for each gene in condition S1.
#' \item \code{truemeans.S2}: the true mean value for each gene in condition S2.
#' \item \code{n.random.outliers.up.S1}: the number of 'random' outliers with extremely high counts for each gene in condition S1.
#' \item \code{n.random.outliers.up.S2}: the number of 'random' outliers with extremely high counts for each gene in condition S2.
#' \item \code{n.random.outliers.down.S1}: the number of 'random' outliers with extremely low counts for each gene in condition S1.
#' \item \code{n.random.outliers.down.S2}: the number of 'random' outliers with extremely low counts for each gene in condition S2.
#' \item \code{n.single.outliers.up.S1}: the number of 'single' outliers with extremely high counts for each gene in condition S1.
#' \item \code{n.single.outliers.up.S2}: the number of 'single' outliers with extremely high counts for each gene in condition S2.
#' \item \code{n.single.outliers.down.S1}: the number of 'single' outliers with extremely low counts for each gene in condition S1.
#' \item \code{n.single.outliers.down.S2}: the number of 'single' outliers with extremely low counts for each gene in condition S2.
#' \item \code{M.value}: the M-value (observed log2 fold change between condition S1 and condition S2) for each gene.
#' \item \code{A.value}: the A-value (observed average expression level across condition S1 and condition S2) for each gene.
#' \item \code{truelog2foldchanges}: the true (simulated) log2 fold changes between condition S1 and condition S2.
#' \item \code{upregulation}: a binary vector indicating which genes are simulated to be upregulated in condition S2 compared to condition S1.
#' \item \code{downregulation}: a binary vector indicating which genes are simulated to be downregulated in condition S2 compared to condition S1.
#' \item \code{differential.expression}: a binary vector indicating which genes are simulated to be differentially expressed in condition S2 compared to condition S1.
#' }
#' @param filtering A character string containing information about the filtering that has been applied to the data set.
#' @param analysis.date If a differential expression analysis has been performed, a character string detailing when it was performed.
#' @param package.version If a differential expression analysis has been performed, a character string giving the version of the differential expression packages that were applied. 
#' @param method.names If a differential expression analysis has been performed, a list with entries \code{full.name} and \code{short.name}, giving the full name of the differential expression method (may including version number and parameter settings) and a short name or abbreviation. 
#' @param code If a differential expression analysis has been performed, a character string containing the code that was run to perform the analysis. The code should be in R markdown format, and can be written to an HTML file using the \code{\link{generateCodeHTMLs}} function.
#' @param result.table If a differential expression analysis has been performed, a data frame containing the results of the analysis. The number of rows should be equal to the number of rows in \code{count.matrix} and if present, the row names should be identical. The only mandatory column is \code{score}, which gives a score for each gene, where a higher score suggests a "more highly differentially expressed" gene. Different comparison functions use different columns of this table, if available. The list below gives the columns that are used by the interfaced methods.
#' \itemize{
#' \item \code{pvalue} nominal p-values
#' \item \code{adjpvalue} p-values adjusted for multiple comparisons
#' \item \code{logFC} estimated log-fold changes between the two conditions
#' \item \code{score} the score that will be used to rank the genes in order of significance. Note that high scores always signify differential expression, that is, a strong association with the predictor. For example, for methods returning a nominal p-value the score can be defined as 1 - pvalue.
#' \item \code{FDR} false discovery rate estimates
#' \item \code{posterior.DE} posterior probabilities of differential expression
#' \item \code{prob.DE} conditional probabilities of differential expression
#' \item \code{lfdr} local false discovery rates
#' \item \code{statistic}  test statistics from the differential expression analysis
#' \item \code{dispersion.S1} dispersion estimates in condition S1
#' \item \code{dispersion.S2} dispersion estimates in condition S2
#' }
#' 
#' @export
#' @author Charlotte Soneson
#' @return A \code{compData} object.
#' @examples
#' count.matrix <- round(matrix(1000*runif(4000), 1000))
#' sample.annotations <- data.frame(condition = c(1, 1, 2, 2))
#' info.parameters <- list(dataset = "mydata", uID = "123456")
#' cpd <- compData(count.matrix, sample.annotations, info.parameters)
compData <- function(count.matrix, sample.annotations, 
                     info.parameters, variable.annotations = data.frame(), 
                     filtering = "no info", analysis.date = "", 
                     package.version = "", method.names = list(), 
                     code = "", result.table = data.frame()) {
  
  count.matrix <- as.matrix(count.matrix)
  mode(count.matrix) <- "numeric"
  
  sample.annotations <- data.frame(sample.annotations, 
                                   stringsAsFactors = FALSE)
  
  cmpd <- new("compData", count.matrix = count.matrix, 
              sample.annotations = sample.annotations, 
              info.parameters = info.parameters, 
              variable.annotations = variable.annotations, 
              filtering = filtering, 
              analysis.date = analysis.date, 
              package.version = package.version, 
              method.names = method.names,
              code = code, 
              result.table = result.table)
  
  validObject(cmpd)
  
  cmpd
}

#' Convert a list with data and results to a \code{compData} object
#' 
#' Given a list with data and results (resulting e.g. from \code{compcodeR} version 0.1.0), convert it to a \code{compData} object.
#' 
#' @export
#' @param inp.list A list with data and results, e.g. generated by \code{compcodeR} version 0.1.0.
#' @author Charlotte Soneson
#' @examples
#' convertListTocompData(list(count.matrix = matrix(round(1000*runif(4000)), 1000),
#'                            sample.annotations = data.frame(condition = c(1,1,2,2)),
#'                            info.parameters = list(dataset = "mydata", 
#'                            uID = "123456")))
convertListTocompData <- function(inp.list) {
  if (is.null(inp.list$count.matrix) || is.null(inp.list$sample.annotations) || 
        is.null(inp.list$info.parameters) || !is.matrix(inp.list$count.matrix) || 
        !is.data.frame(inp.list$sample.annotations) || !is.list(inp.list$info.parameters)) {
    message("Cannot convert list to compData object (not all required slots are available).")
    return(NULL)
  } else {
    compData(count.matrix = inp.list$count.matrix, 
             sample.annotations = inp.list$sample.annotations, 
             info.parameters = inp.list$info.parameters, 
             variable.annotations = as.data.frame(inp.list$variable.annotations,
                                                  stringsAsFactors = FALSE), 
             filtering = as.character(inp.list$filtering), 
             analysis.date = as.character(inp.list$analysis.date), 
             package.version = as.character(inp.list$package.version), 
             method.names = as.list(inp.list$method.names), 
             code = as.character(inp.list$code), 
             result.table = as.data.frame(inp.list$result.table,
                                          stringsAsFactors = FALSE))
  }
}

#' Convert a \code{compData} object to a list
#' 
#' Given a \code{compData} object, convert it to a list.
#' 
#' @export
#' @param cpd A \code{compData} object
#' @author Charlotte Soneson
#' @examples
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 12500,
#'                                     samples.per.cond = 5, n.diffexp = 1250)
#' mydata.list <- convertcompDataToList(mydata.obj)
convertcompDataToList <- function(cpd) {
  list(count.matrix = count.matrix(cpd),
       sample.annotations = sample.annotations(cpd),
       info.parameters = info.parameters(cpd),
       variable.annotations = variable.annotations(cpd),
       filtering = filtering(cpd),
       analysis.data = analysis.date(cpd),
       package.version = package.version(cpd),
       method.names = method.names(cpd),
       code = code(cpd),
       result.table = result.table(cpd))
}

################################################################################
## PhyloCompData
################################################################################

#' Check the validity of a \code{phyloCompData} object
#' 
#' Check the validity of a \code{phyloCompData} object.
#' An object that passes the check can be used as the input for the differential expression analysis methods interfaced by \code{compcodeR}.
#' 
#' @export 
#' @author Charlotte Soneson, Paul Bastide
#' @param object A \code{phyloCompData} object
#' @examples
#' mydata <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
#'                                 samples.per.cond = 5, n.diffexp = 100,
#'                                 id.species = factor(1:10),
#'                                 tree = ape::rphylo(10, 1, 0),
#'                                 lengths.relmeans = "auto", lengths.dispersions = "auto")
#' check_phyloCompData(mydata)
#' 
check_phyloCompData <- function(object) {
  cc <- check_compData(object)
  if (!isTRUE(cc)) return(cc)
  
  ## Length matrix
  if (length(count.matrix(object)) != 0 && 
      length(length.matrix(object)) != 0) {
    if (ncol(count.matrix(object)) != ncol(length.matrix(object)) || nrow(count.matrix(object)) != nrow(length.matrix(object))) {
      return("The dimension of count.matrix is different from the dimension of length.matrix.")
    }
    if (length(colnames(count.matrix(object))) != 0 && 
        length(colnames(length.matrix(object))) != 0 && 
        any(colnames(count.matrix(object)) != 
            colnames(length.matrix(object)))) {
      return("The colnames of count.matrix are different from the colnames of length.matrix.")
    }
    if (length(rownames(count.matrix(object))) != 0 && 
        length(rownames(length.matrix(object))) != 0 && 
        any(rownames(count.matrix(object)) != 
            rownames(length.matrix(object)))) {
      return("The rownames of count.matrix are different from the rownames of length.matrix.")
    }
  }
  
  ## Tree
  if (!is.null(phylo.tree(object)) && length(phylo.tree(object)) != 0) {
    if (!inherits(phylo.tree(object), "phylo")) stop("object 'tree' is not of class 'phylo'.")
    if (length(count.matrix(object)) != 0) {
      tmp <- try(checkParamMatrix(count.matrix(object), "count.matrix", phylo.tree(object)),
                 silent = TRUE)
      if (inherits(tmp, 'try-error')) return(geterrmessage())
    }
    if (length(length.matrix(object)) != 0) {
      tmp <- try(checkParamMatrix(length.matrix(object), "length.matrix", phylo.tree(object)),
                 silent = TRUE)
      if (inherits(tmp, 'try-error')) return(geterrmessage())
    }
    if (length(sample.annotations(object)) != 0) {
      tmp <- try(checkParamMatrix(sample.annotations(object), "sample.annotations", phylo.tree(object), transpose = TRUE),
                 silent = TRUE)
      if (inherits(tmp, 'try-error')) return(geterrmessage())
    }
    if (length(sample.annotations(object)$id.species) != 0) {
      ids <- sample.annotations(object)$id.species
      names(ids) <- rownames(sample.annotations(object))
      tmp <- try(checkSpecies(ids, "id.species", phylo.tree(object), tol = 1e-10, TRUE),
                 silent = TRUE)
      if (inherits(tmp, 'try-error')) return(geterrmessage())
    }
    
    ## Sample Annotation
    if (length(sample.annotations(object)$id.species) == 0) {
      return("The sample.annotations must contain a column named id.species.")
    }
  }
  
  return(TRUE)
}

#' Class phyloCompData
#'
#' The \code{phyloCompData} class extends the \code{\linkS4class{compData}} class
#' with sequence length and phylogeny related information.
#' 
#' @section Slots:
#' \describe{
#' \item{\code{tree}:}{The phylogenetic tree describing the relationships between samples. The taxa names of the \code{tree} should be the same as the column names of the \code{count.matrix}. Class \code{phylo}.}
#' \item{\code{length.matrix}:}{The length matrix, with genes as rows and samples as columns. The column names of the \code{length.matrix} should be the same as the column names of the \code{count.matrix}. Class \code{matrix}.}
#' \item{\code{sample.annotations}:}{In addition to the columns described in the \code{\linkS4class{compData}} class, if the tree is specified, it should contain an extra column named \code{id.species} of factors giving the species for each sample. The row names should be the same as the column names of count.matrix. Class \code{data.frame}.}
#' }
#' 
#' @section Methods:
#' \describe{
#' \item{phylo.tree}{\code{signature(x="phyloCompData")}}
#' \item{phylo.tree<-}{\code{signature(x="phyloCompData",value="phylo")}:
#' Get or set the tree in a \code{phyloCompData} object. \code{value} should be a phylo object.}
#' \item{length.matrix}{\code{signature(x="phyloCompData")}}
#' \item{length.matrix<-}{\code{signature(x="phyloCompData",value="matrix")}:
#' Get or set the length matrix in a \code{phyloCompData} object. \code{value} should be a numeric matrix.}
#' }
#' @section Construction:
#' An object of the class \code{phyloCompData} can be constructed using the \code{\link{phyloCompData}} function. 
#' @name phyloCompData-class
#' @rdname phyloCompData-class
#' @exportClass phyloCompData
#' @author Charlotte Soneson, Paul Bastide
setClass(
  Class = "phyloCompData", 
  contains = "compData",
  slots = c(tree = "ANY",
            length.matrix = "matrix"),
  validity = check_phyloCompData
)

setGeneric("phylo.tree", function(x) standardGeneric("phylo.tree"))
setMethod("phylo.tree", "compData", function(x) NULL)
setMethod("phylo.tree", "phyloCompData", function(x) x@tree)

setGeneric("phylo.tree<-", function(x, value) standardGeneric("phylo.tree<-"))
setReplaceMethod("phylo.tree", "compData",
                 function(x, value) {stop("There is no 'phylo.tree' slot in a 'compData' object. Please use a 'phyloCompData' object.")})
setReplaceMethod("phylo.tree", "phyloCompData",
                 function(x, value) {x@tree <- value; check_phyloCompData(x); x})

setGeneric("length.matrix", function(x) standardGeneric("length.matrix"))
setMethod("length.matrix", "compData", function(x) NULL)
setMethod("length.matrix", "phyloCompData", function(x) x@length.matrix)

setGeneric("length.matrix<-", function(x, value) standardGeneric("length.matrix<-"))
setReplaceMethod("length.matrix", "compData",
                 function(x, value) {stop("There is no 'lenght.matrix' slot in a 'compData' object. Please use a 'phyloCompData' object.")})
setReplaceMethod("length.matrix", "phyloCompData",
                 function(x, value) {x@length.matrix <- value; check_phyloCompData(x); x})


#' Create a \code{phyloCompData} object
#' 
#' The \code{\linkS4class{phyloCompData}} class extends the \code{\link{compData}} class
#' with sequence length and phylogeny related information.
#' 
#' @inheritParams compData
#' @param tree The phylogenetic tree describing the relationships between samples. The taxa names of the \code{tree} should be the same as the column names of the \code{count.matrix}.
#' @param length.matrix The length matrix, with genes as rows and samples as columns. The column names of the \code{length.matrix} should be the same as the column names of the \code{count.matrix}.
#' @param sample.annotations A data frame, containing at least one column named 'condition', encoding the grouping of the observations into two groups, and one column named \code{id.species} of factors giving the species for each sample if the tree is specified. The row names should be the same as the column names of count.matrix. \code{Class data.frame}.
#' 
#' @export
#' @author Charlotte Soneson, Paul Bastide
#' @return A \code{phyloCompData} object.
#' @examples
#' tree <- ape::read.tree(
#'   text = "(((A1:0,A2:0,A3:0):1,B1:1):1,((C1:0,C2:0):1.5,(D1:0,D2:0):1.5):0.5);"
#'   )
#' count.matrix <- round(matrix(1000*runif(8000), 1000))
#' sample.annotations <- data.frame(condition = c(1, 1, 1, 1, 2, 2, 2, 2),
#'                                  id.species = c("A", "A", "A", "B", "C", "C", "D", "D"))
#' info.parameters <- list(dataset = "mydata", uID = "123456")
#' length.matrix <- round(matrix(1000*runif(8000), 1000))
#' colnames(count.matrix) <- colnames(length.matrix) <- rownames(sample.annotations) <- tree$tip.label
#' cpd <- phyloCompData(count.matrix, sample.annotations, info.parameters,
#'                      tree = tree, length.matrix = length.matrix)
#'
phyloCompData <- function(count.matrix, sample.annotations, 
                          info.parameters, variable.annotations = data.frame(), 
                          filtering = "no info", analysis.date = "", 
                          package.version = "", method.names = list(), 
                          code = "", result.table = data.frame(),
                          tree = list(),
                          length.matrix = matrix(NA_integer_, 0, 0)) {
  
  if (is.null(tree)) tree <- list()
  if (length(tree) == 0) attr(tree, "class") <- "phylo"
  
  cmpd <- new("phyloCompData", 
              compData(count.matrix, sample.annotations, info.parameters,
                       variable.annotations, filtering, analysis.date, 
                       package.version, method.names, code, result.table),
              tree = tree,
              length.matrix = length.matrix)
  
  validObject(cmpd)
  
  cmpd
}

#' Create a \code{phyloCompData} object
#' 
#' The \code{phyloCompData} class extends the \code{\link{compData}} class
#' with sequence length and phylogeny related information.
#' 
#' @param compDataObject An object of class \code{\link{compData}}.
#' @param tree A phylogenetic tree describing the relationships between samples.
#' @param length.matrix A length matrix, with genes as rows and observations as columns.
#' 
#' @keywords internal
#' @author Charlotte Soneson, Paul Bastide
#' @return A \code{phyloCompData} object.
#' 
phyloCompDataFromCompData <- function(compDataObject,
                                      tree = list(),
                                      length.matrix = matrix(NA_integer_, 0, 0)) {
  
  if (is.null(tree)) tree <- list()
  if (length(tree) == 0) attr(tree, "class") <- "phylo"
  
  cmpd <- new("phyloCompData", 
              compDataObject,
              tree = tree,
              length.matrix = length.matrix)
  
  validObject(cmpd)
  
  cmpd
}

#' Convert a list with data and results to a \code{phyloCompData} object
#' 
#' Given a list with data and results (resulting e.g. from \code{compcodeR} version 0.1.0), convert it to a \code{phyloCompData} object.
#' 
#' @export
#' @param inp.list A list with data and results, e.g. generated by \code{compcodeR} version 0.1.0.
#' @author Charlotte Soneson, Paul Bastide
#' @examples
#' tree <- ape::read.tree(
#'   text = "(((A1:0,A2:0,A3:0):1,B1:1):1,((C1:0,C2:0):1.5,(D1:0,D2:0):1.5):0.5);"
#'   )
#' count.matrix <- round(matrix(1000*runif(8000), 1000))
#' sample.annotations <- data.frame(condition = c(1, 1, 1, 1, 2, 2, 2, 2),
#'                                  id.species = c("A", "A", "A", "B", "C", "C", "D", "D"))
#' info.parameters <- list(dataset = "mydata", uID = "123456")
#' length.matrix <- round(matrix(1000*runif(8000), 1000))
#' colnames(count.matrix) <- colnames(length.matrix) <- rownames(sample.annotations) <- tree$tip.label
#' convertListTophyloCompData(list(count.matrix = count.matrix,
#'                                 sample.annotations = sample.annotations,
#'                                 info.parameters = list(dataset = "mydata", 
#'                                                        uID = "123456"),
#'                                 tree = tree,
#'                                 length.matrix = length.matrix))
#'                                  
convertListTophyloCompData <- function(inp.list) {
  compData <- convertListTocompData(inp.list)
  if (is.null(inp.list$tree)) inp.list$tree <- list()
  if (length(inp.list$tree) == 0) attr(inp.list$tree, "class") <- "phylo"
  phyloCompDataFromCompData(compData,
                            tree = inp.list$tree,
                            length.matrix = if (is.null(inp.list$length.matrix)) matrix(NA, 0, 0) else inp.list$length.matrix)
}

#' Convert a \code{phyloCompData} object to a list
#' 
#' Given a \code{phyloCompData} object, convert it to a list.
#' 
#' @export
#' @param cpd A \code{phyloCompData} object
#' @author Charlotte Soneson, Paul Bastide
#' @examples
#' tree <- ape::read.tree(
#'   text = "(((A1:0,A2:0,A3:0):1,B1:1):1,((C1:0,C2:0):1.5,(D1:0,D2:0):1.5):0.5);"
#'   )
#' id.species <- factor(c("A", "A", "A", "B", "C", "C", "D", "D"))
#' names(id.species) <- tree$tip.label
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
#'                                     samples.per.cond = 4, n.diffexp = 100,
#'                                     tree = tree,
#'                                     id.species = id.species)
#' mydata.list <- convertcompDataToList(mydata.obj)
#' 
convertphyloCompDataToList <- function(cpd) {
  list(count.matrix = count.matrix(cpd),
       sample.annotations = sample.annotations(cpd),
       info.parameters = info.parameters(cpd),
       variable.annotations = variable.annotations(cpd),
       filtering = filtering(cpd),
       analysis.data = analysis.date(cpd),
       package.version = package.version(cpd),
       method.names = method.names(cpd),
       code = code(cpd),
       result.table = result.table(cpd),
       tree = phylo.tree(cpd),
       length.matrix = length.matrix(cpd))
}