#' Show method for \code{compData} object
#' 
#' Show method for \code{compData} object. 
#' @param object A \code{compData} object
#' @author Charlotte Soneson
#' @examples 
#' mydata <- generateSyntheticData(dataset = "mydata", n.vars = 12500, 
#'                                 samples.per.cond = 5, n.diffexp = 1250)
#' mydata
setMethod(
  f = "show",  # name of function
  signature = "compData",  # class of each argument. other ex: signature = c(x = , y = )
  definition = function(object) {
    nbr.random.outliers <- 
      sum(object@variable.annotations$n.random.outliers.up.S1 + 
            object@variable.annotations$n.random.outliers.down.S1 + 
            object@variable.annotations$n.random.outliers.up.S2 + 
            object@variable.annotations$n.random.outliers.down.S2)
    nbr.single.outliers <- 
      sum(object@variable.annotations$n.single.outliers.up.S1 + 
            object@variable.annotations$n.single.outliers.down.S1 + 
            object@variable.annotations$n.single.outliers.up.S2 + 
            object@variable.annotations$n.single.outliers.down.S2)
    nbr.diffexp <- sum(object@variable.annotations$differential.expression)
    if (length(object@method.names) == 0) {
      diffexp.statement <- "No differential expression analysis has been performed."
    } else {
      diffexp.statement <- 
        paste("Differential expression analysis was performed by the method", 
              object@method.names$full.name, "on", 
              object@analysis.date, ".")
    }
    cat(paste("An object of class", class(object), "\n"))
    cat(paste("Dataset name:", object@info.parameters$dataset, "\n"))
    cat(paste("Number of samples:", ncol(object@count.matrix), "\n"))
    cat(paste("Number of variables:", nrow(object@count.matrix), "\n"))
    cat(paste("Number of random outliers:", nbr.random.outliers, "\n"))
    cat(paste("Number of single outliers:", nbr.single.outliers, "\n"))
    cat(paste("Number of known truly differentially expressed genes:", nbr.diffexp, "\n"))
    cat(paste(diffexp.statement, "\n\n"))
    cat("count.matrix:\n")
    print(head(object@count.matrix[, 1:min(ncol(object@count.matrix), 6)]))
    nbrleftoutrows <- max(c(0, nrow(object@count.matrix) - 6))
    nbrleftoutcols <- max(c(0, ncol(object@count.matrix) - 6))
    if (nbrleftoutrows > 0 & nbrleftoutcols > 0) {
      cat(paste("+", nbrleftoutrows, "rows and", nbrleftoutcols, "cols...\n\n"))
    } else if (nbrleftoutrows > 0) {
      cat(paste("+", nbrleftoutrows, "rows...\n\n"))
    } else if (nbrleftoutcols > 0) {
      cat(paste("+", nbrleftoutrows, "cols...\n\n"))
    }
    cat("sample.annotations:\n")
    print(head(object@sample.annotations))
    nbrleftoutrows <- max(c(0, nrow(object@sample.annotations) - 6))
    if (nbrleftoutrows > 0) {
      cat(paste("+", nbrleftoutrows, "rows...\n\n"))
    } else {
      cat("\n")
    }
    cat("variable.annotations:\n")
    print(head(object@variable.annotations))
    nbrleftoutrows <- max(c(0, nrow(object@variable.annotations) - 6))
    if (nbrleftoutrows > 0) {
      cat(paste("+", nbrleftoutrows, "rows...\n\n"))
    } else {
      cat("\n")
    }
    
    if (length(object@method.names) != 0) {
      cat("Differential expression results:\n\n")
      print(head(object@result.table))
      nbrleftoutrows <- max(c(0, nrow(object@result.table) - 6))
      if (nbrleftoutrows > 0) {
        cat(paste("+", nbrleftoutrows, "rows...\n\n"))
      } else {
        cat("\n")
      }
    }
    
    if (length(object@length.matrix != 0)) {
      cat("length.matrix:\n")
      print(head(object@length.matrix[, 1:min(ncol(object@length.matrix), 6)]))
      cat("\n")
    } else {
      cat("No length matrix. \n")
    }
  }
)

#setGeneric("show")
