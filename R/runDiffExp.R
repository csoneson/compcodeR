checkClass <- function(object, objname, trueclass) {
  if (!(class(object) == trueclass)) {
    stop(paste("The object", objname, "should be of class", trueclass))
  }
}

#' The main function to run differential expression analysis
#' 
#' The main function for running differential expression analysis (comparing two conditions), using one of the methods interfaced through \code{compcodeR} or a user-defined method.
#' 
#' @param data.file The path to a \code{.rds} file containing the data on which the differential expression analysis will be performed, for example a \code{compData} object returned from \code{\link{generateSyntheticData}}.
#' @param result.extent The extension that will be added to the data file name in order to construct the result file name. This can be for example the differential expression method together with a version number. 
#' @param Rmdfunction A function that creates an Rmd file containing the code that should be run to perform the differential expression analysis. All functions available through \code{compcodeR} can be listed using the \code{\link{listcreateRmd}} function. 
#' @param output.directory The directory in which the result object will be saved.
#' @param norm.path Logical, whether to include the full (absolute) path to the output object in the saved code. 
#' @param ... Additional arguments that will be passed to the \code{Rmdfunction}, such as parameter choices for the differential expression method.
#' @export
#' 
#' @author Charlotte Soneson
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
#'                                     samples.per.cond = 5, n.diffexp = 100, 
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' listcreateRmd()
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "voom.limma", 
#'            Rmdfunction = "voom.limma.createRmd", 
#'            output.directory = tmpdir, norm.method = "TMM")
#' 
#' \dontrun{
#' ## The following list covers the currently available 
#' differential expression methods:
#' runDiffExp(data.file = "mydata.rds", result.extent = "baySeq", 
#'            Rmdfunction = "baySeq.createRmd",
#'            output.directory = ".", norm.method = "edgeR", 
#'            distr.choice = "NB", equaldisp = TRUE)
#' runDiffExp(data.file = "mydata.rds", result.extent = "DESeq.GLM", 
#'            Rmdfunction = "DESeq.GLM.createRmd",
#'            output.directory = ".", sharing.mode = "maximum", 
#'            disp.method = "pooled", fit.type = "parametric")
#' runDiffExp(data.file = "mydata.rds", result.extent = "DESeq.nbinom", 
#'            Rmdfunction = "DESeq.nbinom.createRmd",
#'            output.directory = ".", sharing.mode = "maximum", 
#'            disp.method = "pooled", fit.type = "parametric")
#' runDiffExp(data.file = "mydata.rds", result.extent = "DESeq2", 
#'            Rmdfunction = "DESeq2.createRmd",
#'            output.directory = ".", fit.type = "parametric", 
#'            test = "Wald", beta.prior = TRUE, 
#'            independent.filtering = TRUE, cooks.cutoff = TRUE, 
#'            impute.outliers = TRUE)
#' runDiffExp(data.file = "mydata.rds", result.extent = "DSS", 
#'            Rmdfunction = "DSS.createRmd",
#'            output.directory = ".", norm.method = "quantile", 
#'            disp.trend = TRUE)
#' runDiffExp(data.file = "mydata.rds", result.extent = "EBSeq", 
#'            Rmdfunction = "EBSeq.createRmd",
#'            output.directory = ".", norm.method = "median")
#' runDiffExp(data.file = "mydata.rds", result.extent = "edgeR.exact", 
#'            Rmdfunction = "edgeR.exact.createRmd",
#'            output.directory = ".", norm.method = "TMM", 
#'            trend.method = "movingave", disp.type = "tagwise")
#' runDiffExp(data.file = "mydata.rds", result.extent = "edgeR.GLM", 
#'            Rmdfunction = "edgeR.GLM.createRmd",
#'            output.directory = ".", norm.method = "TMM", 
#'            disp.type = "tagwise", disp.method = "CoxReid", 
#'            trended = TRUE)
#' runDiffExp(data.file = "mydata.rds", result.extent = "logcpm.limma", 
#'            Rmdfunction = "logcpm.limma.createRmd",
#'            output.directory = ".", norm.method = "TMM")
#' runDiffExp(data.file = "mydata.rds", result.extent = "NBPSeq", 
#'            Rmdfunction = "NBPSeq.createRmd",
#'            output.directory = ".", norm.method = "TMM", 
#'            disp.method = "NBP")
#' runDiffExp(data.file = "mydata.rds", result.extent = "NOISeq", 
#'            Rmdfunction = "NOISeq.prenorm.createRmd",
#'            output.directory = ".", norm.method = "TMM")
#' runDiffExp(data.file = "mydata.rds", result.extent = "SAMseq", 
#'            Rmdfunction = "SAMseq.createRmd",
#'            output.directory = ".")
#' runDiffExp(data.file = "mydata.rds", result.extent = "sqrtcpm.limma", 
#'            Rmdfunction = "sqrtcpm.limma.createRmd",
#'            output.directory = ".", norm.method = "TMM")
#' runDiffExp(data.file = "mydata.rds", result.extent = "TCC", 
#'            Rmdfunction = "TCC.createRmd",
#'            output.directory = ".", norm.method = "tmm", 
#'            test.method = "edger", iteration = 3, 
#'            normFDR = 0.1, floorPDEG = 0.05)
#' runDiffExp(data.file = "mydata.rds", result.extent = "ttest", 
#'            Rmdfunction = "ttest.createRmd",
#'            output.directory = ".", norm.method = "TMM")
#' runDiffExp(data.file = "mydata.rds", result.extent = "voom.limma", 
#'            Rmdfunction = "voom.limma.createRmd", 
#'            output.directory = ".", norm.method = "TMM")
#' runDiffExp(data.file = "mydata.rds", result.extent = "voom.ttest", 
#'            Rmdfunction = "voom.ttest.createRmd",
#'            output.directory = ".", norm.method = "TMM")
#' runDiffExp(data.file = "mydata.rds", result.extent = "vst.limma", 
#'            Rmdfunction = "vst.limma.createRmd",
#'            output.directory = ".", fit.type = "parametric")
#' runDiffExp(data.file = "mydata.rds", result.extent = "vst.ttest", 
#'            Rmdfunction = "vst.ttest.createRmd",
#'            output.directory = ".", fit.type = "parametric")
#' }
runDiffExp <- function(data.file, result.extent, Rmdfunction, output.directory = ".", norm.path = TRUE, ...) {
  code.file <- "tempcode.Rmd"
  checkClass(data.file, "data.file", "character")
  checkClass(result.extent, "result.extent", "character")
  checkClass(Rmdfunction, "Rmdfunction", "character")
  checkClass(output.directory, "output.directory", "character")
  
  if (norm.path) {
    output.directory <- normalizePath(output.directory, winslash = "/")
  } else {
    ## Just remove any trailing "/"
    if (substr(output.directory, nchar(output.directory), nchar(output.directory)) == "/") {
      output.directory <- substr(output.directory, 1, nchar(output.directory) - 1)
    }
  }

  ## Extract the file name of the data file (after the last /), and extend it with the 
  ## result.extent (before the .rds) to create the result file name.
  data.file.name <- basename(data.file)
  result.file <- file.path(output.directory, 
                           gsub(".rds", paste("_", result.extent, ".rds", sep = ""), 
                                data.file.name))
  
  ## Create the .Rmd file with the differential expression code
  eval(parse(text = Rmdfunction))(data.file, result.file, code.file, ...)
  input.Rmd <- code.file
  
  ## Run the differential expression analysis
  knit(input = input.Rmd, output = str_replace(code.file, ".Rmd", ".md"))
  
  ## Add the code to the result object
  temp.file <- readRDS(result.file)
  code(temp.file) <- readLines(str_replace(code.file, ".Rmd", ".md"))
  saveRDS(temp.file, file = result.file)

  ## Remove the temporary code file
  file.remove(code.file)
  file.remove(str_replace(code.file, ".Rmd", ".md")) 
}
