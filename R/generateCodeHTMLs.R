#' Generate HTML file(s) containing code used to run differential expression analysis.
#' 
#' A function to extract the code used to generate differential expression results from saved \code{compData} result objects (typically obtained by \code{\link{runDiffExp}}), and to write the code to HTML files. This requires that the code was saved as a character string in R markdown format in the \code{code} slot of the result object, which is done automatically by \code{\link{runDiffExp}}. If the differential expression analysis was performed with functions outside \code{compcodeR}, the code has to be added manually to the result object. 
#' 
#' @param input.files A vector with paths to one or several \code{.rds} files containing \code{compData} objects with the results from differential expression analysis. One code HTML file is generated for each file in the vector. 
#' @param output.directory The path to the directory where the code HTML files will be saved.
#' @export
#' @author Charlotte Soneson
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
#'                                     samples.per.cond = 5, n.diffexp = 100, 
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "voom.limma", 
#'            Rmdfunction = "voom.limma.createRmd", output.directory = tmpdir, 
#'            norm.method = "TMM")
#' generateCodeHTMLs(file.path(tmpdir, "mydata_voom.limma.rds"), tmpdir)
generateCodeHTMLs <- function(input.files, output.directory) {
  input.files <- normalizePath(input.files, winslash = "/")
  output.directory <- normalizePath(output.directory, winslash = "/")
  
  for (i in 1:length(input.files)) {
    temp <- as.character(input.files[i])
    if (!(substr(basename(temp), nchar(basename(temp)) - 3, nchar(basename(temp))) == ".rds")) {
      message(paste(temp, "is not an .rds file, no code file will be generated."))
    } else {
      ## Extract only the filename (that will be used to name the code file)
      filename <- sub(".rds", "", basename(temp))
      rdsobj <- readRDS(temp)
      ## If it is a list, convert to a compData object
      if (is.list(rdsobj)) {
        rdsobj <- convertListTocompData(rdsobj)
      }
      if (!is.logical(check_compData(rdsobj))) {
        message(paste(temp, "is not a valid compData object, no code file will be generated."))
      } else {
        if (length(code(rdsobj)) == 0) {
          message(paste("Code not available in", input.files[i]))
        } else {
          ## Create the code file
          codefile <- file.path(output.directory, paste(filename,  "_code.md", sep = ""))
          fileConn <- file(codefile)
          
          ## Write the code to the file
          writeLines(text = code(rdsobj), con = fileConn)
          close(fileConn)
          
          ## Generate HTML file
          markdownToHTML(file = codefile,
                         output = str_replace(codefile, pattern = ".md", replacement = ".html"),
                         title = 'Code')
          
          ## Remove the temporary markdown file
          file.remove(codefile)
        }
      }
    }
  }
}