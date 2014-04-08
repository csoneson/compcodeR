#' Check a list or a \code{compData} object for compatibility with the differential expression functions interfaced by \code{compcodeR}
#' 
#' Check if a list or a \code{compData} object contains the necessary slots for applying the differential expression functions interfaced by the \code{compcodeR} package. This function is provided for backward compatibility, see also \code{\link{check_compData}} and \code{\link{check_compData_results}}.
#' 
#' @param data.obj A list containing data and condition information, or a \code{compData} object.
#' @author Charlotte Soneson
#' @export
#' @examples
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
#'                                     samples.per.cond = 5, n.diffexp = 100)
#' checkDataObject(mydata.obj)
checkDataObject <- function(data.obj) {
  if (is.list(data.obj)) data.obj <- convertListTocompData(data.obj)
  if (check_compData(data.obj) == TRUE) {
    "Data object looks ok."
  }
}
