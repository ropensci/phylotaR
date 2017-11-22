
#' @name genTDObj
#' @title Generate R object from NCBI taxonomy dump
#' @description Generates an interrogable R object
#' from the NCBI taxonomy dump using the \code{CHNOSZ}
#' library. The returned object if a list of taxonomic nodes
#' and names.
#' @export
genTDObj <- function() {
  # TODO: add cache with parameters
  nodes <- CHNOSZ::getnodes()
  names <- CHNOSZ::getnames()
  tdobj <- list('nodes'=nodes,
                'names'=names)
  tdobj
}