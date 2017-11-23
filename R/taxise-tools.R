
#' @name genTDObj
#' @title Generate R object from NCBI taxonomy dump
#' @description Generates an interrogable R object
#' from the NCBI taxonomy dump using the \code{CHNOSZ}
#' library. The returned object is a list of taxonomic nodes
#' and names.
#' @param wd Working directory
#' @details Object will be cached.
#' @export
genTDObj <- function(wd) {
  if(!chkObj(wd, 'tdobj')) {
    nodes <- CHNOSZ::getnodes()
    names <- CHNOSZ::getnames()
    tdobj <- list('nodes'=nodes,
                  'names'=names)
    svObj(wd=wd, obj=tdobj, nm='tdobj')
  } else {
    tdobj <- ldObj(wd, 'tdobj')
  }
  tdobj
}