# This file holds the pipeline functions

#' @name runPhylota
#' @title Run the entire PhyLoTa pipeline
#' @description Run the entire PhyLoTa pipeline.
#' All generated files will be stored in the wd.
#' The process can be stopped at anytime and 
#' restarted with \code{rstrtPhylota}.
#' @param wd Working directory
#' @param txid Root taxonomic ID
#' @details All objects and data will be cached.
#' @export
#' @seealso \code{\link{rstrtPhylota}}
runPhylota <- function(wd, txid) {
  # TODO
  # Taxonomy steps
  #  - genTxNds
  # Cluser steps
}

#' @name rstrtPhylota
#' @title Restart a PhyLoTa pipeline run
#' @description Restarts the running of a pipeline
#' as started with \code{runPhylota}.
#' @param wd Working directory
#' @export
#' @seealso \code{\link{runPhylota}}
rstrtPhylota <- function(wd) {
  # TODO
}