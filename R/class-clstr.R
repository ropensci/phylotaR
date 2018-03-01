chckClstr <- function(object) {
  # TODO
  TRUE
}

#' @name Clstr-class
#' @aliases Clstr-method
#' @param x \code{Clstr} object
#' @param object \code{Clstr} object
#' @title Clstr-class
#' @description Cluster record contains all information on a cluster.
#' @slot id Cluster ID, integer
#' @slot sids Sequence IDs
#' @slot txids Source txids for sequences
#' @slot typ Cluster type: direct, subtree or merged
#' @slot seed Seed sequence ID
#' @slot prnt Parent taxonomic ID
#' @exportClass Clstr
#' @seealso 
#' \code{\link{genClstr}}
setClass('Clstr', representation=representation(
  id='integer',
  sids='vector',
  txids='vector',
  typ='character',
  prnt='character',
  seed='character'),
  validity=chckClstr)

#' @rdname Clstr-class
#' @exportMethod as.character
setMethod('as.character', c('x'='Clstr'),
          function(x) {
            msg <- paste0('Cluster Record [', x@id,']\n')
            msg <- paste0(msg, ' - [', x@typ,
                          '] type\n')
            msg <- paste0(msg, ' - [', x@seed,
                          '] seed sequence\n')
            msg <- paste0(msg, ' - [', length(x@sids),
                          '] sequences\n')
            msg <- paste0(msg, ' - [', length(unique(x@txids)),
                          '] unique txids\n')
          })
#' @rdname Clstr-class
#' @exportMethod show
setMethod('show', 'Clstr',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname Clstr-class
#' @exportMethod print
setMethod('print', 'Clstr',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname Clstr-class
#' @exportMethod str
setMethod('str', c('object'='Clstr'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
#' @rdname Clstr-class
#' @exportMethod summary
setMethod('summary', c('object'='Clstr'),
          function(object){
            msg <- as.character(x)
            cat(msg)
          })