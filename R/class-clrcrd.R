chckClRcrd <- function(object) {
  # TODO
  TRUE
}

#' @name ClRcrd-class
#' @aliases ClRcrd-method
#' @param x \code{ClRcrd} object
#' @param object \code{ClRcrd} object
#' @title ClRcrd-class
#' @description Cluster record contains all information on a cluster.
#' @slot id Cluster ID, integer
#' @slot sids Sequence IDs
#' @slot txids Source txids for sequences
#' @slot typ Cluster type: direct, subtree or merged
#' @slot seed Seed sequence ID
#' @slot prnt Parent taxonomic ID
#' @exportClass ClRcrd
#' @seealso 
#' \code{\link{genClRcrd}}
setClass('ClRcrd', representation=representation(
  id='integer',
  sids='vector',
  txids='vector',
  typ='character',
  prnt='character',
  seed='character'),
  validity=chckClRcrd)

#' @rdname ClRcrd-class
#' @exportMethod as.character
setMethod('as.character', c('x'='ClRcrd'),
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
#' @rdname ClRcrd-class
#' @exportMethod show
setMethod('show', 'ClRcrd',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname ClRcrd-class
#' @exportMethod print
setMethod('print', 'ClRcrd',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname ClRcrd-class
#' @exportMethod str
setMethod('str', c('object'='ClRcrd'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
#' @rdname ClRcrd-class
#' @exportMethod summary
setMethod('summary', c('object'='ClRcrd'),
          function(object){
            msg <- as.character(x)
            cat(msg)
          })