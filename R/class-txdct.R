chckTxDct <- function(object) {
  # TODO
  TRUE
}

#' @name TxDct-class
#' @aliases TxDct-method
#' @param x \code{TxDct} object
#' @param object \code{TxDct} object
#' @title TxDct-class
#' @description Cluster record contains all information on a cluster.
#' @slot txids Taxonomic IDs of taxon records
#' @slot rcrds List of records
#' @slot prnt Parent taxonomic ID
#' @slot txtr Taxonomic tree
#' @exportClass TxDct
#' @seealso 
#' \code{\link{genTxDct}}
setClass('TxDct', representation=representation(
  txids='vector',
  rcrds='environment',
  txtr='TreeMan',
  prnt='character'),
  validity=chckTxDct)

#' @rdname TxDct-class
#' @exportMethod as.character
setMethod('as.character', c('x'='TxDct'),
          function(x) {
            msg <- paste0('TxDct [', length(x@txids),
                          '] rcrds, parent [id ', x@prnt,']\n')
          })
#' @rdname TxDct-class
#' @exportMethod show
setMethod('show', 'TxDct',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname TxDct-class
#' @exportMethod print
setMethod('print', 'TxDct',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname TxDct-class
#' @exportMethod str
setMethod('str', c('object'='TxDct'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
#' @rdname TxDct-class
#' @exportMethod summary
setMethod('summary', c('object'='TxDct'),
          function(object){
            msg <- as.character(x)
            cat(msg)
          })