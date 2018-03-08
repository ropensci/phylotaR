chckTxRcrd <- function(object) {
  # TODO
  TRUE
}

#' @name TxRcrd-class
#' @aliases TxRcrd-method
#' @param x \code{TxRcrd} object
#' @param object \code{TxRcrd} object
#' @title TxRcrd-class
#' @description Taxonomic dictionary contains a taxonomic
#' tree and NCBI taxonomy data for all taxonomic IDs.
#' @slot id Taxonomic ID
#' @slot scnm Scientific name
#' @slot cmnm Common name
#' @slot rnk Rank
#' @slot lng Lineage
#' @slot prnt Parent
#' @exportClass TxRcrd
#' @seealso 
#' \code{\link{genTxRcrd}}
setClass('TxRcrd', representation=representation(
  id='character',
  scnm='character',
  cmnm='character',
  rnk='character',
  lng='list',
  prnt='character'),
  validity=chckTxRcrd)

#' @rdname TxRcrd-class
#' @exportMethod as.character
setMethod('as.character', c('x'='TxRcrd'),
          function(x) {
            msg <- paste0('TxRcrd [id ', x@id,
                          ' (', x@scnm, ')]\n')
          })
#' @rdname TxRcrd-class
#' @exportMethod show
setMethod('show', 'TxRcrd',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname TxRcrd-class
#' @exportMethod print
setMethod('print', 'TxRcrd',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname TxRcrd-class
#' @exportMethod str
setMethod('str', c('object'='TxRcrd'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
#' @rdname TxRcrd-class
#' @exportMethod summary
setMethod('summary', c('object'='TxRcrd'),
          function(object){
            msg <- as.character(x)
            cat(msg)
          })