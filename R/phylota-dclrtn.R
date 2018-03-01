
chckPhyLoTa <- function(object) {
  # TODO
  TRUE
}

#' @name PhyLoTa-class
#' @aliases PhyLoTa-method
#' @param x \code{PhyLoTa} object
#' @param object \code{PhyLoTa} object
#' @title PhyLoTa-class
#' @description PhyLoTa table contains sequences and clusters information.
#' @slot clstr_ids IDs of all clusters
#' @slot txids IDs of all taxa
#' @slot sids IDs of all sequences
#' @slot txdct Taxonomic dictionary
#' @slot sqs All sequence records
#' @slot clstrs All cluster records
#' @exportClass PhyLoTa
#' @seealso 
#' \code{\link{genPhyLoTa}}
setClass('PhyLoTa', representation=representation(
  clstr_ids='vector',
  txids='vector',
  sids='vector',
  txdct='list',
  sqs='SqsRcrd',
  clstrs='list'),
  validity=chckPhyLoTa)

#' @rdname PhyLoTa-class
#' @exportMethod as.character
setMethod('as.character', c('x'='PhyLoTa'),
          function(x) {
            paste0('PhyLoTa Table of [', length(x@clstrs),'] clusters')
          })
#' @rdname PhyLoTa-class
#' @exportMethod show
setMethod('show', 'PhyLoTa',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname PhyLoTa-class
#' @exportMethod print
setMethod('print', 'PhyLoTa',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname PhyLoTa-class
#' @exportMethod str
setMethod('str', c('object'='PhyLoTa'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
#' @rdname PhyLoTa-class
#' @exportMethod summary
setMethod('summary', c('object'='PhyLoTa'),
          function(object){
            msg <- as.character(x)
            cat(msg)
          })
# setMethod('plot', c('object'='PhyLoTa'),
#           function(object){
#             msg <- as.character(x)
#             cat(msg)
#           })
