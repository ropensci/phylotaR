
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
#' @slot cids IDs of all clusters
#' @slot sids IDs of all sequences
#' @slot txids IDs of all taxa
#' @slot sqs All sequence records as SqRcrdBx
#' @slot cls All cluster records as ClRcrdBx
#' @slot txdct Taxonomic dictionary
#' @exportClass PhyLoTa
#' @seealso 
#' \code{\link{genPhyLoTa}}
setClass('PhyLoTa', representation=representation(
  cids='vector',
  txids='vector',
  sids='vector',
  txdct='list',
  sqs='SqRcrdBx',
  cls='ClRcrdBx'),
  validity=chckPhyLoTa)

#' @rdname PhyLoTa-class
#' @exportMethod as.character
setMethod('as.character', c('x'='PhyLoTa'),
          function(x) {
            msg <- 'PhyLoTa Table\n'
            msg <- paste0(msg, '- [', length(x@cids),
                          '] clusters\n')
            msg <- paste0(msg, '- [', length(x@sids),
                          '] sequences\n')
            msg <- paste0(msg, '- [', length(x@txids),
                          '] source taxa\n')
            msg
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
