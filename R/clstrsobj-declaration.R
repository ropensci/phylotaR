
chckClstrsObj <- function(object) {
  # TODO
  TRUE
}

#' @name ClstrsObj-class
#' @aliases ClstrsObj-method
#' @param x \code{ClstrsObj} object
#' @param object \code{ClstrsObj} object
#' @title ClstrsObj-class
#' @description Clusters Object contains sequences and clusters information.
#' @slot clstr_ids IDs of all clusters in object
#' @slot txids IDs of all taxa in object
#' @slot sqids IDs of all sequences in object
#' @slot txdct Taxonomic dictionary
#' @slot sqs All sequences
#' @slot clstrs All clstrs
#' @exportClass ClstrsObj
#' @seealso 
#' \code{\link{genClstrsObj}}
setClass('ClstrsObj', representation=representation(
  clstr_ids='vector',
  txids='vector',
  sqids='vector',
  txdct='list',
  sqs='list',
  clstrs='list'),
  validity=chckClstrsObj)

#' @rdname ClstrsObj-class
#' @exportMethod as.character
setMethod('as.character', c('x'='ClstrsObj'),
          function(x) {
            paste0('ClstrsObj of [', length(x@clstrs),'] clusters')
          })
#' @rdname ClstrsObj-class
#' @exportMethod show
setMethod('show', 'ClstrsObj',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname ClstrsObj-class
#' @exportMethod print
setMethod('print', 'ClstrsObj',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname ClstrsObj-class
#' @exportMethod str
setMethod('str', c('object'='ClstrsObj'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
#' @rdname ClstrsObj-class
#' @exportMethod summary
setMethod('summary', c('object'='ClstrsObj'),
          function(object){
            msg <- as.character(x)
            cat(msg)
          })
# setMethod('plot', c('object'='ClstrsObj'),
#           function(object){
#             msg <- as.character(x)
#             cat(msg)
#           })
