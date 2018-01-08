
chckClstrsObj <- function(object) {
  # TODO
  TRUE
}

#' @name ClstrsObj-class
#' @title ClstrsObj-class
#' @description Cluster holder
#' @exportClass ClstrsObj
setClass('ClstrsObj', representation=representation(
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
