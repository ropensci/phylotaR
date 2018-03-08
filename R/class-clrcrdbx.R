chckClRcrdBx <- function(object) {
  # TODO
  TRUE
}

#' @name ClRcrdBx-class
#' @aliases ClRcrdBx-method
#' @param x \code{ClRcrdBx} object
#' @param object \code{ClRcrdBx} object
#' @title ClRcrdBx-class
#' @description Multiple cluster records.
#' @slot ids Vector of Cluster Record IDs
#' @slot cls List of ClRcrds named by ID
#' @exportClass ClRcrdBx
#' @seealso 
#' \code{\link{genClRcrdBx}}
setClass('ClRcrdBx', representation=representation(
  ids='vector',
  cls='list'),
  validity=chckClRcrdBx)

#' @rdname ClRcrdBx-class
#' @exportMethod as.character
setMethod('as.character', c('x'='ClRcrdBx'),
          function(x) {
            msg <- 'Multiple ClRcrd(s)\n'
            msg <- paste0(msg, ' - [', length(x@ids),
                          '] clusters\n')
            msg
          })
#' @rdname ClRcrdBx-class
#' @exportMethod show
setMethod('show', 'ClRcrdBx',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname ClRcrdBx-class
#' @exportMethod print
setMethod('print', 'ClRcrdBx',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname ClRcrdBx-class
#' @exportMethod str
setMethod('str', c('object'='ClRcrdBx'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
#' @rdname ClRcrdBx-class
#' @exportMethod summary
setMethod('summary', c('object'='ClRcrdBx'),
          function(object){
            msg <- as.character(x)
            cat(msg)
          })

# Accessor methods
setMethod('[[', c('ClRcrdBx', 'character'),
          function(x, i) {
            pull <- which(x@ids %in% i)
            if(length(pull) == 1) {
              return(x@cls[[pull[1]]])
            }
            stop(paste0('[', i , '] not in records'))
          })
setMethod('[', c('ClRcrdBx', 'character', 'missing', 'missing'),
          function(x, i, j, ..., drop=TRUE) {
            pull <- i %in% x@ids
            if(all(pull)) {
              x <- genClRcrdBx(x@cls[x@ids %in% i])
              x@ids <- i
              return(x)
            }
            mssng <- paste0(i[!pull], collapse=', ')
            stop(paste0('[', mssng , '] not in records'))
          })
