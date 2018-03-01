chckSqRcrdBx <- function(object) {
  # TODO
  TRUE
}

#' @name SqRcrdBx-class
#' @aliases SqRcrdBx-method
#' @param x \code{SqRcrdBx} object
#' @param object \code{SqRcrdBx} object
#' @title SqRcrdBx-class
#' @description Multiple sequence records containing sequence data.
#' @details Sequences are stored as raw. Use rawToChar().
#' @slot ids Vector of Sequence Record IDs
#' @slot nncltds Vector of sequence lengths
#' @slot nambgs Vector of number of ambiguous nucleotides
#' @slot txids Vector source txid associated with each sequence
#' @slot sqs List of SqRcrds named by ID
#' @exportClass SqRcrdBx
#' @seealso 
#' \code{\link{genSqRcrdBx}}
setClass('SqRcrdBx', representation=representation(
  ids='vector',
  nncltds='vector',
  nambgs='vector',
  txids='vector',
  sqs='list'),
  validity=chckSqRcrdBx)

#' @rdname SqRcrdBx-class
#' @exportMethod as.character
setMethod('as.character', c('x'='SqRcrdBx'),
          function(x) {
            msg <- 'Multiple SqRcrd(s)\n'
            msg <- paste0(msg, ' - [', length(x@ids),
                          '] sequences\n')
            msg <- paste0(msg, ' - [', length(unique(x@txids)),
                          '] unique txids\n')
            msg <- paste0(msg, ' - [', median(x@nncltds),
                          '] median sequence length\n')
            msg <- paste0(msg, ' - [', median(x@nambgs),
                          '] median ambiguous nucleotides\n')
          })
#' @rdname SqRcrdBx-class
#' @exportMethod show
setMethod('show', 'SqRcrdBx',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname SqRcrdBx-class
#' @exportMethod print
setMethod('print', 'SqRcrdBx',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname SqRcrdBx-class
#' @exportMethod str
setMethod('str', c('object'='SqRcrdBx'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
#' @rdname SqRcrdBx-class
#' @exportMethod summary
setMethod('summary', c('object'='SqRcrdBx'),
          function(object){
            msg <- as.character(x)
            cat(msg)
          })

# Accessor methods
setMethod('[[', c('SqRcrdBx', 'character'),
          function(x, i) {
            pull <- which(x@ids %in% i)
            if(length(pull) == 1) {
              return(x@sqs[[pull[1]]])
            }
            stop(paste0('[', i , '] not in records'))
          })
setMethod('[', c('SqRcrdBx', 'character', 'missing', 'missing'),
          function(x, i, j, ..., drop=TRUE) {
            pull <- i %in% x@ids
            if(all(pull)) {
              return(genSqRcrdBx(x@sqs[x@ids %in% i]))
            }
            mssng <- paste0(i[!pull], collapse=', ')
            stop(paste0('[', mssng , '] not in records'))
          })
