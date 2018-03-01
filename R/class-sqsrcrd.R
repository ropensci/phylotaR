chckSqsRcrd <- function(object) {
  # TODO
  TRUE
}

#' @name SqsRcrd-class
#' @aliases SqsRcrd-method
#' @param x \code{SqsRcrd} object
#' @param object \code{SqsRcrd} object
#' @title SqsRcrd-class
#' @description Multiple sequence records containing sequence data.
#' @details Sequences are stored as raw. Use rawToChar().
#' @slot ids Vector of SqsRcrd IDs
#' @slot nncltds Vector of sequence lengths
#' @slot nambgs Vector of number of ambiguous nucleotides
#' @slot txids Vector source txid associated with each sequence
#' @slot sqs List of SqsRcrds named by ID
#' @exportClass SqsRcrd
#' @seealso 
#' \code{\link{genSqsRcrd}}
setClass('SqsRcrd', representation=representation(
  ids='vector',
  nncltds='vector',
  nambgs='vector',
  txids='vector',
  sqs='list'),
  validity=chckSqsRcrd)

#' @rdname SqsRcrd-class
#' @exportMethod as.character
setMethod('as.character', c('x'='SqsRcrd'),
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
#' @rdname SqsRcrd-class
#' @exportMethod show
setMethod('show', 'SqsRcrd',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname SqsRcrd-class
#' @exportMethod print
setMethod('print', 'SqsRcrd',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname SqsRcrd-class
#' @exportMethod str
setMethod('str', c('object'='SqsRcrd'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
#' @rdname SqsRcrd-class
#' @exportMethod summary
setMethod('summary', c('object'='SqsRcrd'),
          function(object){
            msg <- as.character(x)
            cat(msg)
          })

# Accessor methods
setMethod('[[', c('SqsRcrd', 'character'),
          function(x, i) {
            pull <- which(x@ids %in% i)
            if(length(pull) == 1) {
              return(x@sqs[[pull[1]]])
            }
            stop(paste0('[', i , '] not in records'))
          })
setMethod('[', c('SqsRcrd', 'character', 'missing', 'missing'),
          function(x, i, j, ..., drop=TRUE) {
            pull <- i %in% x@ids
            if(all(pull)) {
              return(genSqsRcrd(x@sqs[x@ids %in% i]))
            }
            mssng <- paste0(i[!pull], collapse=', ')
            stop(paste0('[', mssng , '] not in records'))
          })
