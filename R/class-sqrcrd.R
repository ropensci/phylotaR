chckSqRcrd <- function(object) {
  # TODO
  TRUE
}

#' @name SqRcrd-class
#' @aliases SqRcrd-method
#' @param x \code{SqRcrd} object
#' @param object \code{SqRcrd} object
#' @title SqRcrd-class
#' @description Sequence record contains sequence data.
#' @details Sequence is stored as raw. Use rawToChar().
#' @slot id Unique ID
#' @slot nm Best-guess sequence name
#' @slot accssn Accession
#' @slot vrsn Accession version
#' @slot url URL
#' @slot gi GI
#' @slot txid Taxonomic ID
#' @slot sq Sequence
#' @slot dfln Definition line
#' @slot ml_typ Molecule type, e.g. DNA
#' @slot rcrd_typ Record type: Whole or feature
#' @slot nncltds Number of nucleotides
#' @slot nambgs Number of ambiguous nucleotides
#' @slot pambgs Proportion of ambiguous nucleotides
#' @slot gcr GC ratio
#' @slot age Number of days between sequence upload and running pipeline 
#' @exportClass SqRcrd
#' @seealso 
#' \code{\link{genSqRcrd}}
setClass('SqRcrd', representation=representation(
  id='character',
  nm='character',
  accssn='character',
  vrsn='character',
  gi='character',
  url='character',
  txid='character',
  orgnsm='character',
  sq='raw',
  dfln='character',
  ml_typ='character',
  rcrd_typ='character',
  nncltds='integer',
  nambgs='integer',
  pambgs='numeric',
  gcr='numeric',
  age='integer'),
  validity=chckSqRcrd)

#' @rdname SqRcrd-class
#' @exportMethod as.character
setMethod('as.character', c('x'='SqRcrd'),
          function(x) {
            paste0('SqRcrd [ID: ', x@id,']')
          })
#' @rdname SqRcrd-class
#' @exportMethod show
setMethod('show', 'SqRcrd',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname SqRcrd-class
#' @exportMethod print
setMethod('print', 'SqRcrd',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname SqRcrd-class
#' @exportMethod str
setMethod('str', c('object'='SqRcrd'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
#' @rdname SqRcrd-class
#' @exportMethod summary
setMethod('summary', c('object'='SqRcrd'),
          function(object){
            msg <- as.character(x)
            cat(msg)
          })