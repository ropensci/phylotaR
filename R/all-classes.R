#' @name ClRcrd-class
#' @aliases ClRcrd-method
#' @param x \code{ClRcrd} object
#' @param object \code{ClRcrd} object
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title ClRcrd-class
#' @description Cluster record contains all information on a cluster.
#' @slot id Cluster ID, integer
#' @slot sids Sequence IDs
#' @slot nsqs Number of sequences
#' @slot txids Source txids for sequences
#' @slot ntx Number of taxa
#' @slot typ Cluster type: direct, subtree or merged
#' @slot seed Seed sequence ID
#' @slot prnt Parent taxonomic ID
#' @exportClass ClRcrd
setClass('ClRcrd', representation=representation(
  id='integer',
  sids='vector',
  nsqs='integer',
  txids='vector',
  ntx='integer',
  typ='character',
  prnt='character',
  seed='character'))

#' @rdname ClRcrd-class
#' @exportMethod as.character
setMethod('as.character', c('x'='ClRcrd'),
          function(x) {
            msg <- paste0('Cluster Record [id ', x@id,']\n')
            msg <- paste0(msg, ' - [', x@typ,
                          '] type\n')
            msg <- paste0(msg, ' - [', x@seed,
                          '] seed sequence\n')
            msg <- paste0(msg, ' - [', x@nsqs,
                          '] sequences\n')
            msg <- paste0(msg, ' - [', x@ntx,
                          '] taxa\n')
          })
#' @rdname ClRcrd-class
#' @exportMethod show
setMethod('show', 'ClRcrd',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname ClRcrd-class
#' @exportMethod print
setMethod('print', 'ClRcrd',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname ClRcrd-class
#' @exportMethod str
setMethod('str', c('object'='ClRcrd'),
          function(object, max.level=2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level=max.level, ...)
          })
#' @rdname ClRcrd-class
#' @exportMethod summary
setMethod('summary', c('object'='ClRcrd'),
          function(object){
            msg <- as.character(object)
            cat(msg)
          })

chckClRcrdBx <- function(object) {
  length(object@ids) == length(object@cls)
}

#' @name ClRcrdBx-class
#' @aliases ClRcrdBx-method
#' @param x \code{ClRcrdBx} object
#' @param object \code{ClRcrdBx} object
#' @param i cid(s)
#' @param j Unused
#' @param drop Unused
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title ClRcrdBx-class
#' @description Multiple cluster records.
#' @slot ids Vector of Cluster Record IDs
#' @slot cls List of ClRcrds named by ID
#' @exportClass ClRcrdBx
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
            msg <- as.character(object)
            cat(msg)
          })

# Accessor methods
#' @rdname ClRcrdBx-class
#' @exportMethod [[
setMethod('[[', c('ClRcrdBx', 'character'),
          function(x, i) {
            pull <- which(x@ids %in% i)
            if(length(pull) == 1) {
              return(x@cls[[pull[1]]])
            }
            stop(paste0('[', i , '] not in records'))
          })
#' @rdname ClRcrdBx-class
#' @exportMethod [
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

#' @name SqRcrd-class
#' @aliases SqRcrd-method
#' @param x \code{SqRcrd} object
#' @param object \code{SqRcrd} object
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title SqRcrd-class
#' @description Sequence record contains sequence data.
#' @details Sequence is stored as raw. Use rawToChar().
#' @slot id Unique ID
#' @slot nm Best-guess sequence name
#' @slot accssn Accession
#' @slot vrsn Accession version
#' @slot url URL
#' @slot gi GI
#' @slot txid Taxonomic ID of source taxon
#' @slot orgnsm Scientific name of source taxon
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
  age='integer'))

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
            msg <- as.character(object)
            cat(msg)
          })

chckSqRcrdBx <- function(object) {
  length(object@ids) == length(object@sqs) &
    length(object@ids) == length(object@txids)  
}

#' @name SqRcrdBx-class
#' @aliases SqRcrdBx-method
#' @param x \code{SqRcrdBx} object
#' @param object \code{SqRcrdBx} object
#' @param i sid(s)
#' @param j Unused
#' @param drop Unused
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title SqRcrdBx-class
#' @description Multiple sequence records containing sequence data.
#' @details Sequences are stored as raw. Use rawToChar().
#' @slot ids Vector of Sequence Record IDs
#' @slot nncltds Vector of sequence lengths
#' @slot nambgs Vector of number of ambiguous nucleotides
#' @slot txids Vector source txid associated with each sequence
#' @slot sqs List of SqRcrds named by ID
#' @exportClass SqRcrdBx
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
            msg <- as.character(object)
            cat(msg)
          })

# Accessor methods
#' @rdname SqRcrdBx-class
#' @exportMethod [[
setMethod('[[', c('SqRcrdBx', 'character'),
          function(x, i) {
            pull <- which(x@ids %in% i)
            if(length(pull) == 1) {
              return(x@sqs[[pull[1]]])
            }
            stop(paste0('[', i , '] not in records'))
          })
#' @rdname SqRcrdBx-class
#' @exportMethod [
setMethod('[', c('SqRcrdBx', 'character', 'missing', 'missing'),
          function(x, i, j, ..., drop=TRUE) {
            pull <- i %in% x@ids
            if(all(pull)) {
              return(genSqRcrdBx(x@sqs[x@ids %in% i]))
            }
            mssng <- paste0(i[!pull], collapse=', ')
            stop(paste0('[', mssng , '] not in records'))
          })

chckTxRcrd <- function(object) {
  length(object@lng[['rnks']]) ==
    length(object@lng[['ids']])
}

#' @name TxRcrd-class
#' @aliases TxRcrd-method
#' @param x \code{TxRcrd} object
#' @param object \code{TxRcrd} object
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
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
            msg <- as.character(object)
            cat(msg)
          })

chckTxDct <- function(object) {
  length(object@txids) ==
    length(ls(object@rcrds)) &
    length(object@txids) ==
    length(object@indx)
}

#' @name TxDct-class
#' @aliases TxDct-method
#' @param x \code{TxDct} object
#' @param object \code{TxDct} object
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title TxDct-class
#' @description Taxonomic dictionary contains a taxonomic
#' tree and NCBI taxonomy data for all taxonomic IDs.
#' @slot txids Taxonomic IDs of taxon records
#' @slot rcrds Environment of records
#' @slot prnt Parent taxonomic ID
#' @slot txtr Taxonomic tree
#' @slot indx Taxonomic ID index for tree IDs
#' @exportClass TxDct
setClass('TxDct', representation=representation(
  txids='vector',
  indx='vector',
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
            msg <- as.character(object)
            cat(msg)
          })


chckPhyLoTa <- function(object) {
  length(object@cids) == length(object@cls@cls) &
    length(object@sids) == length(object@sqs@sqs) &
    all(object@txids %in% object@txdct@txids)
}

#' @name PhyLoTa-class
#' @aliases PhyLoTa-method
#' @param x \code{PhyLoTa} object
#' @param object \code{PhyLoTa} object
#' @param y Not used in plot for this class
#' @param i Either sid or cid
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title PhyLoTa-class
#' @description PhyLoTa table contains sequences and clusters information.
#' @slot cids IDs of all clusters
#' @slot sids IDs of all sequences
#' @slot txids IDs of all taxa
#' @slot sqs All sequence records as SqRcrdBx
#' @slot cls All cluster records as ClRcrdBx
#' @slot txdct Taxonomic dictionary as TxDct
#' @slot prnt_id Parent taxonomic ID
#' @slot prnt_nm Parent taxonomic name
#' @exportClass PhyLoTa
setClass('PhyLoTa', representation=representation(
  cids='vector',
  txids='vector',
  sids='vector',
  txdct='TxDct',
  sqs='SqRcrdBx',
  cls='ClRcrdBx',
  prnt_id='character',
  prnt_nm='character'),
  validity=chckPhyLoTa)

#' @rdname PhyLoTa-class
#' @exportMethod as.character
setMethod('as.character', c('x'='PhyLoTa'),
          function(x) {
            msg <- paste0('PhyLoTa Table (',
                          x@prnt_nm, ')\n')
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
            summary_phylota(object)
          })
#' @rdname PhyLoTa-class
#' @exportMethod plot
setMethod('plot', c('x'='PhyLoTa', y='missing'),
          function(x, y, ...){
            plot_phylota(x)
          })

# Accessor methods
#' @rdname PhyLoTa-class
#' @exportMethod [[
setMethod('[[', c('PhyLoTa', 'character'),
          function(x, i) {
            pull <- which(x@cids %in% i)
            if(length(pull) == 1) {
              return(x@cls@cls[[pull[1]]])
            }
            pull <- which(x@sids %in% i)
            if(length(pull) == 1) {
              return(x@sqs@sqs[[pull[1]]])
            }
            stop(paste0('[', i , '] not in table'))
          })
