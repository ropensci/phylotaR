# Style note ----
#  - lots of style guides recommend avoiding S4
#  - those that do recommend, suggest CamelCase
#  - for now I'm following 'TypeData'
#  - functions that interact are typedata_verb
# ClstrRec ----
#' @name ClstrRec-class
#' @aliases ClstrRec-method
#' @param x \code{ClstrRec} object
#' @param object \code{ClstrRec} object
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title Cluster record
#' @description Cluster record contains all information on a cluster.
#' @slot id Cluster ID, integer
#' @slot sids Sequence IDs
#' @slot nsqs Number of sequences
#' @slot txids Source txids for sequences
#' @slot ntx Number of taxa
#' @slot typ Cluster type: direct, subtree or merged
#' @slot seed Seed sequence ID
#' @slot prnt Parent taxonomic ID
#' @family run-public
#' @exportClass ClstrRec
#' @example examples/clstrrec-class.R
setClass('ClstrRec', representation = representation(
  id = 'integer',
  sids = 'vector',
  nsqs = 'integer',
  txids = 'vector',
  ntx = 'integer',
  typ = 'character',
  prnt = 'character',
  seed = 'character'))

#' @rdname ClstrRec-class
#' @exportMethod as.character
setMethod('as.character', c('x' = 'ClstrRec'),
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
#' @rdname ClstrRec-class
#' @exportMethod show
setMethod('show', 'ClstrRec',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname ClstrRec-class
#' @exportMethod print
setMethod('print', 'ClstrRec',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname ClstrRec-class
#' @exportMethod str
setMethod('str', c('object' = 'ClstrRec'),
          function(object, max.level = 2L, ...) {
            if (is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level = max.level, ...)
          })
#' @rdname ClstrRec-class
#' @exportMethod summary
setMethod('summary', c('object' = 'ClstrRec'),
          function(object){
            msg <- as.character(object)
            cat(msg)
          })

# ClstrArc ----
clstrarc_check <- function(object) {
  length(object@ids) == length(object@clstrs)
}
#' @name ClstrArc-class
#' @aliases ClstrArc-method
#' @param x \code{ClstrArc} object
#' @param object \code{ClstrArc} object
#' @param i cid(s)
#' @param j Unused
#' @param drop Unused
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title Cluster record archive
#' @description Multiple cluster records.
#' @slot ids Vector of cluster record IDs
#' @slot clstrs List of ClstrArc named by ID
#' @family run-public
#' @exportClass ClstrArc
#' @example examples/clstrarc-class.R
setClass('ClstrArc', representation = representation(
  ids = 'vector',
  clstrs = 'list'),
  validity = clstrarc_check)

#' @rdname ClstrArc-class
#' @exportMethod as.character
setMethod('as.character', c('x' = 'ClstrArc'),
          function(x) {
            msg <- 'Archive of cluster record(s)\n'
            msg <- paste0(msg, ' - [', length(x@ids),
                          '] clusters\n')
            msg
          })
#' @rdname ClstrArc-class
#' @exportMethod show
setMethod('show', 'ClstrArc',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname ClstrArc-class
#' @exportMethod print
setMethod('print', 'ClstrArc',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname ClstrArc-class
#' @exportMethod str
setMethod('str', c('object' = 'ClstrArc'),
          function(object, max.level = 2L, ...) {
            if (is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level = max.level, ...)
          })
#' @rdname ClstrArc-class
#' @exportMethod summary
setMethod('summary', c('object' = 'ClstrArc'),
          function(object){
            msg <- as.character(object)
            cat(msg)
          })

# Accessor methods
#' @rdname ClstrArc-class
#' @exportMethod [[
setMethod('[[', c('ClstrArc', 'character'),
          function(x, i) {
            pull <- which(x@ids %in% i)
            if (length(pull) == 1) {
              return(x@clstrs[[pull[1]]])
            }
            stop(paste0('[', i , '] not in records'))
          })
#' @rdname ClstrArc-class
#' @exportMethod [
setMethod('[', c('ClstrArc', 'character', 'missing', 'missing'),
          function(x, i, j, ..., drop = TRUE) {
            pull <- i %in% x@ids
            if (all(pull)) {
              clstrrecs <- x@clstrs[x@ids %in% i]
              x <- new('ClstrArc', ids = i, clstrs = clstrrecs)
              return(x)
            }
            mssng <- paste0(i[!pull], collapse = ', ')
            stop(paste0('[', mssng , '] not in records'))
          })

# SeqRec ----
#' @name SeqRec-class
#' @family run-public
#' @aliases SeqRec-method
#' @param x \code{SeqRec} object
#' @param object \code{SeqRec} object
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title Sequence record
#' @description Sequence record contains sequence data.
#' @details Sequence is stored as raw. Use rawToChar().
#' @slot id Unique ID
#' @slot nm Best-guess sequence name
#' @slot accssn Accession
#' @slot vrsn Accession version
#' @slot url URL
#' @slot txid Taxonomic ID of source taxon
#' @slot orgnsm Scientific name of source taxon
#' @slot sq Sequence
#' @slot dfln Definition line
#' @slot ml_typ Molecule type, e.g. DNA
#' @slot rec_typ Record type: Whole or feature
#' @slot nncltds Number of nucleotides
#' @slot nambgs Number of ambiguous nucleotides
#' @slot pambgs Proportion of ambiguous nucleotides
#' @slot gcr GC ratio
#' @slot age Number of days between sequence upload and running pipeline 
#' @exportClass SeqRec
#' @example examples/seqrec-class.R
setClass('SeqRec', representation = representation(
  id = 'character',
  nm = 'character',
  accssn = 'character',
  vrsn = 'character',
  url = 'character',
  txid = 'character',
  orgnsm = 'character',
  sq = 'raw',
  dfln = 'character',
  ml_typ = 'character',
  rec_typ = 'character',
  nncltds = 'integer',
  nambgs = 'integer',
  pambgs = 'numeric',
  gcr = 'numeric',
  age = 'integer'))

#' @rdname SeqRec-class
#' @exportMethod as.character
setMethod('as.character', c('x' = 'SeqRec'),
          function(x) {
            paste0('SeqRec [ID: ', x@id,']')
          })
#' @rdname SeqRec-class
#' @exportMethod show
setMethod('show', 'SeqRec',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname SeqRec-class
#' @exportMethod print
setMethod('print', 'SeqRec',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname SeqRec-class
#' @exportMethod str
setMethod('str', c('object' = 'SeqRec'),
          function(object, max.level = 2L, ...) {
            if (is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level = max.level, ...)
          })
#' @rdname SeqRec-class
#' @exportMethod summary
setMethod('summary', c('object' = 'SeqRec'),
          function(object){
            msg <- as.character(object)
            cat(msg)
          })

# SeqArc ----
seqarc_check <- function(object) {
  length(object@ids) == length(object@sqs) &
    length(object@ids) == length(object@txids)  
}
#' @name SeqArc-class
#' @family run-public
#' @aliases SeqArc-method
#' @param x \code{SeqArc} object
#' @param object \code{SeqArc} object
#' @param i sid(s)
#' @param j Unused
#' @param drop Unused
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title Sequence record archive
#' @description Multiple sequence records containing sequence data.
#' @details Sequences are stored as raw. Use rawToChar().
#' @slot ids Vector of Sequence Record IDs
#' @slot nncltds Vector of sequence lengths
#' @slot nambgs Vector of number of ambiguous nucleotides
#' @slot txids Vector source txid associated with each sequence
#' @slot sqs List of SeqRecs named by ID
#' @exportClass SeqArc
#' @example examples/seqarc-class.R
setClass('SeqArc', representation = representation(
  ids = 'vector',
  nncltds = 'vector',
  nambgs = 'vector',
  txids = 'vector',
  sqs = 'list'),
  validity = seqarc_check)

#' @rdname SeqArc-class
#' @exportMethod as.character
setMethod('as.character', c('x' = 'SeqArc'),
          function(x) {
            msg <- 'Archive of sequence record(s)\n'
            msg <- paste0(msg, ' - [', length(x@ids), '] sequences\n')
            msg <- paste0(msg, ' - [', length(unique(x@txids)),
                          '] unique txids\n')
            msg <- paste0(msg, ' - [', median(x@nncltds),
                          '] median sequence length\n')
            msg <- paste0(msg, ' - [', median(x@nambgs),
                          '] median ambiguous nucleotides\n')
            msg
          })
#' @rdname SeqArc-class
#' @exportMethod show
setMethod('show', 'SeqArc',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname SeqArc-class
#' @exportMethod print
setMethod('print', 'SeqArc',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname SeqArc-class
#' @exportMethod str
setMethod('str', c('object' = 'SeqArc'),
          function(object, max.level = 2L, ...) {
            if (is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level = max.level, ...)
          })
#' @rdname SeqArc-class
#' @exportMethod summary
setMethod('summary', c('object' = 'SeqArc'),
          function(object){
            msg <- as.character(object)
            cat(msg)
          })

# Accessor methods
#' @rdname SeqArc-class
#' @exportMethod [[
setMethod('[[', c('SeqArc', 'character'),
          function(x, i) {
            pull <- which(x@ids %in% i)
            if (length(pull) == 1) {
              return(x@sqs[[pull[1]]])
            }
            stop(paste0('[', i , '] not in records'))
          })
#' @rdname SeqArc-class
#' @exportMethod [
setMethod('[', c('SeqArc', 'character', 'missing', 'missing'),
          function(x, i, j, ..., drop = TRUE) {
            pull <- i %in% x@ids
            if (all(pull)) {
              return(seqarc_gen(x@sqs[x@ids %in% i]))
            }
            mssng <- paste0(i[!pull], collapse = ', ')
            stop(paste0('[', mssng , '] not in records'))
          })

# TaxRec ----
taxrec_check <- function(object) {
  length(object@lng[['rnks']]) ==
    length(object@lng[['ids']])
}

#' @name TaxRec-class
#' @family run-public
#' @aliases TaxRec-method
#' @param x \code{TaxRec} object
#' @param object \code{TaxRec} object
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title Taxonomic record
#' @description Taxonomic dictionary contains a taxonomic
#' tree and NCBI taxonomy data for all taxonomic IDs.
#' @slot id Taxonomic ID
#' @slot scnm Scientific name
#' @slot cmnm Common name
#' @slot rnk Rank
#' @slot lng Lineage
#' @slot prnt Parent
#' @exportClass TaxRec
#' @example examples/taxrec-class.R
setClass('TaxRec', representation = representation(
  id = 'character',
  scnm = 'character',
  cmnm = 'character',
  rnk = 'character',
  lng = 'list',
  prnt = 'character'),
  validity = taxrec_check)

#' @rdname TaxRec-class
#' @exportMethod as.character
setMethod('as.character', c('x' = 'TaxRec'),
          function(x) {
            paste0('TaxRec [id ', x@id, ' (', x@scnm, ')]\n')
          })
#' @rdname TaxRec-class
#' @exportMethod show
setMethod('show', 'TaxRec',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname TaxRec-class
#' @exportMethod print
setMethod('print', 'TaxRec',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname TaxRec-class
#' @exportMethod str
setMethod('str', c('object' = 'TaxRec'),
          function(object, max.level = 2L, ...) {
            if (is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level = max.level, ...)
          })
#' @rdname TaxRec-class
#' @exportMethod summary
setMethod('summary', c('object' = 'TaxRec'),
          function(object){
            msg <- as.character(object)
            cat(msg)
          })

# TaxDict ----
taxdict_check <- function(object) {
  length(object@txids) ==
    length(ls(object@recs)) &
    length(object@txids)
}

#' @name TaxDict-class
#' @family run-public
#' @aliases TaxDict-method
#' @param x \code{TaxDict} object
#' @param object \code{TaxDict} object
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title Taxonomic record dictionary
#' @description Taxonomic dictionary contains a taxonomic
#' tree and NCBI taxonomy data for all taxonomic IDs.
#' @slot txids Taxonomic IDs of taxon records
#' @slot recs Environment of records
#' @slot prnt Parent taxonomic ID
#' @slot txtr Taxonomic tree
#' @exportClass TaxDict
#' @example examples/taxdict-class.R
setClass('TaxDict', representation = representation(
  txids = 'vector',
  recs = 'environment',
  txtr = 'TreeMan',
  prnt = 'character'),
  validity = taxdict_check)

#' @rdname TaxDict-class
#' @exportMethod as.character
setMethod('as.character', c('x' = 'TaxDict'),
          function(x) {
            paste0('Taxonomic dictionary [', length(x@txids),
                   '] recs, parent [id ', x@prnt,']\n')
          })
#' @rdname TaxDict-class
#' @exportMethod show
setMethod('show', 'TaxDict',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname TaxDict-class
#' @exportMethod print
setMethod('print', 'TaxDict',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname TaxDict-class
#' @exportMethod str
setMethod('str', c('object' = 'TaxDict'),
          function(object, max.level = 2L, ...) {
            if (is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level = max.level, ...)
          })
#' @rdname TaxDict-class
#' @exportMethod summary
setMethod('summary', c('object' = 'TaxDict'),
          function(object){
            msg <- as.character(object)
            cat(msg)
          })

# Phylota ----
phylota_check <- function(object) {
  length(object@cids) == length(object@clstrs@clstrs) &
    length(object@sids) == length(object@sqs@sqs) &
    all(object@txids %in% object@txdct@txids)
}
#' @name Phylota-class
#' @family run-public
#' @aliases Phylota-method
#' @param x \code{Phylota} object
#' @param object \code{Phylota} object
#' @param i Either sid or cid
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title Phylota object
#' @description Phylota table contains all sequence, cluster and taxonomic
#' information from a phylotaR pipeline run.
#' @slot cids IDs of all clusters
#' @slot sids IDs of all sequences
#' @slot txids IDs of all taxa
#' @slot sqs All sequence records as SeqArc
#' @slot clstrs All cluster records as ClstrArc
#' @slot txdct Taxonomic dictionary as TaxDict
#' @slot prnt_id Parent taxonomic ID
#' @slot prnt_nm Parent taxonomic name
#' @exportClass Phylota
#' @example examples/phylota-class.R
setClass('Phylota', representation = representation(
  cids = 'vector',
  txids = 'vector',
  sids = 'vector',
  txdct = 'TaxDict',
  sqs = 'SeqArc',
  clstrs = 'ClstrArc',
  prnt_id = 'character',
  prnt_nm = 'character'),
  validity = phylota_check)

#' @rdname Phylota-class
#' @exportMethod as.character
setMethod('as.character', c('x' = 'Phylota'),
          function(x) {
            msg <- paste0('Phylota Table (',
                          x@prnt_nm, ')\n')
            msg <- paste0(msg, '- [', length(x@cids),
                          '] clusters\n')
            msg <- paste0(msg, '- [', length(x@sids),
                          '] sequences\n')
            msg <- paste0(msg, '- [', length(x@txids),
                          '] source taxa\n')
            msg
          })
#' @rdname Phylota-class
#' @exportMethod show
setMethod('show', 'Phylota',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname Phylota-class
#' @exportMethod print
setMethod('print', 'Phylota',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname Phylota-class
#' @exportMethod str
setMethod('str', c('object' = 'Phylota'),
          function(object, max.level = 2L, ...) {
            if (is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level = max.level, ...)
          })
#' @rdname Phylota-class
#' @exportMethod summary
setMethod('summary', c('object' = 'Phylota'),
          function(object){
            summary_phylota(object)
          })

# Accessor methods
#' @rdname Phylota-class
#' @exportMethod [[
setMethod('[[', c('Phylota', 'character'),
          function(x, i) {
            pull <- which(x@cids %in% i)
            if (length(pull) == 1) {
              return(x@clstrs@clstrs[[pull[1]]])
            }
            pull <- which(x@sids %in% i)
            if (length(pull) == 1) {
              return(x@sqs@sqs[[pull[1]]])
            }
            stop(paste0('[', i , '] not in table'))
          })
