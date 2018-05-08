# Style note:
#  - lots of style guides recommend avoiding S4
#  - those that do recommend, suggest camelcase
#  - for now I'm following 'DataType'
#  - functions that interact are data_type_verb
#' @name RecordCluster-class
#' @aliases RecordCluster-method
#' @param x \code{RecordCluster} object
#' @param object \code{RecordCluster} object
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title RecordCluster-class
#' @description Cluster record contains all information on a cluster.
#' @slot id Cluster ID, integer
#' @slot sids Sequence IDs
#' @slot nsqs Number of sequences
#' @slot txids Source txids for sequences
#' @slot ntx Number of taxa
#' @slot typ Cluster type: direct, subtree or merged
#' @slot seed Seed sequence ID
#' @slot prnt Parent taxonomic ID
#' @exportClass RecordCluster
setClass('RecordCluster', representation = representation(
  id = 'integer',
  sids = 'vector',
  nsqs = 'integer',
  txids = 'vector',
  ntx = 'integer',
  typ = 'character',
  prnt = 'character',
  seed = 'character'))

#' @rdname RecordCluster-class
#' @exportMethod as.character
setMethod('as.character', c('x' = 'RecordCluster'),
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
#' @rdname RecordCluster-class
#' @exportMethod show
setMethod('show', 'RecordCluster',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname RecordCluster-class
#' @exportMethod print
setMethod('print', 'RecordCluster',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname RecordCluster-class
#' @exportMethod str
setMethod('str', c('object' = 'RecordCluster'),
          function(object, max.level = 2L, ...) {
            if (is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level = max.level, ...)
          })
#' @rdname RecordCluster-class
#' @exportMethod summary
setMethod('summary', c('object' = 'RecordCluster'),
          function(object){
            msg <- as.character(object)
            cat(msg)
          })

archive_cluster_check <- function(object) {
  length(object@ids) == length(object@cls)
}

#' @name ArchiveCluster-class
#' @aliases ArchiveCluster-method
#' @param x \code{ArchiveCluster} object
#' @param object \code{ArchiveCluster} object
#' @param i cid(s)
#' @param j Unused
#' @param drop Unused
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title ArchiveCluster-class
#' @description Multiple cluster records.
#' @slot ids Vector of Cluster Record IDs
#' @slot cls List of ArchiveCluster named by ID
#' @exportClass ArchiveCluster
setClass('ArchiveCluster', representation = representation(
  ids = 'vector',
  cls = 'list'),
  validity = archive_cluster_check)

#' @rdname ArchiveCluster-class
#' @exportMethod as.character
setMethod('as.character', c('x' = 'ArchiveCluster'),
          function(x) {
            msg <- 'Archive of cluster record(s)\n'
            msg <- paste0(msg, ' - [', length(x@ids),
                          '] clusters\n')
            msg
          })
#' @rdname ArchiveCluster-class
#' @exportMethod show
setMethod('show', 'ArchiveCluster',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname ArchiveCluster-class
#' @exportMethod print
setMethod('print', 'ArchiveCluster',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname ArchiveCluster-class
#' @exportMethod str
setMethod('str', c('object' = 'ArchiveCluster'),
          function(object, max.level = 2L, ...) {
            if (is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level = max.level, ...)
          })
#' @rdname ArchiveCluster-class
#' @exportMethod summary
setMethod('summary', c('object' = 'ArchiveCluster'),
          function(object){
            msg <- as.character(object)
            cat(msg)
          })

# Accessor methods
#' @rdname ArchiveCluster-class
#' @exportMethod [[
setMethod('[[', c('ArchiveCluster', 'character'),
          function(x, i) {
            pull <- which(x@ids %in% i)
            if (length(pull) == 1) {
              return(x@cls[[pull[1]]])
            }
            stop(paste0('[', i , '] not in records'))
          })
#' @rdname ArchiveCluster-class
#' @exportMethod [
setMethod('[', c('ArchiveCluster', 'character', 'missing', 'missing'),
          function(x, i, j, ..., drop = TRUE) {
            pull <- i %in% x@ids
            if (all(pull)) {
              x <- archive_cluster_generate(x@cls[x@ids %in% i])
              x@ids <- i
              return(x)
            }
            mssng <- paste0(i[!pull], collapse = ', ')
            stop(paste0('[', mssng , '] not in records'))
          })

#' @name RecordSequence-class
#' @aliases RecordSequence-method
#' @param x \code{RecordSequence} object
#' @param object \code{RecordSequence} object
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title RecordSequence-class
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
#' @exportClass RecordSequence
setClass('RecordSequence', representation = representation(
  id = 'character',
  nm = 'character',
  accssn = 'character',
  vrsn = 'character',
  gi = 'character',
  url = 'character',
  txid = 'character',
  orgnsm = 'character',
  sq = 'raw',
  dfln = 'character',
  ml_typ = 'character',
  rcrd_typ = 'character',
  nncltds = 'integer',
  nambgs = 'integer',
  pambgs = 'numeric',
  gcr = 'numeric',
  age = 'integer'))

#' @rdname RecordSequence-class
#' @exportMethod as.character
setMethod('as.character', c('x' = 'RecordSequence'),
          function(x) {
            paste0('RecordSequence [ID: ', x@id,']')
          })
#' @rdname RecordSequence-class
#' @exportMethod show
setMethod('show', 'RecordSequence',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname RecordSequence-class
#' @exportMethod print
setMethod('print', 'RecordSequence',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname RecordSequence-class
#' @exportMethod str
setMethod('str', c('object' = 'RecordSequence'),
          function(object, max.level = 2L, ...) {
            if (is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level = max.level, ...)
          })
#' @rdname RecordSequence-class
#' @exportMethod summary
setMethod('summary', c('object' = 'RecordSequence'),
          function(object){
            msg <- as.character(object)
            cat(msg)
          })

archive_sequence_check <- function(object) {
  length(object@ids) == length(object@sqs) &
    length(object@ids) == length(object@txids)  
}

#' @name ArchiveSequence-class
#' @aliases ArchiveSequence-method
#' @param x \code{ArchiveSequence} object
#' @param object \code{ArchiveSequence} object
#' @param i sid(s)
#' @param j Unused
#' @param drop Unused
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title ArchiveSequence-class
#' @description Multiple sequence records containing sequence data.
#' @details Sequences are stored as raw. Use rawToChar().
#' @slot ids Vector of Sequence Record IDs
#' @slot nncltds Vector of sequence lengths
#' @slot nambgs Vector of number of ambiguous nucleotides
#' @slot txids Vector source txid associated with each sequence
#' @slot sqs List of RecordSequences named by ID
#' @exportClass ArchiveSequence
setClass('ArchiveSequence', representation = representation(
  ids = 'vector',
  nncltds = 'vector',
  nambgs = 'vector',
  txids = 'vector',
  sqs = 'list'),
  validity = archive_sequence_check)

#' @rdname ArchiveSequence-class
#' @exportMethod as.character
setMethod('as.character', c('x' = 'ArchiveSequence'),
          function(x) {
            msg <- 'Multiple RecordSequence(s)\n'
            msg <- paste0(msg, ' - [', length(x@ids),
                          '] sequences\n')
            msg <- paste0(msg, ' - [', length(unique(x@txids)),
                          '] unique txids\n')
            msg <- paste0(msg, ' - [', median(x@nncltds),
                          '] median sequence length\n')
            msg <- paste0(msg, ' - [', median(x@nambgs),
                          '] median ambiguous nucleotides\n')
            msg
          })
#' @rdname ArchiveSequence-class
#' @exportMethod show
setMethod('show', 'ArchiveSequence',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname ArchiveSequence-class
#' @exportMethod print
setMethod('print', 'ArchiveSequence',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname ArchiveSequence-class
#' @exportMethod str
setMethod('str', c('object' = 'ArchiveSequence'),
          function(object, max.level = 2L, ...) {
            if (is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level = max.level, ...)
          })
#' @rdname ArchiveSequence-class
#' @exportMethod summary
setMethod('summary', c('object' = 'ArchiveSequence'),
          function(object){
            msg <- as.character(object)
            cat(msg)
          })

# Accessor methods
#' @rdname ArchiveSequence-class
#' @exportMethod [[
setMethod('[[', c('ArchiveSequence', 'character'),
          function(x, i) {
            pull <- which(x@ids %in% i)
            if (length(pull) == 1) {
              return(x@sqs[[pull[1]]])
            }
            stop(paste0('[', i , '] not in records'))
          })
#' @rdname ArchiveSequence-class
#' @exportMethod [
setMethod('[', c('ArchiveSequence', 'character', 'missing', 'missing'),
          function(x, i, j, ..., drop = TRUE) {
            pull <- i %in% x@ids
            if (all(pull)) {
              return(archive_sequence_generate(x@sqs[x@ids %in% i]))
            }
            mssng <- paste0(i[!pull], collapse = ', ')
            stop(paste0('[', mssng , '] not in records'))
          })

record_taxon_check <- function(object) {
  length(object@lng[['rnks']]) ==
    length(object@lng[['ids']])
}

#' @name RecordTaxon-class
#' @aliases RecordTaxon-method
#' @param x \code{RecordTaxon} object
#' @param object \code{RecordTaxon} object
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title RecordTaxon-class
#' @description Taxonomic dictionary contains a taxonomic
#' tree and NCBI taxonomy data for all taxonomic IDs.
#' @slot id Taxonomic ID
#' @slot scnm Scientific name
#' @slot cmnm Common name
#' @slot rnk Rank
#' @slot lng Lineage
#' @slot prnt Parent
#' @exportClass RecordTaxon
setClass('RecordTaxon', representation = representation(
  id = 'character',
  scnm = 'character',
  cmnm = 'character',
  rnk = 'character',
  lng = 'list',
  prnt = 'character'),
  validity = record_taxon_check)

#' @rdname RecordTaxon-class
#' @exportMethod as.character
setMethod('as.character', c('x' = 'RecordTaxon'),
          function(x) {
            paste0('RecordTaxon [id ', x@id, ' (', x@scnm, ')]\n')
          })
#' @rdname RecordTaxon-class
#' @exportMethod show
setMethod('show', 'RecordTaxon',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname RecordTaxon-class
#' @exportMethod print
setMethod('print', 'RecordTaxon',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname RecordTaxon-class
#' @exportMethod str
setMethod('str', c('object' = 'RecordTaxon'),
          function(object, max.level = 2L, ...) {
            if (is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level = max.level, ...)
          })
#' @rdname RecordTaxon-class
#' @exportMethod summary
setMethod('summary', c('object' = 'RecordTaxon'),
          function(object){
            msg <- as.character(object)
            cat(msg)
          })

dictionary_taxon_check <- function(object) {
  length(object@txids) ==
    length(ls(object@rcrds)) &
    length(object@txids) ==
    length(object@indx)
}

#' @name DictionaryTaxon-class
#' @aliases DictionaryTaxon-method
#' @param x \code{DictionaryTaxon} object
#' @param object \code{DictionaryTaxon} object
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title DictionaryTaxon-class
#' @description Taxonomic dictionary contains a taxonomic
#' tree and NCBI taxonomy data for all taxonomic IDs.
#' @slot txids Taxonomic IDs of taxon records
#' @slot rcrds Environment of records
#' @slot prnt Parent taxonomic ID
#' @slot txtr Taxonomic tree
#' @slot indx Taxonomic ID index for tree IDs
#' @exportClass DictionaryTaxon
setClass('DictionaryTaxon', representation = representation(
  txids = 'vector',
  indx = 'vector',
  rcrds = 'environment',
  txtr = 'TreeMan',
  prnt = 'character'),
  validity = dictionary_taxon_check)

#' @rdname DictionaryTaxon-class
#' @exportMethod as.character
setMethod('as.character', c('x' = 'DictionaryTaxon'),
          function(x) {
            msg <- paste0('Taxonomic dictionary [', length(x@txids),
                          '] rcrds, parent [id ', x@prnt,']\n')
          })
#' @rdname DictionaryTaxon-class
#' @exportMethod show
setMethod('show', 'DictionaryTaxon',
          function(object){
            msg <- as.character(object)
            cat(msg)
          })
#' @rdname DictionaryTaxon-class
#' @exportMethod print
setMethod('print', 'DictionaryTaxon',
          function(x){
            msg <- as.character(x)
            print(msg)
          })
#' @rdname DictionaryTaxon-class
#' @exportMethod str
setMethod('str', c('object' = 'DictionaryTaxon'),
          function(object, max.level = 2L, ...) {
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level = max.level, ...)
          })
#' @rdname DictionaryTaxon-class
#' @exportMethod summary
setMethod('summary', c('object' = 'DictionaryTaxon'),
          function(object){
            msg <- as.character(object)
            cat(msg)
          })


phylota_check <- function(object) {
  length(object@cids) == length(object@cls@cls) &
    length(object@sids) == length(object@sqs@sqs) &
    all(object@txids %in% object@DictionaryTaxon@txids)
}

#' @name Phylota-class
#' @aliases Phylota-method
#' @param x \code{Phylota} object
#' @param object \code{Phylota} object
#' @param i Either sid or cid
#' @param max.level Maximum level of nesting for str()
#' @param ... Further arguments for str()
#' @title Phylota-class
#' @description Phylota table contains sequences and clusters information.
#' @slot cids IDs of all clusters
#' @slot sids IDs of all sequences
#' @slot txids IDs of all taxa
#' @slot sqs All sequence records as ArchiveSequence
#' @slot cls All cluster records as ArchiveCluster
#' @slot DictionaryTaxon Taxonomic dictionary as DictionaryTaxon
#' @slot prnt_id Parent taxonomic ID
#' @slot prnt_nm Parent taxonomic name
#' @exportClass Phylota
setClass('Phylota', representation = representation(
  cids = 'vector',
  txids = 'vector',
  sids = 'vector',
  DictionaryTaxon = 'DictionaryTaxon',
  sqs = 'ArchiveSequence',
  cls = 'ArchiveCluster',
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
            if(is.na(max.level)) {
              stop('max.level must be numeric')
            }
            str@default(object, max.level = max.level, ...)
          })
#' @rdname Phylota-class
#' @exportMethod summary
setMethod('summary', c('object' = 'Phylota'),
          function(object){
            summary_Phylota(object)
          })

# Accessor methods
#' @rdname Phylota-class
#' @exportMethod [[
setMethod('[[', c('Phylota', 'character'),
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
