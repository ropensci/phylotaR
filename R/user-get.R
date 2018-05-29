#' @name get_ntaxa
#' @title Count number of unique taxa
#' @description Count the number of unique taxa
#' represented by cluster(s) or sequences in phylota table
#' Use rnk to specify a taxonomic level to count. If NULL
#' counts will be made to the lowest level reported on NCBI.
#' @param phylota Phylota object
#' @param cid CLuster ID(s)
#' @param sid Sequence ID(s)
#' @param rnk Taxonomic rank
#' @param keep_higher Keep higher taxonomic ranks?
#' @return vector
#' @export
get_ntaxa <- function(phylota, cid = NULL, sid = NULL,
                      rnk = NULL, keep_higher = FALSE) {
  count <- function(sids) {
    txids <- get_txids(phylota = phylota, sid = sids,
                       rnk = rnk, keep_higher = keep_higher)
    length(unique(txids))
  }
  get_sids_and_count <- function(cid) {
    cl <- phylota@cls[[cid]]
    count(cl@sids)
  }
  if (!is.null(sid)) {
    return(count(sids = sid))
  }
  vapply(cid, get_sids_and_count, integer(1))
}

#' @name get_txids
#' @title Get taxonomic IDs by rank
#' @description Return taxonomic IDs for
#' a vector of sequence IDs or all sequences in a cluster.
#' User can specify what rank the IDs should be
#' returned. If NULL, the lowest level is returned.
#' @param phylota Phylota object
#' @param cid CLuster ID
#' @param sid Sequence ID(s)
#' @param txids Vector of txids
#' @param rnk Taxonomic rank
#' @param keep_higher Keep higher taxonomic IDs?
#' @details txids can either be provided by user or
#' they can be determined for a vector of sids or for a
#' cid.
#' If keep_higher is TRUE, any sequence that has
#' a identity that is higher than the given rank will be returned.
#' If FALSE, these sequences will return ''.
#' @return vector
#' @export
get_txids <- function(phylota, cid = NULL, sid = NULL,
                      txids = NULL, rnk = NULL,
                      keep_higher = FALSE) {
  get <- function(txid) {
    tx <- phylota@txdct@rcrds[[txid]]
    rnks <- tx@lng[['rnks']]
    ids <- tx@lng[['ids']]
    mtch <- which(rnk == rnks)
    if (length(mtch) == 0) {
      if (keep_higher) {
        # if no match, use lowest available rank
        mtch <- length(rnks)
      } else {
        return('')
      }
    }
    ids[[mtch[[1]]]]
  }
  if (is.null(txids)) {
    if (!is.null(cid)) {
      cl <- phylota@cls[[cid]]
      sid <- cl@sids
    }
    txids <- get_sq_slot(phylota = phylota,
                         sid = sid,
                         slt_nm = 'txid')
  }
  if (is.null(rnk)) {
    return(txids)
  }
  vapply(txids, get, '')
}

#' @name get_nsqs
#' @title Count number of sequences
#' @description Count the number of sequences in a cluster(s)
#' @param phylota Phylota object
#' @param cid CLuster ID(s)
#' @return vector
#' @export
get_nsqs <- function(phylota, cid) {
  count <- function(cid) {
    cl <- phylota@cls[[cid]]
    length(cl@sids)
  }
  vapply(cid, count, 1L)
}

#' @name get_sq_slot
#' @title Get slot data for each sequence
#' @description Get slot data for either or sequences
#' in a cluster of a vector of sequence IDs.
#' Use list_sqrcrd_slots() for a list of available
#' slots.
#' @param phylota Phylota object
#' @param cid Cluster ID
#' @param sid Sequence ID(s)
#' @param slt_nm Slot name
#' @return vector
#' @export
# TODO: GCR
get_sq_slot <- function(phylota, cid = NULL, sid = NULL,
                        slt_nm = list_sqrcrd_slots()) {
  get <- function(sid) {
    i <- which(sid == phylota@sqs@ids)
    sq <- phylota@sqs@sqs[[i]]
    slot(sq, slt_nm)
  }
  if (!is.null(cid)) {
    cl <- phylota@cls[[cid]]
    sid <- cl@sids
  }
  slt_nm <- match.arg(slt_nm)
  expctd <- new(getSlots('SqRcrd')[[slt_nm]], 1)
  vapply(sid, get, expctd)
}

#' @name get_cl_slot
#' @title Get slot data for each cluster record
#' @description Get slot data for cluster(s)
#' @param phylota Phylota object
#' @param cid Cluster ID
#' @param slt_nm Slot name
#' @return vector
#' @export
get_cl_slot <- function(phylota, cid,
                        slt_nm = list_clrcrd_slots()) {
  get <- function(cid) {
    i <- which(cid == phylota@cls@ids)
    cl <- phylota@cls@cls[[i]]
    slot(cl, slt_nm)
  }
  slt_nm <- match.arg(slt_nm)
  expctd <- new(getSlots('ClRcrd')[[slt_nm]], 1)
  vapply(cid, get, expctd)
}

#' @name get_tx_slot
#' @title Get slot data for each taxon record
#' @description Get slot data for taxa(s)
#' @param phylota Phylota object
#' @param txid Taxonnomic ID
#' @param slt_nm Slot name
#' @return vector or list
#' @export
get_tx_slot <- function(phylota, txid,
                        slt_nm = list_txrcrd_slots()) {
  get <- function(txid) {
    tx <- phylota@txdct@rcrds[[txid]]
    slot(tx, slt_nm)
  }
  slt_nm <- match.arg(slt_nm)
  expctd <- new(getSlots('TxRcrd')[[slt_nm]], 1)
  vapply(txid, get, expctd)
}
