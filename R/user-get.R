#' @name get_ntaxa
#' @title Count number of unique taxa
#' @description Count the number of unique taxa represented by cluster(s) or
#' sequences in phylota table Use rnk to specify a taxonomic level to count. If
#' NULL counts will be made to the lowest level reported on NCBI.
#' @param phylota Phylota object
#' @param cid Cluster ID(s)
#' @param sid Sequence ID(s)
#' @param rnk Taxonomic rank
#' @param keep_higher Keep higher taxonomic ranks?
#' @return vector
#' @example examples/get_ntaxa.R
#' @export
#' @family tools-public
get_ntaxa <- function(phylota, cid = NULL, sid = NULL, rnk = NULL,
                      keep_higher = FALSE) {
  count <- function(sids) {
    txids <- get_txids(phylota = phylota, sid = sids, rnk = rnk,
                       keep_higher = keep_higher)
    length(unique(txids))
  }
  get_sids_and_count <- function(cid) {
    clstr <- phylota@clstrs[[cid]]
    count(clstr@sids)
  }
  if (!is.null(sid)) {
    return(count(sids = sid))
  }
  vapply(cid, get_sids_and_count, integer(1))
}

#' @name get_txids
#' @title Get taxonomic IDs by rank
#' @description Return taxonomic IDs for a vector of sequence IDs or all
#' sequences in a cluster. User can specify what rank the IDs should be
#' returned. If NULL, the lowest level is returned.
#' @param phylota Phylota object
#' @param cid Cluster ID
#' @param sid Sequence ID(s)
#' @param txids Vector of txids
#' @param rnk Taxonomic rank
#' @param keep_higher Keep higher taxonomic IDs?
#' @details txids can either be provided by user or they can be determined for
#' a vector of sids or for a cid. If keep_higher is TRUE, any sequence that has
#' a identity that is higher than the given rank will be returned. If FALSE,
#' these sequences will return ''.
#' @return vector
#' @export
#' @example examples/get_txids.R
#' @family tools-public
get_txids <- function(phylota, cid = NULL, sid = NULL, txids = NULL, rnk = NULL,
                      keep_higher = FALSE) {
  get <- function(txid) {
    tx <- phylota@txdct@recs[[txid]]
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
      clstr <- phylota@clstrs[[cid]]
      sid <- clstr@sids
    }
    txids <- get_sq_slot(phylota = phylota, sid = sid, 
                         slt_nm = 'txid')
  }
  if (is.null(rnk)) {
    return(txids)
  }
  vapply(txids, get, '')
}

#' @name get_nsqs
#' @title Count number of sequences
#' @description Count the number of sequences in a cluster(s).
#' @param phylota Phylota object
#' @param cid Cluster ID(s)
#' @return vector
#' @export
#' @family tools-public
#' @example examples/get_nsqs.R
get_nsqs <- function(phylota, cid) {
  count <- function(cid) {
    clstr <- phylota@clstrs[[cid]]
    length(clstr@sids)
  }
  vapply(cid, count, 1L)
}

#' @name get_sq_slot
#' @title Get slot data for each sequence
#' @description Get slot data for either or sequences in a cluster of
#' a vector of sequence IDs. Use list_seqrec_slots() for a list of
#' available slots.
#' @param phylota Phylota object
#' @param cid Cluster ID
#' @param sid Sequence ID(s)
#' @param slt_nm Slot name
#' @return vector
#' @export
#' @family tools-public
#' @example examples/get_sq_slot.R
get_sq_slot <- function(phylota, cid = NULL, sid = NULL,
                        slt_nm = list_seqrec_slots()) {
  get <- function(sid) {
    i <- which(sid == phylota@sqs@ids)
    sq <- phylota@sqs@sqs[[i]]
    slot(sq, slt_nm)
  }
  if (!is.null(cid)) {
    clstr <- phylota@clstrs[[cid]]
    sid <- clstr@sids
  }
  slt_nm <- match.arg(slt_nm)
  expctd <- new(getSlots('SeqRec')[[slt_nm]], 1)
  vapply(sid, get, expctd)
}

#' @name get_clstr_slot
#' @title Get slot data for each cluster record
#' @description Get slot data for cluster(s)
#' @param phylota Phylota object
#' @param cid Cluster ID
#' @param slt_nm Slot name
#' @return vector
#' @family tools-public
#' @export
#' @example examples/get_clstr_slot.R
get_clstr_slot <- function(phylota, cid,
                           slt_nm = list_clstrrec_slots()) {
  get <- function(cid) {
    i <- which(cid == phylota@clstrs@ids)
    clstr <- phylota@clstrs@clstrs[[i]]
    slot(clstr, slt_nm)
  }
  slt_nm <- match.arg(slt_nm)
  expctd <- new(getSlots('ClstrRec')[[slt_nm]], 1)
  vapply(cid, get, expctd)
}

#' @name get_tx_slot
#' @title Get slot data for each taxon record
#' @description Get slot data for taxa(s)
#' @param phylota Phylota object
#' @param txid Taxonomic ID
#' @param slt_nm Slot name
#' @return vector or list
#' @export
#' @family tools-public
#' @example examples/get_tx_slot.R
get_tx_slot <- function(phylota, txid, slt_nm = list_taxrec_slots()) {
  get <- function(txid) {
    tx <- phylota@txdct@recs[[txid]]
    slot(tx, slt_nm)
  }
  slt_nm <- match.arg(slt_nm)
  expctd <- new(getSlots('TaxRec')[[slt_nm]], 1)
  vapply(txid, get, expctd)
}

#' @name get_stage_times
#' @title Get run times for different stages
#' @description Get slot data for taxa(s)
#' @param wd Working directory
#' @return list of runtimes in minutes
#' @export
#' @family tools-public
#' @example examples/get_stage_times.R
get_stage_times <- function(wd) {
  lgfl <- file.path(wd, 'log.txt')
  lines <- readLines(con = lgfl)
  stage_tms <- stage_ends <- stage_starts <- c('taxise' = NA,
                                               'download' = NA,
                                               'cluster' = NA,
                                               'cluster\\^2' = NA)
  stage_nms <- names(stage_tms)
  for (ln in lines) {
    for (stgnm in stage_nms) {
      pttrn <- paste0('Starting stage ', stgnm, ': ')
      if (grepl(pattern = pttrn, x = ln, ignore.case = TRUE)) {
        ln <- sub(pattern = pttrn,
                  replacement = '', x = ln, ignore.case = TRUE)
        ln <- gsub(pattern = '(\\[|\\])', replacement = '', x = ln)
        stage_starts[[stgnm]] <- ln
      }
      pttrn <- paste0('Completed stage ', stgnm, ': ')
      if (grepl(pttrn, ln, ignore.case = TRUE)) {
        ln <- sub(pattern = pttrn,
                  replacement = '', x = ln, ignore.case = TRUE)
        ln <- gsub(pattern = '(\\[|\\])', replacement = '', x = ln)
        stage_ends[[stgnm]] <- ln
      }
    }
  }
  for (stgnm in stage_nms) {
    stage_tms[[stgnm]] <- difftime(as.POSIXct(stage_ends[[stgnm]]),
                                   as.POSIXct(stage_starts[[stgnm]]),
                                   units = 'mins')
  }
  stage_tms
}
