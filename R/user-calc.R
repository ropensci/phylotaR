#' @name calc_mad
#' @title Calculate MAD score
#' @description For all sequences in a cluster(s) the MAD score.
#' @param phylota Phylota object
#' @param cid Cluster ID(s)
#' @details MAD is a measure of the deviation in sequence length of a cluster.
#' Values range from 0 to 1. Clusters with values close to 1 have sequences with
#' similar lengths.
#' @return vector
#' @example examples/calc_mad.R
#' @export
#' @family tools-public
calc_mad <- function(phylota, cid) {
  calc <- function(cid) {
    sids <- phylota@clstrs[[cid]]@sids
    sqlns <- get_sq_slot(phylota = phylota, sid = sids, slt_nm = 'nncltds')
    sum(sqlns/(length(sqlns)*max(sqlns)))
  }
  vapply(cid, calc, numeric(1))
}

#' @name calc_wrdfrq
#' @title Calculate word frequencies
#' @description For all sequences in a cluster(s) calculate the frequency of
#' separate words in either the sequence definitions or the reported feature
#' name.
#' @param phylota Phylota object
#' @param cid Cluster ID(s)
#' @param min_frq Minimum frequency
#' @param min_nchar Minimum number of characters for a word
#' @param type Definitions (dfln) or features (nm)
#' @param ignr_pttrn Ignore pattern, REGEX for text to ignore.
#' @details By default, anything that is not alphanumeric is  ignored. 'dfln'
#' and 'nm' match the slot names in a SeqRec, see list_seqrec_slots().
#' @return list
#' @example examples/calc_wrdfrq.R
#' @export
#' @family tools-public
calc_wrdfrq <- function(phylota, cid, min_frq = 0.1, min_nchar = 1,
                        type = c('dfln', 'nm'), ignr_pttrn = "[^a-z0-9]") {
  calc <- function(cid) {
    sids <- phylota@clstrs[[cid]]@sids
    wrds <- get_sq_slot(phylota = phylota, sid = sids, slt_nm = type)
    wrds <- tolower(wrds)
    wrds <- unlist(strsplit(wrds, '\\s+'))
    wrds <- gsub(pattern = ignr_pttrn, replacement = '', x = wrds)
    pull <- vapply(wrds, nchar, 1) > min_nchar
    wrds <- wrds[pull]
    counts <- table(wrds)
    prps <- counts/length(wrds)
    prps <- sort(prps, decreasing = TRUE)
    prps[prps > min_frq]
  }
  type <- match.arg(type)
  res <- lapply(cid, calc)
  names(res) <- cid
  res
}
