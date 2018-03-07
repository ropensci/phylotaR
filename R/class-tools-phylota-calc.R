
#' @name calc_mad
#' @title Calculate MAD score
#' @description For all sequences in a cluster(s) the
#' MAD score.
#' @param phylota Phylota object
#' @param cid CLuster ID(s)
#' @param sid Sequence ID(s)
#' @details MAD is a measure of the deviation in 
#' sequence length of a cluster. Values range from 0 to
#' 1. Clusters with values close to 1 have sequences with
#' similar lengths.
#' @return vector
#' @export
calc_mad <- function(phylota, cid) {
  calc <- function(cid) {
    sids <- phylota@cls[[cid]]
    sqlns <- get_sq_slot(phylota=phylota,
                         sids=sids,
                         slt_nm='nncltds')
    sum(sqlns/(length(sqlns)*max(sqlns)))
  }
  vapply(cid, calc, integer(1))
}

#' @name calc_wrdfrq
#' @title Calculate word frequencies
#' @description For all sequences in a cluster(s)
#' calculate the frequency of separate words in
#' either the sequnece definitions or the reported
#' feature name.
#' @param phylota Phylota object
#' @param cid CLuster ID(s)
#' @param min_frq Minimum frequency
#' @param type Definitions (dfln) or features (nm)
#' @param ignr_pttrn Ignore pattern, REGEX for text to ignore.
#' @details By default, anything that is not alphanumeric is 
#' ignored.
#' 'dfln' and 'nm' match the slot names in a SqRcrd, see
#' list_sqrcrd_slots().
#' @return list
#' @export
calc_wrdfrq <- function(phylota, cid, min_frq=0.1,
                        type=c('dfln', 'nm'),
                        ignr_pttrn="[^a-z0-9]") {
  calc <- function(cid) {
    sids <- phylota@cls[[cid]]
    wrds <- get_sq_slot(phylota=phylota,
                        sids=sids,
                        slt_nm=type)
    wrds <- tolower(wrds)
    wrds <- strsplit(wrds, '\\s+')[[1]]
    wrds <- gsub(pattern=ignr_pttrn,
                 replacement='', x=wrds)
    wrds <- wrds[wrds != '']
    counts <- table(wrds)
    prps <- counts/length(ftr_nms)
    prps[prps > .1]
  }
  type <- match.arg(type)
  res <- lapply(sqs@ids, calc)
  names(res) <- cid
  res
}