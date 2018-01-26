#' @name getSqDfs
#' @title Return sequence deflines
#' @description Return all deflines for a cluster's sequences
#' @param clstrs_obj clstrs_obj
#' @param id cluster ID(s)
#' @export
getSqDfs <- function(clstrs_obj, id, prse=0.1) {
  gis <- clstrs_obj@clstrs[[id]][['gis']]
  dflns <- sapply(gis, function(x) clstrs_obj@sqs[[x]][['def']])
  if(!is.null(prse)) {
    wrds <- unlist(strsplit(dflns, ' '))
    wrds <- gsub(pattern="[^a-zA-Z0-9]", replacement='',
                 x=wrds)
    wrdfqs <- table(wrds)/length(dflns)
    sort(wrdfqs, decreasing = TRUE)[1:10]
    cmmn <- names(wrdfqs)[wrdfqs > prse]
    return(paste0(cmmn, collapse=' | '))
  }
  dflns
}

#' @name getSqLns
#' @title Return sequence deflines
#' @description Return all deflines for a cluster's sequences
#' @param clstrs_obj clstrs_obj
#' @param id cluster ID(s)
#' @export
getSqLns <- function(clstrs_obj, cid=NULL, sid=NULL) {
  if(!is.null(cid)) {
    sid <- clstrs_obj@clstrs[[cid]][['gis']]
  }
  sapply(sid, function(x) clstrs_obj@sqs[[x]][['length']])
}

#' @name getSqAmbs
#' @title Return sequence deflines
#' @description Return all deflines for a cluster's sequences
#' @param clstrs_obj clstrs_obj
#' @param id cluster ID(s)
#' @export
getSqAmbs <- function(clstrs_obj, id) {
  .calc <- function(x) {
    sq <- clstrs_obj@sqs[[x]][['seq']]
    res <- gregexpr(pattern='[^atcgATCG]',
                    text=sq)[[1]]
    length(res)/nchar(sq)
  }
  gis <- clstrs_obj@clstrs[[id]][['gis']]
  sapply(gis, .calc)
}

#' @name getSqGCR
#' @title Return GC-ratio by sequences
#' @description Return the proportion of G or C in
#' the unambiguous sequence of all sequences in cluster
#' @param clstrs_obj clstrs_obj
#' @param id cluster ID(s)
#' @export
getSqGCR <- function(clstrs_obj, id) {
  .calc <- function(x) {
    sq <- clstrs_obj@sqs[[x]][['seq']]
    unambsq <- gsub(pattern='[^atcgATCG]',
                    replacement='', x=sq)[[1]]
    res <- gregexpr(pattern='[^cgCG]',
                    text=unambsq)[[1]]
    length(res)/nchar(unambsq)
  }
  gis <- clstrs_obj@clstrs[[id]][['gis']]
  sapply(gis, .calc)
}

#' @name getClMAD
#' @title Return Median Alignment Density
#' @description 
#' @param clstrs_obj clstrs_obj
#' @param id cluster ID(s)
#' @export
getClMAD <- function(clstrs_obj, cids) {
  calc <- function(cid) {
    sqlns <- getSqLns(clstrs_obj, cid=cid)
    sum(sqlns/(length(sqlns)*max(sqlns)))
  }
  sapply(cids, calc)
}

#' @name getSqTx
#' @title Return taxonomic ID per sequence
#' @description Get the taxonomic ID by sequence in
#' a cluster.
#' @param clstrs_obj clstrs_obj
#' @param id cluster ID(s)
#' @export
getSqTx <- function(clstrs_obj, cid=NULL, sid=NULL,
                    rank=FALSE) {
  if(!is.null(cid)) {
    txids <- clstrs_obj@clstrs[[cid]][['tis']]
  } else {
    txids <- sapply(sid,
                    function(x) clstrs_obj@sqs[[x]][['ti']])
  }
  if(rank != FALSE) {
    txids <- getIDFrmTxdct(txdct=clstrs_obj@txdct,
                         id=txids, rank=rank)
  }
  txids
}
