#' @name hrrchcDwnld
#' @title Hierarchically get sequences for a txid
#' @description Looks up and downloads sequences for a
#' taxonomic ID.
#' @param txid Taxonomic node ID, numeric
#' @param txdct Taxonomic dictionary
#' @param ps parameters
#' @param lvl Log level
hrrchcDwnld <- function(txid, txdct, ps, lvl=0) {
  # get subtree counts if that is smaller than ps[['mdlthrs']]
  # or if there are no direct descendants
  dds <- getDDs(id=txid, txdct=txdct)
  subtree_count <- nSqs(txid, drct=FALSE, ps=ps)
  if(subtree_count <= ps[['mdlthrs']] | length(dds) == 0) {
    info(lvl=2+lvl, ps=ps, "+ whole subtree ...")
    sqs <- btchDwnldSqRcrds(txid=txid, drct=FALSE, ps=ps, lvl=lvl)
    return(sqs)
  }
  # 1st direct sqs from focal taxon, then from DDs
  info(lvl=2+lvl, ps=ps, "+ direct ...")
  sqs <- btchDwnldSqRcrds(txid=txid, drct=TRUE, ps=ps, lvl=lvl)
  info(lvl=2+lvl, ps=ps, "+ by child ...")
  for(dd in dds) {
    lvl <- lvl + 1
    info(lvl=2+lvl, ps=ps, "Working on child [id ", dd,"]")
    sqs <- c(sqs, hrrchcDwnld(txid=dd, txdct=txdct,
                              ps=ps, lvl=lvl))
  }
  sqs
}

#' @name agmntSqRcrds
#' @title Augment sequence records list
#' @description Convert to sqsrcrds and add taxids
#' @param sqs List of SqRcrds
#' @param txdct Taxonomic Dictionary
agmntSqRcrds <- function(sqs, txdct) {
  # txids are not downloaded as part of sequence, added here
  txdct_nms <- vapply(txdct@rcrds, function(x) x@scnm, '')
  txdct_ids <- vapply(txdct@rcrds, function(x) x@id, '')
  sqs_nms <- vapply(sqs, function(x) x@orgnsm, '')
  sqs_ids <- txdct_ids[match(sqs_nms, txdct_nms)]
  for(i in seq_along(sqs)) {
    sqs[[i]]@txid <- as.character(sqs_ids[[i]])
  }
  genSqRcrdBx(sqs)
}

#' @title prtDwnldSqRcrds
#' @description Download batch of sequences.
#' @details Given a set of IDs, downloads and breaks
#' up by feature information.
#' @param gis Sequence GI IDs
#' @param ps Parameter list
#' @return Vector of sequence objects
prtDwnldSqRcrds <- function(gis, ps) {
  ftch_args <- list(db="nucleotide",
                    rettype='gbwithparts',
                    retmode='xml', id=gis)
  rw_rcrds <- srchNCch(func=rentrez::entrez_fetch,
                       args=ftch_args, fnm='fetch',
                       ps=ps)
  rwRcrd2SqRcrd(rw_rcrds=rw_rcrds, gis=gis, ps=ps)
}

#' @title btchDwnldSqRcrds
#' @description Downloads sequences from GenBank in 500 ID batches.
#' @param txid NCBI taxonomic ID
#' @param drct Node-level only or subtree as well? Default FALSE.
#' @param ps parameters
#' @param lvl Log level
#' @return Vector of sequence records
btchDwnldSqRcrds <- function(txid, ps, drct=FALSE, lvl=0) {
  # Searches for GIs, returns accessions
  # test w/ golden moles 9389
  allsqs <- numeric()
  sqcnt <- nSqs(txid=txid, ps=ps, drct=drct)
  if(sqcnt < 1) {
    return(list())
  }
  gis <- getGIs(txid=txid, drct=drct, sqcnt=sqcnt,
                ps=ps)
  if(length(gis) > ps[['mdlthrs']]) {
    info(lvl=lvl+3, ps=ps, "More than [", ps[['mdlthrs']],
         ' sqs] available. Choosing at random.')
    gis <- gis[sample(1:length(gis), size=ps[['mdlthrs']],
                      replace=FALSE)]
  }
  info(lvl=lvl+3, ps=ps, "Getting [", length(gis), " sqs] ...")
  # Fetch sequences in batches
  btch <- ps[['btchsz']]
  for(i in seq(0, length(gis)-1, btch)) {
    lower <- i+1
    upper <- ifelse(i+btch<length(gis), i+btch, length(gis))
    crrnt_ids <- gis[lower:upper]
    info(lvl=lvl+4, ps=ps, "[", lower, "-", upper, "]");
    sqs <- prtDwnldSqRcrds(gis=crrnt_ids, ps=ps)
    allsqs <- c(allsqs, sqs)
  }
  allsqs
}


