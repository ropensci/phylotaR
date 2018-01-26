
#' @name genClstrsObj
#' @title Create a clstrsObj
#' @description TODO
#' @param wd Working directory
#' @export
genClstrsObj <- function(wd) {
  # TODO check wd has completed cluster stage
  clstrs_sqs <- ldObj(wd, nm='clstrs_sqs')
  clstrs <- clstrs_sqs[['clstrs']]
  names(clstrs) <- sapply(clstrs, function(x) x[['unique_id']])
  txids <- unique(unlist(sapply(clstrs, function(x) x[['tis']])))
  sqids <- unique(unlist(sapply(clstrs, function(x) x[['gis']])))
  clstr_ids <- names(clstrs)
  txdct <- ldObj(wd, nm='txdct')
  new('ClstrsObj', sqs=clstrs_sqs[['sqs']],
      clstrs=clstrs, txids=txids, sqids=sqids,
      clstr_ids=clstr_ids, txdct=txdct)
}

# CLUSTER FUNCTIONS
#' @name getClGCR
#' @title Return mean GC-ratio of sequences in cluster
#' @description Return the mean proportion of G or C in
#' the unambiguous sequence of all sequences in cluster(s)
#' @param clstrs_obj clstrs_obj
#' @param id cluster ID(s)
#' @export
getClGCR <- function(clstrs_obj, id) {
  sapply(id, function(x) mean(getSqGCR(clstrs_obj, x)))
}

#' @name getClNTx
#' @title Get n. taxa per cluster
#' @description TODO
#' @param clstrs_obj Clusters Object
#' @export
getClNTx <- function(clstrs_obj, id) {
  sapply(clstrs_obj@clstrs[id], function(x) x[['n_ti']])
}

#' @name getClMdLn
#' @title Get sequence length per cluster
#' @description TODO
#' @param clstrs_obj Clusters Object
#' @export
getClMdLn <- function(clstrs_obj, id) {
  sapply(id, function(x) median(getSqLns(clstrs_obj, x)))
}

#' @name byLineage
#' @title Find clusters with a given taxonomic group
#' @description Returns a boolean vector of all the
#' clusters that have representatives of a given
#' taxonomic group.
#' @param clstrs_obj clstrs_obj
#' @param nm Taxonomic group name
#' @export
byLineage <- function(clstrs_obj, nm) {
  .calc <- function(clstr) {
    txids <- unique(clstr[['tis']])
    vals <- getTxData(clstrs_obj=clstrs_obj,
                      txid=txids, sltnm='Lineage')
    res <- sapply(vals, function(x) value %in% x)
    any(res)
  }
  sapply(clstrs_obj@clstrs, .calc)
}

# OTHER

#' @name getTxData
#' @title Return taxonomic information
#' @description Return vector of taxonomic data values
#' using taxonomic IDs from a clusters object.
#' @param clstrs_obj clstrs_obj
#' @param txid Taxonomic ID(s)
#' @param sltnm Taxonomic data slot name
#' @export
getTxData <- function(clstrs_obj, txid,
                      sltnm=c("ScientificName", "CommonName",
                              "OtherNames", "ParentTaxId", "Rank",
                              "Lineage", "LineageEx")) {
  sltnm <- match.arg(sltnm)
  if(sltnm %in% c("OtherNames", "Lineage", "LineageEx")) {
    res <- lapply(clstrs_obj@txdct[as.character(txid)],
                  function(x) x[[sltnm]])
  } else if(sltnm == "CommonName") {
    res <- sapply(clstrs_obj@txdct[as.character(txid)],
                  function(x) x[["OtherNames"]][["CommonName"]])
    names(res) <- NULL
  } else {
    res <- sapply(clstrs_obj@txdct[as.character(txid)],
                  function(x) x[[sltnm]])
    names(res) <- NULL
  }
  res
}

#' @name plotTable
#' @title Plot table of txid presence by clusters
#' @description TODO
#' @param clstrs_obj clstrs_obj
#' @param clstr_ids IDs of clusters in table
#' @export
plotTable <- function(clstrs_obj, clstr_ids=clstrs_obj@clstr_ids,
                      txids=clstrs_obj@txids,
                      clstr_names=clstrs_obj@clstr_ids,
                      txid_names=clstrs_obj@txids) {
  pData <- function(clstr_id) {
    clstr <- clstrs_obj@clstrs[[clstr_id]]
    value <- factor(as.numeric(txids %in% clstr[['tis']]))
    data.frame(txid=as.character(txids),
               clstrid=clstr_id,
               value=value)
  }
  p_data <- plyr::mdply(.data=clstr_ids, .fun=pData)[ ,-1]
  p_data[['clstrnm']] <- 
    clstr_names[match(p_data[['clstrid']], clstr_ids)]
  p_data[['txidnm']] <- 
    txid_names[match(p_data[['txid']], txids)]
  # lvls <- clstr_names[match(clstr_ids, names(sort(clstr_id_ord)))]
  # p_data[['clstrnm']] <- factor(p_data[['clstrnm']],
  #                               levels=lvls,
  #                               ordered=TRUE)
  # names(txid_ord) <- txids
  # lvls <- txid_names[match(txids, names(sort(txid_ord)))]
  # p_data[['txidnm']] <- factor(p_data[['txidnm']],
  #                              levels=lvls,
  #                              ordered=TRUE)
  ggplot2::ggplot(p_data, ggplot2::aes(clstrnm, txidnm)) +
    ggplot2::geom_tile(ggplot2::aes(fill=value)) +
    ggplot2::xlab('') + ggplot2::ylab('') +
    ggplot2::scale_fill_manual(values=c('white', 'black')) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position='none') 
}

#' @name getClstr
#' @title Get cluster
#' @description TODO
#' @param clstrs_obj Clusters Object
#' @param id Unique cluster ID(s)
#' @export
getClstr <- function(clstrs_obj, id) {
  clstrs_obj@clstrs[id]
}

#' @name getSqs
#' @title Get seuqneces for a cluster
#' @description TODO
#' @param clstrs_obj Clusters Object
#' @param id Sequence ID(s)
#' @export
getSqs <- function(clstrs_obj, id) {
  clstrs_obj@sqs[id]
}

#' @name writeSqs
#' @title Write sequences as fasta
#' @description TODO
#' @param sids Sequence IDs
#' @param dflns Definition lines
#' @param flpth File path and name
#' @export
writeSqs <- function(sids, dflns, flpth) {
  sqs <- clstrs_obj@sqs[sids]
  fasta <- ''
  for(i in seq_along(sqs)) {
    sq <- sqs[[i]]
    dfln <- dflns[[i]]
    fasta <- paste0(fasta, '>', dfln, '\n',
                    sq[['seq']], '\n\n')
  }
  cat(fasta, file=flpth)
}

#' @name getBestClstrs
#' @title Get best clusters
#' @description TODO
#' @param clstrs_obj Clusters Object
#' @param n Number of cluster IDs to return
#' @export
getBestClstrs <- function(clstrs_obj, n) {
  # return IDs of best clusters
  # TODO: first order by most taxa, then sequence length
  # TODO: shouldn't there be an element of selecting a range of gene regions?
  mst_taxa <- sort(nTaxa(clstrs_obj), decreasing=TRUE)
  mst_sq <- sort(seqLength(clstrs_obj), decreasing=TRUE)
  names(mst_taxa[1:n])
}

#' @name drpClstrs
#' @title Drop clusters from clusters object
#' @description Drop clusters from clusters object. Only
#' the ids provided are kept.
#' @param clstrs_obj Clusters Object
#' @param id Clusters IDs to be kept
#' @export
drpClstrs <- function(clstrs_obj, id) {
  clstrs_obj@clstrs <- clstrs_obj@clstrs[id]
  clstrs_obj@clstr_ids <- names(clstrs_obj@clstrs)
  clstrs_obj
}

#' @name drpSqs
#' @title Drop sequences from cluster
#' @description Drop sequences from a cluster
#' @param clstrs_obj Clusters Object
#' @param cid Clusters ID
#' @param sid Sequence IDs to be kept
#' @export
drpSqs <- function(clstrs_obj, cid, sid) {
  clstr <- clstrs_obj@clstrs[[cid]]
  txids <- getSqTx(clstrs_obj=clstrs_obj,
                   sid=sid, rank=FALSE)
  names(sid) <- names(txids) <- NULL
  clstr$gis <- sid
  clstr$tis <- txids
  clstr$n_ti <- length(unique(txids))
  clstrs_obj@clstrs[[cid]] <- clstr
  clstrs_obj
}


#' @name fltrClstrSqs
#' @title Filter out sequences for a cluster
#' @description Identify all sequences in a cluster
#'  keep only those sequences representing rank of choice (e.g. species)
#'  identify all sequences of the same taxa
#'  select 'best' sequence for each taxon
#'   - drop all with too much ambiguity
#'   - select longest sequence per taxon
#' @param clstrs_obj Clusters Object
#' @param id Clusters ID
#' @param rank Taxonomic rank
#' @param mn_pambg Min. ambiguous bases
#' @export
fltrClstrSqs <- function(clstrs_obj, id, rank='species',
                         mn_pambg=0.1) {
  clstr <- clstrs_obj@clstrs[[id]]
  pambgs <- getSqAmbs(clstrs_obj=clstrs_obj, id=id)
  gis <- names(pambgs)[pambgs < mn_pambg]
  txids <- getSqTx(clstrs_obj, sid=gis, rank=rank)
  sqlns <- getSqLns(clstrs_obj, sid=gis)
  names(sqlns) <- gis
  keep <- tapply(sqlns, factor(txids),
                 function(x) names(x)[which.max(x)])
  drpSqs(clstrs_obj=clstrs_obj, cid=id, sid=keep)
}
