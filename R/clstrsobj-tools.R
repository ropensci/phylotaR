
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

# SEQ FUNCTIONS

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
getSqLns <- function(clstrs_obj, id) {
  gis <- clstrs_obj@clstrs[[id]][['gis']]
  sapply(gis, function(x) clstrs_obj@sqs[[x]][['length']])
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

#' @name getSqTx
#' @title Return taxonomic ID per sequence
#' @description Get the taxonomic ID by sequence in
#' a cluster.
#' @param clstrs_obj clstrs_obj
#' @param id cluster ID(s)
#' @export
getSqTx <- function(clstrs_obj, id) {
  clstrs_obj@clstrs[[id]][['tis']]
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

#' @name writeSeqs
#' @title Write sequences as fasta
#' @description TODO
#' @param clstrs_obj Clusters Object
#' @param id unique cluster ID
#' @param wd Working directory, where files will be saved
#' @export
writeSeqs <- function(clstrs_obj, id, wd) {
  # TODO: allow user defined sequence description
  for(ech in id) {
    clstr <- clstrs_obj@clstrs[[ech]]
    sqs <- clstrs_obj@sqs[clstr$gis]
    filename <- paste0(clstr[['unique_id']], '.fasta')
    fasta <- ''
    for(i in seq_along(sqs)) {
      sq <- sqs[[i]]
      fasta <- paste0(fasta, '>', sq[['gi']],
                      '|', sq[['ti']], '\n',
                      sq[['seq']], '\n\n')
    }
    cat(fasta, file=file.path(wd, filename))
  }
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
  clstr <- clstrs_obj@clstrs[[id]]
  
}