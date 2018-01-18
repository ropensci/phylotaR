
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

#' @name getTxData
#' @title Return taxonomic information
#' @description Return vector of taxonomic data values
#' using taxonomic IDs from a clsuters object.
#' @param clstrs_obj clstrs_obj
#' @param txid Taxonomic ID(s)
#' @param sltnm Taxonomic data slot name
#' @export
getTxData <- function(clstrs_obj, txid,
                      sltnm='scientificname') {
  res <- sapply(clstrs_obj@txdct[as.character(txid)],
         function(x) x[[sltnm]])
  names(res) <- NULL
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
  .pData <- function(clstr_id) {
    clstr <- clstrs_obj@clstrs[[clstr_id]]
    value <- factor(as.numeric(txids %in% clstr[['tis']]))
    data.frame(txid=as.character(txids),
               clstrid=clstr_id,
               value=value)
  }
  p_data <- plyr::mdply(.data=clstr_ids, .fun=.pData)[ ,-1]
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

#' @name nTaxa
#' @title Get n. taxa per cluster
#' @description TODO
#' @param clstrs_obj Clusters Object
#' @export
nTaxa <- function(clstrs_obj) {
  sapply(clstrs_obj@clstrs, function(x) x[['n_ti']])
}

#' @name seqLength
#' @title Get sequence length per cluster
#' @description TODO
#' @param clstrs_obj Clusters Object
#' @export
seqLength <- function(clstrs_obj) {
  # not best measure of sequence length
  # TODO: calc median sequence length in pipeline
  maxsl <- sapply(clstrs_obj@clstrs, function(x) x[['MaxLength']])
  minsl <- sapply(clstrs_obj@clstrs, function(x) x[['MinLength']])
  maxsl - minsl
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