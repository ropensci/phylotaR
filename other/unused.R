#' @name getSqDfs
#' @title Return sequence deflines
#' @description Return all deflines for a cluster's sequences
#' @param clstrs_obj clstrs_obj
#' @param id cluster ID(s)
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


#' @name readPhyLoTa
#' @title Generate a PhyLoTa object in R
#' @description Creates a PhyLoTa object containing
#' information on clusters, sequences and taxonomy
#' from the working directory of a completed pipeline.
#' @param wd Working directory
#' @export
readPhyLoTa <- function(wd) {
  # TODO check wd has completed cluster stage
  clstrs_sqs <- ldObj(wd, nm='clstrs_sqs')
  cls <- clstrs_sqs[['clstrs']]
  sqs <- clstrs_sqs[['sqs']]
  txids <- sort(unique(sqs@txids))
  sids <- sort(unique(sqs@ids))
  cids <- cls@ids
  txdct <- ldObj(wd, nm='txdct')
  new('PhyLoTa', sqs=sqs, cls=cls,
      txids=txids, sids=sids,
      cids=cids, txdct=txdct)
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
  count <- function(id) {
    cl <- clstrs@cls[[id]]
    length(unique(cl@txids))
  }
  n_taxa <- vapply(id, count, 1L)
  hist(n_taxa)
  
  id <- cid
  cl
  cid <- '27'
  cl <- clstrs@cls[[cid]]
  sqs <- clstrs@sqs[cl@sids]
  
  getNm <- function(id) {
    sq <- sqs[[id]]
    tolower(sq@nm)
  }
  ftr_nms <- vapply(sqs@ids, getNm, '')
  ftr_nms <- ftr_nms[ftr_nms != '']
  counts <- table(ftr_nms)
  prps <- counts/length(ftr_nms)
  prps[prps > .1]
  
  which(n_taxa > 250)
  clstrs@cls[['5']]
  
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


#' @name getIDFrmTxdct
#' @title Get taxonomic IDs by rank
#' @description
#' @param txdct Taxonomic dictionary
#' @param id Taxon IDs
#' @param ret Return ID or Name?
#' @param rank Rank of output
getIDFrmTxdct <- function(txdct, id, ret=c('TaxId', 'ScientificName'),
                          rank=get_ncbi_ranks()) {
  calc <- function(id) {
    rcrd <- txdct[[id]]
    rnks <- sapply(rcrd[['LineageEx']],
                   function(x) x[['Rank']])
    ids <- sapply(rcrd[['LineageEx']],
                  function(x) x[[ret]])
    mtch <- which(rank == rnks)
    if(length(mtch) == 0) {
      # if no match, use lowest available rank
      mtch <- length(rnks)
    }
    ids[[mtch[[1]]]]
  }
  rank <- match.arg(rank)
  ret <- match.arg(ret)
  sapply(as.character(id), calc)
}

#' @name dwnldTD
#' @title Download NCBI taxonomy dump
#' @description Download and unpack NCBI
#' taxonomy dump.
#' @param ps Parameters
#' @details If \code{tdpth} is left NULL, downloads
#' the taxdump.tar.gz. Otherwise, manually download
#' the file and provide path.
#' \url{ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz}
dwnldTD <- function(ps) {
  # unpack parameters
  wd <- ps[['wd']]
  v <- ps[['v']]
  tdpth <- ps[['tdpth']]
  txdr <- file.path(wd, 'taxonomy')
  if(!file.exists(txdr)) {
    dir.create(txdr)
  }
  if(is.null(tdpth)) {
    info(lvl=1, ps=ps,
         'Downloading NCBI taxdump.tar.gz ...')
    flpth <- file.path(txdr, 'taxdump.tar.gz')
    url <- 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
    res <- curl::curl_fetch_disk(url=url, path=flpth)
    if(!file.exists(res[['content']])) {
      stop('Download failed. Curl status code [',
           res[['status_code']],']')
    }
    info(lvl=1, ps=ps, 'Done.')
  } else {
    flpth <- tdpth
    if(!file.exists(flpth)) {
      stop('[', flpth, '] does not exist.')
    }
  }
  .untar(tarfile=flpth, files=c('nodes.dmp', 'names.dmp'),
         exdir=txdr)
  expflpths <- c(file.path(txdr, 'nodes.dmp'),
                 file.path(txdr, 'names.dmp'))
  if(!all(file.exists(expflpths))) {
    stop('Failed to unpack [', flpth,']')
  }
}

#' @name getStats
#' @title Generate PhyLoTa statistics for txid
#' @description 
#' @param txid Taxonomic ID
#' @param phylt_nds PhyLoTa data.frame
#' @param td_nds 'nds' element from the \code{tdobj}
#' @param td_nms 'nms' element from the \code{tdobj}
#' @param rcrsv Add stats for kid nodes recursively? T/F
#' @param ps Parameters
#' @details Returns an updated \code{phylt_nds}.
getStats <- function(txid, phylt_nds, td_nds, td_nms,
                     ps, rcrsv=FALSE) {
  info(lvl=3, ps=ps, "Adding species counts for txid [",
       txid, "]")
  rank <- CHNOSZ::getrank(txid, nodes=td_nds)
  parent <- CHNOSZ::parent(txid, nodes=td_nds)
  genus <- getGenus(txid, td_nds)
  stats <- data.frame(ti=txid,
                      ti_anc=parent,
                      rank=rank,
                      n_gi_node=0,
                      n_gi_sub_model=0,
                      n_gi_sub_nonmodel=0,
                      n_leaf_desc=0, # counts itself, for a leaf it is 1
                      n_sp_desc=0, # counts itself, for a species it is 1
                      n_sp_model=0,
                      ti_genus=genus,
                      n_genera=0,
                      n_otu_desc=1)
  if(rank == 'species') {
    stats['n_sp_desc'] <- 1
  }
  n_drctsqs <- nSqs(txid, drct=TRUE, ps=ps)
  stats['n_gi_node'] <- n_drctsqs
  kids <- getKids(txid, td_nds=td_nds)
  if(length(kids) == 0) {
    stats['n_leaf_desc'] <- stats['n_leaf_desc'] + 1
  }
  if(n_drctsqs > ps[['mdlthrs']]) {
    stats['n_gi_sub_model'] <- n_drctsqs
    stats['n_sp_model'] <- stats['n_sp_model'] + 1
  } else {
    stats['n_gi_sub_nonmodel'] <- n_drctsqs
  }
  for(kid in kids) {
    if(rcrsv) {
      phylt_nds <- getStats(txid=kid,
                            phylt_nds=phylt_nds,
                            td_nds=td_nds, td_nms=td_nms,
                            rcrsv=rcrsv, ps=ps)
    }
    # TODO: @Hannes, this section worries me. What if a user
    #  supplies phylt_nds without rows for any kids?
    #  Function crashes.
    # I kept getting kid duplicates, using [1] with `which`
    if(CHNOSZ::getrank(kid, nodes=td_nds) == 'genus') {
      stats['n_genera'] <- stats['n_genera'] + 1
    }
    cols <- c('n_leaf_desc', 'n_otu_desc', 'n_sp_desc',
              'n_sp_model', 'n_genera')
    kid_nd <- phylt_nds[which(phylt_nds[,'ti']==kid)[1],]
    stats['n_gi_sub_model'] <- stats['n_gi_sub_model'] +
      kid_nd['n_gi_sub_model']
    stats['n_gi_sub_nonmodel'] <- stats['n_gi_sub_nonmodel'] +
      kid_nd["n_gi_sub_nonmodel"]
    stats[cols] <- stats[cols] + phylt_nds[
      which(phylt_nds$ti==kid)[1], cols]
    stopifnot(dim(stats)[1]==1)
  }
  info(lvl=3, ps=ps, "Added stats for txid [", txid, "]")
  rbind(phylt_nds, stats[names(phylt_nds)])
}

# TODO: rename -- 'kids' implies tips.
#' @name getKids
#' @title Return descendent IDs
#' @description Return vector of descendent IDs
#' from NCBI taxonomic node
#' @export
# TODO: is this really 'children' and not just sub nodes?
getKids <- function(id, td_nds) {
  td_nds$id[which(td_nds$parent==id)]
}

#' @name getGenus
#' @title Return genus ID
#' @description Return vector of descendent IDs
#' from NCBI taxonomic node
#' @export
# @Hannes: eqv of .get.genus()
getGenus <- function(txid, td_nds) {
  lwr_rnks <- c('genus', 'subgenus', 'species group',
                'species subgroup', 'species',
                'subspecies','varietas', 'forma')
  hghr_rnks <- c('superkingdom', 'kingdom', 'subkingdom',
                 'superphylum', 'phylum',
                 'subphylum', 'superclass', 'class',
                 'subclass', 'infraclass', 'superorder',
                 'order', 'suborder', 'infraorder',
                 'parvorder', 'superfamily', 'family',
                 'subfamily', 'tribe', 'subtribe')
  # if rank is already higher than genus,
  #  do not search for it
  crrnt <- CHNOSZ::getrank(txid, nodes=td_nds)
  if(!crrnt  %in% lwr_rnks) {
    return(NA)
  }
  while(crrnt != 'genus') {
    txid <- CHNOSZ::parent(txid, nodes=td_nds)
    crrnt <- CHNOSZ::getrank(txid, nodes=td_nds)
    # some hybrids between genera can have a
    #  subfamily as a parent...
    if(crrnt %in% hghr_rnks) {
      break
    }
  }
  return(txid)
}


#' @name genTDObj
#' @title Generate R object from NCBI taxonomy dump
#' @description Generates an interrogable R object
#' from the NCBI taxonomy dump using the \code{CHNOSZ}
#' library. The returned object is a list of taxonomic nodes
#' and names.
#' @param ps Parameters
#' @details Object will be cached.
genTDObj <- function(ps) {
  # unpack
  wd <- ps[['wd']]
  v <- ps[['v']]
  if(!chkObj(wd, 'tdobj')) {
    info(lvl=1, ps=ps,
         'Reading from taxonomy dump ...')
    if(!file.exists(file.path(wd, 'taxonomy'))) {
      stop('No taxonomy folder in `wd`.')
    }
    nds <- CHNOSZ::getnodes(file.path(wd, 'taxonomy'))
    nms <- CHNOSZ::getnames(file.path(wd, 'taxonomy'))
    tdobj <- list('nds'=nds,
                  'nms'=nms)
    svObj(wd=wd, obj=tdobj, nm='tdobj')
    info(lvl=1, ps=ps, 'Done.')
  } else {
    tdobj <- ldObj(wd, 'tdobj')
  }
  tdobj
}

#' @name getMngblIds
#' @title Identify manageable taxonomic node IDs
#' @description Given a root \code{txid}, return a set of
#' taxonomic node IDs that only have up to a maximum number
#' of descendants.
#' @param txid Taxomomic IDs
#' @param td_nds 'nds' element from the \code{tdobj}
#' @param ps Parameters
#' @details Returns a \code{list}.
#' If there are more descendants than mx_dscndnts or
#' retrieveing the number takes longer than \code{tmout}
#' seconds, the taxonomic ID is discarded and its children
#' are added to the queue.
getMngblIds <- function(txid, td_nds, ps) {
  queue <- txid
  mngbl_ids <- vector()
  ndscndnts <- vector()
  rjctd_ids <- vector()
  tot <- 0
  while(length(queue) > 0) {
    id <- head(queue, 1)
    queue <- tail(queue, length(queue)-1)
    n <- nNcbiNds(txid=id, ps=ps)
    if (n > ps[['mxnds']]) {
      queue <- c(queue, getKids(id, td_nds))
      rjctd_ids <- c(rjctd_ids, id)
      info(lvl=1, ps=ps, "Taxon [", id,
           "] has too many descendants. Processing child taxa.")
    } else {
      mngbl_ids <- c(mngbl_ids, id)
      ndscndnts <- c(ndscndnts, n)
      info(lvl=1, ps=ps, "Taxon [", id,
           "] has maneagable number of descendants [",
           n, '].')
      tot <- tot + n
      info(lvl=1, ps=ps,
           "Current number of nodes to be processed [", tot, "]")
    }
  }
  list('mngbl_ids'=mngbl_ids, 'rjctd_ids'=rjctd_ids,
       'ndscndnts'=ndscndnts)
}

#' @name genPhylotaNds
#' @title Generate taxonomic nodes in PhyLoTa format
#' @description
#' @param nid_sets Taxonomic processable and unprocessable
#' node IDs, list as generated by \code{getMngblIds}
#' @param v v? T/F
#' @details
genPhylotaNds <- function(nid_sets, td_nds, td_nms, ps) {
  phylt_nds <- data.frame('ti'=NA,
                          'ti_anc'=NA,
                          'rank'=NA,
                          'n_gi_node'=NA,
                          'n_gi_sub_nonmodel'=NA,
                          'n_gi_sub_model'=NA,
                          'n_sp_desc'=NA,
                          'n_sp_model'=NA,
                          'n_leaf_desc'=NA,
                          'n_otu_desc'=NA,
                          'ti_genus'=NA,
                          'n_genera'=NA)
  nids <- nid_sets[['mngbl_ids']]
  info(lvl=1, ps=ps, "Adding [", length(nids),
       "] nodes recursively")
  for(i in seq_along(nids)) {
    txid <- nids[i]
    info(lvl=2, ps=ps, "Recursively processing txid [",
         txid, "] [", i, "/", length(nids), "]")
    phylt_nds <- getStats(txid=txid, phylt_nds=phylt_nds,
                          td_nds=td_nds, td_nms=td_nms,
                          ps=ps, rcrsv=TRUE)
    info(lvl=1, ps=ps, "Finished processing txid [",
         txid, "] [", i, "/", length(nids), "]")
  }
  ## Add the top nodes non-recursively. This has to be in reversed
  # order to make sure a parent is not
  ## inserted before it's children, since we need the info from
  # the children. Therefore, this must not
  ## happen in parallel!
  nids <- rev(nid_sets[['rjctd_ids']])
  info(lvl=1, ps=ps, "Adding [", length(nids),
       "] top-level nodes non-recursively, sequentially")
  for (i in seq_along(nids)) {
    txid <- nids[i]
    info(lvl=2, ps=ps, "Processing txid [",
         txid, "] [", i, "/", length(nids), "]")
    phylt_nds <- getStats(txid=txid,
                          phylt_nds=phylt_nds,
                          td_nds=td_nds,
                          td_nms=td_nms,
                          rcrsv=FALSE, ps=ps)
    info(lvl=1, ps=ps, "Finished processing txid [",
         txid, "] [", i, "/", length(nids), "]")
  }
  # rm first row
  phylt_nds[-1, ]
}

#' @name genTxdct
#' @title Generate taxonomic dictionary
#' @description Look-up all taxonomic data available in NCBI for
#' IDs in PhyLoTa table. Returns a list of esummary objects.
#' @param phylt_nds PhyLoTa nodes as generated by \code{genPhylotaNds}
genTxdct <- function(phylt_nds, ps) {
  txdct <- vector(mode='list',
                  length=nrow(phylt_nds))
  names(txdct) <- phylt_nds[['ti']]
  txids <- names(txdct)
  for(txid in txids) {
    info(lvl=2, ps=ps, "getting data for [", txid,
         "]")
    args <- list(db="taxonomy", id=txid, rettype='xml')
    rcrd <- safeSrch(func=rentrez::entrez_fetch,
                     fnm='fetch', args=args, ps=ps)
    rcrd <- XML::xmlToList(rcrd)[['Taxon']]
    # augment record for txdct tools
    rcrd[['Lineage']] <- strsplit(rcrd[['Lineage']],
                                  split='; ')[[1]]
    itslf <- list('TaxId'=rcrd[['TaxId']],
                  'ScientificName'=rcrd[['ScientificName']],
                  'Rank'=rcrd[['Rank']])
    rcrd[['LineageEx']] <- c(rcrd[['LineageEx']],
                             list('Taxon'=itslf))
    txdct[[txid]] <- rcrd
  }
  txdct
}


#' @name getADs
#' @title Get all descendants
#' @description Return all the taxonomic node IDs descending
#' from given taxonomic ID
#' @param txid Taxonomic node ID, numeric
#' @param phylt_nds PhyLoTa nodes data.frame
# TODO: create separate taxonomy look-up tools
getADs <- function(txid, phylt_nds) {
  dds <- getDDFrmPhyltNds(txid=txid, phylt_nds=phylt_nds)
  res <- dds
  for(dd in dds) {
    res <- c(res, getADs(txid=dd, phylt_nds=phylt_nds))
  }
  res
}

#' @name getDDFrmPhyltNds
#' @title Get direct descendants from PhyLoTa nodes
#' @description Find next node IDs using the PhyLoTa
#' nodes data.frame.
#' @param txid parent
#' @param phylt_nds PhyLoTa nodes data.frame
#' @details Returns vector of node IDs
getDDFrmPhyltNds <- function(txid, phylt_nds) {
  phylt_nds[phylt_nds[,'ti_anc']==txid,'ti']
}


#' @name clstrPhylt
#' @title Reformat cluster obj for PhyLoTa table
#' @description Returns a data frame with columns matchine the fields in
#' PhyLoTA's cluster table
#' @param clstrs Clusters
clstrPhylt <- function(clstrs) {
  # remove gi and id fields which won't be part of cluster table
  l <- lapply(clstrs, function(cl)cl[which(!names(cl)%in%c(
    'gis', 'tis', 'unique_id'))])
  # create data frame
  do.call(rbind, lapply(l, as.data.frame))
}

#' @name clstrCiGi
#' @title Reformat cluster obj into CI GI
#' @description TODO
#' @param clstrs Clusters
clstrCiGi <- function(clstrs) {
  ## select relevant columns from cluster objects and make data frame    
  do.call(rbind, lapply(clstrs, function(c) {
    data.frame(ti=rep(c[['ti_root']], c[['n_gi']]),
               clustid=rep(c[['ci']], c[['n_gi']]),
               cl_type=rep(c[['cl_type']], c[['n_gi']]),
               gi=c[['gis']],
               ti_of_gi=c[['tis']])
  }))
}

#' @name getGnsFrmPhyltNds
#' @title Get genus for txid using PhyLoTa nodes
#' @description Returns genus for a taxonomic node ID
#' @param txid Taxonomic node ID, numeric
#' @param phylt_nds PhyLoTa nodes data.frame
getGnsFrmPhyltNds <- function(txid, phylt_nds) {
  # @hannes can you check that this is doing what we want?
  # this always be 1 (below genus) or none (above genus)
  res <- unique(phylt_nds[phylt_nds[['ti_anc']]==txid, 'ti_genus'])
  if(length(res) > 1 || length(res) == 0) {
    res <- NA
  }
  res
}

#' @name writeClstr
#' @title Write out PhyLoTa cluster
#' @description Takes PhyLoTa data.frame and writes
#' as .tsv
#' @param phylt_nds PhyLoTa nodes data.frame
#' @details PhyLoTa data.frame must be informed by
#' clustering functions before writing out.
writeClstr <- function(txid, cldf, cigidf, sqs, ps) {
  dbpth <- file.path(ps[['wd']], 'dbfiles')
  if(!file.exists(dbpth)) {
    dir.create(dbpth)
  }
  clstrs_fl <- file.path(dbpth, paste0('dbfiles-', txid,
                                       '-clusters.tsv'))
  sqs_fl <- file.path(dbpth, paste0('dbfiles-', txid,
                                    '-seqs.tsv'))
  ci_gi_fl <- file.path(dbpth, paste0('dbfiles-', txid,
                                      '-ci_gi.tsv'))
  sqdf <- do.call(rbind, lapply(sqs, as.data.frame))
  info(lvl=1, ps=ps, "Taxid", ps[['txid']], ": writing",
       nrow(cldf), "clusters,", nrow(sqdf), "sequences,",
       nrow(cigidf), "ci_gi entries to file")
  write.table(cldf, file=clstrs_fl, append=file.exists(clstrs_fl),
              quote=FALSE, sep="\t", row.names=FALSE,
              col.names=!file.exists(clstrs_fl))
  write.table(sqdf, file=sqs_fl, append=file.exists(sqs_fl),
              quote=FALSE, sep="\t", row.names=FALSE,
              col.names=!file.exists(sqs_fl))
  write.table(cigidf, file=ci_gi_fl, append=file.exists(ci_gi_fl),
              quote=FALSE, sep="\t", row.names=FALSE,
              col.names=!file.exists(ci_gi_fl))
}

#' @name writeTax
#' @title Write out the taxonomic PhyLoTa nodes
#' @description Writes a .tsv file to disk in the 
#' PhyLoTa format for all taxonomic information.
#' @param phylt_nds PhyLoTa nodes as generated by \code{genPhylotaNds}
#' @param td_nms Taxonomic dump names
#' @param ps Parameters
writeTax <- function(phylt_nds, td_nms, ps) {
  dbpth <- file.path(ps[['wd']], 'dbfiles')
  if(!file.exists(dbpth)) {
    dir.create(dbpth)
  }
  fl <- file.path(dbpth, paste0('dbfiles-taxonomy-', ps[['txid']], '.tsv'))
  taxon_name<- as.character(td_nms[match(phylt_nds[['ti']], td_nms[['id']]),
                                   "name"])
  commons <- td_nms[grep('common', td_nms[['type']]),]
  common_name <- as.character(commons[['name']][match(phylt_nds[['ti']],
                                                      commons[['id']])])
  if(length(common_name) == 0) {
    common_name <- NA
  }
  res <- data.frame(ti=as.integer(phylt_nds[['ti']]),
                    ti_anc=as.integer(phylt_nds[['ti_anc']]),
                    terminal_flag=NA,
                    rank_flag=as.integer(!phylt_nds[['rank']]=='no rank'),
                    model=NA,
                    taxon_name=taxon_name,
                    common_name=common_name,
                    rank=as.character(phylt_nds[['rank']]),
                    n_gi_node=as.integer(phylt_nds[['n_gi_node']]),
                    n_gi_sub_nonmodel=as.integer(phylt_nds[['n_gi_sub_nonmodel']]),
                    n_gi_sub_model=as.integer(phylt_nds[['n_gi_sub_model']]),
                    n_clust_node=NA,
                    n_clust_sub=NA,
                    n_PIclust_sub=NA,
                    n_sp_desc=as.integer(phylt_nds[['n_sp_desc']]),
                    n_sp_model=as.integer(phylt_nds[['n_sp_model']]),
                    n_leaf_desc=as.integer(phylt_nds[['n_leaf_desc']]),
                    n_otu_desc=as.integer(phylt_nds[['n_otu_desc']]),
                    ti_genus=as.integer(phylt_nds[['ti_genus']]),
                    n_genera=as.integer(phylt_nds[['n_genera']])
  )
  write.table(res, file=fl, sep="\t", row.names=FALSE,
              col.names=TRUE)
  info(lvl=1, ps=ps, "Wrote [", nrow(res), "] nodes to file")
}