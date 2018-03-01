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
