## Taxa can have "direct" sequence links up to some certain level, e.g. genus.
# This means that no sequnces of the children
## of the taxon are included. These clusters will be stored separately
## as "node" clusters. In genbank, one searches for direct links using keyword
# 'noexp' and 'exp' for subtree links.
## Below, we also calulate the direct ('node') clusters for the taxon of
# interest
## In NCBI, subtree search on terminals (taxa without children)
# gives the direct sequences.
## Therefore, we have to exclude terminals from direct cluster calculation,
# because their
## clusters are already calculated in the 'subtree' (direct=false) mode!
## two numbers can only be calculated if we have cluster
# info on multiple taxonomic levels:
##  n_child, the number of child clusters, and ci_anc,
# the parent cluster.
##  Since we are dealing with single-likeage clusters,
# we can identify the parent cluster of a
##  cluster if at least one gi of both clusters is the same.
## Get the indices of the parent clusters that contain
#a gi of the child clusters
## If a gi is in two clusters (e.g. parent and grandparent),
# take the lower one, by taking
## the higher index

#' @name calcClstrs
#' @title Calculate clusters for all sequences in WD
#' @description TODO
#' @param phylt_nds PhyLoTa nodes data.frame
#' @export
calcClstrs <- function(txid, phylt_nds, ps) {
  root_txid <- txid
  # load sequences
  sq_fls <- list.files(file.path(ps[['wd']], 'cache', 'sqs'))
  for(i in seq_along(sq_fls)) {
    sq_fl <- sq_fls[i]
    # TODO: use the cache tool
    sqs <- readRDS(file=file.path(file.path(ps[['wd']], 'cache',
                                            'sqs', sq_fl)))
    txid <- as.numeric(sub('\\.RData', '', sq_fl))
    # cluster
    clstrs <- clstrSqs(txid=txid, sqs=sqs,
                       phylt_nds=phylt_nds,
                       drct=FALSE, infrmtv=FALSE,
                       blst_rs=NULL, ps=ps)
    cldf <- clstrPhylt(clstrs=clstrs)
    cigidf <- clstrCiGi(clstrs=clstrs)
    # output
    info(lvl=1, ps=ps, "Taxid [", txid, "]: writing [", nrow(cldf),
        "] clusters, [", length(sqs), "] sequences, [",
        nrow(cigidf), "] ci_gi entries to file")
    # TODO move these to sep. func
    clstrs_fl <- file.path(ps[['wd']], paste0('dbfiles-', root_txid,
                                      '-clusters.tsv'))
    sqs_fl <- file.path(ps[['wd']], paste0('dbfiles-', root_txid,
                                   '-seqs.tsv'))
    ci_gi_fl <- file.path(ps[['wd']], paste0('dbfiles-', root_txid,
                                     '-ci_gi.tsv'))
    write.table(cldf, file=clstrs_fl, append=file.exists(clstrs_fl),
                quote=FALSE, sep="\t", row.names=FALSE,
                col.names=!file.exists(clstrs_fl))
    # write.table(seqdf, file=sqs_fl, append=file.exists(sqs_fl),
    #             quote=FALSE, sep="\t", row.names=FALSE,
    #             col.names=!file.exists(sqs_fl))
    write.table(cigidf, file=ci_gi_fl, append=file.exists(ci_gi_fl),
                quote=FALSE, sep="\t", row.names=FALSE,
                col.names=!file.exists(ci_gi_fl))
    info(lvl=1, ps=ps, "Finished processing taxid [", txid, "] # [",
        i, "/", length(sq_fls), "]")
  }
}

#' @name clstrSqs
#' @title Cluster sequences for sequences
#' @description Identify sequence clusters using BLAST
#' @param phylt_nds PhyLoTa nodes data.frame
#' @export
# TODO: too big, break down
clstrSqs <- function(txid, sqs, phylt_nds, ps,
                     drct=FALSE, infrmtv=FALSE,
                     blst_rs=NULL) {
  all_clstrs <- list()
  clstr_typ <- ifelse(drct, "direct", "subtree")
  rnks <- as.character(phylt_nds[['rank']])
  rnk <- as.character(rnks[match(txid, phylt_nds[['ti']])])
  info(lvl=1, ps=ps, "Processing taxid [", txid, "] of rank [",
      rnk, "] attempting to make [", clstr_typ, "] clusters")
  if(!drct) {
    dds <- getDDFrmPhyltNds(txid=txid, phylt_nds=phylt_nds)
    if(length(dds) > 0) {
      all_clstrs <- c(all_clstrs, clstrSqs(txid, ps=ps,
                                           phylt_nds=phylt_nds,
                                           sqs=sqs, drct=TRUE,
                                           infrmtv=infrmtv,
                                           blst_rs=blst_rs))
    }
  }
  if(drct) {
    txids <- txid
  } else {
    txids <- getADs(txid=txid, phylt_nds=phylt_nds)
  }
  all_sq_txids <- sapply(sqs, '[[', 'ti')
  sq_txids <- names(sqs[which(all_sq_txids %in% txids)])
  info(lvl=1, ps=ps, "Found [", length(sq_txids), '] [',
      clstr_typ, "] GIs for taxid [", txid, "]")
  # @Hannes: you have a max seqs, why not a min?
  if(length(sq_txids) < 3) {
    info(lvl=1, ps=ps, "Insufficient [", clstr_typ,
        "] sequences for taxid [",
        txid, "], cannot make clusters")
    return(all_clstrs)
  }
  sqs_prt <- sqs[sq_txids[sq_txids %in% names(sqs)]]
  if(is.null(blst_rs)) {
    info(lvl=1, ps=ps, "Current BLAST result is NULL from parent, performing BLAST")
    txids_blst_rs <- blstSqs(txid=txid, typ=clstr_typ, sqs=sqs_prt, ps=ps)
    # TODO: Can we do this more elegantly?
    if(is.null(txids_blst_rs)) {
      info(lvl=1, ps=ps, "Current BLAST result is NULL after blasting")
      return(all_clstrs)
    }
  } else {
    txids_blst_rs <- list()
    info(lvl=1, ps=ps, "Using BLAST results from parent cluster")
    txids_blst_rs <- blst_rs[which(blst_rs$subject.id %in% sq_txids),]
    txids_blst_rs <- txids_blst_rs[which(txids_blst_rs$query.id %in% sq_txids),]
  }
  # Sq clusters are stored in a list of lists, named by gi
  # Get top-level clusters
  info(lvl=1, ps=ps, "Clustering BLAST results for taxid [", txid, "]")
  # TODO
  raw_clstrs <- clstrBlstRs(blst_rs=txids_blst_rs, infrmtv=infrmtv)
  info(lvl=1, ps=ps, "Finished clustering BLAST results for taxid [",
      txid, "]")
  # make dataframe with fields as in PhyLoTa database
  # TODO
  clstrs <- addClstrInf(clstrs=raw_clstrs, phylt_nds=phylt_nds,
                        txid=txid, sqs=sqs, drct=drct)
  info(lvl=1, ps=ps, "Generated [", length(clstrs),
      "] clusters")
  all_clstrs <- c(all_clstrs, clstrs)
  # If we do not calculate a direct cluster here, iterate over
  # txid's children to retrieve clusters
  if(!drct) {
    dds <- getDDFrmPhyltNds(txid=txid, phylt_nds=phylt_nds)
    for(dd in dds) {
      info(lvl=1, ps=ps, "Processing child taxon of [", txid,
          "] [", dd, "]")
      dd_clstrs <- clstrSqs(txid=dd, phylt_nds=phylt_nds, 
                            sqs=sqs, blst_rs=txids_blst_rs,
                            infrmtv=infrmtv, ps=ps,
                            drct=FALSE)
      # TODO: break this up. Too nested.
      dd_clstrs <- lapply(dd_clstrs, function(cc) {
        idx <- max(which(sapply(clstrs, function(c) any(
          cc[['gis']] %in% c[['gis']]))))
        cc[['ci_anc']] <- clstrs[[idx]][['ci']]
        clstrs[[idx]][['n_child']] <-
          clstrs[[idx]][['n_child']] + 1
        cc
      })
      all_clstrs <- c(all_clstrs, dd_clstrs)
    }
  }
  all_clstrs
}

#' @name writeClstr
#' @title Write out PhyLoTa cluster
#' @description Takes PhyLoTa data.frame and writes
#' as .tsv
#' @param phylt_nds PhyLoTa nodes data.frame
#' @details PhyLoTa data.frame must be informed by
#' clustering functions before writing out.
#' @export
# TODO
writeClstr <- function(phylt_nds, ps) {
  ## Write all data to file
  info(lvl=1, ps=ps, "Taxid", txid, ": writing",
       nrow(cldf), "clusters,", nrow(seqdf), "sequences,",
       nrow(cigidf), "ci_gi entries to file")
  write.table(cldf, file=clusters.file, append=file.exists(clusters.file),
              quote=FALSE, sep="\t", row.names=FALSE,
              col.names=!file.exists(clusters.file))
  write.table(seqdf, file=seqs.file, append=file.exists(seqs.file),
              quote=FALSE, sep="\t", row.names=FALSE,
              col.names=!file.exists(seqs.file))
  write.table(cigidf, file=ci_gi.file, append=file.exists(ci_gi.file),
              quote=FALSE, sep="\t", row.names=FALSE,
              col.names=!file.exists(ci_gi.file))
  info(lvl=1, ps=ps, "Finished processing taxid ", txid,
       " # ", i, " / ", length(txids), "]")
}

#' @name getADs
#' @title Get all descendants
#' @description Return all the taxonomic node IDs descending
#' from given taxonomic ID
#' @param txid Taxonomic node ID, numeric
#' @param phylt_nds PhyLoTa nodes data.frame
#' @export
# TODO: create separate taxonomy look-up tools
getADs <- function(txid, phylt_nds) {
  dds <- getDDFrmPhyltNds(txid=txid, phylt_nds=phylt_nds)
  res <- dds
  for(dd in dds) {
    res <- c(res, getADs(txid=dd, phylt_nds=phylt_nds))
  }
  res
}

#' @name getGnsFrmPhyltNds
#' @title Get genus for txid using PhyLoTa nodes
#' @description Returns genus for a taxonomic node ID
#' @param txid Taxonomic node ID, numeric
#' @param phylt_nds PhyLoTa nodes data.frame
#' @export
getGnsFrmPhyltNds <- function(txid, phylt_nds) {
  # @hannes can you check that this is doing what we want?
  # this always be 1 (below genus) or none (above genus)
  res <- unique(phylt_nds[phylt_nds[['ti_anc']]==txid, 'ti_genus'])
  if(length(res) > 1 || length(res) == 0) {
    res <- NA
  }
  res
}

#' @name blstSqs
#' @title BLAST All vs All
#' @description Return blast results from
#' BLASTing all vs all for given sequences
#' @param txid Taxonomic node ID, numeric
#' @param typ Cluster type, 'direct' or 'subtree'
#' @param sqs Sequences
#' @param wd Working directory
#' @param verbose Verbose? T/F
#' @export
blstSqs <- function(txid, typ, sqs, ps) {
  info(lvl=1, ps=ps, "BLAST all vs all for [",
      length(sqs), "] sequences")
  dbfl <- paste0('taxon-', txid, '-typ-', typ,
                 '-db.fa')
  outfl <- paste0('taxon-', txid, '-typ-', typ,
                  '-blastout.txt')
  mkBlstDB(sqs, dbfl=dbfl, ps=ps)
  blst_rs <- blstN(dbfl=dbfl, outfl=outfl, ps=ps)
  # TODO: Not so elegant
  if(is.null(blst_rs)) {
    return(NULL)
  }
  info(lvl=1, ps=ps, "Number of BLAST results [",
       nrow(blst_rs), "]")
  info(lvl=1, ps=ps, "Filtering BLAST results")
  fltrBlstRs(blst_rs=blst_rs, ps=ps)
}


#' @name clstrBlstRs
#' @title Cluster BLAST Results
#' @description TODO
#' @param blst_rs BLAST results
#' @param infrmtv Informative? T/F
#' @export
clstrBlstRs <- function(blst_rs, infrmtv=TRUE) {
  g <- igraph::graph.data.frame(blst_rs[ ,c("query.id",
                                            "subject.id")],
                                directed=FALSE)
  clstrs <- igraph::clusters(g)
  # filter for phylogenetically informative clusters
  if(infrmtv) {
    clstrs <- clstrs[['membership']][which(
      clstrs[['membership']] %in% which(clstrs[['csize']] > 2))]
  } else {
    clstrs <- clstrs[['membership']]
  }
  # we will return a list, one entry with sequence IDs
  #  for each cluster
  clstr_lst <- lapply(unique(clstrs), function(x) {
    list(gis=sort(names(clstrs)[which(clstrs==x)]))
  })
  # Get the seed gi, we will chose it to be the sequence
  # in the cluster that has
  # the most hits with the other members in the cluster;
  # i.e. the most connected
  # node in the graph
  degrees <- igraph::degree(g)
  # get seed gis and as field to clusters
  clstr_lst <- lapply(clstr_lst, function(cl){
    idx <- order(degrees[cl[['gis']]], decreasing=TRUE)[1]
    # index of most connected component
    cl[['seed_gi']] <- cl[['gis']][idx]
    cl
  })
  clstr_lst
}

#' @name addClstrInf
#' @title Add taxonomic and sequence info to clusters
#' @description input: List of clstrs
#' @param clstrs Clusters
#' @param phylt_nds PhyLoTa data.frame
#' @param txid Taxonomic node ID
#' @param sqs Sequneces
#' @param drct Direct? T/F
#' @export
addClstrInf <- function(clstrs, phylt_nds, txid, sqs, drct=FALSE) {
  # we will take the column names in the database as names in the list
  for(i in 1:length(clstrs)) {
    cl <- clstrs[[i]]
    cl_sqs <- lapply(cl[['gis']], function(gi)sqs[[as.character(gi)]])
    cl[['ti_root']] <- txid
    cl[['ci']] <- i-1
    cl[['cl_type']] <- ifelse(drct, 'node', 'subtree')
    cl[['n_gi']] <- length(cl[['gis']])
    # all taxon ids for gis in cluster
    cl[['tis']] <- sapply(cl_sqs, '[[', 'ti')
    cl[['n_ti']] <- length(unique(cl[['tis']]))
    # sequence lengths
    l <- sapply(cl_sqs, '[[', 'length')
    cl[['MinLength']] <- min(l)
    cl[['MaxLength']] <- max(l)
    # get n. genera
    cl[['n_gen']] <- length(unique(sapply(cl[['tis']],
                                          getGnsFrmPhyltNds, phylt_nds)))
    # n_child and ci_anc will be set (or incremented) later,
    # when multiple hierarchies are calculated
    cl[['n_child']] <- 0
    cl[['ci_anc']] <- NA
    # make a unique cluster id consisting of seed gi, taxon id,
    # cluster id, cluster type
    unique_id <- paste0(cl[['seed_gi']], '-', txid, '-',
                        cl[['ci']], '-', cl[['cl_type']])
    cl[['unique_id']] <- unique_id
    clstrs[[i]] <- cl
  }
  clstrs
}

#' @name clstrPhylt
#' @title Reformat cluster obj for PhyLoTa table
#' @description Returns a data frame with columns matchine the fields in
#' PhyLoTA's cluster table
#' @param clstrs Clusters
#' @export
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
#' @export
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