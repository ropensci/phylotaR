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

#' @name clstr
#' @title Cluster sequences
#' @description Identify sequence clusters using BLAST
#' @param phylt_nds PhyLoTa nodes data.frame
#' @export
# TODO: a bit big, break down
clstr <- function(wd, txid, sqs, phylt_nds,
                  drct=FALSE, infrmtv=FALSE,
                  blst_rs=NULL, verbose=FALSE) {
  clstr_typ <- ifelse(drct, "direct", "subtree")
  rnks <- as.character(phylt_nds[['rank']])
  rnk <- as.character(rnks[match(txid, phylt_nds[['ti']])])
  .cp(v=verbose, "Processing taxid [", txid, "] of rank [",
      rnk, "] attempting to make [", clstr_typ, "] clusters")
  if(!drct) {
    dds <- getDDFrmPhyltNds(txid=txid, phylt_nds=phylt_nds)
    if(length(dds) > 0) {
      all_clstrs <- clstr(wd, txid, phylt_nds=phylt_nds,
                          drct=TRUE, infrmtv=infrmtv,
                          blst_rs=blst_rs, verbose=verbose)
    }
  }
  if(drct) {
    txids <- txid
  } else {
    txids <- getADs(txid=txid, phylt_nds=phylt_nds)
  }
  all_sq_txids <- sapply(sqs, '[[', 'ti')
  sq_txids <- names(sqs[which(all_sq_txids %in% txids)])
  .cp(v=verbose, "Found [", length(sq_txids), '] [',
      clstr_typ, "] GIs for taxid [", txid, "]")
  if(length(sq_txids) == 0) {
    .cp(v=verbose, "No [", clstr_typ, "] sequences for taxid [",
        txid, "], cannot make clusters")
    return(all_clstrs)
  }
  sqs_prt <- sqs[sq_txids[sq_txids %in% names(sqs)]]
  txids_blst_rs <- list()
  if(is.null(blst_rs)) {
    .cp(v=verbose, "Current BLAST result is NULL from parent, performing BLAST")
    txids_blst_rs <- blstSqs(wd=wd, txid=txid, typ=clstr_typ, sqs=sqs_prt,
                             verbose=verbose)
    # TODO: Can we do this more elegantly?
    if(is.null(txids_blst_rs)) {
      .cp(v=verbose, "Current BLAST result is NULL after blasting")
      return(all_clstrs)
    }
  } else {
    .cp(v=verbose, "Using BLAST results from parent cluster")
    txids_blst_rs <- blst_rs[which(blst_rs$subject.id %in% gis),]
    txids_blst_rs <- txids_blst_rs[which(txids_blst_rs$query.id %in% gis),]
  }
  # Sq clusters are stored in a list of lists, named by gi
  # Get top-level clusters
  .cp(v=verbose, "Clustering BLAST results for taxid [", txid, "]")
  # TODO
  raw_clstrs <- cluster.blast.results(txids_blst_rs,
                                      informative=informative)
  .cp(v=verbose, "Finished clustering BLAST results for taxid [",
      txid, "]")
  # make dataframe with fields as in PhyLoTa database
  # TODO
  clstrs <- .add.cluster.info(raw_clstrs, nodes, txid, seqs,
                              direct=direct)
  .cp(v=verbose, "Generated [", length(clstrs),
      "] clusters")
  all_clstrs <- c(all_clstrs, clstrs)
  # If we do not calculate a direct cluster here, iterate over
  # txid's children to retrieve clusters
  if(!drct) {
    dds <- getDDFrmPhyltNds(txid=txid, phylt_nds=phylt_nds)
    for(dd in dds) {
      .cp(v=verbose, "Processing child taxon of [", txid,
          "] [", dd, "]")
      # TODO
      dd_clstrs <- clstr(dd, nodes, seqs, txids_blst_rs)
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
writeClstr <- function(phylt_nds) {
  ## Write all data to file
  cat("Taxid", txid, ": writing", nrow(cldf), "clusters,",
      nrow(seqdf), "sequences,", nrow(cigidf), "ci_gi entries to file\n")
  write.table(cldf, file=clusters.file, append=file.exists(clusters.file),
              quote=FALSE, sep="\t", row.names=FALSE,
              col.names=!file.exists(clusters.file))
  write.table(seqdf, file=seqs.file, append=file.exists(seqs.file),
              quote=FALSE, sep="\t", row.names=FALSE,
              col.names=!file.exists(seqs.file))
  write.table(cigidf, file=ci_gi.file, append=file.exists(ci_gi.file),
              quote=FALSE, sep="\t", row.names=FALSE,
              col.names=!file.exists(ci_gi.file))
  cat("Finished processing taxid ", txid, " # ", i, " / ",
      length(txids), "\n")
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
blstSqs <- function(txid, typ, sqs, wd, verbose=FALSE) {
  .cp(v=verbose, "BLAST all vs all for [",
      length(sqs), "] sequences")
  dbfl <- paste0('taxon-', txid, '-typ-', typ,
                 '-db.fa')
  outfl <- paste0('taxon-', txid, '-typ-', typ,
                  '-blastout.txt')
  mkBlstDB(sqs, dbfl=dbfl, wd=wd, verbose=verbose)
  blst_rs <- blstN(dbfl=dbfl, outfl=outfl, wd=wd,
                   verbose=verbose)
  # TODO: Not so elegant
  if(is.null(blst_rs)) {
    return(NULL)
  }
  .cp(v=verbose, "Number of BLAST results [", nrow(blst_rs), "]")
  .cp(v=verbose, "Filtering BLAST results")
  fltrBlstRs(blst_rs=blst_rs, min_cvrg=0.51, verbose=verbose)
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
    clstrs <- clstrs$membership[which(
      clstrs$membership %in% which(clstrs$csize > 2))]
  } else {
    clstrs <- clstrs$membership
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
