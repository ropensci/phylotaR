#' @name read_phylota
#' @title Generate a Phylota object in R
#' @description Creates a Phylota object containing information on
#' clusters, sequences and taxonomy from the working directory of a
#' completed pipeline.
#' @param wd Working directory
#' @return Phylota
#' @export
#' @family tools-public
read_phylota <- function(wd) {
  if (!file.exists(file.path(wd, 'cache',
                            'clusters_sqs.RData'))) {
    msg <- paste0('Data file not found in [', wd,
                  '].\nAre you sure pipeline has completed?')
    stop(msg)
  }
  clstrs_sqs <- obj_load(wd, nm = 'clstrs_sqs')
  clstrs <- clstrs_sqs[['clstrs']]
  sqs <- clstrs_sqs[['sqs']]
  # TODO: how does this occur?
  # prevent dups
  non_dups <- which(!duplicated(sqs@ids))
  sqs <- seqarc_gen(sqs@sqs[non_dups])
  txids <- sort(unique(sqs@txids))
  sids <- sort(unique(sqs@ids))
  cids <- clstrs@ids
  txdct <- obj_load(wd, nm = 'txdct')
  prnt_id <- txdct@prnt
  prnt_nm <- txdct@recs[[prnt_id]]@scnm
  phylota <- new('Phylota', sqs = sqs, clstrs = clstrs, txids = txids,
                 sids = sids, cids = cids, txdct = txdct,
                 prnt_id = prnt_id, prnt_nm = prnt_nm)
  update_phylota(phylota)
}

#' @name write_phylota
#' @title Write out Phylota-like Table
#' @description Create a Phylota-like table from phylota object.
#' @param phylota Phylota
#' @param outfile Output file
#' @return NULL
#' @export
#' @family tools-public
write_table <- function(phylota, outfile) {
  # TODO
}

#' @name write_sqs
#' @title Write out sequences
#' @description Write out sequences, as .fasta, for a given vector of
#' IDs.
#' @param phylota Phylota
#' @param outfile Output file
#' @param sid Sequence ID(s)
#' @param sq_nm Sequence name(s)
#' @details 
#' The user can control the output definition lines of the sequences
#' using the sq_nm. By default sequences IDs are used. Note, ensure
#' the sq_nm are in the same order as sid.
#' @return NULL
#' @export
#' @family tools-public
write_sqs <- function(phylota, outfile,
                      sid, sq_nm = sid) {
  get <- function(i) {
    sq <- phylota@sqs[[sid[i]]]
    paste0('>', sq_nm[i], '\n',
           rawToChar(sq@sq), '\n\n')
  }
  fasta <- vapply(seq_along(sid), get, '')
  fasta <- paste0(fasta, collapse = '')
  write(x = fasta, file = outfile)
}

#' @name plot_phylota_pa
#' @title Plot presence/absence matrix
#' @description Plot presence/absence of taxa by each
#' cluster in phylota object.
#' @param phylota Phylota object
#' @param cids Vector of cluster IDs
#' @param txids Vector of taxonomic IDs
#' @param cnms Cluster names
#' @param txnms Taxonomic names
#' @details Cluster names and taxonomic names
#' can be given to the function, by default IDs are used.
#' @return geom_object
#' @examples
#' library(phylotaR)
#' data(cycads)
#' # drop all but first ten
#' cycads <- drop_cls(cycads, cycads@cids[1:10])
#' # plot all
#' p <- plot_phylota_pa(phylota = cycads, cids = cycads@cids,
#'                      txids = cycads@txids)
#' print(p)  # lots of information, difficult to interpret
#' # get genus-level taxonomic names
#' genus_txids <- get_txids(cycads,
#'                          txids = cycads@txids,
#'                          rnk = 'genus')
#' genus_txids <- unique(genus_txids)
#' # dropping missing
#' genus_txids <- genus_txids[genus_txids ! =  '']
#' genus_nms <- get_tx_slot(cycads, genus_txids, slt_nm = 'scnm')
#' # make alphabetical for plotting
#' genus_nms <- sort(genus_nms, decreasing = TRUE)
#' # generate geom_object
#' p <- plot_phylota_pa(phylota = cycads, cids = cycads@cids,
#'                      txids = genus_txids, txnms = genus_nms)
#' # plot
#' print(p)  # easier to interpret
#' @export
#' @family tools-public
plot_phylota_pa <- function(phylota, cids, txids, cnms = cids,
                            txnms = txids) {
  mkdata <- function(cid) {
    clstr <- phylota@clstrs@clstrs[[which(cid == phylota@cids)]]
    value <- apply(X = tis_mtrx[ ,cl@sids], 1, any)
    value <- as.numeric(value)
    data.frame(txid = as.character(txids), cid = cid, value = value)
  }
  value <- NULL
  # gen p_data
  tis_mtrx <- mk_txid_in_sq_mtrx(phylota = phylota,
                                 txids = txids)
  p_data <- lapply(cids, mkdata)
  p_data <- do.call('rbind', p_data)
  p_data[['cnm']] <- cnms[match(p_data[['cid']], cids)]
  p_data[['txnm']] <- txnms[match(p_data[['txid']], txids)]
  # reorder
  p_data[['cnm']] <- factor(p_data[['cnm']], levels = cnms,
                            ordered = TRUE)
  p_data[['txnm']] <- factor(p_data[['txnm']], levels = txnms,
                             ordered = TRUE)
  # plot
  # https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
  cnm <- txnm <- NULL
  ggplot2::ggplot(p_data, ggplot2::aes(cnm, txnm)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value)) +
    ggplot2::xlab('') + ggplot2::ylab('') +
    ggplot2::scale_fill_gradient(low = '#e0e0e0',
                                 high = '#303030') +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = 'none')
}

#' @name plot_phylota_treemap
#' @title Plot treemap of Phylota object
#' @description Treemaps show relative size with boxes. The user can
#' explore which taxa or clusters are most represented either by
#' sequence or cluster number. If cluster IDs are provided, the plot
#' is made for clusters. If taxonomic IDs are provided, the plot is
#' made for taxa.
#' @param phylota Phylota object
#' @param cids Cluster IDs
#' @param txids Taxonomic IDs
#' @param cnms Cluster names
#' @param txnms Taxonomic names
#' @param with_labels Show names per box?
#' @param area What determines the size per box?
#' @param fill What determines the coloured fill per box?
#' @details The function can take a long time to run for large Phylota
#' objects over many taxonomic IDs because searches are made across
#' lineages. The idea of the function is to assess the data dominance
#' of specific clusters and taxa.
#' @return geom_object
#' @examples 
#' data("tinamous")
#' # Plot clusters, size by n. sq, fill by n. tx
#' p <- plot_phylota_treemap(phylota = tinamous, cids = tinamous@cids,
#'                           area = 'nsq', fill = 'ntx')
#' print(p)
#' # Plot taxa, size by n. sq, fill by ncl
#' txids <- get_txids(tinamous, txids = tinamous@txids,
#'                    rnk = 'genus')
#' txids <- txids[txids ! =  '']
#' txids <- unique(txids)
#' txnms <- get_tx_slot(tinamous, txids, slt_nm = 'scnm')
#' p <- plot_phylota_treemap(phylota = tinamous, txids = txids,
#'                           txnms = txnms, area = 'nsq',
#'                           fill = 'ncl')
#' print(p)
#' @export
#' @family tools-public
plot_phylota_treemap <- function(phylota, cids = NULL, txids = NULL,
                                 cnms = cids, txnms = txids,
                                 with_labels = TRUE,
                                 area = c('ntx', 'nsq', 'ncl'),
                                 fill = c('NULL', 'typ', 'ntx',
                                        'nsq', 'ncl')) {
  mkcldata <- function(i) {
    clstr <- phylota@clstrs@clstrs[[which(cids[i] == phylota@cids)]]
    data.frame(id = cids[i], label = cnms[i], nsq = clstr@nsqs,
               ntx = clstr@ntx, typ = clstr@typ)
  }
  getcltxids <- function(i) {
    clstr <- phylota@clstrs@clstrs[[i]]
    res <- clstr@txids
    names(res) <- NULL
    res
  }
  mktxdata <- function(i) {
    anycltxids <- function(ids) {
      any(ids %in% all_ids)
    }
    ads <- descendants_get(id = txids[i], txdct = phylota@txdct,
                           direct = FALSE)
    all_ids <- c(ads, txids[i])
    nsq <- sum(sids_txids %in% all_ids)
    nclstr <- sum(vapply(cids_txids, anycltxids, logical(1)))
    data.frame(id = txids[i], label = txnms[i], nsq = nsq,
               ncl = nclstr)
  }
  area <- match.arg(area)
  fill <- match.arg(fill)
  if (!is.null(cids)) {
    p_data <- lapply(seq_along(cids), mkcldata)
  } else if (!is.null(txids)) {
    sids_txids <- get_sq_slot(phylota, sid = phylota@sids,
                              slt_nm = 'txid')
    cids_txids <- lapply(seq_along(phylota@cids), getcltxids)
    p_data <- lapply(seq_along(txids), mktxdata)
  } else {
    stop('Either cids or txids not provided.')
  }
  p_data <- do.call('rbind', p_data)
  p <- ggplot2::ggplot(p_data, ggplot2::aes_string(area = area,
                                                   fill = fill,
                                                   subgroup = 'label',
                                                   label = 'label')) +
    treemapify::geom_treemap()
  if (with_labels) {
    p <- p + treemapify::geom_treemap_subgroup_text()
  }
  p
}

#' @name mk_txid_in_sq_mtrx
#' @title Return matrix of txid in sequence
#' @description Searches through lineages of sequences' source
#' organisms to determine whether each txid is represented by the
#' sequence.
#' @param phylota Phylota
#' @param txids Taxonomic IDs
#' @param sids Sequence IDs
#' @return matrix
#' @family tools-private
mk_txid_in_sq_mtrx <- function(phylota, txids,
                               sids = phylota@sids) {
  is_txid_in_sqs <- function(txid) {
    ads <- descendants_get(id = txid, txdct = phylota@txdct,
                           direct = FALSE)
    all_ids <- c(ads, txid)
    sids_txids %in% all_ids
  }
  sids_txids <- get_sq_slot(phylota, sid = sids, slt_nm = 'txid')
  tis_mtrx <- do.call('rbind', lapply(txids, is_txid_in_sqs))
  colnames(tis_mtrx) <- sids
  rownames(tis_mtrx) <- txids
  tis_mtrx
}

#' @name is_txid_in_sq
#' @title Is txid in sequence?
#' @description Checks if given txid is represented by sequence by
#' looking at sequence source organism's lineage.
#' @param phylota Phylota
#' @param txid Taxonomic ID
#' @param sid Sequence ID
#' @return boolean
#' @examples
#' data(tinamous)
#' sid <- tinamous@sids[[1]]
#' sq <- tinamous[[sid]]
#' txid <- sq@txid
#' # expect true
#' is_txid_in_sq(phylota = tinamous, txid = txid, sid = sid)
#' @export
#' @family tools-public
is_txid_in_sq <- function(phylota, txid, sid) {
  sq <- phylota@sqs@sqs[[which(sid == phylota@sids)]]
  sq_tx <-  phylota@txdct@recs[[sq@txid]]
  txid %in% sq_tx@lng[['ids']]
}

#' @name is_txid_in_clstr
#' @title Is txid in cluster?
#' @description Checks if given txid is represented by any of the
#' sequences of a cluster by searching through all the sequence search
#' organism lineages.
#' @param phylota Phylota
#' @param txid Taxonomic ID
#' @param cid Cluster ID
#' @return boolean
#' @examples
#' data(tinamous)
#' cid <- tinamous@cids[[1]]
#' clstr <- tinamous[[cid]]
#' sq <- tinamous[[cl@sids[[1]]]]
#' txid <- sq@txid
#' # expect true
#' is_txid_in_cl(phylota = tinamous, txid = txid, cid = cid)
#' @export
#' @family tools-public
is_txid_in_clstr <- function(phylota, txid, cid) {
  clstr <- phylota@clstrs@clstrs[[which(cid == phylota@cids)]]
  bool <- vapply(clstr@sids, is_txid_in_sq, logical(1),
                 txid = txid, phylota = phylota)
  any(bool)
}

#' @name summary_phylota
#' @title Summarise clusters in Phylota Table
#' @description Generates a summary data.frame from all clusters in
#' Phylota object.
#' @param phylota Phylota object
#' @family tools-private
summary_phylota <- function(phylota) {
  print_frq_wrds <- function(wrd_prps, max_n = 2) {
    if (length(wrd_prps) == 0) {
      return('-')
    }
    n <- ifelse(length(wrd_prps) > max_n, max_n,
                length(wrd_prps))
    wrd_prps <- wrd_prps[1:n]
    wrd_prps <- signif(wrd_prps, digits = 1)
    res <-  paste0(names(wrd_prps), ' (', wrd_prps, ')')
    paste0(res, collapse = ', ')
  }
  get_row <- function(cid) {
    clstr <- phylota@clstrs[[cid]]
    mad_scr <- calc_mad(phylota = phylota, cid = cid)
    sqlngth <- median(get_sq_slot(phylota = phylota,
                                  cid = cid,
                                  slt_nm = 'nncltds'))
    dflns <- calc_wrdfrq(phylota = phylota, cid = cid,
                         type = 'dfln', min_frq = 0)[[1]]
    dflns <- print_frq_wrds(dflns)
    ftr_nms <- calc_wrdfrq(phylota = phylota, cid = cid,
                           type = 'nm', min_frq = 0)[[1]]
    ftr_nms <- print_frq_wrds(ftr_nms)
    res <- c(clstr@id, clstr@typ, clstr@seed, clstr@prnt,
             length(unique(clstr@txids)), length(clstr@sids),
             sqlngth, mad_scr, dflns, ftr_nms)
    res
  }
  res <- lapply(phylota@cids, get_row)
  res <- matrix(unlist(res), nrow = length(phylota@cids),
                byrow = TRUE)
  colnames(res) <- c('ID', 'Type', 'Seed', 'Parent',
                     'N_taxa', 'N_seqs', 'Med_sql', 'MAD',
                     'Definition', 'Feature')
  res <- data.frame(res, stringsAsFactors = FALSE)
  res[['N_taxa']] <- as.integer(res[['N_taxa']])
  res[['N_seqs']] <- as.integer(res[['N_seqs']])
  res[['Med_sql']] <- as.numeric(res[['Med_sql']])
  res[['MAD']] <- as.numeric(res[['MAD']])
  res
}

#' @name update_phylota
#' @title Update slots
#' @description After change, run to update slots.
#' @param phylota Phylota
#' @family tools-private
#' @return Phylota
update_phylota <- function(phylota) {
  get_sids <- function(i) {
    clstr <- clstrs@clstrs[[i]]
    clstr@sids
  }
  get_seeds <- function(i) {
    clstr <- clstrs@clstrs[[i]]
    clstr@seed
  }
  clstrs <- phylota@clstrs # clstrs informs all other slots
  all_sids <- lapply(seq_along(clstrs@ids), get_sids)
  all_seeds <- vapply(seq_along(clstrs@ids), get_seeds, '')
  all_sids <- unique(c(unlist(all_sids), all_seeds))
  sqs <- phylota@sqs
  sqs <- sqs[all_sids]
  phylota@txids <- unique(sqs@txids)
  phylota@sids <- sqs@ids
  phylota@sqs <- sqs
  initialize(phylota)
}
