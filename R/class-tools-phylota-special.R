
#' @name read_phylota
#' @title Generate a PhyLoTa object in R
#' @description Creates a PhyLoTa object containing
#' information on clusters, sequences and taxonomy
#' from the working directory of a completed pipeline.
#' @param wd Working directory
#' @return phylota
#' @export
read_phylota <- function(wd) {
  if(!file.exists(file.path(wd, 'cache',
                            'clstrs_sqs.RData'))) {
    msg <- paste0('Data file not found in [', wd,
                  '].\nAre you sure pipeline has completed?')
    stop(msg)
  }
  clstrs_sqs <- ldObj(wd, nm='clstrs_sqs')
  cls <- clstrs_sqs[['clstrs']]
  sqs <- clstrs_sqs[['sqs']]
  # TODO: how does this occur?
  # prevent dups
  non_dups <- which(!duplicated(sqs@ids))
  sqs <- genSqRcrdBx(sqs@sqs[non_dups])
  txids <- sort(unique(sqs@txids))
  sids <- sort(unique(sqs@ids))
  cids <- cls@ids
  txdct <- ldObj(wd, nm='txdct')
  phylota <- new('PhyLoTa', sqs=sqs, cls=cls,
                 txids=txids, sids=sids,
                 cids=cids, txdct=txdct)
  update_phylota(phylota)
}

#' @name write_phylota
#' @title Write out phylota
#' @description Write out phylota object either
#' as .fasta or PhyLoTa table.
#' @param phylota Phylota object
#' @param drcty Output directory
#' @param type PhyLoTa or fasta?
#' @param cls_names Cluster names
#' @param sqs_names Sequence names
#' @details User can either output all the clusters
#' in the phylota object as a PhyLoTa-style table,
#' or as each of the clusters' sequences as separate
#' .fasta files.
#' The user can control the output names of the sequences
#' and clusters using the cls_names and sqs_names arguments.
#' By default, the IDs for the clusters and sequences
#' are used as names.
#' If PhyLoTa, the output file will be called 'phylota.csv'.
#' If fasta, the output files will each be named after the
#' cls_names argument.
#' Note, ensure the cls_names and sqs_names are in the same
#' order as cids and sids.
#' @return NULL
#' @export
write_phylota <- function(phylota, drctry='.',
                          type=c('fasta', 'phylota'),
                          cls_names=phylota@cids,
                          sqs_names=phylota@sids) {
  get_fasta <- function(sid) {
    sq <- sqs[[i]]
    dfln <- dflns[[i]]
    
  }
  type <- match.arg(type)
  if(type == 'phylota') {
    stop('PhyLoTa table output not yet available.')
  }
  fasta <- ''
  for(i in seq_along(phylota@cids)) {
    cid <- phylota@cids[i]
    cl <- phylota@cls[[cid]]
    sqs <- cl@sqs
    for(j in seq_along(sqs)) {
      sid <- cl@sids[i]
      sq <- sqs[[sid]]
      fasta <- paste0(fasta, '>', sqs_names[i],
                      '\n', rawToChar(sq@sq),
                      '\n\n')
    }
    otpt_fl <- paste0(cid, '.fasta')
    write(x=fasta, file=file.path(drctry, cid))
  }
}

plot_phylota <- function(phylota){
  # TODO
}

#' @name summary_phylota
#' @title Summarise clusters in PhyLoTa Table
#' @description Generates a summary data.frame
#' from all clusters in PhyLoTa object.
#' @param phylota PhyLoTa object
summary_phylota <- function(phylota) {
  print_frq_wrds <- function(wrd_prps, max_n=2) {
    if(length(wrd_prps) == 0) {
      return('-')
    }
    n <- ifelse(length(wrd_prps) > max_n, max_n,
                length(wrd_prps))
    wrd_prps <- wrd_prps[1:n]
    wrd_prps <- signif(wrd_prps, digits=1)
    res <-  paste0(names(wrd_prps), ' (',
                   wrd_prps, ')')
    paste0(res, collapse=', ')
  }
  get_row <- function(cid) {
    cl <- phylota@cls[[cid]]
    dflns <- calc_wrdfrq(phylota, cid, type='dfln',
                         min_frq=0)[[1]]
    dflns <- print_frq_wrds(dflns)
    ftr_nms <- calc_wrdfrq(phylota, cid, type='nm',
                           min_frq=0)[[1]]
    ftr_nms <- print_frq_wrds(ftr_nms)
    res <- c(cl@id, cl@typ, cl@seed, #cl@prnt,
             length(unique(cl@txids)), length(cl@sids),
             dflns, ftr_nms)
    res
  }
  res <- lapply(phylota@cids, get_row)
  res <- matrix(unlist(res), nrow=length(phylota@cids),
                byrow=TRUE)
  colnames(res) <- c('ID', 'Type', 'Seed', #'Parent',
                     'N_taxa', 'N_seqs', 'Definition',
                     'Feature')
  res <- data.frame(res, stringsAsFactors=FALSE)
  res[['N_taxa']] <- as.integer(res[['N_taxa']])
  res[['N_seqs']] <- as.integer(res[['N_seqs']])
  res
}

#' @name update_phylota
#' @title Update slots
#' @description After change, run to update
#' slots.
#' @param phylota Phylota object
update_phylota <- function(phylota) {
  get_sids <- function(i) {
    cl <- cls@cls[[i]]
    cl@sids
  }
  get_seeds <- function(i) {
    cl <- cls@cls[[i]]
    cl@seed
  }
  cls <- phylota@cls # cls informs all other slots
  all_sids <- lapply(seq_along(cls@ids), get_sids)
  all_seeds <- vapply(seq_along(cls@ids), get_seeds, '')
  all_sids <- unique(c(unlist(all_sids), all_seeds))
  sqs <- phylota@sqs
  sqs <- sqs[all_sids]
  # TODO: update TxDct
  # txids <- sqs@txids
  phylota@sids <- sqs@ids
  phylota@sqs <- sqs
  initialize(phylota)
}