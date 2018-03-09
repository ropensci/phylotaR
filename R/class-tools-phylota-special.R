
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
#' @title Write out PhyLoTa-like Table
#' @description Create a PhyLoTa-like table
#' from phylota object.
#' @param phylota Phylota object
#' @param outfile Output file
#' @return NULL
#' @export
write_table <- function(phylota, outfile) {
  # TODO
}

#' @name write_sqs
#' @title Write out sequences
#' @description Write out sequences, as .fasta,
#' for a given vector of IDs.
#' @param phylota Phylota object
#' @param outfile Output file
#' @param sid Sequence ID(s)
#' @param sq_nm Sequence name(s)
#' @details 
#' The user can control the output definition
#' lines of the sequences using the sq_nm.
#' By default sequences IDs are used.
#' Note, ensure the sq_nm are in the same
#' order as sid.
#' @return NULL
#' @export
write_sqs <- function(phylota, outfile,
                          sid, sq_nm=sid) {
  get <- function(i) {
    sq <- phylota@sqs[[sid[i]]]
    paste0('>', sq_nm[i], '\n',
           rawToChar(sq@sq), '\n\n')
  }
  fasta <- vapply(seq_along(sid), get, '')
  fasta <- paste0(fasta, collapse='')
  write(x=fasta, file=outfile)
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
    mad_scr <- calc_mad(phylota=phylota, cid=cid)
    sqlngth <- median(get_sq_slot(phylota=phylota,
                                  cid=cid,
                                  slt_nm='nncltds'))
    dflns <- calc_wrdfrq(phylota=phylota, cid=cid,
                         type='dfln', min_frq=0)[[1]]
    dflns <- print_frq_wrds(dflns)
    ftr_nms <- calc_wrdfrq(phylota=phylota, cid=cid,
                           type='nm', min_frq=0)[[1]]
    ftr_nms <- print_frq_wrds(ftr_nms)
    res <- c(cl@id, cl@typ, cl@seed, cl@prnt,
             length(unique(cl@txids)), length(cl@sids),
             sqlngth, mad_scr, dflns, ftr_nms)
    res
  }
  res <- lapply(phylota@cids, get_row)
  res <- matrix(unlist(res), nrow=length(phylota@cids),
                byrow=TRUE)
  colnames(res) <- c('ID', 'Type', 'Seed', 'Parent',
                     'N_taxa', 'N_seqs', 'Med_sql', 'MAD',
                     'Definition', 'Feature')
  res <- data.frame(res, stringsAsFactors=FALSE)
  res[['N_taxa']] <- as.integer(res[['N_taxa']])
  res[['N_seqs']] <- as.integer(res[['N_seqs']])
  res[['Med_sql']] <- as.numeric(res[['Med_sql']])
  res[['MAD']] <- as.numeric(res[['MAD']])
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
  phylota@txids <- unique(sqs@txids)
  phylota@sids <- sqs@ids
  phylota@sqs <- sqs
  initialize(phylota)
}