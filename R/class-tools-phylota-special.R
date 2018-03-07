
#' @name read_phylota
#' @title Generate a PhyLoTa object in R
#' @description Creates a PhyLoTa object containing
#' information on clusters, sequences and taxonomy
#' from the working directory of a completed pipeline.
#' @param wd Working directory
#' @return phylota
#' @export
read_phylota <- function(wd) {
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

#' @name write_phylota
#' @title Write out phylota
#' @description Write out phylota object either
#' as .fasta or PhyLoTa table.
#' @param phylota Phylota object
#' @param type PhyLoTa or fasta?
#' @param cls_names Cluster names
#' @param sqs_names Sequence names
#' @return NULL
#' @export
write_phylota <- function(phylota, type=c('fasta', 'phylota'),
                          cls_names, sqs_names) {
  
}

plot_phylota <- function(phylota){
  # TODO
}

summary_phylota <- function(phylota) {
  # TODO
}

update_phylota <- function(phylota) {
  # TODO
  # Make sure seed sqs are not lost
  phylota
}