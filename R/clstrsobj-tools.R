
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
  new('ClstrsObj', sqs=clstrs_sqs[['sqs']],
      clstrs=clstrs)
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