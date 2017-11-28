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
  cat("Taxid", taxid, ": writing", nrow(cldf), "clusters,",
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
  cat("Finished processing taxid ", taxid, " # ", i, " / ",
      length(txids), "\n")
}