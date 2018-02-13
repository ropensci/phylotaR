#' @importFrom methods is new
#' @importFrom stats median 
#' @importFrom utils head packageVersion read.table tail write.table

#' @name parameters
#' @title Default parameters
#' @description Returns a parameter list with default
#' parameter values.
#' @param wd The working directory where all output files are saved.
#' @param txid Taxonomic group of interest, allows vectors.
#' @param mkblstn File path to mkblastdb
#' @param blstn File path to blastn
#' @param tdpth File path to taxdump.tar.gz.
#' If NULL, file will be downloaded automatically.
#' @param v Print progress statements to console?
#' Statements will always be printed to log.txt.
#' @param ncps The number of CPUs to run clustering in parallel.
#' If 1, no parallelisation performed.
#' @param mxnds The maximum number of nodes descending from a taxonomic group.
#' If there are more than this number, nodes at the lower taxonomic level are analysed.
#' @param mdlthrs 'Model organism threshold'. Taxa with more sequences than this number
#' will be considered model organisms and a random mdlthrs subset of their sequences will
#' be downloaded.
#' @param mnsql The minimum length of sequence in nucleotide base pairs to download.
#' @param mxsql The maximum length of sequence in nucleotide base pairs to download.
#' Any longer sequences will be ignored.
#' @param mxrtry The maximum number of attempts to make when downloading.
#' @param mxsqs The maximum number of sequences to BLAST in all-vs-all searches.
#' If there are more sequences for a node, BLAST is performed at the lower taxonomic level.
#' @param mxevl The maximum E-value for a successful BLAST.
#' @param mncvrg The maximum percentile coverage defining an overlapping BLAST hit.
#' Sequences with BLAST matches with lower values are not considered orthologous.
#' @export
parameters <-function(wd='.', txid=numeric(),
                      mkblstdb='', blstn='',
                      v=FALSE, ncps=1, tdpth=NULL,
                      mxnds=100000, mdlthrs=3000,
                      mnsql=250, mxsql=7500, 
                      mxrtry=100, mxsqs=10000,
                      mxevl=1.0e-10, mncvrg=51) {
  ps <- as.list(environment())
  ps[['wt_tms']] <- c(1, 3, 6, 10, 60, 300)
  ps
}
