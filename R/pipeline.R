
#' @name setUp
#' @title Set-up parameters
#' @description 
#' @param wd Working directory
#' @param txid Root taxonomic ID(s), vector or numeric
#' @param ncbi_dr Directory to NCBI BLAST tools, default '.'
#' @param verbose Verbose, T/F
#' @param ... Additional parameters
#' @export
#' @seealso \code{\link{setUp}}
setUp <- function(wd, txid, ncbi_dr='.', verbose=FALSE,
                  ...) {
  # header log
  msg <- paste0('phylotaR: Implementation of PhyLoTa in R [v',
                packageVersion('phylotaR'), ']')
  brdr <- paste0(rep('-', nchar(msg)), collapse='')
  msg <- paste0(brdr, '\n', msg, '\n', brdr)
  .log(v=verbose, wd=wd, lvl=1, msg)
  # set up
  ncbi_execs <- setUpNcbiTools(d=ncbi_dr,
                               verbose=verbose,
                               wd=wd)
  setUpPrmtrs(wd=wd, txid=txid,
              ncbi_execs=ncbi_execs, ...)
  # end
  .log(v=verbose, wd=wd, lvl=1, brdr)
}

#' @name run
#' @title Run PhyLoTa pipeline
#' @description Run the entire PhyLoTa pipeline.
#' All generated files will be stored in the wd.
#' The process can be stopped at anytime and 
#' restarted with \code{restart}.
#' \code{nstages} must be a numeric value representing
#' the number of stages that will be run. Stages are run
#' in the following order:  1 - taxise, 2 - download,
#' 3 - cluster and 4 - align.
#' For example, specifying \code{nstages} = 3, will run
#' taxise, download and cluster.
#' Stages can also be run individually, see linked
#' functions below.
#' @param wd Working directory
#' @param nstages number of stages to run
#' @export
#' @seealso \code{\link{restart}}, \code{\link{runTaxise}},
#' \code{\link{runDownload}}, \code{\link{runClusters}},
#' \code{\link{runAlign}}
run <- function(wd, nstages=4) {
  # TODO: save progress
  # header log
  verbose <- ldPrmtrs(wd)[['verbose']]
  msg <- paste0('Running pipeline on [', .Platform$OS.type, '] at [',
                Sys.time(), ']')
  brdr <- paste0(rep('-', nchar(msg)), collapse='')
  msg <- paste0(brdr, '\n', msg, '\n', brdr)
  info(v=verbose, wd=wd, lvl=1, msg)
  if(nstages < 1) {
    stop('`nstages` is less than 1.')
  }
  if(nstages > 4) {
    stop('`nstages` is greater than 4.')
  }
  if(nstages >= 1) {
    # Generate taxonomic 'nodes'
    runTaxise(wd)
  }
  if(nstages >= 2) {
    # Download sequences
    runDownload(wd)
  }
  if(nstages >= 3) {
    # Generate clusters
    runClusters(wd)
  }
  if(nstages == 4) {
    # Generate alignments
    runAlign(wd)
  }
  # footer log
  msg <- paste0('\nCompleted pipeline on at [', Sys.time(), ']')
  brdr <- paste0(rep('-', nchar(msg)), collapse='')
  msg <- paste0(brdr, '\n', msg, '\n', brdr)
  info(v=verbose, wd=wd, lvl=1, msg)
}

#' @name restart
#' @title Restart a PhyLoTa pipeline run
#' @description Restarts the running of a pipeline
#' as started with \code{run}.
#' @param wd Working directory
#' @export
#' @seealso \code{\link{run}}
restart <- function(wd) {
  # TODO
}

#' @name reset
#' @title Reset a PhyLoTa pipeline run
#' @description Resets the pipeline to a specified stage.
#' @param wd Working directory
#' @export
#' @seealso \code{\link{run}}
reset <- function(wd) {
  # TODO
}