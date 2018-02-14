
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
setUp <- function(wd, txid, ncbi_dr='.', v=FALSE,
                  ...) {
  # header log
  msg <- paste0('phylotaR: Implementation of PhyLoTa in R [v',
                packageVersion('phylotaR'), ']')
  brdr <- paste0(rep('-', nchar(msg)), collapse='')
  msg <- paste0(brdr, '\n', msg, '\n', brdr, '\n')
  .log(v=v, wd=wd, msg)
  # set up
  ncbi_execs <- setUpNcbiTools(d=ncbi_dr, v=v,
                               wd=wd)
  setUpPrmtrs(wd=wd, txid=txid,
              ncbi_execs=ncbi_execs, v=v, ...)
  # end
  .log(v=v, wd=wd, paste0(brdr, '\n'))
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
#' @param nstages Number of total stages to run, max 4.
#' @export
#' @seealso \code{\link{restart}}, \code{\link{runTaxise}},
#' \code{\link{runDownload}}, \code{\link{runClusters}},
#' \code{\link{runAlign}}
run <- function(wd, nstages=4) {
  stgs_msg <- chckStgs(frm=1, to=nstages)
  runStgs(wd=wd, frm=1, to=nstages, stgs_msg=stgs_msg)
}

#' @name runStgs
#' @title Sequentially run each stage
#' @description Runs stages from \code{frm} to \code{to}.
#' Records stage progress in cache.
#' @param wd Working directory
#' @param to Total number of stages to run
#' @param frm Starting stage to run from
#' @param stgs_msg Printout stage message for log
#' @param rstrt Restarting, T/F
#' @export
runStgs <- function(wd, to, frm, stgs_msg, rstrt=FALSE) {
  .run <- function() {
    if(frm <= 1 & to >= 1) {
      if(!ps[['v']]) {
        cat('... Taxise\n')
      }
      runTaxise(wd)
      svPrgrss(wd, 'taxise')
    }
    if(frm <= 2 & to >= 2) {
      if(!ps[['v']]) {
        cat('... Download\n')
      }
      runDownload(wd)
      svPrgrss(wd, 'download')
    }
    if(frm <= 3 & to >= 3) {
      if(!ps[['v']]) {
        cat('... Cluster\n')
      }
      runClusters(wd)
      svPrgrss(wd, 'cluster')
    }
    if(frm <= 4 & to >= 4) {
      if(!ps[['v']]) {
        cat('... Cluster2\n')
      }
      runClusters2(wd)
      svPrgrss(wd, 'align')
    }
  }
  ps <- ldPrmtrs(wd)
  # header log
  if(rstrt) {
    msg <- paste0('Restarting pipeline on [', .Platform$OS.type,
                  '] at [', Sys.time(), ']')
  } else {
    msg <- paste0('Running pipeline on [', .Platform$OS.type,
                  '] at [', Sys.time(), ']')
  }
  brdr <- paste0(rep('-', nchar(msg)), collapse='')
  msg <- paste0(brdr, '\n', msg, '\n', brdr)
  info(ps=ps, lvl=1, msg)
  info(ps=ps, lvl=1, stgs_msg)
  errmsg <- try(.run(), silent=TRUE)
  if('try-error' %in% is(errmsg)) {
    # ctrl+c
    if(grepl('Operation was aborted by an application callback',
             errmsg[[1]])) {
      msg <- paste0('---- Halted by user [', Sys.time(),
                    '] ----')
      .log(v=ps[['v']], wd=ps[['wd']], msg)
      stop(msg)
    }
    # unexpected pipeline error
    msg <- paste0('Unexpected ', errmsg[[1]], '\n',
                  'Occurred [', Sys.time(), ']\n',
                  'Contact package maintainer for help.\n')
    .log(v=ps[['v']], wd=ps[['wd']], msg)
    stop(msg)
  }
  # footer log
  msg <- paste0('Completed pipeline at [', Sys.time(), ']')
  brdr <- paste0(rep('-', nchar(msg)), collapse='')
  msg <- paste0(brdr, '\n', msg, '\n', brdr)
  info(ps=ps, lvl=1, msg)
}

#' @name restart
#' @title Restart a PhyLoTa pipeline run
#' @description Restarts the running of a pipeline
#' as started with \code{run}.
#' @param wd Working directory
#' @param nstages Number of total stages to run, max 4.
#' @export
#' @seealso \code{\link{run}}
restart <- function(wd, nstages=4) {
  stg <- rdPrgrss(wd)
  if(is.na(stg)) {
    stop('Pipeline already complete. Use `reset()` to re-run pipeline.')
  }
  frm <- which(c('taxise', 'download', 'cluster', 'cluster2')
                %in% stg)
  if(frm > nstages) {
    stop('Pipeline has already completed [', nstages,
         '] stages. Increase `nstages`.')
  }
  stgs_msg <- chckStgs(frm=frm, to=nstages)
  runStgs(wd=wd, frm=frm, to=nstages, stgs_msg=stgs_msg, rstrt=TRUE)
}

#' @name chckStgs
#' @title Check to and frm stage numbers
#' @description Ensure stage numbers, return stage print message.
#' @param frm Starting stage
#' @param to Ending stage
#' @export
#' @seealso \code{\link{run}}
chckStgs <- function(frm, to) {
  if(to < 1 | frm < 1) {
    stop('Total stages to run cannot be less than 1.')
  }
  if(to > 4 | frm > 4) {
    stop('Total stages to run cannot be more than 4.')
  }
  if(frm > to) {
    stop('Starting stage must always come before ending stage.')
  }
  stgs <- c('taxise', 'download', 'cluster', 'cluster2')
  stgs_msg <- paste0(stgs[frm:to], collapse = ', ')
  paste0('Running stages: ', stgs_msg)
}

#' @name reset
#' @title Reset a PhyLoTa pipeline run
#' @description Resets the pipeline to a specified stage.
#' @param wd Working directory
#' @param stage Name of stage to which the pipeline will be reset
#' @param hard T/F, delete all cached data?
#' @export
#' @seealso \code{\link{restart}}, \code{\link{setParameters}}
reset <- function(wd, stage, hard=FALSE) {
  if(!stage %in% c('taxise', 'download', 'cluster', 'cluster2')) {
    stop('Invalid stage name.')
  }
  ps <- ldPrmtrs(wd)
  # TODO: hard/soft
  rstPrgrss(wd=wd, stg=stage)
  msg <- paste0('Reset pipeline to [', stage, ']')
  brdr <- paste0(rep('-', nchar(msg)), collapse='')
  msg <- paste0(brdr, '\n', msg, '\n', brdr)
  info(ps=ps, lvl=1, msg)
}

#' @name setParameters
#' @title Change parameters in a working directory
#' @description Reset parameters afater running \code{setUp()}.
#' @param wd Working directory
#' @param parameters Parameters to be changed
#' @param values New values for each parameter
#' @export
setParameters <- function(wd, parameters, values) {
  # TODO: make parameters an object with pre-defined parameter types
  ps <- ldPrmtrs(wd)
  for(i in 1:length(parameters)) {
    ps[[parameters[i]]] <- values[[i]]
  }
  svObj(wd=wd, obj=ps, nm='prmtrs')
  msg <- paste0('The following parameters have been reset:')
  brdr <- paste0(rep('-', nchar(msg)), collapse='')
  info(lvl=1, ps=ps, paste0(brdr, '\n', msg))
  mxnchrs <- max(sapply(parameters, nchar)) + 3
  for(prmtr in parameters) {
    spcr <- paste0(rep(' ', mxnchrs - nchar(prmtr)), collapse='')
    prmtr_msg <- paste0(prmtr, spcr, '[', ps[[prmtr]], ']')
    info(lvl=2, ps=ps, prmtr_msg)
  }
  info(lvl=1, ps=ps, brdr)
}