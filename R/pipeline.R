#' @name setup
#' @title Set-up parameters
#' @description Set up working directory with parameters.
#' @details See \code{\link{parameters}}() for a description of all parameters
#' and their defaults. You can change parameters after a folder has been set up
#' with \code{\link{parameters_reset}}().
#' @param wd Working directory
#' @param txid Root taxonomic ID(s), vector or numeric
#' @param ncbi_dr Directory to NCBI BLAST tools, default '.'
#' @param v Verbose, T/F
#' @param ... Additional parameters
#' @export
#' @return NULL
#' @family run-public
#' @example examples/setup.R
setup <- function(wd, txid, ncbi_dr='.', v=FALSE, ...) {
  # no ~
  checks <- vapply(X = c(wd, ncbi_dr), FUN = grepl, FUN.VALUE = logical(1),
                   pattern = '~')
  if (any(checks)) {
    stop(paste0('Do not use `~` in filepaths.'))
  }
  # header log
  msg <- paste0('phylotaR: Implementation of PhyLoTa in R [v',
                packageVersion('phylotaR'), ']')
  brdr <- paste0(rep('-', nchar(msg)), collapse = '')
  msg <- paste0(brdr, '\n', msg, '\n', brdr, '\n')
  .log(v = v, wd = wd, msg)
  # set up
  ncbi_execs <- blast_setup(d = ncbi_dr, v = v, wd = wd)
  parameters_setup(wd = wd, txid = txid, ncbi_execs = ncbi_execs, v = v, ...)
  # record session info
  writeLines(capture.output(sessionInfo()), file.path(wd, "session_info.txt"))
  # end
  .log(v = v, wd = wd, paste0(brdr, '\n'))
}

#' @name run
#' @title Run phylotaR pipeline
#' @description Run the entire phylotaR pipeline. All generated files will be
#' stored in the wd. The process can be stopped at anytime and  restarted with
#' \code{restart}. \code{nstages} must be a numeric value representing the
#' number of stages that will be run. Stages are run in the following order:
#' 1 - taxise, 2 - download, 3 - cluster and 4 - cluster2.
#' 
#' For example, specifying \code{nstages} = 3, will run taxise, download and
#' cluster. Stages can also be run individually, see linked functions below.
#' @param wd Working directory
#' @param nstages Number of total stages to run, max 4.
#' @family run-public
#' @export
#' @example examples/run.R
#' @return NULL
run <- function(wd, nstages=4) {
  stgs_msg <- stage_args_check(frm = 1, to = nstages)
  stages_run(wd = wd, frm = 1, to = nstages, stgs_msg = stgs_msg)
}

#' @name restart
#' @title Restart a phylotaR pipeline run
#' @description Restarts the running of a pipeline
#' as started with \code{run}.
#' @param wd Working directory
#' @param nstages Number of total stages to run, max 4.
#' @export
#' @family run-public
#' @example examples/restart.R
#' @return NULL
restart <- function(wd, nstages=4) {
  stg <- progress_read(wd)
  if (is.na(stg)) {
    stop('Pipeline already complete. Use `reset()` to re-run pipeline.')
  }
  frm <- which(c('taxise', 'download', 'cluster', 'cluster2')
                %in% stg)
  if (frm > nstages) {
    stop('Pipeline has already completed [', nstages,
         '] stages. Increase `nstages`.')
  }
  stgs_msg <- stage_args_check(frm = frm, to = nstages)
  stages_run(wd = wd, frm = frm, to = nstages, stgs_msg = stgs_msg,
             rstrt = TRUE)
}

#' @name reset
#' @title Reset a phylotaR pipeline run
#' @description Resets the pipeline to a specified stage.
#' @param wd Working directory
#' @param stage Name of stage to which the pipeline will be reset
#' @param hard T/F, delete all cached data?
#' @family run-public
#' @export
#' @example examples/reset.R
#' @return NULL
reset <- function(wd, stage, hard=FALSE) {
  if (!stage %in% c('taxise', 'download', 'cluster', 'cluster2')) {
    stop('Invalid stage name.')
  }
  ps <- parameters_load(wd)
  if (hard) {
    flpth <- file.path(ps[['wd']], 'cache')
    rdata_fls <- list.files(flpth, '.RData')
    rdata <- vector(mode = 'list', length = length(rdata_fls))
    names(rdata) <- rdata_fls
    for (rdata_fl in rdata_fls) {
      rdata[[rdata_fl]] <- readRDS(file = file.path(flpth, rdata_fl))
    }
    cache_setup(ps = ps, ovrwrt = TRUE)
    for (rdata_fl in rdata_fls) {
      saveRDS(object = rdata[[rdata_fl]], file = file.path(flpth, rdata_fl))
    }
  }
  progress_reset(wd = wd, stg = stage)
  if (hard) {
    msg <- paste0('Reset (hard) pipeline to [', stage, ']')
  } else {
    msg <- paste0('Reset (soft) pipeline to [', stage, ']')
  }
  brdr <- paste0(rep('-', nchar(msg)), collapse = '')
  msg <- paste0(brdr, '\n', msg, '\n', brdr)
  info(ps = ps, lvl = 1, msg)
}

#' @name parameters_reset
#' @title Change parameters in a working directory
#' @description Reset parameters after running \code{setup()}.
#' @param wd Working directory
#' @param parameters Parameters to be changed, vector.
#' @param values New values for each parameter, vector.
#' @family run-public
#' @export
#' @example examples/parameters_reset.R
#' @return NULL
parameters_reset <- function(wd, parameters, values) {
  # TODO: make parameters an object with pre-defined parameter types
  ps <- parameters_load(wd)
  for (i in seq_along(parameters)) {
    ps[[parameters[i]]] <- values[[i]]
  }
  obj_save(wd = wd, obj = ps, nm = 'prmtrs')
  msg <- paste0('The following parameters have been reset:')
  brdr <- paste0(rep('-', nchar(msg)), collapse = '')
  info(lvl = 1, ps = ps, paste0(brdr, '\n', msg))
  mxnchrs <- max(vapply(parameters, nchar, integer(1))) + 3
  for (prmtr in parameters) {
    spcr <- paste0(rep(' ', mxnchrs - nchar(prmtr)), collapse = '')
    prmtr_msg <- paste0(prmtr, spcr, '[', ps[[prmtr]], ']')
    info(lvl = 2, ps = ps, prmtr_msg)
  }
  info(lvl = 1, ps = ps, brdr)
}
