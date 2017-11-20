#' @name mkPrmtrs
#' @title Return run parameters
#' @description Generates a paramters list that is passed on to
#' phylotaR functions.
#' @param ncbi_execs File directories for NCBI tools, see \code{setUpNcbiTools()}
#' @param mdl_thrshld Maximum number of sequences per species
#' @param mx_blst_sqs Maximum number of sequences to blast in a single run;
#' @param mx_sq_lngth Maximum characters in one sequence
#' @param sq_cch_dr Directory for sequence cache; will be created if does not exist
#' @param cores Number of cores for parellisation 
#' @details
#' ncbi_execs: users must run \code{setUpNcbiTools()} before running pipeline
#' mx_blst_sqs: if taxon has more subtree sequences than that,
#' its children will get clustered.
#' cores: set to 1 to run sequentially
#' @export
#' @seealso
#' \link{setUpNcbiTools}
mkPrmtrs <- function(ncbi_execs,
                     mdl_thrshld=3000,
                     mx_blst_sqs=10000,
                     mx_sq_lngth=25000,
                     sq_cch_dr="sequences",
                     cores=1) {
  prmtrs <- list(mdl_thrshld=mdl_thrshld,
                 mx_blst_sqs=mx_blst_sqs,
                 mx_sq_lngth=mx_sq_lngth,
                 sq_cch_dr=sq_cch_dr,
                 cores=cores)
  prmtrs <- c(prmtrs, ncbi_execs)
  if(cores > 1) {
    doMC::registerDoMC(cores)
  }
  prmtrs
}

#' @name setUpNcbiTools
#' @title Ensures NCBI BLAST tools are installed
#' @description Ensures NCBI BLAST executables are installed on the 
#' system. Tests version number of BLAST tools.
#' @param dr Directory to NCBI BLAST tools
#' @details BLAST tools must be version >2.7.1.
#' @export
#' @seealso
#' \link{mkPrmtrs}
setUpNcbiTools <- function(dr) {
  sccdd <- TRUE
  mkblstdb <- file.path(dr, 'makeblastdb')
  blstn <- file.path(dr, 'blastn')
  for(ech in c(mkblstdb, blstn)) {
    cmd <- paste0(ech, ' -version')
    res <- try(system(cmd, intern = TRUE, ignore.stderr=TRUE),
               silent=TRUE)
    if(grepl('error', res[[1]], ignore.case=TRUE)) {
      tst <- FALSE
      cat('Invalid path: [', ech, '] \n', sep='')
      sccdd <- FALSE
    } else {
      # test version
      vrsn <- gsub('[a-zA-Z:+]', '', res[1])
      vrsn <- gsub('\\s', '', vrsn)
      vrsn <- as.numeric(strsplit(vrsn, '\\.')[[1]])
      tst <- vrsn[1] >= 2 & vrsn[2] >= 7 & vrsn[3] >= 1
      if(tst) {
        cat('Found: [', res[1], '] \n', sep='')
      } else {
        cat('Incorrect version: [', res[1], '] \n', sep='')
        sccdd <- FALSE
      }
    }
  }
  if(!sccdd) {
    stop('Unable to find correct versions of NCBI BLAST tools\n')
  }
  ncbi_execs <- list('mkblstdb'=mkblstdb,
                     'blstn'=blstn)
}