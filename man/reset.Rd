% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline.R
\name{reset}
\alias{reset}
\title{Reset a phylotaR pipeline run}
\usage{
reset(wd, stage, hard = FALSE)
}
\arguments{
\item{wd}{Working directory}

\item{stage}{Name of stage to which the pipeline will be reset}

\item{hard}{T/F, delete all cached data?}
}
\description{
Resets the pipeline to a specified stage.
}
\examples{
\dontrun{
  
  # Note: this example requires BLAST and internet to run.
  
  # example with temp folder
  wd <- file.path(tempdir(), 'aotus')
  # setup for aotus, make sure aotus/ folder already exists
  if (!dir.exists(wd)) {
    dir.create(wd)
  }
  ncbi_dr <- '[SET BLAST+ BIN PATH HERE]'
  setup(wd = wd, txid = 9504, ncbi_dr = ncbi_dr)  # txid for Aotus primate genus
  # individually run taxise
  taxise_run(wd = wd)
  # reset back to taxise as if it has not been run
  reset(wd = 'aotus', stage = 'taxise')
  # run taxise again ....
  taxise_run(wd = wd)
}
}
\seealso{
Other run-public: 
\code{\link{ClstrArc-class}},
\code{\link{ClstrRec-class}},
\code{\link{Phylota-class}},
\code{\link{SeqArc-class}},
\code{\link{SeqRec-class}},
\code{\link{TaxDict-class}},
\code{\link{TaxRec-class}},
\code{\link{clusters2_run}()},
\code{\link{clusters_run}()},
\code{\link{parameters_reset}()},
\code{\link{restart}()},
\code{\link{run}()},
\code{\link{setup}()},
\code{\link{taxise_run}()}
}
\concept{run-public}
