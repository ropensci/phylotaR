% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stage4.R
\name{clusters2_run}
\alias{clusters2_run}
\title{Run the cluster2 stage}
\usage{
clusters2_run(wd)
}
\arguments{
\item{wd}{Working directory}
}
\description{
Run the fourth stage of the phylotaR pipeline,
cluster2. Identify clusters at higher taxonomic levels by
merging sister clusters.
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
  # individually run stages
  taxise_run(wd = wd)
  download_run(wd = wd)
  clusters_run(wd = wd)
  clusters2_run(wd = wd)
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
\code{\link{clusters_run}()},
\code{\link{parameters_reset}()},
\code{\link{reset}()},
\code{\link{restart}()},
\code{\link{run}()},
\code{\link{setup}()},
\code{\link{taxise_run}()}
}
\concept{run-public}
