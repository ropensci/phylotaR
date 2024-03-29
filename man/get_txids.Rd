% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-get.R
\name{get_txids}
\alias{get_txids}
\title{Get taxonomic IDs by rank}
\usage{
get_txids(
  phylota,
  cid = NULL,
  sid = NULL,
  txids = NULL,
  rnk = NULL,
  keep_higher = FALSE
)
}
\arguments{
\item{phylota}{Phylota object}

\item{cid}{Cluster ID}

\item{sid}{Sequence ID(s)}

\item{txids}{Vector of txids}

\item{rnk}{Taxonomic rank}

\item{keep_higher}{Keep higher taxonomic IDs?}
}
\value{
vector
}
\description{
Return taxonomic IDs for a vector of sequence IDs or all
sequences in a cluster. User can specify what rank the IDs should be
returned. If NULL, the lowest level is returned.
}
\details{
txids can either be provided by user or they can be determined for
a vector of sids or for a cid. If keep_higher is TRUE, any sequence that has
a identity that is higher than the given rank will be returned. If FALSE,
these sequences will return ''.
}
\examples{
data('bromeliads')
# get all the genus IDs and names
genus_ids <- get_txids(phylota = bromeliads, txids = bromeliads@txids,
                       rnk = 'genus')
genus_ids <- unique(genus_ids)
# drop empty IDs -- this happens if a given lineage has no ID for specified rank
genus_ids <- genus_ids[genus_ids != '']
# get names
(get_tx_slot(phylota = bromeliads, txid = genus_ids, slt_nm = 'scnm'))
}
\seealso{
Other tools-public: 
\code{\link{calc_mad}()},
\code{\link{calc_wrdfrq}()},
\code{\link{drop_by_rank}()},
\code{\link{drop_clstrs}()},
\code{\link{drop_sqs}()},
\code{\link{get_clstr_slot}()},
\code{\link{get_nsqs}()},
\code{\link{get_ntaxa}()},
\code{\link{get_sq_slot}()},
\code{\link{get_stage_times}()},
\code{\link{get_tx_slot}()},
\code{\link{is_txid_in_clstr}()},
\code{\link{is_txid_in_sq}()},
\code{\link{list_clstrrec_slots}()},
\code{\link{list_ncbi_ranks}()},
\code{\link{list_seqrec_slots}()},
\code{\link{list_taxrec_slots}()},
\code{\link{plot_phylota_pa}()},
\code{\link{plot_phylota_treemap}()},
\code{\link{read_phylota}()},
\code{\link{write_sqs}()}
}
\concept{tools-public}
