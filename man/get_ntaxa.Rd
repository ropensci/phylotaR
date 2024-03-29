% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user-get.R
\name{get_ntaxa}
\alias{get_ntaxa}
\title{Count number of unique taxa}
\usage{
get_ntaxa(phylota, cid = NULL, sid = NULL, rnk = NULL, keep_higher = FALSE)
}
\arguments{
\item{phylota}{Phylota object}

\item{cid}{Cluster ID(s)}

\item{sid}{Sequence ID(s)}

\item{rnk}{Taxonomic rank}

\item{keep_higher}{Keep higher taxonomic ranks?}
}
\value{
vector
}
\description{
Count the number of unique taxa represented by cluster(s) or
sequences in phylota table Use rnk to specify a taxonomic level to count. If
NULL counts will be made to the lowest level reported on NCBI.
}
\examples{
data('bromeliads')
# how many species are there?
(get_ntaxa(phylota = bromeliads, cid = '0', rnk = 'species'))
# how many genera are there?
(get_ntaxa(phylota = bromeliads, cid = '0', rnk = 'genus'))
# how many families are there?
(get_ntaxa(phylota = bromeliads, cid = '0', rnk = 'family'))
# use list_ncbi_ranks() to see available rank names
(list_ncbi_ranks())
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
\code{\link{get_sq_slot}()},
\code{\link{get_stage_times}()},
\code{\link{get_tx_slot}()},
\code{\link{get_txids}()},
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
