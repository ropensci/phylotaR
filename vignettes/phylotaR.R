## ----init, eval=TRUE, include=FALSE, warning=FALSE, error=FALSE---------------
ncbi_dr = "/workspaces/blast/bin/"
eval_pipe = dir.exists(ncbi_dr)
if (eval_pipe) {
  library(phylotaR)
  wd = tempdir()
  unlink(wd)
}

## ----setup, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  library(phylotaR)
#  wd <- '[YOUR PATH TO AOTUS FOLDER]'
#  ncbi_dr <- '[YOUR PATH TO NCBI BLAST TOOLS]'

## ----setup2, eval=eval_pipe, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
txid <- 9504
setup(wd = wd, txid = txid, ncbi_dr = ncbi_dr, v = TRUE)

## ----running, eval=eval_pipe, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
run(wd = wd)

## ----restarting, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  reset(wd = wd, stage = 'cluster')
#  restart(wd = wd)

## ----parameters reset, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  # use ctrl+c or Esc to halt R
#  # increase the btchsz from the default to 300
#  parameters_reset(wd = wd, parameters = 'btchz', values = 300)
#  restart(wd = wd)
#  # ^ restart from whatever point it was halted

## ----selection1, eval=TRUE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
library(phylotaR)
# pre-load already run aotus from package data
data('aotus')
all_clusters <- aotus
print(all_clusters)
# otherwise, run:
# all_clusters <- read_phylota(wd)
cids <- all_clusters@cids
n_taxa <- get_ntaxa(phylota = all_clusters, cid = cids)

## ----selection2, eval=TRUE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
keep <- cids[n_taxa > 6]
selected <- drop_clstrs(phylota = all_clusters, cid = keep)
smmry <- summary(selected)
print(smmry)

## ----selection3, eval=TRUE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
cid <- smmry[2, 'ID']
# get the cluster record
cluster_record <- selected@clstrs[[cid]]
# use the seq. IDs to get the sequence records
seq_records <- selected@sqs[cluster_record@sids]
# extract a single record
seq_record <- seq_records[[seq_records@ids[[1]]]]
summary(seq_record)
# get the sequence
seq <- rawToChar(seq_record@sq)
print(substr(x = seq, start = 1, stop = 80))

## ----selection4, eval=TRUE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
# choose best sequence per species
reduced <- drop_by_rank(phylota = selected, rnk = 'species', n = 1)
# get txids at the species level for each sequence
txids <- get_txids(phylota = reduced, cid = cid, rnk = 'species')
# look up name for txids
scientific_names <- get_tx_slot(phylota = reduced, txid = txids, slt_nm = 'scnm')
# clean the names
scientific_names <- gsub('\\.', '', scientific_names)
scientific_names <- gsub('\\s+', '_', scientific_names)
print(scientific_names)
# look up sequence IDs for our chosen cluster
sids <- reduced@clstrs[[cid]]@sids
# write out
write_sqs(phylota = reduced, sid = sids, sq_nm = scientific_names,
          outfile = file.path(tempdir(), 'cytb.fasta'))
# ^ to avoid clutter, we're writing to a temporary folder

