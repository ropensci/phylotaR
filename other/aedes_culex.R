# Libs ----
# devtools::install_github('ropensci/phylotaR', ref = 'multiple_ids')
library(phylotaR)
library(taxize)

# Txids ----
# look-up random selection of 60 aedes and culex spp
# 7158 -- Aedes genus
aedes_spp <- downstream(x = '7158', db = 'ncbi', downto = 'species')
# 7174 -- Culex genus
culex_spp <- downstream(x = '7174', db = 'ncbi', downto = 'species')

txids <- c(aedes_spp$`7158`$childtaxa_id, culex_spp$`7174`$childtaxa_id)
txids <- sample(txids, 60)
# save(txids, file = 'sample_txids.RData')
# load('sample_txids.RData')

# phylotaR ----
ncbi_dr <- '/usr/bin/'
wd <- 'aedes_culex'
if (file.exists(wd)) {
  unlink(wd, recursive = TRUE)
}
dir.create(wd)
setup(wd = wd, txid = txids, ncbi_dr = ncbi_dr, v = TRUE)
run(wd)
