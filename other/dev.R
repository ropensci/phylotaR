
devtools::load_all('~/Coding/phylotaR')
library(phylotaR)

library(taxize)
# 7158 -- Aedes genus
aedes_spp <- downstream(x = '7158', db = 'ncbi', downto = 'species')
# 7174 -- Culex genus
culex_spp <- downstream(x = '7174', db = 'ncbi', downto = 'species')

txids <- culex_spp$`7174`$childtaxa_id
txids <- c(aedes_spp$`7158`$childtaxa_id, culex_spp$`7174`$childtaxa_id)
txids <- sample(txids, 60)
# save(txids, 'sample_txids.RData')
# load('sample_txids.RData')

# test lots of txids
ncbi_dr <- '/usr/bin/'
wd <- 'aedes_culex'
if (file.exists(wd)) {
  unlink(wd, recursive = TRUE)
}
dir.create(wd)
library(phylotaR)
setup(wd = wd, txid = txids, ncbi_dr = ncbi_dr, v = TRUE)
run(wd)
taxise_run(wd)
download_run(wd)


# test multiple txids
ncbi_dr <- '/usr/bin/'
wd <- 'aotus_aloutta'
if (file.exists(wd)) {
  unlink(wd, recursive = TRUE)
}
dir.create(wd)
setup(wd = wd, txid = c('9504', '9499'), ncbi_dr = ncbi_dr, v = TRUE)
run(wd)
taxise_run(wd)
download_run(wd)


# Polyporaltes error
# problem: KY475568 does not have the correct organism name, cannot be found in txdct
# 'Polyporaceae sp. W-2017a'
