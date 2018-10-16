
devtools::load_all('~/Coding/phylotaR')
library(phylotaR)

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
