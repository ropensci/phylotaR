



devtools::load_all('~/Coding/phylotaR_restez')

wd <- '~/Coding/restez/hystricomorpha'
ps <- parameters_load(wd)
ps[['wd']] <- '~/Coding/restez/hystricomorpha'


ncbi_dr <- file.path('NCBI', 'bin')
wd <- 'aotus'
if (file.exists(wd)) {
  unlink(wd, recursive = TRUE)
}
dir.create(wd)
setup(wd = wd, txid = '9504', ncbi_dr = ncbi_dr, v = TRUE)
run(wd)
taxise_run(wd)
download_run(wd)


# Polyporaltes error
# problem: KY475568 does not have the correct organism name, cannot be found in txdct
# 'Polyporaceae sp. W-2017a'
