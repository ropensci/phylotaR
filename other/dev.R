
devtools::load_all('~/Coding/phylotaR')


ncbi_dr <- file.path('NCBI', 'bin')
wd <- 'aotus'
if (file.exists(wd)) {
  unlink(wd, recursive = TRUE)
}
dir.create(wd)
setup(wd = wd, txid = '9479', ncbi_dr = ncbi_dr, v = TRUE)
taxise_run(wd)
