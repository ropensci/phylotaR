



devtools::load_all('~/Coding/phylotaR_restez')
library(phylotaR)

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

