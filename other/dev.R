
cmd_path <- file.path('/usr/local/ncbi/blast/bin', 'makeblastdb')
res <- phylotaR:::cmdln(cmd = cmd_path, args = '-version')
# Is the status TRUE?
(res[['status']] == 0)
# What are the details?
stdout <- rawToChar(res[['stdout']])
(stdout <- strsplit(x = stdout, split = '\n')[[1]])
# Can we extract the version number?
vrsn <- gsub('[a-zA-Z:+]', '', stdout[[1]])
vrsn <- gsub('\\s', '', vrsn)
(vrsn <- as.numeric(strsplit(vrsn, '\\.')[[1]]))
# Version should be 2 +, expect TRUE
(vrsn[1] >= 2 & vrsn[2] >= 0)

cmd_path <- '/usr/local/ncbi/blast/bin'
phylotaR:::blast_setup(d = cmd_path, v = TRUE, wd = getwd())

stdout <- "MAKEBL~1: 2.7.1+\n Package: blast 2.7.1, build Oct 18 2017 19:55:35\n"


devtools::load_all('~/Coding/phylotaR_restez')

wd <- '~/Coding/restez/hystricomorpha'
ps <- parameters_load(wd)
ps[['wd']] <- '~/Coding/restez/hystricomorpha'


library(taxize)
# 7158 -- Aedes genus
aedes_spp <- downstream(x = '7158', db = 'ncbi', downto = 'species')
# 7174 -- Culex genus
culex_spp <- downstream(x = '7174', db = 'ncbi', downto = 'species')

txids <- culex_spp$`7174`$childtaxa_id
txids <- c(aedes_spp$`7158`$childtaxa_id, culex_spp$`7174`$childtaxa_id)
txids <- sample(txids, 60)
# save(txids, file = 'sample_txids.RData')
# load('sample_txids.RData')

# test lots of txids
ncbi_dr <- '/usr/local/ncbi/blast/bin/'
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
clusters_run(wd)
clusters2_run(wd)

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
