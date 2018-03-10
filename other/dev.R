
devtools::load_all('~/Coding/phylotaR')




# INPUT
ptable <- read_phylota(wd)

# CLTYPE STATS
# drop clusters of 10
# n_taxa <- get_ntaxa(ptable, cid=ptable@cids)
# n_taxa <- get_cl_slot(ptable, cid=ptable@cids, slt_nm='ntx')
# keep <- names(n_taxa)[n_taxa > 10]
keep <- ptable@cids[1:100]
ptable <- drop_cls(ptable, keep)
table(sapply(ptable@cls@cls, function(x) x@typ))

# REDUCE
# get n taxa per cluster
n_taxa <- get_ntaxa(ptable, cid=ptable@cids)
# drop all clusters with fewer than 100 taxa
keep <- names(n_taxa)[n_taxa > 100]
ptable <- drop_cls(ptable, keep)

# FILTER
fltrd <- drop_by_rank(ptable, rnk='genus', n=2,
                      choose_by=c('nncltds'),
                      greatest=c(TRUE))

# SUMMARISE
smmry <- summary_phylota(fltrd)
smmry[grepl('rna', smmry$Feature),]
smmry <- smmry[order(smmry$N_taxa, decreasing=TRUE), ]
slctd_smmry <- smmry[smmry[['MAD']] > 0.6, ]
slctd_smmry <- slctd_smmry[1:10, ]
slctd_smmry$ID <- as.numeric(slctd_smmry$ID)

# KEEP ONLY TOP TEN
slctd <- drop_cls(fltrd, as.character(slctd_smmry$ID))

# GET GENUS NAMES
txids <- get_txids(slctd, cid='0', rnk='genus')
genus_names <- get_tx_slot(slctd, txids, 'scnm')
table(genus_names)
