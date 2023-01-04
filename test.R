# Test example from issue https://github.com/ropensci/phylotaR/issues/48
ncbi_dr = "/workspaces/blast/bin/"
devtools::load_all()
wd = "/workspaces/phylotaR-test"
#unlink(wd)

# library(phylotaR)
# wd <- '/media/...
# ncbi_dr <- '/home/...
# txid <- 4747
# txid <- 9504 # Small
txid <- 1912919
setup(wd = wd, txid = txid, ncbi_dr = ncbi_dr, v=TRUE, overwrite = TRUE)
run(wd = wd)

all_clusters <- read_phylota(wd)
all_clusters
cids <- all_clusters@cids
length(cids)
n_taxa <- get_ntaxa(phylota = all_clusters, cid = cids)
summary(n_taxa)
keep <- cids[n_taxa >= 10]
length(keep)
selected <- drop_clstrs(phylota = all_clusters, cid = keep)
selected

debug(drop_by_rank)
reduced <- drop_by_rank(selected, rnk = 'species', n=1)

slct(unqids[1])
slct(unqids[20])
slct(unqids[21])
keep <- unlist(lapply(unqids, slct))

lapply(seq_along(unqids), function(i){
  print(i)
  slct(unqids[i])
})

slct(unqids[21])
txid = unqids[21]
slct(txid)
# choose_by
# [1] "pambgs"  "age"     "nncltds"

get_sq_slot(
  phylota = phylota, sid = pssbls,
  slt_nm = "nncltds"
)

debug(slct)
slct(txid)

is.na(vals)

slct = function(txid) {
    pssbls <- sids[txid == txids]
    for (i in seq_along(choose_by)) {
      vals <- get_sq_slot(
        phylota = phylota, sid = pssbls,
        slt_nm = choose_by[[i]]
      )
      if (anyNA(vals)) {
        next()
      }
      names(vals) <- pssbls
      mx_n <- ifelse(length(vals) > n, n, length(vals))
      vals <- sort(x = vals, decreasing = greatest[i])[1:mx_n]
      pssbls <- names(vals)
    }
    
    pssbls
  }
