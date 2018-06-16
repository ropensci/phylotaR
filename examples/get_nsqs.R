data("cycads")
# count seqs for a random 10 clusters
random_cids <- sample(cycads@cids, 10)
nsqs <- get_nsqs(phylota = cycads, cid = random_cids)
