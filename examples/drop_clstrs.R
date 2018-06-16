data("dragonflies")
# specify cids to *keep*
random_cids <- sample(dragonflies@cids, 100)
# drop an entire cluster
nbefore <- length(dragonflies@cids)
dragonflies <- drop_clstrs(phylota = dragonflies, cid = random_cids)
nafter <- length(dragonflies@cids)
# now there are only 100 clusters
(nafter < nbefore)
