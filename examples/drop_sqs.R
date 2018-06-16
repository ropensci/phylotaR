data("dragonflies")
# drop random sequences from cluster 0
clstr <- dragonflies[['0']]
# specify the sids to *keep*
sids <- sample(clstr@sids, 100)
(dragonflies <- drop_sqs(phylota = dragonflies, cid = '0', sid = sids))
# Note, sequences dropped may be represented in other clusters
