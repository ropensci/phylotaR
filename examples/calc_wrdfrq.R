data('dragonflies')
# work out what gene region the cluster is likely representing with word freqs.
random_cids <- sample(dragonflies@cids, 10)
# most frequent words in definition line
(calc_wrdfrq(phylota = dragonflies, cid = random_cids, type = 'dfln'))
# most frequent words in feature name
(calc_wrdfrq(phylota = dragonflies, cid = random_cids, type = 'nm'))
