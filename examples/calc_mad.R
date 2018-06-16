data("bromeliads")
random_cids <- sample(bromeliads@cids, 10)
(calc_mad(phylota = bromeliads, cid = random_cids))
