data('bromeliads')
# how many species are there?
(get_ntaxa(phylota = bromeliads, cid = '0', rnk = 'species'))
# how many genera are there?
(get_ntaxa(phylota = bromeliads, cid = '0', rnk = 'genus'))
# how many families are there?
(get_ntaxa(phylota = bromeliads, cid = '0', rnk = 'family'))
# use list_ncbi_ranks() to see available rank names
(list_ncbi_ranks())
