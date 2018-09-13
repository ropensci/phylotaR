data("dragonflies")
# For faster computations, let's only work with the 5 clusters.
dragonflies <- drop_clstrs(phylota = dragonflies, cid = dragonflies@cids[10:15])


# We can use drop_by_rank() to reduce to 10 sequences per genus for each cluster
(reduced_1 <- drop_by_rank(phylota = dragonflies, rnk = 'genus', n = 10,
                           choose_by = c('pambgs', 'age', 'nncltds'),
                           greatest = c(FALSE, FALSE, TRUE)))

# We can specify what aspects of the sequences we would like to select per genus
# By default we select the sequences with fewest ambiguous nucleotides (e.g.
# we avoid Ns), the youngest age and then longest sequence.
# We can reverse the 'greatest' to get the opposite.
(reduced_2 <- drop_by_rank(phylota = dragonflies, rnk = 'genus', n = 10,
                           choose_by = c('pambgs', 'age', 'nncltds'),
                           greatest = c(TRUE, TRUE, FALSE)))


# Leading to smaller sequnces ...
r1_sqlngth <- mean(get_sq_slot(phylota = reduced_1,
                                sid = reduced_1@sids, slt_nm = 'nncltds'))
r2_sqlngth <- mean(get_sq_slot(phylota = reduced_2,
                                sid = reduced_2@sids, slt_nm = 'nncltds'))
(r1_sqlngth > r2_sqlngth)
# ... with more ambigous characters ....
r1_pambgs <- mean(get_sq_slot(phylota = reduced_1, sid = reduced_1@sids,
                              slt_nm = 'pambgs'))
r2_pambgs <- mean(get_sq_slot(phylota = reduced_2, sid = reduced_2@sids,
                              slt_nm = 'pambgs'))
(r1_pambgs < r2_pambgs)
# .... and older ages (measured in days since being added to GenBank).
r1_age <- mean(get_sq_slot(phylota = reduced_1, sid = reduced_1@sids,
                           slt_nm = 'age'))
r2_age <- mean(get_sq_slot(phylota = reduced_2, sid = reduced_2@sids,
                           slt_nm = 'age'))
(r1_age < r2_age)


# Or... we can simply reduce the clusters to just one sequence per genus
(dragonflies <- drop_by_rank(phylota = dragonflies, rnk = 'genus', n = 1))
