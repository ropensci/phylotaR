data('aotus')
random_sid <- sample(aotus@sids, 1)
(get_sq_slot(phylota = aotus, sid = random_sid, slt_nm = 'dfln'))
# see list_seqrec_slots() for available slots
(list_seqrec_slots())
