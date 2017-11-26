# For testing taxonomy tools
# nodes.dmp and names.dmp with complete
#  lineages are required

library(phylotaR)
# put 'NCBI' taxonomy folder in wd
tdobj <- genTDObj(wd='.')
td_nds <- tdobj[['nds']]
td_nms <- tdobj[['nms']]
rm(tdobj)
platyrrhini <- 9479
queue <- ids <- vector()
id <- platyrrhini
while(TRUE) {
  ids <- c(id, ids)
  queue <- c(getKids(id, td_nds=td_nds), queue)
  if(length(queue) == 0) {
    break
  }
  id <- queue[1]
  queue <- queue[-1]
}
# create new td_nds and td_nms
new_nds <- td_nds[td_nds[['id']] %in% ids, ]
new_nms <- td_nms[td_nms[['id']] %in% ids, ]
# reduce factors
new_nds[['rank']] <- factor(new_nds[['rank']])
new_nms[['name']] <- factor(new_nms[['name']])
new_nms[['type']] <- factor(new_nms[['type']])
tdobj <- list('nds'=new_nds, 'nms'=new_nms)
save(tdobj, file="data/tdobj.rda", compress="xz")
