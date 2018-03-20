
devtools::load_all('~/Coding/phylotaR')

aotus_ps <- ldPrmtrs('demos/aotus')
wd <- '/home/dom/Coding/phylotaR_demo/palms'
ps <- ldPrmtrs(wd)
ps[['wd']] <- wd
ps[['v']] <- TRUE
ps$mkblstdb <- aotus_ps$mkblstdb
ps$blstn <- aotus_ps$blstn
clstrClstrs(ps)
