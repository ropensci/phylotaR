library(ips)
load(file='dev_data.RData')
wd <- '/home/dom/Desktop/testing_phylotaR/aotus'
phylt_nds <- ldObj(wd=wd, nm='phylt_nds')

names(clstrs[[1]])

clstrs[[1]]
pull <- sapply(clstrs, function(x) length(unique(x[['tis']]))) > 10
pull <- (sapply(clstrs, function(x) x[['MinLength']]) > 500) & pull
sapply(clstrs[pull], function(x) length(x[['tis']]))
clstrs[[]]

which(pull)[1]

clstr_sqs <- sqs[clstrs[[20]][['gis']]]

sq <- clstr_sqs[[1]]
fastasqs <- rep(NA, length(clstr_sqs))
for(i in 1:length(clstr_sqs)) {
  sq <- clstr_sqs[[i]]
  fastasq <- paste0('>', sq[['gi']], '|',
                    sq[['ti']], '\n',
                    sq[['seq']], '\n\n')
  fastasqs[[i]] <- fastasq
}
x <- read.fas(text=fastasqs)
align <- mafft(x=x, path='/usr/bin/mafft')
tree <- raxml(DNAbin=align, m="GTRGAMMA",
              f="a", N=10, p=1234, x=1234, 
              exec='/usr/bin/raxmlHPC')
tree[[1]]
is(tree)
plot(tree[[2]])

phylo <- tree[[2]]
library(treeman)
tree <- as(phylo, 'TreeMan')

calcDstMtrx(tree, tree['tips'])
