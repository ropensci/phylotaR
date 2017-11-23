library(phylotaR)

set.seed(111)

#options(error=recover)

# Please don't use global variables
# They're very hard to track and I can't incorporate them in an R package.
ncbi_execs <- setUpNcbiTools(dr='/home/dom/Programs/ncbi-blast-2.7.1+/bin')
prmtrs <- mkPrmtrs(ncbi_execs=ncbi_execs, cores=2)
## Download file ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz,
## unzip and specify directory where file 'nodes.dmp' is located
# always best to use file.path -- OS independent
taxdir <- '/home/dom/Desktop/testing_phylotaR'

## Do analysis for Bromeliaceae family
taxid <- 4613
nodes.create(root.taxa=taxid, taxdir=taxdir,
             file.name='dbfiles-bromeliaceae-nodes.tsv')

clusters.ci_gi.seqs.create(15123, 'dbfiles-bromeliaceae-nodes.tsv',
                           files=list(clusters='dbfiles-bromeliaceae-clusters.tsv',
                                      ci_gi='dbfiles-bromeliaceae-ci_gi.tsv',
                                      seqs='dbfiles-bromeliaceae-seqs.tsv'),
                           informative=FALSE, prmtrs=prmtrs)
