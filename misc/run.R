library(phylotaR)

set.seed(111)

options(error=recover)

# Global variables are 'evil'
# https://stackoverflow.com/questions/12598242/global-variables-in-packages-in-r
# Let's use an adjustable parameters list instead
ncbi_execs <- setUpNcbiTools(dr='/home/dom/Programs/ncbi-blast-2.7.1+/bin')
prmtrs <- mkPrmtrs(ncbi_execs=ncbi_execs)
## Download file ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz,
## unzip and specify directory where file 'nodes.dmp' is located
taxdir <- '/Users/hettling/taxdump'

## Do analysis for Bromeliaceae family
taxid <- 4613
## nodes.create(taxid, taxdir=taxdir, file.name='dbfiles-bromeliaceae-nodes.tsv')

clusters.ci_gi.seqs.create(15123, 'dbfiles-bromeliaceae-nodes.tsv',
                           files=list(clusters='dbfiles-bromeliaceae-clusters.tsv',
                                      ci_gi='dbfiles-bromeliaceae-ci_gi.tsv',
                                      seqs='dbfiles-bromeliaceae-seqs.tsv'))
