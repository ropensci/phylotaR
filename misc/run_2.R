# Set working directory
wd <- '/home/dom/Desktop/testing_phylota'
# Set up execs and parameters
ncbi_execs <- setUpNcbiTools(d='/home/dom/Programs/ncbi-blast-2.7.1+/bin')
setUpPrmtrs(wd=wd,
            ncbi_execs=ncbi_execs,
            mdl_thrshld=3000,
            mx_blst_sqs=10000,
            mx_sq_lngth=25000,
            cores=2)


