
require('parallel')
require('data.table')
cl <- makeCluster(getOption("cl.cores", 4))
##f = file('test_file')
#inc = 1000
#max = 10000
#skip <- seq(0, max-inc, inc)
#skip[1] <- 1


max <- 110920890
skip <- seq(0, max, length.out=100)
skip[1] <- 1

l <- parLapply(cl, skip, function(sk) {
    require('data.table')
    inc <- 1120413 
    n <- scan(file='nucl_gb.accession2taxid', what=list('character', 'character', 'numeric', 'numeric'), nlines=inc, skip=sk, sep='\t')
    ##n <- fread('nucl_gb.accession2taxid', colClasses=c('character', 'character', 'numeric', 'numeric'), nrows=inc, skip=sk, sep='\t',showProgress=T, data.table=F )
    acc <- n[[1]]
    acc.v <- n[[2]]
    taxid <- n[[3]]
    gi <- n[[4]]
    data.frame(acc, acc.v, taxid, gi)        
})

require('doParallel')
require('doSNOW')
cl <- makeCluster(4, outfile="")
registerDoSNOW(cl)
##registerDoParallel(cl)
max <- 110920890
skip <- seq.int(0, max, 100000)
skip[1] <- 1


r <- foreach(i=seq_along(skip), .combine=rbind) %dopar% {
##    require('data.table')
##    require('ff')
    require('sqldf')
    inc <- 100000
    sk <- skip[i]
    cat("Processing chunk ", i, "\n")
    ##n <- scan(file='nucl_gb.accession2taxid', what=list('character', 'character', 'numeric', 'numeric'), nlines=inc, skip=sk, sep='\t')
    ##n <- fread('nucl_gb.accession2taxid', colClasses=c('character', 'character', 'numeric', 'numeric'),
    ##           nrows=inc, skip=sk, sep='\t',showProgress=T, data.table=F )
    #n <- read.table.ffdf(file='nucl_gb.accession2taxid', nrows=inc, skip=sk)
    f <- file('nucl_gb.accession2taxid')    
    n <- sqldf("select * from f limit 100000", dbname=tempfile(), file.format=list(header=T, row.names=F, skip=sk))
                                        #                                    #    acc <- n[[1]]
#    acc.v <- n[[2]]
#    taxid <- n[[3]]
#    gi <- n[[4]]
#    data.frame(acc, acc.v, taxid, gi)        
#    n <- read.csv.ffdf(nrows=inc, )
    cat("Finished processing chunk ", i, "\n")
    n
}
