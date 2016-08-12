## script to populate the 'accession2taxid' mapping table,
## needed to create the other phylota tables
## The table is created from the mapping file
## accession2taxid/nucl_gb.accession2taxid

accession2taxid.create <- function(dbloc, taxdir, droptable=F) {

    ## specify location of NCBI mapping file
    subdir <- 'accession2taxid'
    file.loc <- paste0(taxdir, '/', subdir, '/nucl_gb.accession2taxid')

    ## print sqlite statements to temporary file
    tmp <- tempfile()
    cat(file=tmp, '.separator "\\t"\n')
    ## table might already exist
    if (droptable) {
        cat(file=tmp, append=T, 'DROP TABLE accession2taxid;\n')
    }
    ## create table to get correct types
    cat(file=tmp, append=T, 'CREATE TABLE accession2taxid("accession" TEXT, "accession.version" TEXT, "taxid" INTEGER, "gi" INTEGER);\n')
    cat(file=tmp, append=T, paste('.import', file.loc, 'accession2taxid\n'))

    ## create an index on taxid
    cat(file=tmp, append=T, 'CREATE INDEX taxid on accession2taxid(taxid);\n')

    ## prepare command and run
    cmd <- paste('sqlite3', dbloc, '<', tmp)
    ret <- system(cmd)

    ## remove temp file
    unlink(tmp)
    return(ret==0)
}
