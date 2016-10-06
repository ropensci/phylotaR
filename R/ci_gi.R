## script to populate the phylota 'clusters' table
## Dependencies: dataframe of clusters, as produced by .make.clusters from clusters.R

## CREATE TABLE "ci_gi_194" (
##   "ti" int(10)  DEFAULT NULL,
##   "clustid" int(10)  DEFAULT NULL,
##   "cl_type" text  DEFAULT NULL,
##   "gi" bigint(20)  DEFAULT NULL,
##   "ti_of_gi" int(10)  DEFAULT NULL
## );
## CREATE INDEX "ci_gi_194_cl_type" ON "ci_gi_194" ("cl_type");
## CREATE INDEX "ci_gi_194_gi" ON "ci_gi_194" ("gi");
## CREATE INDEX "ci_gi_194_ti" ON "ci_gi_194" ("ti","clustid","cl_type");
## CREATE INDEX ti_of_gi ON ci_gi_194(ti_of_gi);

## Input: list of clusters, as produced by .make.clusters
ci_gi.create <- function(clusters) {
    ## select relevant columns and make data frame
    df <- do.call(rbind, lapply(cl, function(x)as.data.frame(x)[ , c('ti_root', 'ci', 'cl_type', 'gis', 'tis') ]))
    names(df)[names(df) == 'gis'] <- 'gi'
    names(df)[names(df) == 'tis'] <- 'ti_of_gi'
    names(df)[names(df) == 'ci'] <- 'clustid'
    names(df)[names(df) == 'ti_root'] <- 'ti'

    return(df)
}
