# LIBS
library(phylotaR)

# VARS ----
demos <- list('anisoptera' = '6962', 'acipenseridae' = '7900',
              'tinamiformes' = '8802', 'aotus' = '9504',
              'bromeliaceae' = '4613', 'cycadidae' = '1445963',
              'eutardigrada' = '42242', 'kazachstania' = '71245',
              'platyrrhini' = '9479')
# set these paths for your own system
ncbi_dr <- readLines(file.path('demos', 'ncbi_dr.txt'))

# FOR LOOP ----
for (i in seq_along(demos)) {
  txid <- demos[[i]]
  wd <- file.path(getwd(), 'demos', names(demos)[[i]])
  # create folder, delete if already exists
  if (file.exists(wd)) {
    unlink(wd, recursive = TRUE)
  }
  dir.create(wd)
  # run
  setup(wd = wd, txid = txid, ncbi_dr = ncbi_dr, v = TRUE)
  run(wd = wd)
}

# TIMINGS ----
timings <- vector(mode = 'list', length = length(demos))
names(timings) <- names(demos)
for (i in seq_along(demos)) {
  wd <- file.path(getwd(), 'demos', names(demos)[[i]])
  timings[[i]] <- get_stage_times(wd)
}

# MARKDOWN ----
mrkdwn <- 'Taxon|Taxise|Download|Cluster|Cluster2|Total|\n'
mrkdwn <- paste0(mrkdwn, '|:--|--:|--:|--:|--:|--:|\n')
for (i in seq_along(timings)) {
  mrkdwn <- paste0(mrkdwn, Hmisc::capitalize(names(timings)[[i]]), '|',
                   signif(timings[[i]][['taxise']], 2), '|',
                   signif(timings[[i]][['download']], 2), '|',
                   signif(timings[[i]][['cluster']], 2), '|',
                   signif(timings[[i]][['cluster\\^2']], 2), '|',
                   signif(sum(timings[[i]]), 2), '|\n')
}
cat(mrkdwn)