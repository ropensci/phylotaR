#' @name search_and_cache
#' @title Run rentrez function and cache results
#' @description Safely run a rentrez function. If the query fails, the
#' function will retry. All query results are cached. To remove cached
#' data use hard reset.
#' @param func rentrez function
#' @param args rentrez function arguments, list
#' @param fnm rentrez function name
#' @template ps
#' @return rentrez function results
#' @family run-private
search_and_cache <- function(func, args, fnm, ps) {
  res <- ncbicache_load(fnm = fnm, args = args, wd = ps[['wd']])
  if (!is.null(res)) {
    return(res)
  }
  res <- safely_connect(func = func, args = args, fnm = fnm, ps = ps)
  ncbicache_save(fnm = fnm, args = args, wd = ps[['wd']], obj = res)
  res
}

#' @name safely_connect
#' @title Safely run rentrez function
#' @description Safely run a rentrez function. If the query fails,
#' the function will retry.
#' @param func rentrez function
#' @param args rentrez function arguments, list
#' @param fnm rentrez function name
#' @template ps
#' @return rentrez function results
#' @family run-private
safely_connect <- function(func, args, fnm, ps) {
  res <- NULL
  for (wt_tm in ps[['wt_tms']]) {
    # limit query to 1 hour
    # TODO: allow user interruption with tryCatch?
    query <- try(R.utils::withTimeout(do.call(func, args), timeout = 3600),
                 silent = TRUE)
    #query <- try(do.call(func, args), silent = TRUE)
    if (download_obj_check(query)) {
      res <- query
      break
    } else {
      # ctrl+c
      if (grepl('Operation was aborted by an application callback',
                query[[1]])) {
        stop(query[[1]])
      }
      # too large a request
      if (grepl('the request is too large', query[[1]])) {
        error(ps = ps, 'NCBI is limiting the size of your request. ',
              'Consider reducing btchsz with parameters_reset().')
      }
      info(lvl = 1, ps = ps, "Retrying in [", wt_tm, "s] for [", fnm, ']')
      Sys.sleep(wt_tm)
    }
  }
  res
}

#' @name download_obj_check
#' @title Check an object returned from rentrez function
#' @description Returns T/F. Checks if object returned from rentrez
#' function is as expected.
#' @param obj Object returned from rentrez function
#' @return T/F
#' @family run-private
download_obj_check <- function(obj) {
  if (inherits(x = obj, what = 'try-error')) {
    return(FALSE)
  }
  if (length(obj) == 1 & inherits(x = obj, what = 'character')) {
    if (grepl(pattern = 'timeout', x = obj)) {
      return(FALSE)
    }
  }
  TRUE
}

#' @name batcher
#' @title Download in batches
#' @description Run downloader function in batches for sequences or
#' taxonomic records
#' @return Vector of records
#' @param ids Vector of record ids
#' @param func Downloader function
#' @template ps
#' @template lvl
#' @family run-private
#' @return vector of rentrez function results
batcher <- function(ids, func, ps, lvl = 0) {
  res <- NULL
  n <- length(ids)
  btch <- ps[['btchsz']]
  for (i in seq(0, n - 1, btch)) {
    lower <- i + 1
    upper <- ifelse(i + btch < n, i + btch, n)
    prt_ids <- ids[lower:upper]
    info(lvl = lvl, ps = ps, "[", lower, "-", upper, "]")
    res <- c(res, func(ids = prt_ids, ps = ps))
  }
  res
}
