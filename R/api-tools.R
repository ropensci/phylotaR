#' @name srchNCch
#' @title Run rentrez function and cache results
#' @description Safely run a rentrez function.
#' If the query fails, the function will retry.
#' All query results are cached. To remove cached data
#' use hard reset.
#' @param func rentrez function
#' @param args rentrez function arguments, list
#' @param fnm rentrez function name
#' @param ps Parameters
#' @export
srchNCch <- function(func, args, fnm, ps) {
  res <- ldNcbiCch(fnm=fnm, args=args, wd=ps[['wd']])
  if(!is.null(res)) {
    return(res)
  }
  res <- safeSrch(func=func, args=args, fnm=fnm, ps=ps)
  svNcbiCch(fnm=fnm, args=args, wd=ps[['wd']], obj=res)
  res
}

#' @name safeSrch
#' @title Safely run rentrez function
#' @description Safely run a rentrez function.
#' If the query fails, the function will retry.
#' @param func rentrez function
#' @param args rentrez function arguments, list
#' @param fnm rentrez function name
#' @param ps Parameters
#' @export
safeSrch <- function(func, args, fnm, ps) {
  res <- NULL
  for(wt_tm in ps[['wt_tms']]) {
    query <- try(do.call(func, args), silent=TRUE)
    if(chckSrchObj(query)) {
      res <- query
      break
    } else {
      # ctrl+c
      if(grepl('Operation was aborted by an application callback',
               query[[1]])) {
        stop(query[[1]])
      }
      info(lvl=1, ps=ps, "Retrying in [", wt_tm, "s] for [",
           fnm, ']')
      Sys.sleep(wt_tm)
    }
  }
  res
}


#' @name chckSrchObj
#' @title Check an object returned from rentrez
#' @description Returns T/F. Checks if object
#' returned from rentrez function is as expected.
#' @param obj
#' @export
chckSrchObj <- function(obj) {
  if(inherits(x=obj, what='try-error')) {
    return(FALSE)
  }
  if(length(obj) == 1 & inherits(x=obj, what='character')) {
    if(grepl(pattern='timeout', x=obj)) {
      return(FALSE)
    }
  }
  TRUE
}
