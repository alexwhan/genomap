#' A wrapper around ASMap::mstmap with some error handling.
#' Associates success or error and p.value with the output.
#'
#' @param obj A cross object suitable for mstmap().
#' @param p.value A value between 0 and 1.
asmap_step <- function(obj, p.value = 1e-6) {
  obj_new <- tryCatch({
    ASMap::mstmap(obj, p.value = p.value, bychr = TRUE)
  }, error = function(e){
    "error"
  })
  if(length(obj_new) == 1) {
    obj_out <- obj
    obj_out$geno <- lapply(obj_out$geno, function(x) {
      x$p.value <- p.value
      x$outcome <- "error"
      return(x)
    })
  } else {
    obj_out <- obj_new
    obj_out$geno <- lapply(obj_new$geno, function(x) {
      x$p.value <- p.value
      x$outcome <- "success"
      return(x)
    })
  }
  return(obj_out)
}

#' Returns the maximum distance between markers for a
#' linkage group.
#'
#' @param obj A cross object suitable for mstmap().
get_lg_stats <- function(obj) {
  lg_stats <- lapply(obj$geno, function(x) {
    if(length(x$map) > 1){
      diff_ord <- diff(x$map)[order(diff(x$map))]
      max_diff <- diff_ord[length(diff_ord)]
      return(max_diff = stats::setNames(max_diff, ""))
    } else {
      return()
    }
  }) %>% unlist()
}

#' Runs mstmap for a linkage groups of map objects that
#' satisfy threshold criteria along sequence of p.values.
#' Intermediate results are kept to trace the point at which
#' linkage groups break.
#' @param obj A cross objects suitable for mstmap.
#' @param comp A comparison object - currently unused
#' @param p.values A decreasing vector of p.values between 0 and 1.
#' @param redo_threshold Distance in cM. If a linkage group has a
#' gap large than this threshold, it will be recalculated in the next cycle.
#' @export
asmap_prog <- function(obj, comp, p.values = 10^-(c(5:10)), redo_threshold = 20) {
  out <- vector(mode = "list", length = length(p.values) + 1)
  out[[1]] <- obj
  i <- 1
  continue <- TRUE
  while(continue) {
    print(paste("i =", i, "p.value =", p.values[[i]], "/n"))
    map_step <- out[[i]]

    #Get chr that need reassessment
    if(i == 1) {
      lgs_redo <- names(map_step$geno)
    } else{
      lg_stats <- get_lg_stats(map_step)
      lgs_redo <- names(lg_stats)[lg_stats > redo_threshold]
    }

    if(length(lgs_redo) == 0) {
      continue <- FALSE
    } else {
      lgs_leave <- names(map_step$geno)[!names(map_step$geno) %in% lgs_redo]

      map_redo_list <- vector(mode = "list", length = length(lgs_redo))
      for(j in 1:length(lgs_redo)) {
        map_redo_list[[j]] <- map_step
        map_redo_list[[j]]$geno <- map_step$geno[names(map_step$geno) == lgs_redo[j]]
      }
      map_leave <- map_step
      if(length(lgs_leave) == 0) {
        map_leave$geno <- NULL
      } else {
        map_leave$geno <- map_leave$geno[names(map_leave$geno) %in% lgs_leave]
      }


      map_redone_list <- lapply(map_redo_list, asmap_step, p.value = p.values[[i]])

      map_redone <- map_redone_list[[1]]
      map_redone$geno <- unlist(lapply(map_redone_list, function(x) {
        return(x$geno)
      }), recursive = FALSE)

      map_redone$geno <- c(map_redone$geno, map_leave$geno)
      map_redone$geno <- map_redone$geno[order(names(map_redone$geno))]
      out[[i + 1]] <- map_redone
    }
    if(i < length(p.values)){
      i <- i + 1
    } else {
      continue <- FALSE
    }
  }
  out <- out[1:(i + 1)]
  return(out)
}
