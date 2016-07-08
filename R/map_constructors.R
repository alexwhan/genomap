#' A wrapper around ASMap::mstmap with some error handling,
#' and associating the input p.value with the output.
#'
#' @param obj A cross object suitable for mstmap().
#' @param p.value A value between 0 and 1.
asmap_step <- function(obj, p.value = 1e-6) {
  obj <- tryCatch({
    obj <- ASMap::mstmap(obj, p.value = p.value, bychr = TRUE)
  }, error = function(e){
    "error"
  })
  obj$p.value <- p.value
  return(obj)
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
#' @param p.values A decreasing vector of p.values between 0 and 1.
#' @param redo_threshold Distance in cM. If a linkage group has a
#' gap large than this threshold, it will be recalculated in the next cycle.
#' @export
asmap_prog <- function(obj, comp, p.values = 10^-(c(5:10)), redo_threshold = 20) {
  out <- vector(mode = "list", length = length(p.values + 1))
  out[[1]] <- obj
  k <- 1
  i <- 1
  continue <- TRUE
  while(continue) {
    print(paste("i =", i, "p.value =", p.values[[i]], "/n"))
    map_step <- out[[k]]

    #Get chr that need reassessment
    lg_stats <- get_lg_stats(map_step)
    lgs_redo <- names(lg_stats)[lg_stats > redo_threshold]
    if(length(lgs_redo) == 0) {
      continue <- FALSE
    } else {
      lgs_leave <- names(map_step$geno)[!names(map_step$geno) %in% lgs_redo]

      map_redo <- map_step
      map_redo$geno <- map_redo$geno[names(map_redo$geno) %in% lgs_redo]
      map_leave <- map_step
      map_leave$geno <- map_leave$geno[names(map_leave$geno) %in% lgs_leave]

      map_redone <- asmap_step(map_redo, p.value = p.values[[i]])

      if(names(map_redone)[1] == "geno") {
        k <- i + 1
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
  }
  return(out)
}
