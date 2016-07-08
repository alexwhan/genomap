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
#' @obj A cross object suitable for mstmap().
get_lg_stats <- function(obj) {
  lg_stats <- lapply(obj$geno, function(x) {
    if(length(x$map) > 1){
      diff_ord <- diff(x$map)[order(diff(x$map))]
      max_diff <- diff_ord[length(diff_ord)]
      return(max_diff = setNames(max_diff, ""))
    } else {
      return()
    }
  }) %>% unlist()
}
