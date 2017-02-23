#' Get linkage group extents
#'
#' @param obj An object of class cross, map or tidy_gen_df
#'
#' @return map coordinates as a vector
#' @export
get_lg_lengths <- function(obj) {
  UseMethod("get_lg_lengths")
}

#' @export
get_lg_lengths.cross <- function(obj) {
  obj <- genomap::map2df(obj)
  get_lg_lengths(obj)
}

#' @export
get_lg_lengths.map <- function(obj) {
  obj <- genomap::map2df(obj)
  get_lg_lengths(obj)
}

#' @export
get_lg_lengths.mpcross <- function(obj) {
  obj <- genomap::map2df(obj)
  get_lg_lengths(obj)
}

globalVariables("mapdist")
#' @export
get_lg_lengths.tidy_gen_map <- function(obj) {
  df <- dplyr::group_by_(obj, "lg")
  dplyr::summarise(df, lg_length = max(mapdist))
}

#' Title
#'
#' @param obj An object of class cross, map or tidy_gen_df
#'
#' @return A vector of linkage group names
#' @export
get_lg_names <- function(obj) {
  UseMethod("get_lg_names")
}

#' @export
get_lg_names.cross <- function(obj) {
  return(names(obj$geno))
}

#' @export
get_lg_names.map <- function(obj) {
  return(names(obj))
}

#' @export
get_lg_names.mpcross <- function(obj) {
  return(names(obj$map))
}

#' @export
get_lg_names.tidy_gen_map <- function(obj) {
  return(unique(obj$lg))
}

#' Get offsets for linkage groups
#'
#' @param obj An object of class cross, map or tidy_gen_df
#' @param order A vector of linkage group names to sort map by
#'
#' @return offsets as a vector
#' @export
get_lg_offsets <- function(obj, order = NULL) {
  lg_lengths <- get_lg_lengths(obj)
  
  if(!is.null(order)) {
    stopifnot(class(order) == "character")
    stopifnot(all(order %in% lg_lengths$lg))
    stopifnot(all(lg_lengths$lg %in% order))
    lg_lengths <- lg_lengths[match(order, lg_lengths$lg)]
  }
  
  lg_lengths$lg_offset <- cumsum(c(0, lg_lengths$lg_length[1:nrow(lg_lengths) - 1]))
  lg_lengths
}