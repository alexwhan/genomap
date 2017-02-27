#' A biparental cross object
#'
#' A cross produced by qtl
#'
#' @format A list of class 'cross' produced by the qtl package
#' \describe{
#' \item{geno}{A list of the genotyping data}
#' \item{map}{A list describing the genetic map}
#' }
"bp_cross"

#' A biparental map object
#'
#' A map object produced by qtl. This is the same object as
#' genomap::bp_cross$map
#'
#' @format A list of class 'map'
#' \describe{
#' \item{1}{A vector of class 'A' describing a linkage group}
#' ...
#' }
"bp_map"

#' A multiparental cross object
#'
#' A 4-way mpcross object produced by mpMap
#'
#' @format A list of class 'mpcross'
#' \describe{
#' \item{founders}{The genotyping of the founders}
#' \item{map}{A list describing the genetic map}
#' ...
#' }
"m4_cross"

#' A multiparental cross object
#'
#' An 8-way mpcross object produced by mpMap
#'
#' @format A list of class 'mpcross'
#' \describe{
#' \item{founders}{The genotyping of the founders}
#' \item{map}{A list describing the genetic map}
#' ...
#' }
"m8_cross"
