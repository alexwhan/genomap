#' Checks genotype scores for heterozygosity
#'
#' Takes a vector of genotype scores (each made up of a two character string)
#' and a character class to match against, and returns a logical value depending
#' on whether the two characters match (\code{FALSE}) or not (\code{TRUE}).
#'
#' @param score A vector of genotype scores, each of which needs to be a two
#'   character string. Numeric values are converted to strings
#' @param characterClass A regex character class to match against. If characters
#'   in \code{score} do not exist in the character class, an error will be
#'   returned
#' @export
#' @examples
#' isHet("AA")
#' isHet(c("AA", "AB", "BB"))
#' isHet("GT", "[AB]")
isHet <- function(score, characterClass = "[[:alnum:]]") {
  if(any(nchar(score) != 2)) stop("Expecting a string of length 2")
  if(any(grepl(sub("^(\\[)", "\\1^", characterClass), score))) {
    stop("There are genotype scores that do not match the character class")
  }
  !grepl(paste0("(", characterClass, ")(?=\\1)"), score, perl = T)
}
