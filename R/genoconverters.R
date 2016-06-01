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
#' @param nonMatch c("false", "error") a string defining whether to return false
#'   or throw an error if scores do not match characterClass.
#' @export
#' @examples
#' isHet("AA")
#' isHet(c("AA", "AB", "BB"))
#' isHet("GT", "[AB]")
isHet <- function(score, characterClass = "[[:alnum:]]", nonMatch = "false") {
  if(class(score) == "numeric") warning("Numeric input will be coerced to character")
  if(any(nchar(score) != 2)) stop("Expecting a string of length 2")
  if(!nonMatch %in% c("false", "error")) stop("nonMatch should be either 'false' or 'error'")
  if(nonMatch == "error" & any(grepl(sub("^(\\[)", "\\1^", characterClass), score))) {
    stop("There are genotype scores that do not match the character class")
  }
  grepl(paste0("^(", characterClass, ")(?!\\1)"), score, perl = T)
}

#' Checks genotype scores for homozygosity
#'
#' @param score A vector of genotype scores, each of which needs to be a two
#'   character string. Numeric values are converted to strings
#' @param characterClass A regex character class to match against. If characters
#'   in \code{score} do not exist in the character class, an error will be
#'   returned
#' @param nonMatch c("false", "error") a string defining whether to return false
#'   or throw an error if scores do not match characterClass.
#' @export
isHomo <- function(score, characterClass = "[[:alnum:]]", nonMatch = "false") {
  if(class(score) == "numeric") warning("Numeric input will be coerced to character")
  if(any(nchar(score) != 2)) stop("Expecting a string of length 2")
  if(!nonMatch %in% c("false", "error")) stop("nonMatch should be either 'false' or 'error'")
  if(nonMatch == "error" & any(grepl(sub("^(\\[)", "\\1^", characterClass), score))) {
    stop("There are genotype scores that do not match the character class")
  }
  grepl(paste0("(", characterClass, ")(?=\\1)"), score, perl = T)
}


#' Convert absolute scores to relative
#'
#' Takes absolute genotype scores and converts to relative scores (AA/AB/BB to
#' match maternal/het/paternal).
#'
#' @param maternal The maternal score. A string with nchar(maternal) == 2
#' @param paternal The paternal score. A string with nchar(paternal) == 2
#' @param progeny A vector of progeny scores
#' @param markerName The name of the marker that the scores are for
#' @param missingString A string representing a missing score
#' @export
convertScore <- function(maternal, paternal, progeny, markerName = "unknown", missingString = "--") {
  if(class(maternal) != "character" | class(paternal) != "character") stop("maternal of paternal score is not character class")
  if(length(maternal) != 1) stop("maternal score is not of length == 1")
  if(length(paternal) != 1) stop("paternal score is not of length == 1")
  #Check both parent scores are homozygous
  if(any(isHet(c(maternal, paternal)))) stop("Parental scores are not homozygous")
  #check both parent scores are informative
  if(any(grepl(missingString, c(maternal, paternal)))) stop("Parental scores are not fully informative")
  #Check if any maternal equal paternal
  if(any(maternal == paternal)) stop("Some parental scores are non-segregating")

  #Create possible het codes from paternal
  #This is to avoid a progeny call that doesn't match parentals being called AB
  hetmp <- paste0(sub("^([[:alnum:]]).", "\\1", maternal), sub("^([[:alnum:]]).", "\\1", paternal))
  hetpm <- paste0(sub("^([[:alnum:]]).", "\\1", paternal), sub("^([[:alnum:]]).", "\\1", maternal))

  progeny.out <- vector(mode = "character", length = length(progeny))
  progeny.out[progeny == maternal] <- "AA"
  progeny.out[progeny == paternal] <- "BB"
  progeny.out[progeny == missingString] <- missingString
  progeny.out[progeny %in% c(hetmp, hetpm)] <- "AB"

  progeny.total <- sum(progeny %in% c(maternal, paternal, hetmp, hetpm))
  maternal.prop <- sum(progeny == maternal)/progeny.total
  paternal.prop <- sum(progeny == paternal)/progeny.total
  het.prop <- sum(progeny %in% c(hetmp, hetpm))/progeny.total
  minor.allele <- ifelse(all(!is.na(c(maternal.prop, paternal.prop))),
                         ifelse(maternal.prop < paternal.prop, "AA", 
                                ifelse(maternal.prop == paternal.prop, "equal", "BB")), NA)
  minor.allele.freq <- ifelse(all(!is.na(c(maternal.prop, paternal.prop))), 
                              c(maternal.prop, paternal.prop)[which.min(c(maternal.prop, paternal.prop))], NA)
  missing <- sum(progeny == missingString)/length(progeny)
  other <- sum(is.na(progeny.out))/length(progeny.out)

  maternal.out <- stats::setNames("AA", names(maternal))
  paternal.out <- stats::setNames("BB", names(paternal))
  names(progeny.out) <- names(progeny)
  

  return(data.frame(markerName = markerName, AA = maternal.prop,
           AB = het.prop,
           BB = paternal.prop,
           missing = missing,
           other = other,
           minor_allele = minor.allele,
           minor_allele_freq = minor.allele.freq,
           progeny_scores = I(list(progeny.out))))
}


#' Convert R/qtl map into data.frame
#'
#' Similar to pull.map, but returns a more useful object.
#'
#' @param map An R/qtl cross object
#' @export
#' @importFrom magrittr %>%
map2df <- function(map) {
  mapdf <- lapply(map$geno, function(x) {
    out <- as.data.frame(x$map)
    out$markerName <- rownames(out)
    return(out)
  })
  mapdf <- dplyr::bind_rows(mapdf, .id = "lg") %>%
    dplyr::rename_("mapdist" = quote(`x$map`))
  return(mapdf)
}

#' Convert R/qtl genotyping data into data.frame
#'
#' Similar to map2df, but for the genotying values.
#'
#' @param map An R/qtl cross object
#' @export

geno2df <- function(map) {
  genodf <- lapply(map$geno, function(x) {
    out <- as.data.frame(x$data)
    markers <- colnames(out)
    out$id <- rownames(out)
    out.g <- tidyr::gather_(out, "markerName", "score", colnames)
    return(out.g)
  })
  genodf <- dplyr::bind_rows(genodf, .id = "lg")
  return(genodf)
}

#' Take output from /code{join2maps()} and a map object to be reversed, and
#' returns the map object with appropriate groups reversed.
#'
#' @param maps.comp_df Output from /code{join2maps()}.
#' @param refmapid  The reference map name to identify columns in maps.comp_df.
#' @param revmapid The name of the map to be reversed, used to identify columns
#'   in maps.comp_df.
#' @param revmap a cross object which is to have linkage groups reversed.
#' @export
#' @importFrom magrittr %>%
revmaplgs <- function(maps.comp_df, refmapid, revmapid, revmap) {
  refmapdist_ <- paste0(deparse(substitute(refmapid)), "_mapdist")
  map2dist_ <- paste0(deparse(substitute(revmapid)), "_mapdist")
  refmaplg_ <- paste0(deparse(substitute(refmapid)), "_lg")
  map2lg_ <- paste0(deparse(substitute(revmapid)), "_lg")

  dist2_var <- "dist2"

  revs <- maps.comp_df %>%
    dplyr::group_by_(map2lg_) %>%
    dplyr::arrange_(refmapdist_) %>%
    dplyr::mutate_(dist2 = map2dist_) %>%
    dplyr::summarise(rev = ifelse(stats::cor(1:dplyr::n(),
                                      lazyeval::interp(~var, var = as.name(dist2_var))) < 0,
                                      TRUE, FALSE)) %>%
    dplyr::filter(rev)

  outmap <- revmap(revmap, revs[[1]])
  return(outmap)
}

#' Provide the 'mirror image' of a vector, especially of positions.
#'
#' @param vec A (possibly named) numeric vector.
#' @export
#' @examples
#' distvec <- c(A = 0, B = 3, C = 10)
#' revvec(distvec)
revvec <- function(vec) {
  if(!is.numeric(vec)) stop("The vector needs to be numeric")
  mid <- (min(vec, na.rm = T) + max(vec, na.rm = T))/2
  newvec <- vec + 2 * (mid - vec)
  newvec <- rev(newvec)
  return(newvec)
}

#' Reverse particular linkage groups in a map from a cross object
#'
#' @param map A cross objects
#' @param revgroups The linkage groups to be reversed
revmap <- function(map, revgroups) {
  if(any(!revgroups %in% names(map$geno))) stop("The supplied linkage groups do not all exist")
  for(i in revgroups) {
    map$geno[[i]]$data <- map$geno[[i]]$data[,ncol(map$geno[[i]]$data):1]
    map$geno[[i]]$map <- revvec(map$geno[[i]]$map)
  }
  return(map)
}

#' Join two map dfs
#'
#' Take two maps that have been converted to df (possibly through
#' /code{map2df()}) and join them.
#' @param map1 The data.frame for the first map.
#' @param map2 The data.frame for the second map.
#' @param map1markerName An unquoted string defining the name for the marker
#'   column in map1.
#' @param map2markerName An unquoted string defining the name for the marker
#'   column in map2.
#' @param map1distName An unquoted string defining the name for the distance
#'   column in map1.
#' @param map2distName An unquoted string defining the name for the distance
#'   column in map2.
#' @param map1lgName An unquoted string defining the name for the linkage group
#'   column in map1.
#' @param map2lgName An unquoted string defining the name for the linkage group
#'   column in map2.
#' @export
#' @importFrom magrittr %>%
join2maps <- function(map1, map2, map1markerName, map2markerName, map1distName, map2distName, map1lgName, map2lgName) {
  map1name <- deparse(substitute(map1))
  map1markerName_ <- deparse(substitute(map1markerName))
  map1distName_ <- deparse(substitute(map1distName))
  map1lgName_ <- deparse(substitute(map1lgName))
  map2name <- deparse(substitute(map2))
  map2markerName_ <- deparse(substitute(map2markerName))
  map2distName_ <- deparse(substitute(map2distName))
  map2lgName_ <- deparse(substitute(map2lgName))
  map1.r <- map1 %>%
    dplyr::rename_(.dots = with(map1, stats::setNames(c(map1markerName_, map1distName_, map1lgName_),
                                        c("markerName", paste(map1name, c("mapdist", "lg"), sep = "_")))))
  map2.r <- map2 %>%
    dplyr::rename_(.dots = with(map2, stats::setNames(c(map2markerName_, map2distName_, map2lgName_),
                                        c("markerName", paste(map2name, c("mapdist", "lg"), sep = "_")))))
  out <- map1.r %>%
    dplyr::left_join(map2.r)
  return(out)
}

#' Take output from /code{join2maps()} and put into long form for plotting
#'
#' Refers to a reference linkage group for facetting, and one for joining,
#' this avoids errors when a marker exists twice in one of the maps.
#' @param df A data.frame (output from /code{join2maps()}) of joined map data.
#' @param reflg_facet An unquoted string. The name for the reference linkage group identifier for
#' faceting.
#' @param reflg_join An unquoted string. The name for the reference linkage group identifier for
#' joining.
#' @param markerName An unquoted string. The name of the variable giving marker names.
#' @export
#' @importFrom magrittr %>%
longmaps <- function(df, reflg_facet, reflg_join, markerName = markerName) {
  reflg_join_ <- deparse(substitute(reflg_join))
  markerName_ <- deparse(substitute(markerName))
  df$markerName.lg <- paste0(df[[markerName_]], "_", df[[reflg_join_]])
  lgs.g <- df %>%
    dplyr::select_(names(df)[grepl("_lg", names(df))], "markerName.lg")

  names(lgs.g) <- gsub("_.*", "", names(lgs.g))

  lgs.g <- lgs.g %>%
    tidyr::gather_("source", "lg", names(lgs.g)[!grepl("markerName.lg", names(lgs.g))])

  dist.g <- df %>%
    dplyr::select_(names(df)[grepl("_mapdist", names(df))], "markerName.lg")
  names(dist.g) <- gsub("_.*", "", names(dist.g))
  dist.g <- dist.g %>%
    tidyr::gather_("source", "mapdist",
                   names(dist.g)[!grepl("markerName.lg", names(dist.g))])

  lgs.c <- lgs.g %>%
    dplyr::inner_join(dist.g)

  reflg_facet_ <- deparse(substitute(reflg_facet))
  out.df <- df %>%
    dplyr::select_(reflg_facet_, "markerName.lg") %>%
    dplyr::left_join(lgs.c)
  out.df <- out.df[!is.na(out.df$lg),]

}

#' Take a linkage group from a cross object, convert into long form and sort
#' according to genotype scores for a marker of interest.
#'
#' @param lgObject A linkage group from a cross object. Usually in the form
#'   cross$geno$name.
#' @param markerOfInterest A string identifying the marker by which the data
#'   should be sorted.
#' @importFrom magrittr %>%
genoComp <- function(lgObject, markerOfInterest) {
  df <- as.data.frame(lgObject$data)
  df$Genotype <- rownames(lgObject$data)
  dfg <- df %>%
    tidyr::gather_("markerName", "score",
                   dplyr::select_vars_(names(df),
                                       names(df),
                                       exclude = "Genotype"))
  ord <- dfg %>%
    dplyr::filter_("markerName" == markerOfInterest) %>%
    dplyr::arrange_("score") %>%
    dplyr::mutate_(ord = order(as.name("score"))) %>%
    dplyr::arrange(ord)
  dfg$Gen.sort <- factor(dfg$Genotype, levels = ord$Genotype, ordered = TRUE)
  dfg$Mark.sort <- factor(dfg$markerName, levels = names(lgObject$map), ordered = TRUE)
  return(dfg)
}

#' Converts a character vector of genotype scores into two classes - one of which represents one homozygous state, and the other represents the other homozygous state and all heterozygotes
#' @param vec A character vector of genotype scores, each of which is made of two characters. An na_string may also be present
#' @param level Either 1 or 2, indicating which homozygous state is to be compared to the other state/heterozygotes. Default is 1
#' @param na_string a string representing NA. Default is "u"
#' @export

make_nptest_classes <- function(vec, level = 1, na_string = "u") {
  stopifnot(class(vec) == "character")
  vec[vec == na_string] <- NA
  vec_fac <- as.factor(vec)
  if(any(genomap::isHet(vec)) & length(levels(vec_fac)) > 3) {
    return(rep("too many", length(vec)))
    #stop("Too many levels in vec")
  }
  if(!any(genomap::isHet(vec)) & length(levels(vec_fac)) > 2) {
    return(rep("too many", length(vec)))
    #stop("Too many levels in vec")
  }
  if(length(levels(vec_fac)) < 2) {
    return(rep(NA_character_, length(vec)))
  }

  stopifnot(length(levels(vec_fac)) >= 2)

  if(level == 1) {
    main_level <- levels(vec_fac)[1]
    gen_levels <- c(main_level, paste0("not_", main_level))
    if(any(genomap::isHet(vec))) {
      out_vec <- as.numeric(vec_fac)
      out_vec[out_vec > 1] <- 2
    }

    if(!any(genomap::isHet(vec))) {
      out_vec <- as.numeric(vec_fac)
    }
  }

  if(level == 2) {

    if(any(genomap::isHet(vec))) {
      if(length(levels(vec_fac)) < 3) {
        return(rep(NA_character_, length(vec)))
      } else {
        main_level <- levels(vec_fac)[3]
        gen_levels <- c(main_level, paste0("not_", main_level))
        out_vec <- as.numeric(vec_fac)
        out_vec[out_vec == 1] <- 2
        out_vec[out_vec == 3] <- 1
      }
    }

    if(!any(genomap::isHet(vec))) {
      main_level <- levels(vec_fac)[2]
      gen_levels <- c(main_level, paste0("not_", main_level))
      out_vec <- as.numeric(vec_fac)
      out_vec[out_vec == 1] <- 3
      out_vec[out_vec == 2] <- 1
      out_vec[out_vec == 3] <- 2
    }
  }
  out_vec <- factor(out_vec, labels = gen_levels)
  return(as.character(out_vec))
}
