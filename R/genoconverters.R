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
#' Works one marker at a time. Takes absolute genotype scores and converts to relative scores (AA/AB/BB to
#' match maternal/het/paternal).
#'
#' @param data A data frame
#' @param markerName Specification of the column containing the marker name
#' @param maternal Specification of the column containing the maternal score
#' @param paternal Specification of the column containing the paternal score
#' @param ... Specification of the columns containing progeny scores
#' @param missingString A string representing a missing score
#' @export
convert_rel <- function(data, markerName, maternal, paternal, ..., missingString = "--") {
  markerName_ <- col_name(substitute(markerName))
  maternal_ <- col_name(substitute(maternal))
  paternal_ <- col_name(substitute(paternal))

  if (n_dots(...) == 0) {
    progeny_cols_ <- setdiff(colnames(data), c(markerName_, maternal_, paternal_))
  } else {
    progeny_cols_ <- unname(dplyr::select_vars(colnames(data), ...))
  }

  convert_rel_(data, markerName_, maternal_, paternal_, progeny_cols_, missingString)
}

#' Convert absolute scores to relative (standard-evaluation)
#'
#' This is a S3 generic.
#'
#' @param data A data frame
#' @param markerName_,maternal_,paternal_ Strings giving names of marker
#'    maternal and paternal columns
#' @param progeny_cols_ Character vector giving column names of progeny scores
#' @param missingString A string representing a missing score
#' @keywords internal
#' @export
convert_rel_ <- function(data, markerName_, maternal_, paternal_, progeny_cols_, missingString = "--") {
  UseMethod("convert_rel_")
}

#' @export
convert_rel_.data.frame <- function(data, markerName_, maternal_, paternal_, progeny_cols_,
                                    missingString = "--") {
  if (length(progeny_cols_) == 0) {
    warning("No progeny columns were defined, return original data frame")
    return(data)
  }

  progeny_idx <- match(progeny_cols_, names(data))
  if (anyNA(progeny_idx)) {
    missing_cols <- paste(progeny_cols_[is.na(progeny_idx)], collapse = ",")
    stop("Unknown column names: ", missing_cols, call. = FALSE)
  }

  if (any(unlist(lapply(data[progeny_idx], class)) != "character")) {
    stop("Progeny columns are not all character class")
  }

  maternal_idx <- match(maternal_, names(data))
  if (!inherits(data[[maternal_]], "character")) {
    stop("The maternal column is not of class character")
  }
  if (any(nchar(unlist(data[data[[maternal_]] != missingString, maternal_idx])) != 2)) {
    stop("There are non-missing maternal scores that have more or less than two characters")
  }

  paternal_idx <- match(paternal_, names(data))
  if (!inherits(data[[paternal_]], "character")) {
    stop("The paternal column is not of class character")
  }
  if (any(nchar(unlist(data[data[[paternal_]] != missingString, paternal_idx])) != 2)) {
    stop("There are non-missing paternal scores that have more or less than two characters")
  }

  data$converted <- apply(data, 1, function(x) {
    converted <- convertScore(x[maternal_idx], x[paternal_idx], x[progeny_idx], missingString = missingString)
    return(data.frame(t(converted)))
  })

  # data_out <- data %>%
  #   dplyr::select_(markerName_, "converted") %>%
  #   tidyr::unnest()

  return(data)
}

#' Gets allele stats
#'
#' @param data A data.frame.
#' @param converted Specify the list column containing the converted allele calls (outptu from convert_rel)
#' @export
allele_stats <- function(data, converted) {
  if(!inherits(data, "data.frame")) stop("data needs to be a data.frame")

  converted_ <- col_name(substitute(converted))

  if(!inherits(data[[converted_]], "list")) stop("The converted column should be a list")

  allele_stats_(data, converted_)
}

#' Gets allele stats
#'
#' @inheritParams allele_stats
#' @importFrom magrittr %>%
#' @export
allele_stats_ <- function(data, converted_) {
  if(!inherits(data, "data.frame")) stop("data needs to be a data.frame")

  if(!inherits(data[[converted_]], "list")) stop("The converted column should be a list")

  converted_df <- data %>%
    dplyr::select_(converted_) %>%
    tidyr::unnest()

  data$allele_stats <- apply(converted_df, 1, function(x) {
    scores <- sum(!is.na(x))
    missing <- (length(x) - scores)/length(x)
    if(scores == 0) {
      AA <- AB <- BB <- minor_allele <- minor_allele_freq <- NA
    } else {
      AA <- sum(x == "AA")/scores
      AB <- sum(x == "AB")/scores
      BB <- sum(x == "BB")/scores
      if(AA == BB) {
        minor_allele <- "equal"
        minor_allele_freq <- AA
      } else {
        minor_allele <- c("AA", "BB")[which.min(c(AA, BB))]
        minor_allele_freq <- c(AA, BB)[which.min(c(AA, BB))]
      }
    }
    data_out <- data.frame(AA = AA, AB = AB, BB = BB,
                           minor_allele = minor_allele,
                           minor_allele_freq = minor_allele_freq)
  })

  return(data)

}

#' Checks informativeness of parents
#'
#' @param data A data.frame where each row is a marker score
#' @param maternal,paternal Strings giving names of maternal, paternal columns
#' @param missingString A string defining missing scores
#' @importFrom magrittr %>%

check_parents_ <- function(data, maternal, paternal, missingString = "--") {
  if(length(dplyr::setdiff(c(maternal, paternal), colnames(data))) != 0) stop("check names for parents in check_parents_")
  other_cols <- dplyr::setdiff(colnames(data), c(maternal, paternal))
  mutate_call <- lazyeval::interp(~ ifelse(any(missingString %in% c(maternal, paternal)),
                                           "missing",
                                           ifelse(any(isHet(c(maternal, paternal))),
                                                  "heterozygous",
                                                  ifelse(maternal == paternal,
                                                         "nonseg",
                                                         ifelse(all(isHomo(c(maternal, paternal))),
                                                                "informative", "unknown")))),
                                  # missingString = as.name(missingString),
                                  maternal = as.name(maternal),
                                  paternal = as.name(paternal))


  data_out <- data %>%
    # dplyr::select_(maternal, paternal) %>%
    dplyr::rowwise() %>%
    dplyr::mutate_(.dots = stats::setNames(list(mutate_call), "parent_status")) %>%
    dplyr::select_(maternal, paternal, "parent_status", .dots = other_cols)

  return(data_out)
}

#' Converts a single marker, typically a row in a data.frame.
#'
#' @param maternal,paternal A two character string, or missingString value.
#' @param progeny A (named) character vector of progeny scores.
#' @param missingString A string defining missing scores.
#'
convertScore <- function(maternal, paternal, progeny, missingString = "--") {
  if(is.null(names(progeny))) names(progeny) <- paste0("progeny", 1:length(progeny))
  progeny.out <- vector(mode = "character", length = length(progeny)) %>%
    stats::setNames(names(progeny))
  if(any(isHet(c(maternal, paternal))) |
     any(grepl(missingString, c(maternal, paternal))) |
     maternal == paternal) {
    progeny.out[1:length(progeny)] <- NA
  } else {

    #Create possible het codes from paternal
    #This is to avoid a progeny call that doesn't match parentals being called AB
    hetmp <- paste0(sub("^([[:alnum:]]).", "\\1", maternal), sub("^([[:alnum:]]).", "\\1", paternal))
    hetpm <- paste0(sub("^([[:alnum:]]).", "\\1", paternal), sub("^([[:alnum:]]).", "\\1", maternal))

    progeny.out[progeny == maternal] <- "AA"
    progeny.out[progeny == paternal] <- "BB"
    progeny.out[progeny == missingString] <- missingString
    progeny.out[progeny %in% c(hetmp, hetpm)] <- "AB"
    progeny.out[progeny.out == ""] <- NA
  }
  return(progeny.out)
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
