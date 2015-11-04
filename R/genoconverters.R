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

#' Checks genotype scores for homozygosity
#'
#' @param score A vector of genotype scores, each of which needs to be a two
#'   character string. Numeric values are converted to strings
#' @param characterClass A regex character class to match against. If characters
#'   in \code{score} do not exist in the character class, an error will be
#'   returned
#' @export
isHomo <- function(score) {
  if(any(nchar(score) != 2)) stop("Expecting a string of length 2")
  if(any(grepl(sub("^(\\[)", "\\1^", characterClass), score))) {
    stop("There are genotype scores that do not match the character class")
  }
  grepl(paste0("(", characterClass, ")(?=\\1)"), score, perl = T)
}


#' Convert absolute scores to relative
#'
#' Takes absolute genotype scores and converts to relative scores (AA/AB/BB to
#' match maternal/het/paternal).
#'
#' @param parent1
#' @export
convertScore <- function(maternal, paternal, progeny, markerName = "unknown", missingString = "--") {
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
  missing <- sum(progeny == missingString)/length(progeny)
  other <- sum(is.na(progeny.out))/length(progeny.out)

  maternal.out <- setNames("AA", names(maternal))
  paternal.out <- setNames("BB", names(paternal))
  names(progeny.out) <- names(progeny)

  return(c(markerName, AA = maternal.prop,
           AB = het.prop,
           BB = paternal.prop,
           missing = missing,
           other = other,
           maternal.out, paternal.out, progeny.out))
}


#' Convert R/qtl map into data.frame
#'
#' Similar to pull.map, but returns a more useful object.
#'
#' @param map An R/qtl cross object
#' @export
map2df <- function(map) {
  mapdf <- lapply(map$geno, function(x) {
    out <- as.data.frame(x$map)
    out$markerName <- rownames(out)
    return(out)
  })
  mapdf <- bind_rows(mapdf, .id = "lg") %>%
    rename(mapdist = `x$map`)
  return(mapdf)
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
revmaplgs <- function(maps.comp_df, refmapid, revmapid, revmap) {
  refmapdist_ <- paste0(deparse(substitute(refmapid)), "_mapdist")
  map2dist_ <- paste0(deparse(substitute(revmapid)), "_mapdist")
  refmaplg_ <- paste0(deparse(substitute(refmapid)), "_lg")
  map2lg_ <- paste0(deparse(substitute(revmapid)), "_lg")

  revs <- maps.comp_df %>%
    group_by_(map2lg_) %>%
    arrange_(refmapdist_) %>%
    mutate_(dist2 = map2dist_) %>%
    summarise(rev = ifelse(cor(1:n(), dist2) < 0, TRUE, FALSE)) %>%
    filter(rev)

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
    rename_(.dots = with(map1, setNames(c(map1markerName_, map1distName_, map1lgName_),
                                        c("markerName", paste(map1name, c("mapdist", "lg"), sep = "_")))))
  map2.r <- map2 %>%
    rename_(.dots = with(map2, setNames(c(map2markerName_, map2distName_, map2lgName_),
                                        c("markerName", paste(map2name, c("mapdist", "lg"), sep = "_")))))
  out <- map1.r %>%
    left_join(map2.r)
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
#' @export
longmaps <- function(df, reflg_facet, reflg_join, markerName = markerName) {
  reflg_join_ <- deparse(substitute(reflg_join))
  markerName_ <- deparse(substitute(markerName))
  df$markerName.lg <- paste0(df[[markerName_]], "_", df[[reflg_join_]])
  lgs.g <- df %>%
    select(contains("_lg"), markerName.lg)
  names(lgs.g) <- gsub("_.*", "", names(lgs.g))
  lgs.g <- lgs.g %>%
    gather(source, lg, -contains("markerName.lg"))

  dist.g <- df %>%
    select(contains("_mapdist"), markerName.lg)
  names(dist.g) <- gsub("_.*", "", names(dist.g))
  dist.g <- dist.g %>%
    gather(source, mapdist, -contains("markerName.lg"))

  lgs.c <- lgs.g %>%
    inner_join(dist.g)

  reflg_facet_ <- deparse(substitute(reflg_facet))
  out.df <- df %>%
    select_(reflg_facet_, "markerName.lg") %>%
    left_join(lgs.c) %>%
    filter(!is.na(lg))
}

#' Take a linkage group from a cross object, conver into long form and sort
#' according to genotype scores for a marker of interest.
#'
#' @param lgObject A linkage group from a cross object. Usually in the form
#'   cross$geno$name
#' @param markerOfInterest A string identifying the marker by which the data
#'   should be sorted
genoComp <- function(lgObject, markerOfInterest) {
  df <- as.data.frame(lgObject)
  df$Genotype <- rownames(lgObject)
  dfg <- df %>%
    gather(marker, score, -Genotype)
  ord <- dfg %>%
    filter(marker == markerOfInterest) %>%
    arrange(score) %>%
    mutate(ord = order(score)) %>%
    arrange(ord)
  dfg$Gen.sort <- factor(dfg$Genotype, levels = ord$Genotype, ordered = TRUE)
  return(dfg)
}

#' Plot genotype data for a linkage group
#'
#' @param df Sorted data.frame output from /code{genoComp}
genoPlot <- function(df) {
  df.p <- df %>%
    filter(!is.na(score)) %>%
    ggplot(aes(marker, Gen.sort)) + geom_tile(aes(fill = factor(score))) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45), panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 10))

  return(df.p)
}

#' Compares segregation of markers based on agreement of two markers (population parents)
genoSeg <- function(df, genos) {
  df.o <- df %>%
    mutate(ord = H45 == `SWDH0020-1_90K`) %>%
    arrange(desc(ord), H45)

  dfg <- df %>%
    gather(Genotype, score, -markerName)

  dfg$marker.sort <- factor(dfg$markerName, levels = rev(df.o$markerName), ordered = T)

  df.p <- dfg %>%
    filter(!is.na(score), score != "--") %>%
    ggplot(aes(marker.sort, Genotype)) + geom_tile(aes(fill = factor(score))) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) + xlab("Marker")
}
