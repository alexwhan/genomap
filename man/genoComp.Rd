% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genoconverters.R
\name{genoComp}
\alias{genoComp}
\title{Take a linkage group from a cross object, convert into long form and sort
according to genotype scores for a marker of interest.}
\usage{
genoComp(lgObject, markerOfInterest)
}
\arguments{
\item{lgObject}{A linkage group from a cross object. Usually in the form
cross$geno$name.}

\item{markerOfInterest}{A string identifying the marker by which the data
should be sorted.}
}
\description{
Take a linkage group from a cross object, convert into long form and sort
according to genotype scores for a marker of interest.
}
