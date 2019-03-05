#' LLSR's database
#'
#' @name llsr_data
#' @description A database is a highly structured collection of data generally stored and accessed from a computer system and stores raw data and parameters for all analised ATPS phase diagrams. 
#' @format multi-level data.frame()
#' \describe{
#'   \item{db.ref}{Table from LLSR's database in which the references for all used manuscripts are stored.}
#'   \item{db.sys}{Table from LLSR's database in which parameters for all implemented mathematical descriptors are stored. Statistic data is also available.}
#'   \item{db.cas}{Table from LLSR's database in which information regarded all chemicals (such as name and CAS number) used in the collected data.}
#'   \item{db.data}{Table from LLSR's database in which raw experimental data is tabulated. The data was used to calculate all properties made available in this package.}
#'   \item{db.tielines}{Table from LLSR's database in which raw experimental data, used to calculate tielines compositions and slopes, are stored.}
#'   ...
#' }
#' @source \url{https://github.com/diegofcoelho/LLSR/}
#' @examples
#' #
#' XYdt <- llsr_data[["db.data"]][, 1:2]
#' Xdt <- llsr_data[["db.data"]][, 1]
#' Ydt <- llsr_data[["db.data"]][, 2]
NULL
#'
#' Dataset of experimental binodal data of an ATPS
#' 
#' @name peg4kslt
#' @description A dataset containing the experimental binodal data for a PEG/SALT Aqueous Two-Phases System (ATPS)
#' @format multi-level data.frame()
#' \describe{
#'   \item{XC}{Ammonium Sulphate mass fraction}
#'   \item{YC}{Poly(ethylene glycol) mass fraction}
#'   ...
#' }
#' @examples
#' #
#' XYdt <- peg4kslt[,1:2]
#' #
#' Xdt<-peg4kslt[,1]
#' #
#' Ydt<-peg4kslt[,1]
NULL
