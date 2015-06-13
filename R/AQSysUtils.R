AQSysHR <- function (HR){
  if (HR==TRUE){
    par(mar = c(6,6,6,4) + 0.1)
    cex=2.5
    cexlab=2.5
    cexaxis=2.5
    cexmain=2.5
    cexsub=2.5
  }else{
    par(mar = c(5, 4, 4, 2) + 0.1)
    cex=1
    cexlab=1
    cexaxis=1
    cexmain=1
    cexsub=1
  }
  
}

.onAttach <- function(libname, pkgname){
  packageStartupMessage('
    Be aware that LLSR is a collaborative package that still in
    development and your help is essential.\n
    If you found any bugs or have a suggestion, do not hesitate and
    contact us on https://github.com/eqipehub/LLSR/issues.\n
    You also can fork this project directly from github and commit
    improvements to us (https://github.com/eqipehub/LLSR).\n
    The information used in the database was obtained free of charge
    but it might be copyrighted by third parties and references must
    be included appropriately.\n
    Please use LLSR.info() to read more details about the current
    package version.', domain = NULL, appendLF = TRUE)
}

AQSys.merge <- function(wrbk, sheets){
  #
  sys.nrow <- NULL
  sys.mrow <- NULL
  sys.ncol <- NULL
  sys.data <- NULL
  #
  for (nSh in grep("datasource_", sheets)){
    #
    sys.nrow <- nrow(readWorksheet(wrbk, nSh, header = FALSE))
    sys.ncol <- ncol(readWorksheet(wrbk, nSh, header = FALSE))
    #
    if (is.null(sys.mrow)) sys.mrow <- sys.nrow
    if (sys.nrow > sys.mrow) sys.mrow <- sys.nrow 
    if (is.odd(sys.ncol)) AQSys.err("2")
    }
  #
  for (nSh in grep("datasource_", sheets)){
    #
    sys.temp <- readWorksheet(wrbk, nSh, header = FALSE)
    sys.nrow <- nrow(readWorksheet(wrbk, nSh, header = FALSE))
    #
    if (sys.nrow < sys.mrow) {
      sys.temp[sys.mrow, ] <- NA
    }
    if (is.null(sys.data)){
      sys.data <- sys.temp
    } else {
      sys.data <- as.data.frame(c(sys.data, sys.temp),stringsAsFactors=FALSE)
    }
    
  }
  invisible(sys.data)
}
