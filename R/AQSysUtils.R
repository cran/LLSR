#' @importFrom dplyr count_
####################################################################################################################
# Set plot area to export high resolution pictures
AQSysHR <- function (HR) {
  par.defaults <- par(no.readonly = TRUE)
  par(par.defaults)
  if (HR == TRUE) {
    clwd = 10
    par(
      mar = c(5, 6, 5, 3) + 0.5,
      cex = 2.5,
      cex.lab = 2.5,
      cex.axis = 2.5,
      cex.main = 2.5,
      cex.sub = 2.5,
      lwd = clwd,
      mgp = c(4, 2, 0)
    )
  } else{
    clwd = 1
    par(
      mar = c(5, 4, 4, 2) + 0.1,
      cex = 1,
      cex.lab = 1,
      cex.axis = 1,
      cex.main = 1,
      cex.sub = 1,
      lwd = clwd
    )
  }
  return(clwd)
}
# Display message when package is loaded
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    '
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
    package version.',
    domain = NULL,
    appendLF = TRUE
  )
}
#
getBNDL <- function(workBook, sheets) {
  # initialize variables
  sys.nrow <- NULL
  sys.mrow <- NULL
  sys.ncol <- NULL
  sys.data <- NULL
  # make a loop through all sheets that satisfy the condition
  # the integrity of the worksheets will be evaluated before the merge
  # and determine which system in the workbook have bigger dataset (mrow)
  for (SheetIndex in grep("BINODAL", sheets)) {
    # determine the number of row and columns in the worksheet
    sys.nrow <- nrow(readWorksheet(workBook, SheetIndex, header = FALSE))
    sys.ncol <- ncol(readWorksheet(workBook, SheetIndex, header = FALSE))
    # initialize variables for first run only
    if (is.null(sys.mrow))
      sys.mrow <- sys.nrow
    if (sys.nrow > sys.mrow)
      sys.mrow <- sys.nrow
    # Each system must have two columns. If the total number in a sheet is odd
    #  it triggers an error (check AQSys.err.R for details)
    if (is.odd(sys.ncol)) {AQSys.err("9")}
  }
  #
  for (SheetIndex in grep("BINODAL", sheets)) {
    # Read data from a sheet and workbook provided by the user-invoked function
    sys.temp <- readWorksheet(workBook, SheetIndex, header = FALSE)
    # DISCONSIDER THE FIRST TWO COLUMNS FROM THE STANDARD DATASHEET, WHICH EXIST TO EXEMPLIFY USERS HOW TO FILL THE WORKSHEET
    sys.temp <- sys.temp[, 0:-2] 
    # COUNT NUMBER OF ROWS
    sys.nrow <- nrow(sys.temp)
    # populate all rows are initialized with NA
    if (sys.nrow < sys.mrow) {
      sys.temp[sys.mrow,] <- NA
    }
    # if no system has been added, just add the first sheet
    if (is.null(sys.data)) {
      sys.data <- sys.temp
    } else {
      # but if sys.data have data, convert it to list and concatenate it with data
      # from the current sheet. Then convert it to dataframe and store it.
      sys.data <- bindDATA(list('SET1'=sys.data, 'SET2'=sys.temp))
      #sys.data <- as.data.frame(c(sys.data, sys.temp), stringsAsFactors = FALSE)
    }
  }
  # return all data merged into a single dataframe
  invisible(sys.data)
}

getCAS <- function(workBook, sheets) {
  # Initiate data.frame
  casdb <- data.frame(stringsAsFactors = FALSE)
  # find a sheet with the "CASDB" fragment in its name and load it
  casdb <- readWorksheet(workBook, grep("CASDB", sheets), header = TRUE)
  #define casdb headers
  names(casdb) <- c("CAS.INDEX", "CAS.CODE", "CHEM.NAME", "CHEM.COMMON")
  return(casdb)
}

getREF <- function(workBook, sheets) {
  # initiate data.frame
  refdb <- data.frame(stringsAsFactors = FALSE)
  # find a sheet with the "REFDB" fragment in its name and load it
  refdb <- readWorksheet(workBook, grep("REFDB", sheets), header = TRUE)
  # initiate its second column
  #refdb[, 2] <- NA
  # define refdb headers
  names(refdb) <- c("REF.INDEX", "REF.NAME", "REF.MD5", "REF.URL", "REF.YEAR")
  # encrypt entries found in the file using md5 and store it in refdb
  refdb[, 3] <- sapply(refdb[, 2], digest, algo = "md5")
  # return data.frame containing all read references
  return(refdb)
}

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r + 1, nrow(existingDF) + 1), ] <- existingDF[seq(r, nrow(existingDF)), ]
  existingDF[r, ] <- newrow
  existingDF
}

saveDATA <- function(path, data){
  workBook <- loadWorkbook(path)
  sheets <- getSheets(workBook)
  #
  if (any(grepl("results", sheets))){
    res_idx <- grep("results", sheets)
    clearSheet(workBook, sheet = res_idx) 
  } else {
    createSheet(workBook, name = "results")
  }
  #
  writeWorksheet(workBook, data, sheet = "results", startRow = 1, startCol = 1)
  saveWorkbook(workBook)
  #return()
}

getTL <- function(workBook, sheets) {
  sys.nrow <- 0
  sys.mrow <- NULL
  sys.ncol <- NULL
  sys.data <- NULL
  sys.block <- NULL
  # path must point to a xlsx or xls file
  sys.temp <- data.frame()
  for (SheetIndex in grep("TIELINE", sheets)) {
    sys.temp <- rbind(sys.temp, readWorksheet(workBook, SheetIndex, endCol = 19, header = TRUE))
  }
  #
  names(sys.temp) <- c("REF.MD5", "XY", "PH", "T", "X", "Y", "Y1", "Y2", "X1", 
                       "X2", "G1", "G2", "W", "Y3", "X3", "Z", "Y4", "X4", "TLSlope")
  #
  for (tlr in seq(1, nrow(sys.temp))){
    if (sys.temp[tlr, "Y1"] < sys.temp[tlr, "Y2"]){
      sys.temp[tlr, "TLSlope"] <- ((sys.temp[tlr, "X1"] - sys.temp[tlr, "Y1"]) / (sys.temp[tlr, "X2"] - sys.temp[tlr, "Y2"]))
      sys.temp[tlr, "XY"] <- "XY"
    } else {
      sys.temp[tlr, "TLSlope"] <- ((sys.temp[tlr, "X2"] - sys.temp[tlr, "Y2"])/(sys.temp[tlr, "X1"] - sys.temp[tlr, "Y1"]))
      sys.temp[tlr, "XY"] <- "YX"
    }
  }
  #sys.temp[, 19] <- ((sys.temp[, 9] - sys.temp[, 7]) / (sys.temp[, 10] - sys.temp[, 8]))
  #
  uniqeKeys <- count_(sys.temp, vars = c('PH', 'T', 'X', 'Y'))
  #
  return(list("uniqeKeys" = uniqeKeys, "systems" = sys.temp))
}

XLSCheck <- function(workBook, sheets){
  req_sheets <- c("BINODAL", "TIELINE", "PAR", "REFDB", "CAS")
  cat("Checking XLS file integrity...\n\n")
  for (sheet in req_sheets) {
    cat(sheet, ' worksheet',': ', sep = '')
    if (any(grepl(sheet, sheets))){
      cat('[OK]\n')
    } else {cat('[ERROR] (Worksheet is missing!)\n')}
  }
  cat('\n')
}

bindDATA <- function(XYData) {
  #
  SET1.DATA <- XYData[["SET1"]]
  SET1.NROW <- nrow(SET1.DATA)
  SET2.DATA <- XYData[["SET2"]]
  SET2.NROW <- nrow(SET2.DATA)
  #
  if (SET1.NROW > SET2.NROW){
    SET3.DATA <- as.data.frame(matrix(nrow = (SET1.NROW - SET2.NROW), ncol = ncol(SET2.DATA)), stringsAsFactors = FALSE)
    colnames(SET3.DATA) <- colnames(SET2.DATA)
    OUT.DATA <- rbind(SET2.DATA, SET3.DATA)
    SET2.DATA <- OUT.DATA
  } else { 
    SET3.DATA <- as.data.frame(matrix(nrow = (SET2.NROW - SET1.NROW), ncol = ncol(SET1.DATA)), stringsAsFactors = FALSE)
    colnames(SET3.DATA) <- colnames(SET1.DATA)
    OUT.DATA <- rbind(SET1.DATA, SET3.DATA)
    SET1.DATA <- OUT.DATA
  }
  #
  return(cbind(SET1.DATA, SET2.DATA))
}

toBNDL <- function(workBook, sheets) {
  # initialize variables
  sys.nrow <- 0
  sys.block <- NULL
  # path must point to a xlsx or xls file
  TLData <- getTL(workBook, sheets)
  uniqeKeys <- TLData[['uniqeKeys']]
  sys.temp <- TLData[['systems']]
  sys.nrow <- (max(uniqeKeys$n) - 1) * 2 + 5
  sys.block <- unname(data.frame(matrix(nrow = sys.nrow), stringsAsFactors = FALSE))
  #
  for (sys in 1:nrow(uniqeKeys)) {
    system_data <- sys.temp[which(
          (sys.temp$PH == uniqeKeys[sys,]$PH | is.na(sys.temp$PH)) &
          sys.temp$T == uniqeKeys[sys,]$T &
          sys.temp$X == uniqeKeys[sys,]$X &
          sys.temp$Y == uniqeKeys[sys,]$Y
      ), names(sys.temp)]
    XYdt <- NULL
    #
    for (tieline in 1:nrow(system_data)) {
      if (is.null(XYdt)) {
        XYdt <- data.frame(
          c(
            system_data[tieline, 'PH'], system_data[tieline, 'T'],
            system_data[tieline, 'Y'], system_data[tieline, 'REF.MD5'], NA
          ),
          c(
            system_data[tieline, 'W'], system_data[tieline, 'X3'],
            system_data[tieline, 'X'], system_data[tieline, 'XY'], NA
          ),
          stringsAsFactors = FALSE
        )
        names(XYdt) <- c('COLUMN_1', 'COLUMN_2')
      } else{
        dFF <-data.frame(
            c(system_data[tieline, 9], system_data[tieline, 10]),
            c(system_data[tieline, 7], system_data[tieline, 8]),
            stringsAsFactors = FALSE
          )
        names(dFF) <- c('COLUMN_1', 'COLUMN_2')
        XYdt <- rbind(XYdt, dFF)
      }
    }
    #
    if (nrow(XYdt) < sys.nrow) {
      XYdt <- data.frame(rbind(as.matrix(XYdt), matrix(ncol = 2, nrow = (sys.nrow - nrow(XYdt)))), stringsAsFactors = FALSE)
    }
    #
    sys.block <- cbind(sys.block, XYdt)
  }
  #
  sys.ncol <- ncol(sys.block[, 0:-1])
  if (is.odd(sys.ncol)) {AQSys.err("9")}
  return(sys.block[, 0:-1])
}
# merges data contained in all sheets from the same workbook which are named using
# the pattern datasource_sheetname_YX in its name.
# the order in the end is optional by now but will soon be part of the main function
AQSys.merge <- function(wrbk, sheets) {
  # initialize variables
  sys.nrow <- NULL
  sys.mrow <- NULL
  sys.ncol <- NULL
  sys.data <- NULL
  # make a loop through all sheets that satisfy the condition
  # the integrity of the worksheets will be evaluated before the merge
  # and determine which system in the workbook have bigger dataset (mrow)
  for (nSh in grep("datasource_", sheets)) {
    # determine the number of row and columns in the worksheet
    sys.nrow <- nrow(readWorksheet(wrbk, nSh, header = FALSE))
    sys.ncol <- ncol(readWorksheet(wrbk, nSh, header = FALSE))
    # initialize variables for first run only
    if (is.null(sys.mrow))
      sys.mrow <- sys.nrow
    if (sys.nrow > sys.mrow)
      sys.mrow <- sys.nrow
    # Each system must have two columns. If the total number in a sheet is odd
    #  it triggers an error (check AQSys.err.R for details)
    if (is.odd(sys.ncol))
      AQSys.err("2")
  }
  #
  for (nSh in grep("datasource_", sheets)) {
    #
    sys.temp <- readWorksheet(wrbk, nSh, header = FALSE)
    sys.nrow <- nrow(sys.temp)
    # populate all rows are initialized with NA
    if (sys.nrow < sys.mrow) {
      sys.temp[sys.mrow,] <- NA
    }
    # if no system has been added, just add the first sheet
    if (is.null(sys.data)) {
      sys.data <- sys.temp
      # but if sys.data have data, convert it to list and concatenate it with data
      # from the current sheet. Then convert it to dataframe and store it.
    } else {
      sys.data <-
        as.data.frame(c(sys.data, sys.temp), stringsAsFactors = FALSE)
    }
    
  }
  # return all data merged into a single dataframe
  invisible(sys.data)
}

TLAnalysis <- function(workBook, sheets) {
  # path must point to a xlsx or xls file
  TLData <- getTL(workBook, sheets)
  uniqeKeys <- TLData[['uniqeKeys']]
  sys.temp <- TLData[['systems']]
  #
  sys.slopes <- data.frame(matrix(nrow = nrow(uniqeKeys), ncol = 6))
  names(sys.slopes) <- c("REF.MD5", "PH", "T", "X", "Y", "TLSlope")
  for (sys in 1:nrow(uniqeKeys)) {
    #
    system_data <- sys.temp[which(
      (sys.temp$PH == uniqeKeys[sys,]$PH | is.na(sys.temp$PH)) &
        sys.temp$T == uniqeKeys[sys,]$T &
        sys.temp$X == uniqeKeys[sys,]$X &
        sys.temp$Y == uniqeKeys[sys,]$Y
    ), names(sys.temp)]
    #
    TLSlopes <- system_data[, 19]
    sys.slopes[sys, 1:6] <- system_data[1, -2][1, -6:-17]
    TLSlope <- mean(TLSlopes[which((abs(median(TLSlopes) - TLSlopes) / max(abs(TLSlopes))) < 0.15)])
    if (!is.nan(TLSlope)){sys.slopes[sys, 6] <- TLSlope}
  }
  return(sys.slopes)
}

toNumeric <- function(XYData, ColDis) {
  # convert and name variables accordingly into vectors
  if (tolower(ColDis) == "xy") {
    xc <- as.vector(as.numeric(sub(",", ".", XYData[, 1], fixed = TRUE)))
    yc <- as.vector(as.numeric(sub(",", ".", XYData[, 2], fixed = TRUE)))
  } else{
    xc <- as.vector(as.numeric(sub(",", ".", XYData[, 2], fixed = TRUE)))
    yc <- as.vector(as.numeric(sub(",", ".", XYData[, 1], fixed = TRUE)))
  }
  # and combine them into a dataframe
  XYdt <- data.frame(XC = xc, YC = yc)
  # Remove NA's
  XYdt <- XYdt[complete.cases(XYdt), ]
  #
  if ((nrow(XYdt) > 0) && (max(XYdt)<=1)){
    XYdt <- XYdt * 100
  }
  #return data silently - should it be Visible or hidden?
  invisible(XYdt)
}
