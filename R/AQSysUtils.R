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
    contact us on https://github.com/diegofcoelho/LLSR/issues.\n
    You also can fork this project directly from github and commit
    improvements to us (http://diegofcoelho.github.io/LLSR/).\n
    The information used in the database was obtained free of charge
    but it might be copyrighted by third parties and references must
    be included appropriately.\n',
    domain = NULL,
    appendLF = TRUE
  )
}
#
matchUID <- function(UIDList, UIDMatrix) {
  uniqeUIDs <- as.character(unique(UIDList))
  #
  req_results <- sapply(uniqeUIDs, function(uid) {
    unlist(which(UIDMatrix == uid, TRUE))
  })
  #
  if (is.list(req_results)) {
    db.binodals.index <- do.call(rbind, req_results)[, "row"]
  } else {
    db.binodals.index <- req_results[1:(length(req_results) / 2)]
  }
  #
  return(db.binodals.index)
}
#
# matchBNDL <- function(matchingList, matchingMatrix) {
#   matchingUIDs <- as.character(unique(matchingList))
#   #
#   db.binodals.index <-
#     as.numeric(setNames(na.exclude(sapply(matchingUIDs, function(uid) {
#       which(matchingMatrix == uid, TRUE)[2]
#     })), NULL))
#   #
#   db.binodals.index <- sort(c(db.binodals.index, db.binodals.index + 1))
#   #
#   return(db.binodals.index)
# }
#
matchTpH<- function(TpH, BinodalMatrix, pH) {
  #
  if (pH){RowIdx <- 1} else {RowIdx <- 2}
  #
  db.TpH.results <- which(BinodalMatrix == TpH, TRUE)
  db.TpH.check <- which(db.TpH.results[, "row"] == RowIdx)
  #
  db.binodals.index <- db.TpH.results[db.TpH.check, "col"]
  #
  db.binodals.index <-
    sapply(db.binodals.index, function(idx) {
      if (is.odd(idx)) {
        idx
      } else {
        idx - 1
      }
    })
  #
  db.binodals.index <- sort(c(db.binodals.index, db.binodals.index + 1))
  #
  return(db.binodals.index)
}
#
matchBNDL<- function(compNameList, BinodalMatrix) {
  #
  if (length(compNameList)==0){return(c())}
  #
  req_results <- sapply(compNameList, function(compName) {
    unlist(which(BinodalMatrix == compName, TRUE))
  })
  #
  if (is.list(req_results)){
    db.binodals.index <- as.numeric(do.call(rbind, req_results)[, 2])
  } else {
    db.binodals.index <- req_results[(1 + length(req_results) / 2):length(req_results)]
  }
  #
  if (length(db.binodals.index)==0){return(c())}
  #
  db.binodals.index <-
    sapply(db.binodals.index, function(idx) {
      if (is.odd(idx)) {
        idx
      } else {
        idx - 1
      }
    })
  #
  db.binodals.index <- sort(c(db.binodals.index, db.binodals.index + 1))
  #
  return(db.binodals.index)
}
#
idx2name <- function(idx, casdb) {
  casdb[casdb$CAS.INDEX == idx,]$CAS.NAME
}
#
commaReplacer <- function(str) {
  as.numeric(gsub(",", ".", str))
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
    sys.nrow <-
      nrow(readWorksheet(workBook, SheetIndex, header = FALSE))
    sys.ncol <-
      ncol(readWorksheet(workBook, SheetIndex, header = FALSE))
    # initialize variables for first run only
    if (is.null(sys.mrow))
      sys.mrow <- sys.nrow
    if (sys.nrow > sys.mrow)
      sys.mrow <- sys.nrow
    # Each system must have two columns. If the total number in a sheet is odd
    #  it triggers an error (check AQSys.err.R for details)
    if (is.odd(sys.ncol)) {
      AQSys.err("9")
    }
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
      sys.temp[sys.mrow, ] <- NA
    }
    # if no system has been added, just add the first sheet
    if (is.null(sys.data)) {
      sys.data <- sys.temp
    } else {
      # but if sys.data have data, convert it to list and concatenate it with data
      # from the current sheet. Then convert it to dataframe and store it.
      sys.data <-
        bindDATA(list("SET1" = sys.data, "SET2" = sys.temp))
      #sys.data <- as.data.frame(c(sys.data, sys.temp), stringsAsFactors = FALSE)
    }
  }
  # return all data merged into a single dataframe
  invisible(sys.data)
}
#
getCAS <- function(workBook, sheets) {
  # Initiate data.frame
  casdb <- data.frame(stringsAsFactors = FALSE)
  # find a sheet with the "CASDB" fragment in its name and load it
  casdb <-
    readWorksheet(workBook, grep("CASDB", sheets), header = TRUE)
  #define casdb headers
  names(casdb) <-
    c("CAS.INDEX", "CAS.CODE", "CAS.NAME", "CAS.COMMON")
  return(casdb)
}
#
getREF <- function(workBook, sheets) {
  # initiate data.frame
  refdb <- data.frame(stringsAsFactors = FALSE)
  # find a sheet with the "REFDB" fragment in its name and load it
  refdb <-
    readWorksheet(workBook, grep("REFDB", sheets), header = TRUE)
  # initiate its second column
  #refdb[, 2] <- NA
  # define refdb headers
  names(refdb) <-
    c("REF.INDEX", "REF.NAME", "REF.MD5", "REF.URL", "REF.YEAR")
  # encrypt entries found in the file using md5 and store it in refdb
  refdb[, 3] <- sapply(refdb[, 2], digest, algo = "md5")
  # return data.frame containing all read references
  return(refdb)
}
#
SysIdxToRef <- function(refSheet, db.cas, db.data) {
  sysNum <- ncol(db.data) / 2
  for (i in seq(1, sysNum)) {
    db.data[3, i * 2] <- db.cas[db.data[3, i * 2], "CAS.NAME"]
    db.data[3, i * 2 - 1] <- db.cas[db.data[3, i * 2 - 1], "CAS.NAME"]
    db.data[4, i * 2 - 1] <- refSheet[db.data[4, i * 2 - 1], 3]
  }
  db.data[6:nrow(db.data),] <- sapply(db.data[6:nrow(db.data),], commaReplacer)
  return(UIDGen(db.data))
}
#
IdxToRef <- function(db.ref, db.data, db.cas) {
  # encrypt entries found in the file using md5 and store it in refdb
  db.data[, "REF.MD5"] <- db.ref[db.data[, "REF.MD5"], 3]
  db.data[, "A"] <- sapply(db.data[, "A"], idx2name, casdb = db.cas)
  db.data[, "B"] <- sapply(db.data[, "B"], idx2name, casdb = db.cas)
  return(UIDGen(db.data))
}
#
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r + 1, nrow(existingDF) + 1),] <-
    existingDF[seq(r, nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}
#
saveDATA <- function(path, data) {
  workBook <- loadWorkbook(path)
  sheets <- getSheets(workBook)
  #
  if (any(grepl("results", sheets))) {
    res_idx <- grep("results", sheets)
    clearSheet(workBook, sheet = res_idx)
  } else {
    createSheet(workBook, name = "results")
  }
  #
  writeWorksheet(
    workBook,
    data,
    sheet = "results",
    startRow = 1,
    startCol = 1
  )
  saveWorkbook(workBook)
  #return()
}
#
getTL <- function(workBook, sheets) {
  sys.nrow <- 0
  sys.mrow <- NULL
  sys.ncol <- NULL
  sys.data <- NULL
  sys.block <- NULL
  # path must point to a xlsx or xls file
  sys.temp <- data.frame()
  for (SheetIndex in grep("TIELINE", sheets)) {
    sys.temp <-
      rbind(sys.temp,
            readWorksheet(
              workBook,
              SheetIndex,
              endCol = 17,
              header = TRUE
            ))
  }
  #
  names(sys.temp) <-
    c(
      "REF.MD5",
      "A",
      "B",
      "ORDER",
      "PH",
      "TEMP",
      "TOP.A",
      "TOP.B",
      "BOT.A",
      "BOT.B",
      "GLB.A",
      "GLB.B",
      "C",
      "TOP.C",
      "BOT.C",
      "GLB.C",
      "TLSlope"
    )
  # Account for user mistakes while typing data to the worksheet by detecting which components are the top and bottom ones
  for (tlr in seq(1, nrow(sys.temp))) {
    if (sys.temp[tlr, "TOP.A"] > sys.temp[tlr, "TOP.B"]) {
      sys.temp[tlr, "TLSlope"] <-
        ((sys.temp[tlr, "BOT.A"] - sys.temp[tlr, "TOP.A"]) / (sys.temp[tlr, "BOT.B"] - sys.temp[tlr, "TOP.B"]))
      sys.temp[tlr, "ORDER"] <- "YX"
    } else {
      sys.temp[tlr, "TLSlope"] <-
        ((sys.temp[tlr, "BOT.B"] - sys.temp[tlr, "TOP.B"]) / (sys.temp[tlr, "BOT.A"] - sys.temp[tlr, "TOP.A"]))
      sys.temp[tlr, "ORDER"] <- "XY"
    }
  }
  # sys.temp[, 19] <- ((sys.temp[, 9] - sys.temp[, 7]) / (sys.temp[, 10] - sys.temp[, 8]))
  #
  uniqeKeys <- count_(sys.temp, vars = c("PH", "TEMP", "A", "B"))
  #
  return(list("uniqeKeys" = uniqeKeys, "systems" = sys.temp))
}
#
XLSCheck <- function(workBook, sheets) {
  #req_sheets <- c("BINODAL", "TIELINE", "PAR", "REFDB", "CAS")
  req_sheets <- c("BINODAL", "TIELINE", "REFDB", "CAS")
  cat("Checking XLS file integrity...\n\n")
  for (sheet in req_sheets) {
    cat(sheet, " worksheet", ": ", sep = "")
    if (any(grepl(sheet, sheets))) {
      cat("[OK]\n")
    } else {
      cat("[ERROR] (Worksheet is missing!)\n")
    }
  }
  cat("\n")
}
#
nameIT <- function(XYData) {
  #  NAMING SYSTEMS AND COMPONENT COLUMNS ACCORDINGLY
  DATA.LENGTH <- ncol(XYData)
  DATA.NSYS <- DATA.LENGTH / 2
  DATA.NAMES <-
    paste(
      rep("S", DATA.LENGTH),
      rep(1:DATA.NSYS, each = 2),
      "_",
      "C",
      rep(LETTERS[1:2], each = 1),
      sep = ""
    )
  names(XYData) <- DATA.NAMES
  return(XYData)
}
#
bindDATA <- function(XYData) {
  # ID-ING DATASETS
  SET1.DATA <- XYData[["SET1"]]
  SET1.NROW <- nrow(SET1.DATA)
  SET2.DATA <- XYData[["SET2"]]
  SET2.NROW <- nrow(SET2.DATA)
  # CHECK WHICH DATASET IS BIGGER TO EVEN OUT THE SMALLEST
  if (SET1.NROW > SET2.NROW) {
    SET3.DATA <-
      as.data.frame(matrix(
        nrow = (SET1.NROW - SET2.NROW),
        ncol = ncol(SET2.DATA)
      ), stringsAsFactors = FALSE)
    colnames(SET3.DATA) <- colnames(SET2.DATA)
    OUT.DATA <- rbind(SET2.DATA, SET3.DATA)
    SET2.DATA <- OUT.DATA
  } else {
    SET3.DATA <-
      as.data.frame(matrix(
        nrow = (SET2.NROW - SET1.NROW),
        ncol = ncol(SET1.DATA)
      ), stringsAsFactors = FALSE)
    colnames(SET3.DATA) <- colnames(SET1.DATA)
    OUT.DATA <- rbind(SET1.DATA, SET3.DATA)
    SET1.DATA <- OUT.DATA
  }
  DATA.OUTPUT <- cbind(SET1.DATA, SET2.DATA)
  # RETURN BINDED DATASET
  return(DATA.OUTPUT)
}
#
toBNDL <- function(workBook, sheets) {
  # initialize variables
  sys.nrow <- 0
  sys.block <- NULL
  # path must point to a xlsx or xls file
  TLData <- getTL(workBook, sheets)
  uniqeKeys <- TLData[["uniqeKeys"]]
  sys.temp <- TLData[["systems"]]
  names(sys.temp)
  sys.nrow <- (max(uniqeKeys$n) - 1) * 2 + 5
  sys.block <-
    unname(data.frame(matrix(nrow = sys.nrow), stringsAsFactors = FALSE))
  #
  for (sys in 1:nrow(uniqeKeys)) {
    system_data <- sys.temp[which(
      (sys.temp$PH == uniqeKeys[sys, ]$PH | is.na(sys.temp$PH)) &
        sys.temp$TEMP == uniqeKeys[sys, ]$TEMP &
        sys.temp$A == uniqeKeys[sys, ]$A &
        sys.temp$B == uniqeKeys[sys, ]$B
    ), names(sys.temp)]
    XYdt <- NULL
    #
    for (tieline in 1:nrow(system_data)) {
      if (is.null(XYdt)) {
        #
        XYdt <- data.frame(
          c(
            system_data[tieline, "PH"],
            system_data[tieline, "TEMP"],
            system_data[tieline, "A"],
            system_data[tieline, "REF.MD5"],
            NA
          ),
          c(
            system_data[tieline, "C"],
            system_data[tieline, "GLB.C"],
            system_data[tieline, "B"],
            system_data[tieline, "ORDER"],
            NA
          ),
          stringsAsFactors = FALSE
        )
        names(XYdt) <- c("COLUMN_1", "COLUMN_2")
      } else{
        dFF <- data.frame(
          as.numeric(c(system_data[tieline, "TOP.A"], system_data[tieline, "BOT.A"])),
          as.numeric(c(system_data[tieline, "TOP.B"], system_data[tieline, "BOT.B"])),
          stringsAsFactors = FALSE
        )
        names(dFF) <- c("COLUMN_1", "COLUMN_2")
        XYdt <- rbind(XYdt, dFF)
      }
    }
    # ORDER DATA - POINTS ARE SCATTERED DUE TIE-LINE DATA BINDING
    XYdt.temp <- XYdt[6:nrow(XYdt),]
    XYdt.temp <- data.frame(sapply(XYdt[6:nrow(XYdt),], function(x) as.numeric(x)))
    #
    XYdt.temp <- XYdt.temp[with(XYdt.temp, order(-COLUMN_1)), ]
    XYdt[6:nrow(XYdt),] <- XYdt.temp
    #
    if (nrow(XYdt) < sys.nrow) {
      XYdt <-
        data.frame(rbind(as.matrix(XYdt), matrix(
          ncol = 2, nrow = (sys.nrow - nrow(XYdt))
        )), stringsAsFactors = FALSE)
    }
    #
    sys.block <- cbind(sys.block, XYdt)
  }
  #
  sys.ncol <- ncol(sys.block[, 0:-1])
  if (is.odd(sys.ncol)) {
    AQSys.err("9")
  }
  return(sys.block[, 0:-1])
}
#
AQSys.merge <- function(wrbk, sheets) {
  # merges data contained in all sheets from the same workbook which are named using
  # the pattern datasource_sheetname_YX in its name.
  # the order in the end is optional by now but will soon be part of the main function
  #
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
      sys.temp[sys.mrow, ] <- NA
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
#
TLAnalysis <- function(workBook, sheets) {
  # path must point to a xlsx or xls file
  TLData <- getTL(workBook, sheets)
  uniqeKeys <- TLData[["uniqeKeys"]]
  sys.temp <- TLData[["systems"]]
  #
  sys.slopes <- data.frame(matrix(nrow = nrow(uniqeKeys), ncol = 7))
  # names(sys.slopes) <- c("REF.MD5", "PH", "TEMP", "A", "B", "TLSlope")
  names(sys.slopes) <-
    c("REF.MD5", "A", "B", "ORDER", "PH", "TEMP", "TLSlope")
  for (sys in 1:nrow(uniqeKeys)) {
    #
    system_data <- sys.temp[which(
      (sys.temp$PH == uniqeKeys[sys, ]$PH | is.na(sys.temp$PH)) &
        sys.temp$TEMP == uniqeKeys[sys, ]$TEMP &
        sys.temp$A == uniqeKeys[sys, ]$A &
        sys.temp$B == uniqeKeys[sys, ]$B
    ), names(sys.temp)]
    #
    TLSlopes <- system_data[, ncol(system_data)]
    sys.slopes[sys, 1:7] <-
      system_data[1, -7:-(ncol(system_data) - 1)]
    # sys.slopes[sys, 1:6] <- system_data[1, -2][1, -6:-17]
    # sys.slopes[sys, 1:7] <- system_data[1 , c(1, 5, 6, 2, 3, 4, 19)]
    #
    TLSlope <-
      mean(TLSlopes[which((abs(median(TLSlopes) - TLSlopes) / max(abs(TLSlopes))) < 0.15)])
    #
    if (!is.nan(TLSlope)) {
      sys.slopes[sys, 7] <- TLSlope
    }
  }
  return(sys.slopes)
}
#
toNumeric <- function(XYData, Order) {
  # convert and name variables accordingly into vectors
  if (tolower(Order) == "xy") {
    xc <-
      as.vector(as.numeric(sub(",", ".", XYData[, 1], fixed = TRUE)))
    yc <-
      as.vector(as.numeric(sub(",", ".", XYData[, 2], fixed = TRUE)))
  } else{
    xc <-
      as.vector(as.numeric(sub(",", ".", XYData[, 2], fixed = TRUE)))
    yc <-
      as.vector(as.numeric(sub(",", ".", XYData[, 1], fixed = TRUE)))
  }
  # and combine them into a dataframe
  XYdt <- data.frame(XC = xc, YC = yc)
  # Remove NA's
  XYdt <- XYdt[complete.cases(XYdt),]
  #
  if ((nrow(XYdt) > 0) && (max(XYdt) <= 1)) {
    XYdt <- XYdt * 100
  }
  #return data silently - should it be Visible or hidden?
  invisible(XYdt)
}
#
UIDGen <- function(db.data) {
  OUTPUT.DATA <- NULL
  BY.ROW <- ("REF.MD5" %in% colnames(db.data))
  #
  if (BY.ROW) {
    UID <- sapply(paste(db.data[, "REF.MD5"],
                        db.data[, "A"],
                        db.data[, "B"],
                        db.data[, "PH"],
                        db.data[, "TEMP"],
                        sep = "_"),
                  digest,
                  algo = "md5")
    OUTPUT.DATA <- cbind(UID, db.data)
    rownames(OUTPUT.DATA) <- NULL
  } else {
    DATA.LENGTH <- ncol(db.data)
    DATA.NSYS <- DATA.LENGTH / 2
    for (IDX in seq(1, DATA.NSYS)) {
      db.data[5, 2 * IDX - 1] <-
        digest(paste(db.data[4, 2 * IDX - 1],
                     db.data[3, 2 * IDX - 1],
                     db.data[3, 2 * IDX],
                     # db.data[1, 2 * IDX - 1],
                     db.data[2, 2 * IDX - 1],
                     sep = "_"),
               algo = "md5")
    }
    OUTPUT.DATA <- db.data
  }
  return(OUTPUT.DATA)
}
#
ChckBndrs <- function(BNFn, slope, BNDL, xmax){
  #
  yMin <- BNFn(xmax)
  Gn <- function (yMin, slope, xmax, x) {
    yMin + slope * (x - xmax)
  }
  #
  x_min_b <- min(uniroot.all(function(x) (BNFn(x) - Gn(yMin, slope, xmax, x)), c(0, xmax), tol = 0.001)) 
  y_max_b <- max(BNDL["Y"])
  ymax <- BNFn(x_min_b)
  #
  while (ymax > y_max_b) {
    xmax <- xmax - 0.01
    yMin <- BNFn(xmax)
    x_min_b <- min(uniroot.all(function(x) (BNFn(x) - Gn(yMin, slope, xmax, x)), c(0, xmax), tol = 0.001)) 
    ymax <- BNFn(x_min_b)
  }
  return(xmax)
  
}
