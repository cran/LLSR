###############################################################################
options(digits = 14)
###############################################################################
#' @import openxlsx digest
###############################################################################
#' @rdname AQSysDB
#' @name AQSysDB
#' @title AQSysDB
#' @description Import DB data from an Excel Worksheet and process it through 
#' mathematical descriptors to output a highly structured variable comparable 
#' to a Database and which hold a list of references, chemicals and parameters
#'  for any implemented mathematical descriptors.
#' @export
#' @param path String containing the full path to the XLS or XLSX file.
# ' @param maxiter	- A positive integer specifying the maximum number of 
# ' iterations allowed.
# ' @param Order Defines how the data is organized in the Worksheet. Use "xy" 
# ' whether the first column corresponds to the lower phase fraction and "yx" 
# ' whether the opposite.
#' @examples
#' \dontrun{
#' AQSysDB("C:/data.xlsx")
#'}
###############################################################################
# AQSysDB() is a simple approach that is ready to use any three-parameter 
# equation and thus
#
# AQSysDB <- function(path, Order = "xy", CAS = FALSE) {
AQSysDB <- function(path) {
  # path must point to a xlsx or xls file
  if (grepl(".xlsx", path) | grepl(".xls", path)) {
    # Clean terminal and load the specified file
    cat('\014')
    workBook <- loadWorkbook(path)
    # load sheets
    sheets <- getSheetNames(path)
    # Check file to make sure all required sheets exists
    XLSCheck(workBook, sheets)
    # initiate refdb - dataframe that will get reference data
    refdb <- getREF(workBook, sheets)
    # initiate casdb - dataframe that will get CAS data
    casdb <- getCAS(workBook, sheets)
    #
    tldb <- TLAnalysis(workBook, sheets)
    tldb.raw <- getTL(workBook, sheets)[["systems"]]
    # AQSys.getBNDL and to AQSys.toBNDL collect data from the specified data 
    # and return a datastream ready for processing Check AQSysUtils.R for 
    # details.
    bndlDATA <- getBNDL(workBook, sheets)
    #
    # EXTRAPOLATE TIELINES INTO BINODAL CURVES
    # tlDATA <- toBNDL(workBook, sheets)
    # BINDS EXTRAPOLATED TIE-LINES AND BINODAL DATA
    # DISABLED TEMPORARILY
    # XPData <- nameIT(bindDATA(list('SET1' = bndlDATA, 'SET2' = tlDATA)))
    #
    #
    XPData <- nameIT(bndlDATA)
    # Each system have two columns, thus the total number of columns divided by 
    # two gives the number of systems
    nSys <- ncol(XPData) / 2
    # set llsrb as a dataframe which data are not converted automatically to 
    # factors
    llsr_db <- list()
    # System's info and data location starts in the row below
    db.info <- 1
    db.data <- 6
    # Just giving user an output on R prompt, showing what system is under
    # analysis
    cat(paste("Analysing ", nSys, " systems. \n\n", collapse = NULL))
    # the experimental phase diagram data fetched in the lines above will be 
    # used to calculate the nonlinear parameters for all equations in AQSysList
    for (i in AQSysList()) {
      model_db <- data.frame(stringsAsFactors = FALSE)
      cat(i, ": ", sep = "")
      xcpt <- 0
      for (j in 1:nSys) {
        # COLUMN_1 e COLUMN_2 are the index for the systems unders analysis at 
        # the momment
        COLUMN_1 <- 2 * (j) - 1
        COLUMN_2 <- COLUMN_1 + 1
        # get the data length of system under analysis
        lSys <- length(XPData[, COLUMN_1])
        # select phase diagram's data only
        rawSys <- XPData[db.data:lSys, COLUMN_1:COLUMN_2]
        # naMAtrix<-as.data.frame(matrix(ncol=2,nrow=(lSys-db.data+1)))
        # XPData[db.data:lSys, COLUMN_1:COLUMN_2] <- naMAtrix
        # remove NA entries and convert to dataframe
        db.Sys <- as.data.frame(na.exclude(rawSys), stringsAsFactors = FALSE)
        numData <- LLSRxy(db.Sys[, 1], db.Sys[, 2], XPData[4, COLUMN_2])
        # XPData[db.data:nrow(numData), COLUMN_1:COLUMN_2] <- numData
        # Adjust parameters according to data
        regData <- AQSys(numData, modelName = i)
        #
        if (!is.null(regData)) {
          resSys <- summary(regData)
          # populate sysDATA with the appropriated parameters from the nonlinear
          # regression
          summary_nrow <- nrow(resSys$parameters)
          summary_ncol <- ncol(resSys$parameters)
          sysDATA <- data.frame(matrix(
            nrow = 1, ncol = summary_nrow * summary_ncol + 13), 
            stringsAsFactors = FALSE)
          # add md5 encoded ref to model_db
          sysDATA[1, db.info] <- refdb[which(
            refdb$REF.INDEX == to.numeric(XPData[db.info + 3, COLUMN_1])), 3] 
          # to.numeric(XPData[db.info+3,COLUMN_1])
          # if cas field in sysdb is filled with the cas
          CA.CAS.INDEX <- to.numeric(XPData[db.info + 2, COLUMN_1])
          CB.CAS.INDEX <- to.numeric(XPData[db.info + 2, COLUMN_2])
          #
          #
          if (XPData[4, COLUMN_2] == "XY") {
            TP.IDX <- (db.info + 2)
            BT.IDX <- (db.info + 1)
          } else if (XPData[4, COLUMN_2] == "YX") {
            TP.IDX <- (db.info + 1)
            BT.IDX <- (db.info + 2)
          }
          #
          sysDATA[1, TP.IDX] <- casdb[which(casdb$CAS.INDEX == CA.CAS.INDEX), 3]
          sysDATA[1, BT.IDX] <- casdb[which(casdb$CAS.INDEX == CB.CAS.INDEX), 3]
          # }
          # populate db with system's pH, additive, additive conc and temp
          sysDATA[1, (db.info + 3)] <- XPData[4, COLUMN_2]
          sysDATA[1, (db.info + 4)] <- to.numeric(XPData[db.info, COLUMN_1])
          sysDATA[1, (db.info + 5)] <- to.numeric(XPData[db.info + 1, COLUMN_1])
          sysDATA[1, (db.info + 6)] <- XPData[db.info, COLUMN_2]
          sysDATA[1, (db.info + 7)] <- to.numeric(XPData[db.info + 1, COLUMN_2])
          ParamNames <- c("REF.MD5", "A", "B", "ORDER","PH", "TEMP", "C", 
                          "GLB.C")
          # The loop below accounts for models with different number of params
          idx <- 0
          for (col in 1:summary_ncol) {
            for (row in 1:summary_nrow) {
              idx <- idx + 1
              sysDATA[1, idx + (db.info + 7)] <- resSys$parameters[row, col]
              sysColName <- gsub(' ', '', paste(rownames(resSys$parameters)[row],
                                               colnames(resSys$parameters)[col],
                                               sep = '-'))
              ParamNames <- c(ParamNames, sysColName)
            }
          }
          # add NLS error-related analysis data to the data.frame
          sysDATA[1, summary_nrow * summary_ncol + (db.info + 8)] <-
            resSys$sigma
          sysDATA[1, summary_nrow * summary_ncol + (db.info + 9)] <-
            sum(resSys$residuals ^ 2)
          sysDATA[1, summary_nrow * summary_ncol + (db.info + 10)] <-
            resSys$convInfo$finTol
          sysDATA[1, summary_nrow * summary_ncol + (db.info + 11)] <-
            length(db.Sys[, 1])
          sysDATA[1, summary_nrow * summary_ncol + (db.info + 12)] <-
            i
          # name the above acessed parameters
          ParamNames <- c(ParamNames, "Res.Std.Err", "SSR", "Ach.Conv.Tol", 
                          "n.Points", "math.Desc")
          colnames(sysDATA) <- ParamNames
          # Add the results for a given model to the output data.frame
          model_db <- rbind(model_db, UIDGen(sysDATA))
        } else{
          # account for the exception
          xcpt <- xcpt + 1
        }
      }
      llsr_db[[i]] <- model_db
      # At the end of the analysis for each equation,  return an OK
      if (xcpt > 0) { cat("[", xcpt, "exceptions ] ") }
      cat("[OK] \n")
      #
    }
    # return silently all data obtained from the worksheet in a list of three 
    # data-frames
    #
    invisible(
      to.ascii(
        list(
          "db.ref" = refdb,
          "db.sys" = llsr_db,
          "db.cas" = casdb,
          "tldbl" = tldb,
          "tldb.raw" = tldb.raw,
          "db.data" = SysIdxToRef(refdb, casdb, XPData),
          "db.tielines" = list(
            "data" = IdxToRef(refdb, tldb.raw, casdb),
            "slopes" = IdxToRef(refdb, tldb, casdb)
            )
          )
        )
      )
  } else {
    # if an invalid path is loaded, it triggers an error
    # (check AQSys.err.R for details)
    AQSys.err("1")
  }
}
