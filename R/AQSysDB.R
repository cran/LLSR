####################################################################################################################
options(digits = 14)
####################################################################################################################
#' @import XLConnect digest
####################################################################################################################
#' @rdname AQSysDB
#' @name AQSysDB
#' @title AQSysDB
#' @description Import DB data from an Excel Worksheet.
#' @export
#' @param path String containing the full path to the XLS or XLSX file.
#' @param ColDis Defines how the data is organized in the Worksheet. Use "xy" whether the first column corresponds to the lower phase fraction and "yx" whether the opposite.
#' @param CAS The user has the option to identify the component's cells in the worksheet with the CAS (CAS = TRUE) or with the row number that matches a CAS entry in the CASDB worksheet (CAS = FALSE)
#' @examples
#' \dontrun{
#' AQSysDB("C:/data.xls", ColDis = "xy", CAS = FALSE)
#'}
####################################################################################################################
# AQSysDB() is a simple approach that is ready to use any three-parameter equation
# and thus
#
AQSysDB <- function(path, ColDis = "xy", CAS = FALSE) {
  # path must point to a xlsx or xls file
  if (grepl(".xlsx", path) | grepl(".xls", path)) {
    # Clean terminal and load the specified file
    cat('\014')
    workBook <- loadWorkbook(path)
    # load sheets
    sheets <- getSheets(workBook)
    # Check file to make sure all required sheets exists
    XLSCheck(workBook, sheets)
    # initiate refdb - dataframe that will get reference data
    refdb <- getREF(workBook, sheets)
    # initiate casdb - dataframe that will get CAS data
    casdb <- getCAS(workBook, sheets)
    #
    tldb <- TLAnalysis(workBook, sheets)
    # AQSys.getBNDL and to AQSys.toBNDL collect data from the specified data and return a datastream ready for processing
    # Check AQSysUtils.R for details.
    bndlDATA <- getBNDL(workBook, sheets)
    tlDATA <- toBNDL(workBook, sheets)
    #return(tlDATA)
    XPData <- bindDATA(list('SET1'=bndlDATA, 'SET2'=tlDATA))
    # Each system have two columns, thus the total number of columns divided by two
    # gives the number of systems
    nSys <- ncol(XPData) / 2
    # set llsrb as a dataframe which data are not converted automatically to factors
    llsr_db <- list()
    # System's info and data location starts in the row below
    db.info <- 1
    db.data <- 6
    # Just giving user an output on R prompt, showing what system is under analysis
    cat(paste("Analysing ", nSys, " systems. \n\n", collapse = NULL))
    # the experimental phase diagram data fetched in the lines above will be used
    # to calculate the nonlinear parameters for all equations in AQSysList()
    for (i in AQSysList()) {
      model_db <- data.frame(stringsAsFactors = FALSE)
      cat(i, ": ", sep = "")
      xcpt <- 0
      for (j in 1:nSys) {
        # COLUMN_1 e COLUMN_2 are the index for the systems unders analysis at the momment
        COLUMN_1 <- 2 * (j) - 1
        COLUMN_2 <- COLUMN_1 + 1
        # get the data length of system under analysis
        lSys <- length(XPData[, COLUMN_1])
        # select phase diagram's data only
        rawSys <- XPData[db.data:lSys, COLUMN_1:COLUMN_2]
        # remove NA entries and convert to dataframe
        db.Sys <- as.data.frame(na.exclude(rawSys), stringsAsFactors = FALSE)
        # Adjust parameters according to data
        regData <- AQSys(LLSRxy(db.Sys[, 1], db.Sys[, 2], XPData[4, COLUMN_2]), mathDesc = i)
        #
        if (!is.null(regData)) {
          resSys <- summary(regData)
          # populate sysDATA with the appropriated parameters from the nonlinear regression
          summary_nrow <- nrow(resSys$parameters)
          summary_ncol <- ncol(resSys$parameters)
          sysDATA <- data.frame(matrix(nrow = 1, ncol = summary_nrow * summary_ncol + 12), stringsAsFactors = FALSE)
          # add md5 encoded ref to model_db
          sysDATA[1, 1] <- refdb[which(refdb$REF.INDEX == to.numeric(XPData[db.info + 3, COLUMN_1])), 3] # to.numeric(XPData[db.info + 3, COLUMN_1]) #to.numeric(XPData[db.info + 3, COLUMN_1]) #
          # if cas field in sysdb is filled with the cas
          C1.CAS.INDEX <- to.numeric(XPData[db.info + 2, COLUMN_1])
          C2.CAS.INDEX <- to.numeric(XPData[db.info + 2, COLUMN_2])
          if (CAS == TRUE) {
            # add Component's CAS directly to sysDATA
            sysDATA[1, 2] <- casdb[which(casdb$CAS.INDEX == C1.CAS.INDEX), 2]
            sysDATA[1, 3] <- casdb[which(casdb$CAS.INDEX == C2.CAS.INDEX), 2]
            # if cas field in sysdb in filled with an index refering to casdb
          } else{
            # Cross reference indexes and add only the Component's NAME
            sysDATA[1, 2] <- casdb[which(casdb$CAS.INDEX == C1.CAS.INDEX), 3]
            sysDATA[1, 3] <- casdb[which(casdb$CAS.INDEX == C2.CAS.INDEX), 3]
          }
          # populate db with system's pH, additive, additive conc and temperature
          sysDATA[1, 4] <- to.numeric(XPData[db.info, COLUMN_1])
          sysDATA[1, 5] <- XPData[db.info, COLUMN_2]
          sysDATA[1, 6] <- to.numeric(XPData[db.info + 1, COLUMN_2])
          sysDATA[1, 7] <- to.numeric(XPData[db.info + 1, COLUMN_1])
          ParamNames <- c("REF.MD5", "UP.Rich", "LP.Rich", "PH", "Additive", "Additive.Conc", "Temperature")
          # The loop below accounts for models with different number of parameters
          idx <- 0
          for (col in 1:summary_ncol) {
            for (row in 1:summary_nrow) {
              idx <- idx + 1
              sysDATA[1, idx + 7] <- resSys$parameters[row, col]
              sysColName <-gsub(' ', '', paste(rownames(resSys$parameters)[row], colnames(resSys$parameters)[col], sep = '-'))
              ParamNames <- c(ParamNames, sysColName)
            }
          }
          # add NLS error-related analysis data to the data.frame
          sysDATA[1, summary_nrow * summary_ncol + 8] <- resSys$sigma
          sysDATA[1, summary_nrow * summary_ncol + 9] <- sum(resSys$residuals ^ 2)
          sysDATA[1, summary_nrow * summary_ncol + 10] <- resSys$convInfo$finTol
          sysDATA[1, summary_nrow * summary_ncol + 11] <- length(db.Sys[, 1])
          sysDATA[1, summary_nrow * summary_ncol + 12] <- i
          # name the above acessed parameters
          ParamNames <- c(ParamNames, "Res.Std.Err", "SSR", "Ach.Conv.Tol", "n.Points", "math.Desc")
          colnames(sysDATA) <- ParamNames
          # Add the results for a given model to the output data.frame
          model_db <- rbind(model_db, sysDATA)
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
    # return silently all data obtained from the worksheet in a list of three dataframes
    invisible(list("db.ref" = refdb, "db.sys" = llsr_db, "db.cas" = casdb, "db.data" = XPData, "db.tielines" = tldb))
  } else {
    # if an invalid path is loaded, it triggers an error
    # (check AQSys.err.R for details)
    AQSys.err("1")
  }
}
