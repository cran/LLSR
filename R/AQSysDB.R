####################################################################################################################
options(digits=14)
require(XLConnect)
require(digest)
####################################################################################################################
#' @import XLConnect
#' @import digest
####################################################################################################################
#' @rdname AQSysDB
#' @name AQSysDB
#' @title AQSysDB
#' @export
#' @param path String containing the full path to the XLS or XLSX file.
#' @param order Defines how the data is organized in the Worksheet. Use "xy" whether the first column corresponds to the lower phase fraction and "yx" whether the opposite.
#' @param CAS The user has the option to identify the component's cells in the worksheet with the CAS (CAS = TRUE) or with the row number that matches a CAS entry in the CASDB worksheet (CAS = FALSE)
#' @examples 
#' \dontrun{
#' AQSysDB("C:/data.xls", order = "xy", CAS = FALSE)
#'}
####################################################################################################################
AQSysDB <- function(path, order = "xy", CAS = FALSE){
  if (grepl(".xlsx",path) | grepl(".xls", path)){
    #
    wrbk <- loadWorkbook(path)
    sheets <- getSheets(wrbk)
    #
    refdb <- data.frame()
    refdb <- readWorksheet(wrbk, grep("REFDB", sheets), header = FALSE)
    refdb[,2] <- NA
    names(refdb) <- c("REF.NAME","REF.MD5")
    refdb[,2] <- sapply(refdb[,1], digest, algo="md5")
    #
    casdb <- data.frame()
    casdb <- readWorksheet(wrbk, grep("CASDB", sheets), header = FALSE)
    names(casdb) <- c("CAS.CODE", "CHEM.NAME", "CHEM.COMMON")
    #
    #CONSIDER FIND OTHERS SHEETS FOLLOWING A PATTERN AND EVALUATE ALL OF THEM
    #
    #Loop through worksheet index that have the datasource string in its name
    #and aggregate data under one data prior going got analysis.
    #
    #wsdt <- readWorksheet(wrbk, 3, header = FALSE)
    #
    #if (is.odd(ncol(wsdt))) AQSys.err("2")
    wsdt <- AQSys.merge(wrbk, sheets)
    #
    nSys <- ncol(wsdt)/2
    #
    llsrdb <- data.frame(stringsAsFactors=FALSE)
    db.info <- 1
    db.data <- 6 
    #
  } else {
    AQSys.err("1")
  }
  #
  cat('\014')
  cat(paste("Analysing ",nSys," systems. \n\n", collapse = NULL))
  #
  for (i in AQSys.List()){
    #
    cat(i,": ", sep = "")
    #
    db.n.row <- nrow(llsrdb)
    if(db.n.row == 0){
      db.first.row <- 1
      db.last.row <- nSys
    }else{
      db.first.row <- db.n.row + 1
      db.last.row <- db.n.row + nSys
    }
    #
    for (j in db.first.row:db.last.row){
      #
      c1 <- 2*(j - db.n.row) - 1
      c2 <- c1 + 1
      lSys<-length(wsdt[,c1])
      #
      llsrdb[j,1] <- refdb[wsdt[db.info + 3, c1],2]
      #
      if (CAS == TRUE) {
        llsrdb[j,2] <- wsdt[db.info + 2, c1]
        llsrdb[j,3] <- wsdt[db.info + 2, c2]
      } else{
        llsrdb[j,2] <- casdb[wsdt[db.info + 2, c1], 1]
        llsrdb[j,3] <- casdb[wsdt[db.info + 2, c2], 1]
      }
      #
      llsrdb[j,4] <- wsdt[db.info, c1]
      llsrdb[j,5] <- wsdt[db.info, c2]
      llsrdb[j,6] <- wsdt[db.info + 1, c2]
      llsrdb[j,7] <- wsdt[db.info + 1, c1]
      #
      rawSys <- wsdt[db.data:lSys, c1:c2]
      #
      db.Sys <- as.data.frame(na.exclude(rawSys))
      #
      if (is.numeric(db.Sys[1,1])){
      } else{
        db.first.col <- as.numeric(sub(",", ".", db.Sys[,1], fixed = TRUE))
        db.second.col <- as.numeric(sub(",", ".", db.Sys[,2], fixed = TRUE))
      }
      #
      if (tolower(order)=="xy"){
        resSys<-summary(AQSys(LLSRxy(db.first.col,db.second.col),mathDesc=i))
      }
      else{
        resSys<-summary(AQSys(LLSRxy(db.second.col,db.first.col),mathDesc=i))
      }
      #
      llsrdb[j,8] <- resSys$parameters[1,1]
      llsrdb[j,9] <- resSys$parameters[2,1]
      llsrdb[j,10] <- resSys$parameters[3,1]
      llsrdb[j,11] <- resSys$sigma
      llsrdb[j,12] <- sum(resSys$residuals^2)
      llsrdb[j,13] <- resSys$parameters[1,2]
      llsrdb[j,14] <- resSys$parameters[2,2]
      llsrdb[j,15] <- resSys$parameters[3,2]
      llsrdb[j,16] <- resSys$parameters[1,3]
      llsrdb[j,17] <- resSys$parameters[2,3]
      llsrdb[j,18] <- resSys$parameters[3,3]
      llsrdb[j,19] <- resSys$convInfo$finTol
      llsrdb[j,20] <- length(db.first.col)
      llsrdb[j,21] <- i
    }
    #
    cat("[OK] \n")
    #
  }
  #
  names(llsrdb) <- c("REF.MD5", "UP.Rich", "LP.Rich", "Sys.pH", "Sys.Additive",
                     "Sys.Additive.Conc", "Sys.Temp","P1", "P2", "P3",  "Res.Std.Err",
                     "SSR",	"P1.Std.Err",	"P2.Std.Err",	"P3.Std.Err",
                    "P1.tValue",	"P2.tValue",	"P3.tValue", "Ach.Conv.Tol",
                    "n.Points", "math.Desc")
  invisible(list("db.ref" = refdb, "db.sys" = llsrdb, "db.cas" = casdb))
}



