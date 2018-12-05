####################################################################################################################
#' @rdname AQSearch
#' @name AQSearch
#' @description This function allow the user to search the LLSR database to find any ATPS that matches the used criteria.
#' @export
####################################################################################################################
AQSearch <- function(db = LLSR::llsr_data, db.ph = NULL, db.upper = NULL, db.lower = NULL,
                     db.temp = NULL, db.additive = NULL, ...)
  UseMethod("AQSearch")
####################################################################################################################
#' @rdname AQSearch
#' @title Search function for ATPS Systems data
#' @description This function allow the user to search the package database to find any ATPS that matches the used criteria.
#' @details The function return the systems that matches the criteria submit by the user. 
#' @param db A highly structure db containing data from previously analised data. LLSR database is used by default but user may input his own db if formatted properly.
#' @param db.ph A numeric variable containing the pH to be searched within DB.
#' @param db.upper A String variable containing either the CAS, chemical formula or name of the upper phase enriched component..
#' @param db.lower A String variable containing either the CAS, chemical formula or name of the lower phase component.
#' @param db.temp A numeric variable containing the Temperature (in Kelvin) to be searched within DB.
#' @param db.additive A String variable containing either the CAS, chemical formula or name of the additive component.
#' @param ... Additional optional arguments. None are used at present.
#' @method AQSearch default
#' @export
#' @return Returns a data.frame containing system's parameters which match searched conditions
#' @examples
#' \dontrun{
#' AQSearch(db.upper="85100-78-3")
#'}
####################################################################################################################
AQSearch.default <- function(db = LLSR::llsr_data, db.ph = NULL, db.upper = NULL, db.lower = NULL,
           db.temp = NULL, db.additive = NULL, ...) {
    #
    #db <- LLSR::llsr_data
    #inialize db.ans
    db.ans <- list()
    # create and initialise a list using the function's parameters
    db.params <- c(db.ph, db.upper, db.lower, db.temp, db.additive)
    # if all parameters are null, the search is not valid and it triggers an error (check AQSys.err.R for details)
    if (all(unlist(lapply(db.params, is.null))))
      AQSys.err("6")
    # each conditional will return a result set that excludes no match data.
    # search a system that matchs the upper-phase component, if search parameter is not null.
    for (modelName in AQSysList()){
      # output variable is initialised with data from db.
      db.grep <- db$db.sys[[modelName]]
      if (!is.null(db.upper)) {
          db.upper.altNames <- db$db.cas[grep(db.upper, db$db.cas$CHEM.COMMON, ignore.case = TRUE), "CHEM.NAME"]
          db.upper.cas <- db$db.cas[grep(db.upper, db$db.cas$CAS.CODE, ignore.case = TRUE), "CHEM.NAME"]
          db.grep <- db.grep[unique(c(grep(db.upper, db.grep$UP.Rich, ignore.case = TRUE), 
                                           which(tolower(db.grep$UP.Rich) %in% db.upper.altNames),
                                           which(tolower(db.grep$UP.Rich) %in% db.upper.cas))),]
      }
      # search a system that matchs the lower-phase component, if search parameter is not null.
      if (!is.null(db.lower)) {
        db.lower.altNames <- db$db.cas[grep(db.lower, db$db.cas$CHEM.COMMON, ignore.case = TRUE), "CHEM.NAME"]
        db.lower.cas <- db$db.cas[grep(db.lower, db$db.cas$CAS.CODE, ignore.case = TRUE), "CHEM.NAME"]
        db.grep <- db.grep[unique(c(grep(db.lower, db.grep$LP.Rich, ignore.case = TRUE), 
                                         which(tolower(db.grep$LP.Rich) %in% db.lower.altNames),
                                         which(tolower(db.grep$LP.Rich) %in% db.lower.cas))),]
      }
      # search a system that matchs the system's pH, if search parameter is not null.
      if (!is.null(db.ph)) {
        db.grep <- db.grep[grep(db.ph, db.grep[,4]),]
      }
      # search a system that matchs the additive component, if search parameter is not null.
      if (!is.null(db.additive)) {
        db.additive.altNames <- db$db.cas[grep(db.additive, db$db.cas$CHEM.COMMON, ignore.case = TRUE), "CHEM.NAME"]
        db.additive.cas <- db$db.cas[grep(db.additive, db$db.cas$CAS.CODE, ignore.case = TRUE), "CHEM.NAME"]
        db.grep <- db.grep[unique(c(grep(db.additive, db.grep$Additive, ignore.case = TRUE), 
                                         which(tolower(db.grep$Additive) %in% db.additive.altNames),
                                         which(tolower(db.grep$Additive) %in% db.additive.cas))),]
      }
      # search a system that matchs the system's temperature, if search parameter is not null.
      if (!is.null(db.temp)) {
        db.grep <- db.grep[grep(db.temp, db.grep[,7]),]
      }
      if (length(db.grep) != 0) {
        db.grep[, 1] <- as.data.frame(sapply(db.grep[, 1], undigest))
        names(db.grep)[1] <- "REFERENCE"
        db.ans[[modelName]] <- db.grep
        print(paste("[", modelName, "] ", "Your search had [", nrow(db.grep), "] results.", sep = ""))
      } else {
        print(paste("[", modelName, "] ", "Your search had no results.", sep = ""))
      }
    }
    # If search isn't null, evaluate data
    if (length(db.ans) != 0)  {
     invisible(db.ans)
      #
    }else{
      # Triggers an "no results" error
      AQSys.err("5")
    }

  }
