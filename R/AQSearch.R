####################################################################################################################
#' @rdname AQSearch
#' @name AQSearch
#' @description This function allow the user to search the LLSR database to find any ATPS that matches the used criteria.
#' @export
####################################################################################################################
AQSearch <- function(db = LLSR::llsr_data, ...)
  UseMethod("AQSearch")
####################################################################################################################
#' @rdname AQSearch
#' @title Search function for ATPS Systems data
#' @description This function allow the user to search the package database to find any ATPS that matches the available criteria.
#' @details The function return the systems that matches the criteria submitted by the user. 
#' @param db A highly structure db containing data from previously analised data. LLSR database is used by default but user may input his own db if formatted properly.
#' @param db.CompA A String variable containing either the CAS, chemical formula or name of the upper phase enriched component..
#' @param db.CompB A String variable containing either the CAS, chemical formula or name of the lower phase component.
#' @param db.CompC A String variable containing either the CAS, chemical formula or name of the additive component.
#' @param db.Temp A numeric variable containing the Temperature (in Kelvin) to be searched within DB.
#' @param db.ph A numeric variable containing the pH to be searched within DB.
#' @param db.uid An Unique md5 hash Identification. User can retrieve data for a specific system if in possesion of its UID.
#' @param stacked A boolean variable used to return value as a nested list or a data.frame. Used internally to organize data output. 
#' @param ... Additional optional arguments. None are used at present.
#' @method AQSearch default
#' @export
#' @return Returns a data.frame containing system's parameters which match searched conditions
#' @examples
#' \dontrun{
#' AQSearch(db.CompA="Ammonium")
#'}
####################################################################################################################
AQSearch.default <-
  function(db = LLSR::llsr_data,
           db.CompA = NULL,
           db.CompB = NULL,
           db.CompC = NULL,
           db.Temp = NULL,
           db.ph = NULL,
           db.uid = NULL,
           stacked = FALSE,
           ...) {
    
    # initialize db.ans]
    cat("[Searching DB]\n")
    db.ans <- list()
    db.uids <- NULL
    # create and initialise a list using the function's parameters
    db.params <- c(db.ph, db.CompA, db.CompB, db.Temp, db.CompC, db.uid)
    # if all parameters are null, the search is not valid and it triggers an error (check AQSys.err.R for details)
    if (all(unlist(lapply(db.params, is.null)))) AQSys.err("6")
    #
    parameters <- AQSearch.Parameter(db = db, db.CompA, db.CompB, db.CompC, db.Temp, db.ph, db.uid, stacked = TRUE)
    tielines <- AQSearch.Tieline(db = db, db.CompA, db.CompB, db.CompC, db.Temp, db.ph, db.uid, stacked = TRUE)
    binodals <- AQSearch.Binodal(db = db, db.CompA, db.CompB, db.CompC, db.Temp, db.ph, db.uid, stacked = TRUE)
    slope <- AQSearch.Slope(db = db, db.CompA, db.CompB, db.CompC, db.Temp, db.ph, db.uid, stacked = TRUE)
    #
    db.ans <-c(parameters, tielines, binodals, slope)
    # If search isn't null, evaluate data
    if (length(db.ans) != 0)  {
     invisible(db.ans)
      #
    }else{
      # Triggers an "no results" error
      AQSys.err("5")
    }

  }
