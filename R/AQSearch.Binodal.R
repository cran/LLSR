###############################################################################
#' @importFrom crayon red
###############################################################################
#' @rdname AQSearch.Binodal
#' @title Search function for ATPS Systems data
#' @description This function allow the user to search the package database to 
#' find any ATPS that matches the available criteria.
#' @details The function return the systems that matches the criteria submitted 
#' by the user.
#' @param db A highly structure db containing data from previously analyzed
#'  data. LLSR database is used by default but user may input his own db if 
#'  formatted properly.
#' @param db.CompA A String variable containing either the CAS, chemical
#'  formula or name of the upper phase enriched component..
#' @param db.CompB A String variable containing either the CAS, chemical
#'  formula or name of the lower phase component.
#' @param db.CompC A String variable containing either the CAS, chemical
#'  formula or name of the additive component.
#' @param db.Temp A numeric variable containing the Temperature (in Kelvin) to
#'  be searched within DB.
#' @param db.ph A numeric variable containing the pH to be searched within DB.
#' @param db.uid An Unique md5 hash Identification. User can retrieve data
#'  for a specific system if in possesion of its UID.
#' @param stacked A boolean variable used to return value as a nested list
#'  or a data.frame. Used internally to organize data output. 
#' @param ... Additional optional arguments. None are used at present.
#' @method AQSearch Binodal
#' @export AQSearch.Binodal
#' @export
#' @return Returns a data.frame containing system's parameters which match
#'  searched conditions
#' @examples
#' \dontrun{
#' AQSearch.Binodal(db.CompA="Ammonium")
#'}
###############################################################################
AQSearch.Binodal <-
  function(db = LLSR::llsr_data,
           db.CompA = NULL,
           db.CompB = NULL,
           db.CompC = NULL,
           db.Temp = NULL,
           db.ph = NULL,
           db.uid = NULL,
           stacked = FALSE,
           ...) {
    cat("  [Binodal]\n")
    # initialize db.ans
    db.ans <- list()
    # create and initialise a list using the function's parameters
    db.params <- c(db.ph, db.CompA, db.CompB, db.Temp, db.CompC, db.uid)
    # if all parameters are null, the search is not valid and it triggers an
    # error (check AQSys.err.R for details)
    if (all(unlist(lapply(db.params, is.null)))) AQSys.err("6")
    # output variable is initialised with data from db.
    db.grep <- db$db.data
    #
    if ((nrow(db.grep)==0)&(ncol(db.grep)==0)) 
      {db.ph = db.CompA = db.CompB = db.Temp = db.CompC = db.uid = NULL}
    #
    if (is.valid(db.CompA)) {
      db.CompA.names <- db$db.cas[grep(tolower(db.CompA), 
                                       tolower(db$db.cas$CAS.NAME), 
                                       fixed = TRUE), "CAS.NAME"]
      db.CompA.altNames <- db$db.cas[grep(tolower(db.CompA), 
                                          tolower(db$db.cas$CAS.COMMON),
                                          fixed = TRUE), "CAS.NAME"]
      db.CompA.cas <- db$db.cas[grep(tolower(db.CompA), 
                                     tolower(db$db.cas$CAS.CODE), 
                                     fixed = TRUE), "CAS.NAME"]
      #
      db.chem.names <- list(db.CompA.names, db.CompA.altNames, db.CompA.cas)
      db.chem.names <- ifelse(length(db.chem.names) == 0, "", db.chem.names)
      db.grep <- db.grep[, matchBNDL(db.chem.names, db.grep)]
    }
    # search a system that matchs the lower-phase component, if search parameter
    #  is not null.
    if (is.valid(db.CompB)) {
      db.CompB.names <- db$db.cas[grep(tolower(db.CompB), 
                                       tolower(db$db.cas$CAS.NAME), 
                                       fixed = TRUE), "CAS.NAME"]
      db.CompB.altNames <- db$db.cas[grep(tolower(db.CompB), 
                                          tolower(db$db.cas$CAS.COMMON), 
                                          fixed = TRUE), "CAS.NAME"]
      db.CompB.cas <- db$db.cas[grep(tolower(db.CompB), 
                                     tolower(db$db.cas$CAS.CODE), 
                                     fixed = TRUE), "CAS.NAME"]
      #
      db.chem.names <- list(db.CompB.names, db.CompB.altNames, db.CompB.cas)
      db.chem.names <- ifelse(length(db.chem.names) == 0, "", db.chem.names)
      db.grep <- db.grep[, matchBNDL(db.chem.names, db.grep)]
    }
    # search a system that matchs the additive component, if search parameter 
    # is not null.
    if (is.valid(db.CompC)) {
      db.CompC.names <- db$db.cas[grep(tolower(db.CompC), 
                                       tolower(db$db.cas$CAS.NAME),
                                       fixed = TRUE), "CAS.NAME"]
      db.CompC.altNames <- db$db.cas[grep(tolower(db.CompC), 
                                          tolower(db$db.cas$CAS.COMMON),
                                          fixed = TRUE), "CAS.NAME"]
      db.CompC.cas <- db$db.cas[grep(tolower(db.CompC), 
                                     tolower(db$db.cas$CAS.CODE),
                                     fixed = TRUE), "CAS.NAME"]
      #
      db.chem.names <- list(db.CompC.names, db.CompC.altNames, db.CompC.cas)
      db.chem.names <- ifelse(length(db.chem.names) == 0, "", db.chem.names)
      db.grep <- db.grep[, matchBNDL(db.chem.names, db.grep)]
    }
    # search a system that match the system's temperature, if search parameter
    #  is not null.
    if (is.valid(db.Temp)) {
      db.grep <- matchTpH(db.Temp, db.grep, FALSE)
    }
    # search a system that matchs the system's pH, if search parameter is not
    #  null.
    if (is.valid(db.ph)) {
      db.grep <- matchTpH(db.ph, db.grep, TRUE)
    }
    # search a system that matchs the system's UID, if search parameter is not 
    # null.
    if (is.valid(db.uid)) {
      db.grep <- db.grep[, matchBNDL(db.uid, db.grep)]
    }
    if (ncol(db.grep) != 0) {
      cat(paste("    Your search had [", red(ncol(db.grep) / 2), "] results.",
                "\n",sep = ""))
      #
      if (stacked) {
        db.ans[["Binodal"]] <- db.grep
      } else {
        db.ans <- db.grep
      }
      invisible(RevREF(db.ans, db))
    } else {
      # Triggers an "no results" error
      AQSys.err("5")
      invisible(NULL)
    }
  }
