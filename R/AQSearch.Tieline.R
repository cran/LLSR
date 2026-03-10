###############################################################################
# ' @importFrom crayon red bold
###############################################################################
#' @rdname AQSearch.Tieline
#' @title Search function for ATPS Systems data
#' @description This function allow the user to search the package database to 
#' find any ATPS that matches the available criteria.
#' @details The function return the systems that matches the criteria submitted 
#' by the user.
#' @param db A highly structure db containing data from previously analyzed 
#' data. LLSR database is used by default but user may input his own db if
#'  formatted properly.
#' @param db.CompA A String variable containing either the CAS, chemical 
#' formula or name of the upper phase enriched component..
#' @param db.CompB A String variable containing either the CAS, chemical 
#' formula or name of the lower phase component.
#' @param db.CompC A String variable containing either the CAS, chemical 
#' formula or name of the additive component.
#' @param db.Temp A numeric variable containing the Temperature (in Kelvin) 
#' to be searched within DB.
#' @param db.ph A numeric variable containing the pH to be searched within DB.
#' @param db.uid An Unique md5 hash Identification. User can retrieve data 
#' for a specific system if in possesion of its UID.
#' @param stacked A boolean variable used to return value as a nested list 
#' or a data.frame. Used internally to organize data output.
#' @param ... Additional optional arguments. None are used at present.
#' @method AQSearch Tieline
#' @export AQSearch.Tieline
#' @export 
#' @return Returns a data.frame containing system's parameters which match 
#' searched conditions
#' @examples
#' \dontrun{
#' AQSearch.Tieline(db.CompA="Ammonium")
#'}
###############################################################################
AQSearch.Tieline <-
  function(db = LLSR::llsr_data,
           db.CompA = NULL,
           db.CompB = NULL,
           db.CompC = NULL,
           db.Temp = NULL,
           db.ph = NULL,
           db.uid = NULL,
           stacked = FALSE,
           ...) {
    cat("  [Tieline]\n")
    # initialize db.ans
    db.ans <- list()
    # create and initialize a list using the function's parameters
    db.params <- c(db.ph, db.CompA, db.CompB, db.Temp, db.CompC, db.uid)
    # if all parameters are null, the search is not valid and it triggers an
    #  error (check AQSys.err.R for details)
    if (all(unlist(lapply(db.params, is.null)))) AQSys.err("6")
    # output variable is initialized with data from db.
    db.grep <- db$db.tielines$data
    cas_search_db <- build_cas_search_db(db$db.cas)
    #
    if ((nrow(db.grep) == 0) & (ncol(db.grep) == 0)) 
      {db.ph = db.CompA = db.CompB = db.Temp = db.CompC = db.uid = NULL}
    #
    if (is.valid(db.CompA)) {
      db.CompA.names <- search_component_names(db.CompA, cas_search_db)
      db.chem.names <- ifelse(length(db.CompA.names) == 0, "", db.CompA.names)
      db.grep <- db.grep[matchComp(db.chem.names, db.grep), ]
    }
    # search a system that matchs the lower-phase component, if search parameter
    #  is not null.
    if (is.valid(db.CompB)) {
      db.CompB.names <- search_component_names(db.CompB, cas_search_db)
      db.chem.names <- ifelse(length(db.CompB.names) == 0, "", db.CompB.names)
      db.grep <- db.grep[matchComp(db.chem.names, db.grep), ]
    }
    # search a system that matchs the additive component, if search parameter 
    # is not null.
    if (is.valid(db.CompC)) {
      db.CompC.names <- search_component_names(db.CompC, cas_search_db)
      db.chem.names <- ifelse(length(db.CompC.names) == 0, "", db.CompC.names)
      db.grep <- db.grep[matchComp(db.chem.names, db.grep), ]
    }
    # search a system that matchs the system's temperature, if search parameter
    #  is not null.
    if (is.valid(db.Temp)) {
      db.grep <- db.grep[(!is.na(db.grep$TEMP) & db.grep$TEMP == db.Temp), ]
    }
    # search a system that matchs the system's pH, if search parameter is not
    #  null.
    if (is.valid(db.ph)) {
      db.grep <- db.grep[(!is.na(db.grep$PH) & db.grep$PH == db.ph), ]
    }
    # search a system that matchs the system's UID, if search parameter is not
    #  null.
    if (is.valid(db.uid)) {
      db.grep <- db.grep[matchUID(db.uid, db.grep), ]
    }
    #
    if (ncol(db.grep) != 0) {
      cat(paste("    Your search had [", red(nrow(db.grep)), "] results.",
                "\n",sep = ""))
      if (stacked) {
        db.ans[["TieLines"]] <- db.grep
      } else {
        db.ans <- db.grep
      }
      #
      invisible(RevREF(db.ans, db))
    } else {
      # Triggers an "no results" error
      AQSys.err("5")
      invisible(NULL)
    }
  }
