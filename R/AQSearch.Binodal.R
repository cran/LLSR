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
    db.params <- list(db.ph, db.CompA, db.CompB, db.Temp, db.CompC, db.uid)
    # if all parameters are null, the search is not valid and it triggers an
    # error (check AQSys.err.R for details)
    if (all(sapply(db.params, is.null))) AQSys.err("6")
    
    # output variable is initialised with data from db.
    db.grep <- db$db.data
    
    # IMPROVEMENT: Exit early if there's no data to search. Use && for short-circuiting.
    if (is.null(db.grep) || (nrow(db.grep) == 0 && ncol(db.grep) == 0)) {
        AQSys.err("5") # No data in db$db.data, so no results possible
        return(invisible(NULL))
    }

    # IMPROVEMENT: Create a local, pre-processed CAS database for searching.
    # This is the single most important performance enhancement.
    # The tolower() operation is now performed only ONCE, not for every component search.
    cas_search_db <- db$db.cas
    cas_search_db$CAS.NAME_lower <- tolower(cas_search_db$CAS.NAME)
    cas_search_db$CAS.COMMON_lower <- tolower(cas_search_db$CAS.COMMON)
    cas_search_db$CAS.CODE_lower <- tolower(cas_search_db$CAS.CODE)

    # --------------------------------------------------------------------------
    # IMPROVEMENT: Abstract the repetitive search logic into a helper function.
    # This makes the main function body much cleaner and easier to maintain (DRY principle).
    #
    # @param component_query The search term (e.g., db.CompA).
    # @param local_cas_db The pre-processed CAS database.
    # @return A character vector of matching canonical CAS.NAMEs.
    #
    search_component <- function(component_query, local_cas_db) {
      if (!is.valid(component_query)) return(NULL)
      
      # Convert search query to lowercase once
      query_lower <- tolower(component_query)
      
      # Find rows where the query matches any of the relevant columns
      match_idx <- grepl(query_lower, local_cas_db$CAS.NAME_lower, fixed = TRUE) |
                   grepl(query_lower, local_cas_db$CAS.COMMON_lower, fixed = TRUE) |
                   grepl(query_lower, local_cas_db$CAS.CODE_lower, fixed = TRUE)
      
      # Return the canonical names for all matches
      unique(local_cas_db$CAS.NAME[match_idx])
    }
    # --------------------------------------------------------------------------
    
    # IMPROVEMENT: Use the helper function for cleaner, faster, and non-repetitive searches.
    # The logic is now streamlined and much more readable.
    for (component in list(db.CompA, db.CompB, db.CompC)) {
      if (is.valid(component) && ncol(db.grep) > 0) {
        chem_names <- search_component(component, cas_search_db)
        if (length(chem_names) > 0) {
          db.grep <- db.grep[, matchBNDL(chem_names, db.grep), drop = FALSE]
        } else {
          db.grep <- db.grep[, FALSE, drop = FALSE] # No names found, so no systems will match
        }
      }
    }

    # search a system that match the system's temperature, if search parameter
    #  is not null.
    if (is.valid(db.Temp) && ncol(db.grep) > 0) {
      db.grep <- matchTpH(db.Temp, db.grep, FALSE)
    }
    
    # search a system that matchs the system's pH, if search parameter is not
    #  null.
    if (is.valid(db.ph) && ncol(db.grep) > 0) {
      db.grep <- matchTpH(db.ph, db.grep, TRUE)
    }
    
    # search a system that matchs the system's UID, if search parameter is not
    # null.
    if (is.valid(db.uid) && ncol(db.grep) > 0) {
      db.grep <- db.grep[, matchBNDL(db.uid, db.grep), drop = FALSE]
    }

    # Add a check for ncol to prevent errors with crayon
    if (ncol(db.grep) > 0) {
      cat(paste("    Your search had [", red(ncol(db.grep) / 2), "] results.",
                "\n", sep = ""))
      
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