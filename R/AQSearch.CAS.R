####################################################################################################################
#' @rdname AQSearch.CAS
#' @name AQSearch.CAS
#' @title Search CAS codes in Chemical database.
#' @method AQSearch CAS
#' @export AQSearch.CAS
#' @export
#' @param db A list of three data.frame variables within all systems data are stored.
#' @param ChemString A String with a text fragment of the chemical name the user want to search the correspondent CAS.
#' @param ... Additional optional arguments. None are used at present.
#' @examples 
#' \dontrun{
#' AQSearch.CAS(db, "NaOH")
#'}
####################################################################################################################
AQSearch.CAS <- function(db, ChemString, ...){
  #
  #names(db) <- c("CAS.CODE", "CHEM.NAME", "CHEM.COMMON")
  db.check(db)
  ans <- data.frame()
  #
  if (is.null(ChemString) == FALSE){
    db.grep <- grep(toupper(ChemString), toupper(db$db.cas[,2]))
    #
    if ((length(db.grep) != 0) & (db.grep[1] != 0)){
      db <- db$db.cas[db.grep, ]
    }else{
      db.grep <- grep(toupper(ChemString), toupper(db$db.cas[,3]))
    }
    #
    if ((length(db.grep) != 0) & (db.grep[1] != 0)){
      ans <- db$db.cas[db.grep, ]
    } else{
      AQSys.err("8")
    }
  } else {
    AQSys.err("4")
  }
  invisible(ans)
}
