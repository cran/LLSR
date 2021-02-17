# Each number used as argument triggers an error with its respective message
AQSys.err <- function (err, ...) {
  switch(
    err,
    "0" = {
      errmsg <- "The selected Equation doesn't exist."
      stop(errmsg, call. = FALSE)
    },
    "1" = {
      errmsg <- paste0("The path must point to a '.xls' or '.xlsx' Microsoft",
                       " Excel worksheet.")
      stop(errmsg, call. = FALSE)
    },
    "2" = {
      errmsg <- paste0("There was a problem when calculating the number of ",
                       "systems to be analysed. Please make sure the data ",
                       "is formatted as the standard provided.")
      stop(errmsg, call. = FALSE)
    },
    "3" = {
      errmsg <- " must be a data.frame."
      stop(..., errmsg, call. = FALSE)
    },
    "4" = {
      errmsg <- "Argument 'db' is missing. Parameter can not be NULL."
      stop(errmsg, call. = FALSE)
    },
    "5" = {
      errmsg <- paste0("       Search had no results. Try removing a few",
      " parameters and search again.\n")
      cat(errmsg)
      #stop(errmsg, call. = FALSE)
    },
    "6" = {
      errmsg <- paste0("      At least one of the parameters (UID, pH, Temp,",
                       " CompA, CompB or CompC) must be not NULL.\n")
      cat(errmsg)
      #stop(errmsg, call. = FALSE)
    },
    "7" = {
      errmsg <- paste0("Input variable db must be a list coontaining three",
                       " data.frame variables (db.cas, db.ref and db.sys).")
      stop(errmsg, call. = FALSE)
    },
    "8" = {
      errmsg <- "Your search had no results.\n"
      cat(errmsg)
      #stop(errmsg, call. = FALSE)
    },
    "9" = {
      errmsg <- paste0("A valid data set has multiple two-columns data for ",
      "calculations but an odd number of columns was observed in the input ",
      "file.")
      stop(errmsg, call. = FALSE)
    },
    "10" = {
      errmsg <- paste0("Target-TLL must be BIGGER than the minimum TLL and",
                       " SMALLER than the maximum TLL for the system chosen.")
      stop(errmsg, call. = FALSE)
    },
    "11" = {
      errmsg <- paste0("Slope's array must have the same length of systems",
                       " provided inthe dataset.")
      stop(errmsg, call. = FALSE)
    },
    "12" = {
      errmsg <- paste0("Slope's database does not have entries for the systems",
                       " provided in the dataset. Please input slope data.")
      stop(errmsg, call. = FALSE)
    },
    "13" = {
      errmsg <- paste0("The weight of each phase (Wt and Wb) must be provided",
                       " when byW is TRUE.")
      stop(errmsg, call. = FALSE)
    },
    "14" = {
      errmsg <- paste0("Each Phase's volume and aparent density must be",
                       " provided when byW is FALSE.")
      stop(errmsg, call. = FALSE)
    },
    {
      errmsg <- "An Unknown error ocourred."
      stop(errmsg, call. = FALSE)
    }
  )
}
