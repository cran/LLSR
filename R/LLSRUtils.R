# quiets concerns of R CMD check re: the .'s that appear in pipelines
# if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
#
#' @importFrom dplyr case_when
is.valid <- function(y) {
  caseFun <- function(x) {
    isTRUE(
      case_when(
        (length(x) == 0) && (typeof(x) == "character") ~ FALSE,
        is.null(x) ~ FALSE,
        is.na(x) ~ FALSE,
        is.nan(x) ~ FALSE,
        (trimws(x) == "") ~ FALSE,
        TRUE ~ TRUE
      )
    )
  }
  all(c(unlist(lapply(y, caseFun)),!is.null(y)))
}
#
is.odd <- function(x)
  x %% 2 != 0
#
is.even <- function(x)
  x %% 2 == 0
#
to.numeric <- function(x) {
  if (is.na(x) ||
      is.null(x) || x == "NA") {
    NA
  } else{
    as.double(sub(",", ".", x))
  }
}
#
is.equal <- function(TLData, tol) {
  # Verify if all Xs and Ys are equal among themselves
  return((all((
    max(TLData[, 1]) - min(TLData[, 1])
  ) < tol) == TRUE) &&
    (all((
      max(TLData[, 2]) - min(TLData[, 2])
    ) < tol) == TRUE))
}
#
FindMinTL <- function(SysCP, maxGP, xMax, slope, BLFn, tol, dfr = 0.05){
  # Calculate distance between the critical point (SysCP) and the furthest 
  # viable TL's Global Point (maxGP)
  DMaxTL <- sqrt((maxGP[1] - SysCP[1]) ^ 2 + (maxGP[2] - SysCP[2]) ^ 2)
  # Establishes that the distance between the critical point and the closest 
  # viable tieline is a fraction (dfr) of the distance to MaxTL
  d <- DMaxTL * dfr
  # Calculate the coordinates for a point existing in the minTL 
  # (closest viable TL)
  X <- as.numeric(SysCP[1] + sqrt(d ^ 2 / (1 + (slope) ^ 2)))
  Y <- as.numeric(SysCP[2] + (-1 / slope) * (X - SysCP[1]))
  # Initiate a generic line function but using calculated X, Y and provided 
  # Slope (assuming all tielines are parallel)
  TLFn <- function(x) { Y + slope * (x - X) }
  # find the intersections between the minTL and the binodal
  xRoots <- uniroot.all(function(x) (BLFn(x) - TLFn(x)), c(0, xMax * 1.5), 
                        tol = 1e-3) # REPLACE XMAX TO THE LIMIT OF SOLUBILITY?
  # Creates an array containing all Xs which characterizes minTL
  xTL <- c(min(xRoots), sum(xRoots) / 2, max(xRoots))
  # return a data.frame containing all Xs and Ys for minTL
  return(setNames(data.frame(xTL, TLFn(xTL)), c("X", "Y")))
}
#
dPoints <- function(P1, P2){
  d <- sqrt(((P1$X - P2$X) ^ 2) + ((P1$Y - P2$Y) ^ 2))
  return(d)
}
#
TLL <- function(minTL, maxTL){
  TLL <- list()
  #
  P1Min <- minTL[1, 1:2]
  P2Min <- minTL[3, 1:2]
  #
  TLL$MinTLL <- dPoints(P1Min, P2Min)
  #
  P1Max <- maxTL[1, 1:2]
  P2Max <- maxTL[3, 1:2]
  #
  TLL$MaxTLL <- dPoints(P1Max, P2Max)
  #
  return(TLL)
}
#
findTL <- function(dTLL, SysTLL, BLFn, slope){
  # If the target TLL is smaller than the minimum calculated TLL, 
  # throw an error.
  if ((dTLL < SysTLL$TLL$MinTLL) | (dTLL > SysTLL$TLL$MaxTLL)) {
    #AQSys.err("10")
    print(c(dTLL, SysTLL$TLL$MinTLL, SysTLL$TLL$MaxTLL))
  }
  # Initial guess for the tieline, calculated using the tieline length 
  # proportions
  X <- (dTLL / SysTLL$TLL$MaxTLL)*max(SysTLL$maxTL$X)
  # Initializing variables
  TL <- setNames(as.data.frame(matrix(nrow = 3, ncol = 2)), c("X", "Y"))
  OUTPUT <- list()
  dt <- 1
  #
  whl_idx <- 1
  while ((dt > 1e-4) & (whl_idx < 30000)) {
    Y <- BLFn(X)
    TLFn <- function(x) { Y + slope * (x - X) }
    xRoots <- uniroot.all(function(x) (BLFn(x) - TLFn(x)), c(0, X), tol = 1e-3)
    if (identical(xRoots, numeric(0))) {
      X <- X - dt * (-(TLL - dTLL) / abs(TLL - dTLL)) / 10
      break
    }
    #
    TL[1, 1] <- min(xRoots)
    TL[3, 1] <- max(xRoots)
    TL[1, 2] <- unname(BLFn(TL[1, 1]))
    TL[3, 2] <- unname(BLFn(TL[3, 1]))
    #
    TLL <- dPoints(TL[1, ], TL[3, ])
    dt <- (abs(TLL - dTLL))^2
    #
    X <- X + dt * (-(TLL - dTLL) / abs(TLL - dTLL)) / 100
    whl_idx <- whl_idx + 1
  }
  TL[2, 1] <- (TL[1, 1] + TL[3, 1]) / 2
  TL[2, 2] <- (TL[1, 2] + TL[3, 2]) / 2
  #
  OUTPUT$TL <- TL
  OUTPUT$TLL <- TLL
  #
  return(OUTPUT)
}
#
findSlope <- function(db, dataSET){
  # iterate through multiple columns and return a list of slopes
  slope <- c()
  # check how many systems were provided
  nSys <- (ncol(dataSET) / 2)
  if ((ncol(dataSET) %% 2) == 0) {
    for (sys in seq(1, nSys)) {
      # Get the system characterization variables
      idx_Y <- idx2name(dataSET[3, sys * 2 - 1], db$db.cas)
      idx_X <- idx2name(dataSET[3, sys * 2], db$db.cas)
      idx_PH <- dataSET[1, sys * 2 - 1]
      idx_T <- round(as.numeric(dataSET[2, sys * 2 - 1]),2)
      #
      TL_db <- db[["db.tielines"]][["slopes"]]
      slope[sys] <- TL_db[which(
        (TL_db$PH == idx_PH | is.na(TL_db$PH)) &
          TL_db$TEMP == idx_T &
          (TL_db$A == idx_X | TL_db$B == idx_X) &
          (TL_db$A == idx_Y | TL_db$B == idx_Y) ), "TLSlope"]
      #
    }
    if (length(slope) == 0) {
      AQSys.err("12")
    } else {
      return(as.numeric(slope))
    }
  } else{
    # Return an error if an invalid dataset is provided.
    AQSys.err("9")
  }
}
#
seqTL <- function(minTL, maxTL, slope, BLFn, nTL = 3, nSYS = 3) {
  dataNames <- c("X", "Y", "System", "Point")
  #
  xMin <- max(minTL["X"])
  xMax <- max(maxTL["X"])
  # Bottom Phase Compositions
  xRange <- seq(xMin, xMax, (xMax - xMin) / (nTL - 1))
  oDATA <- setNames(data.frame(xRange, BLFn(xRange), seq(1, nTL), "B", 
                               row.names = NULL), dataNames) 
  #
  for (p in seq(1, nrow(oDATA))) {
    X <- oDATA[p, "X"]
    Y <- oDATA[p, "Y"]
    #
    TLFn <- function(x) { Y + slope * (x - X) }
    xRoots <- uniroot.all(function(x) (BLFn(x) - TLFn(x)), c(0, X*2), 
                          tol = 0.1) # REPLACE XMAX TO THE LIMIT OF SOLUBILITY?
    xTL <- c(min(xRoots), sum(xRoots) / 2)
    #
    temp.TLC <- setNames(data.frame(xTL, TLFn(xTL), rep(oDATA[p, "System"], 2),
                                    c("T", "M")), dataNames)
    #
    xSYS <- seq(min(xRoots), X, (X - min(xRoots)) / (nSYS + 1))
    #
    temp.SYS <- setNames(data.frame(xSYS, 
                                    TLFn(xSYS), rep(oDATA[p, "System"], 
                                                    nSYS + 2), 
                                    rep("S", nSYS + 2)), dataNames)
    #
    oDATA <- rbind(oDATA, temp.TLC, temp.SYS)
    #
  }
  return(oDATA)
}
#
saveConfig <- function(plot_obj, save, HR, filename, wdir, silent) {
  if (save == TRUE) {
    #
    if (HR == TRUE) {
      image_format <- ".svg"
    } else{
      image_format <- ".png"
    }
    #
    if (is.null(filename)) {
      # Get user choice for a filename to save the plot
      filename <- dlgInput(message = "Enter the figure filename:")$res
    }
    # complete filename with the appropriated extension
    filename <- paste(filename, image_format, sep = "")
    # Check if filename is invalid and quite if so
    if (filename == image_format) {
      stop("Filename is NULL or INVALID.", call. = TRUE)
    }
    #
    #
    if (is.null(wdir)) {
      # Get user choice for a directory to save the plot
      wdir <- dlgDir()$res
    }
    # Check if path is invalid and quite if so
    if ((wdir == "") && (silent == FALSE)) {
      #
      stop("Path is NULL or INVALID.", call. = TRUE)
      #
    } else if ((wdir == "") && (silent == TRUE)) {
      #
      wdir <- getwd()
      dir_and_file <- paste(wdir, filename, sep = .Platform$file.sep)
      #
    } else{
      #
      dir_and_file <- paste(wdir, filename, sep = .Platform$file.sep)
      #
    }
    #
    ggsave(
      filename = dir_and_file,
      plot = plot_obj,
      width = 21.14 / 2,
      height = 14.39 / 2
    )
  }
  return(wdir)
}
#
GenPlotSeries <- function(SysData, xMAX, NP, modelFn, i, seriesNames){
  min_x <- min(SysData[, 1])
  x <- seq( min_x, xMAX, ((xMAX - (min_x / 1.5)) / NP))
  BNDL <- setNames(data.frame(x, modelFn(x)), c("X", "Y"))
  BNDL["System"] <- seriesNames[i]
  return(BNDL)
}
#
LLSRxy <- function(FirstCol, SecondCol, Order = 'xy') {
  # convert and name variables accordingly into vectors
  if (tolower(Order) == "xy") {
    xc <- as.vector(as.numeric(sub(",", ".", FirstCol, fixed = TRUE)))
    yc <- as.vector(as.numeric(sub(",", ".", SecondCol, fixed = TRUE)))
  } else{
    xc <- as.vector(as.numeric(sub(",", ".", SecondCol, fixed = TRUE)))
    yc <- as.vector(as.numeric(sub(",", ".", FirstCol, fixed = TRUE)))
  }
  # and combine them into a dataframe
  XYdt <- data.frame(XC = xc, YC = yc)
  # Remove NA's
  XYdt <- XYdt[complete.cases(XYdt),]
  #
  if ((nrow(XYdt) > 0) && (max(XYdt) <= 1)) {
    XYdt <- XYdt * 100
  }
  #return data silently - should it be Visible or hidden?
  invisible(XYdt)
}
#
Name2Index <- function(chem_name) {
  db <- LLSR::llsr_data[["db.cas"]]
  chem_idx <- db[which(db$CAS.NAME == chem_name), "CAS.INDEX"]
  return(chem_idx)
}
#
###############################################################################
#' @rdname ExportTemplate
#' @export 
#' @title LLSR Template Exporter
#' @description The function makes a copy of LLSR's template file and copy it 
#' to the folder pointed by the user.
export_template <- function() {
  llsr_path <- system.file("extdata", package = "LLSR")
  template <- file.path(llsr_path, "template.xlsx")
  #
  output_dir <- dlgDir()$res
  #
  file.copy(
    from = template,
    to = output_dir,
    copy.mode = TRUE,
    copy.date = TRUE
  )
  #
}
#
###############################################################################
#' @rdname ExportData
#' @export 
#' @title LLSR Data Exporter
#' @description The function saves a copy of a specified variable to a file in 
#' the folder pointed by the user.
#' @param localData A variable existing in R environment and that will be saved
#'  locally.
export_data <- function(localData = NULL) {
  #
  if (is.null(localData)) {
    cat("Error: A variable containing data to be saved must be specified.")
  } else{
    output_dir <- dlgDir()$res
    file_local_database <- file.path(output_dir, "local_data.rda")
    #
    save(localData, file = file_local_database)
  }
  #
}
#
###############################################################################
# Recursive function to ensure that data on every data.frame is ASCII compliant 
to.ascii <- function(x) {
  if (class(x) == 'list') {
    for (nm in names(x)) {
      x[[nm]] <- to.ascii(x[[nm]])
    }
  } else
    if (class(x) == "data.frame") {
      for (col in colnames(x)) {
        x[, col] <- iconv(x[, col], from = "UTF-8", to = "ASCII", sub = "")
        } 
      }
  return(x)
}
