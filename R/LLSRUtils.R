## quiets concerns of R CMD check re: the .'s that appear in pipelines
##if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

options(digits = 14)

#'  @rdname LLSR.info
#'  @export LLSR.info
#'  @title .
#'  @description .
LLSR.info <- function() {
  print("info Test")
}
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

is.equal <- function(TLData, tol) {
  # Verify if all Xs and Ys are equal among themselves
  return ((all((
    max(TLData[, 1]) - min(TLData[, 1])
  ) < tol) == TRUE) &&
    (all((
      max(TLData[, 2]) - min(TLData[, 2])
    ) < tol) == TRUE))
}

FindMinTL <- function(SysCP, maxGP, xMax, slope, BLFn, tol, dfr = 0.2){
  # Calculate distance between the critical point (SysCP) and the furthest viable TL's Global Point (maxGP)
  DMaxTL <- sqrt((maxGP[1] - SysCP[1]) ^ 2 + (maxGP[2] - SysCP[2]) ^ 2)
  # Establishes that the distance between the critical point and the closest viable tieline is a fraction (dfr) of the distance to MaxTL
  d <- DMaxTL * dfr
  # Calculate the coordinates for a point existing in the minTL (closest viable TL)
  X <- as.numeric(SysCP[1] + sqrt(d ^ 2 / (1 + (slope) ^ 2)))
  Y <- as.numeric(SysCP[2] + (-1 / slope) * (X - SysCP[1]))
  # Initiate a generic line function but using calculated X, Y and provided Slope (assuming all tielines are parallel)
  TLFn <- function(x) { Y + slope * (x - X) }
  # find the intersections between the minTL and the binodal
  xRoots <- uniroot.all(function(x) (BLFn(x) - TLFn(x)), c(0, xMax * 1.5), tol = 0.1) # REPLACE XMAX TO THE LIMIT OF SOLUBILITY?
  # Creates an array containing all Xs which characterizes minTL
  xTL <- c(min(xRoots), sum(xRoots) / 2, max(xRoots))
  # return a data.frame containing all Xs and Ys for minTL
  return(setNames(data.frame(xTL, TLFn(xTL)), c("X", "Y")))
}

dPoints <- function(P1, P2){
  d <- sqrt(((P1$X - P2$X) ^ 2) + ((P1$Y - P2$Y) ^ 2))
  return(d)
}

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

findTL <- function(dTLL, SysTLL, BLFn, slope){
  # If the target TLL is smaller than the minimum calculated TLL, throw an error.
  if ((dTLL < SysTLL$TLL$MinTLL) | (dTLL > SysTLL$TLL$MaxTLL)){
    #AQSys.err("10")
    print(c(dTLL, SysTLL$TLL$MinTLL, SysTLL$TLL$MaxTLL))
  }
  # Initial guess for the tieline, calculated using the tieline length proportions
  X <- (dTLL / SysTLL$TLL$MaxTLL)*max(SysTLL$maxTL$X)
  # Initializing variables
  TL <- setNames(as.data.frame(matrix(nrow = 3, ncol = 2)), c("X", "Y"))
  OUTPUT <- list()
  dt <- 1
  #
  while (dt > 1e-7){
    Y <- BLFn(X)
    TLFn <- function(x) { Y + slope * (x - X) }
    xRoots <- uniroot.all(function(x) (BLFn(x) - TLFn(x)), c(0, X), tol = 0.1) 
    #
    TL[1, 1] <- min(xRoots)
    TL[3, 1] <- max(xRoots)
    TL[1, 2] <- unname(BLFn(TL[1, 1]))
    TL[3, 2] <- unname(BLFn(TL[3, 1]))
    #
    TLL <- dPoints(TL[1, ], TL[3, ])
    dt <- abs(TLL - dTLL)
    #
    X <- X + dt * ( -(TLL - dTLL) / abs(TLL - dTLL)) / 10
  }
  TL[2, 1] <- (TL[1, 1] + TL[3, 1]) / 2
  TL[2, 2] <- (TL[1, 2] + TL[3, 2]) / 2
  #
  OUTPUT$TL <- TL
  OUTPUT$TLL <- TLL
  #
  return(OUTPUT)
}

findSlope <- function(dataSET){
  # iterate through multiple columns and return a list of slopes
  slope <- c()
  # check how many systems were provided
  nSys <- (ncol(dataSET) / 2)
  if ((ncol(dataSET) %% 2) == 0) {
    for (sys in seq(1, nSys)){
      # Get the system characterization variables
      idx_Y <- dataSET[3, sys * 2 - 1]
      idx_X <- dataSET[3, sys * 2]
      idx_PH <- dataSET[1, sys * 2 - 1]
      idx_T <- dataSET[2, sys * 2 - 1]
      #
      # CAS_db <- llsr_data[["db.cas"]]
      #print(c(idx_Y, idx_X, idx_PH, idx_T) )
      #print(c(TL_db$Y, TL_db$X, TL_db$PH, TL_db$T) )
      TL_db <- LLSR::llsr_data[["db.tielines"]]
      slope[sys] <- TL_db[which(
        (TL_db$PH == idx_PH | is.na(TL_db$PH)) &
          TL_db$T == idx_T &
          TL_db$X == idx_X &
          TL_db$Y == idx_Y ), "TLSlope"]
      #
    }
    if (length(slope)==0) {
      AQSys.err("12")
    } else {
      return(slope)
    }
  } else{
    # Return an error if an invalid dataset is provided.
    AQSys.err("9")
  }
}

seqTL <- function(minTL, maxTL, slope, BLFn, nTL = 3, nSYS = 3) {
  dataNames <- c("X", "Y", "System", "Point")
  #
  xMin <- max(minTL["X"])
  xMax <- max(maxTL["X"])
  # Bottom Phase Compositions
  xRange <- seq(xMin, xMax, (xMax - xMin) / (nTL - 1))
  oDATA <- setNames(data.frame(xRange, BLFn(xRange), seq(1, nTL), "B", row.names = NULL), dataNames) 
  #
  for (p in seq(1, nrow(oDATA))) {
    X <- oDATA[p, "X"]
    Y <- oDATA[p, "Y"]
    #
    TLFn <- function(x) { Y + slope * (x - X) }
    xRoots <- uniroot.all(function(x) (BLFn(x) - TLFn(x)), c(0, X*2), tol = 0.1) # REPLACE XMAX TO THE LIMIT OF SOLUBILITY?
    xTL <- c(min(xRoots), sum(xRoots) / 2)
    #
    temp.TLC <- setNames(data.frame(xTL, TLFn(xTL), rep(oDATA[p, "System"], 2), c("U", "G")), dataNames)
    #
    xSYS <- seq(min(xRoots), X, (X - min(xRoots)) / (nSYS + 1))
    #
    temp.SYS <- setNames(data.frame(xSYS, TLFn(xSYS), rep(oDATA[p, "System"], nSYS + 2), rep("S", nSYS + 2)), dataNames)
    #
    oDATA <- rbind(oDATA, temp.TLC, temp.SYS)
    #
  }
  return(oDATA)
}

saveConfig <- function (plot_obj, save, HR, filename, wdir, silent) {
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


LLSRxy <- function(FirstCol, SecondCol, ColDis = 'xy') {
  # convert and name variables accordingly into vectors
  if (tolower(ColDis) == "xy") {
    xc <- as.vector(as.numeric(sub(",", ".", FirstCol, fixed = TRUE)))
    yc <-
      as.vector(as.numeric(sub(",", ".", SecondCol, fixed = TRUE)))
  } else{
    xc <- as.vector(as.numeric(sub(",", ".", SecondCol, fixed = TRUE)))
    yc <-
      as.vector(as.numeric(sub(",", ".", FirstCol, fixed = TRUE)))
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
