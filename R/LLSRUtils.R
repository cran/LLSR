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
  xRoots <- uniroot.all(function(x) (BLFn(x) - TLFn(x)), c(0, xMax), tol = 0.1) # REPLACE XMAX TO THE LIMIT OF SOLUBILITY?
  # Creates an array containing all Xs which characterizes minTL
  xTL <- c(min(xRoots), sum(xRoots) / 2, max(xRoots))
  # return a data.frame containing all Xs and Ys for minTL
  return(setNames(data.frame(xTL, TLFn(xTL)), c("X", "Y")))
}

seqTL <- function(minTL, maxTL, slope, BLFn, nTL = 3, nSYS = 3) {
  dataNames <- c("X", "Y", "System", "Point")
  #
  xMin <- max(minTL["X"])
  xMax <- max(maxTL["X"])
  # Bottom Phase Compositions
  xRange <- seq(xMin, xMax, (xMax - xMin) / (nTL - 1))
  oDATA <- setNames(data.frame(xRange, BLFn(xRange), seq(1, nTL), "B"), dataNames) 
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
