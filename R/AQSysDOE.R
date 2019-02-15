####################################################################################################################
options(digits = 14)
####################################################################################################################
#' @import rootSolve
####################################################################################################################
#' @rdname AQSysDOE
#' @name AQSysDOE
#' @title AQSysDOE
#' @description XXX.
#' @export AQSysDOE
#' 
#' @param dataSET - Binodal Experimental data that will be used in the nonlinear fit.
#' @param db A highly structure db containing data from previously analised data. LLSR database is used by default but user may input his own db if formatted properly.
#' @param slope The method assumes all tielines for a given ATPS are parallel, thus only one slope is required. [type:double]
#' @param xmax Maximum value for the Horizontal axis' value (bottom-rich component). [type:double]
#' @param modelName Character String specifying the nonlinear empirical equation to fit data.
#' The default method uses Merchuk's equation. Other mathematical descriptors can be listed using AQSysList(). [type:string]
#' @param nTL Number of tielines plotted for a given ATPS. Default is 3. [type:Integer]
#' @param nPoints Number of points chosen for a given tieline. Default is 3. [type:Integer]
#' @param tol limit of tolerance to reach to assume convergence. Default is 1e-5. [type:Integer]
# ' @param maxiter	- A positive integer specifying the maximum number of iterations allowed.
AQSysDOE <- function(dataSET,
                     db = LLSR::llsr_data,
                     slope = NULL,
                     xmax = NULL,
                     modelName = "merchuk",
                     nTL = 3,
                     nPoints = 3,
                     tol = 1e-5) {
  #
  if (is.null(slope)){
    slope = findSlope(db, dataSET)
  } else if (!((ncol(dataSET) / 2) == length(slope))){
    AQSys.err("11")
  }
  #
  rawEvalData <- AQSysEval(
    dataSET,
    tol = tol,
    nTL = nTL,
    nPoints = nPoints,
    modelName = modelName,
    slope = slope
  )
  #
  if (length(rawEvalData$data)==length(rawEvalData$plot)) {
    SysCharData <- rawEvalData$data
  } else {
    SysCharData <- list()
    SysCharData[[1]] <- rawEvalData$data
  }
  #
  dataNames <- c("X", "Y", "System", "TLL", "Point")
  OUTPUT <- setNames(as.data.frame(matrix(ncol = 5)), dataNames)
  TLLs <- as.data.frame(matrix(ncol = 2, nrow = length(SysCharData)))
  for (idx in seq(1, length(SysCharData))){
    TLLs [idx, 1:2] <- as.data.frame(SysCharData[[idx]]$TLL)
  }
  # values are rounded 1% to make sure they are within the method boundaries
  seqTLL <- seq(round(max(TLLs[, 1]), 2), round(min(TLLs[which(TLLs$V2 > max(TLLs$V1)), names(TLLs)][, 2]), 2)*0.99, length.out = nTL)
  #return(list(A=TLLs, B=seqTLL))
  # Select which model will be used to generate the plot. Function return list of plots and respective number of parameters
  models_npars <- AQSysList(TRUE)
  # Select Model based on the user choice or standard value
  Fn <- ifelse(
    modelName %in% names(models_npars),
    AQSys.mathDesc(modelName),
    AQSys.err("0")
  )
  #
  for (idx in seq_along(SysCharData)) {
    for (TLL in seq_along(seqTLL)) {
      modelFn <- function(x) Fn(SysCharData[[idx]]$PARs, x)
      data <- findTL(seqTLL[TLL], SysCharData[[idx]], modelFn, slope[idx])
      #
      temp.TLC <- setNames(data.frame(data$TL, rep(idx, 3), rep(data$TLL, 3), c("T", "M", "B")), dataNames)
      X <- mean(data$TL$X)
      Y <- mean(data$TL$Y)
      #
      TLFn <- function(x) { Y + slope[idx] * (x - X) }
      #
      xSYS <- seq(min(data$TL$X), max(data$TL$X), length.out = (nPoints + 2))
      temp.SYS <- setNames(data.frame(round(xSYS, 4), round(TLFn(xSYS), 4), rep(idx, (nPoints + 2)), rep(data$TLL, (nPoints + 2)), rep("S", (nPoints + 2))), dataNames)
      #
      OUTPUT <- rbind(OUTPUT, temp.TLC, temp.SYS)
    }
  }
  invisible(list("DOE" = OUTPUT[-1,], "data" = SysCharData))
}
