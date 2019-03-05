####################################################################################################################
options(digits = 14)
####################################################################################################################
#' @import ggplot2
#' @import rootSolve
####################################################################################################################
crit_point_seq <- function(dataSET,
                           db = LLSR::llsr_data,
                           modelName,
                           slope,
                           NP,
                           xmax,
                           xlbl,
                           ylbl,
                           Order,
                           ext) {
  #
  #
  #
  tol <- 1e-5
  if (is.null(slope)){
    slope = findSlope(db, dataSET)
  } 
  # Select which model will be used to generate the plot. Function return list of plots and respective number of parameters
  models_npars <- AQSysList(TRUE)
  #
  #
  #
  # Check data.frame validity and make an array of names for the systems if none is provided
  if ((ncol(dataSET) %% 2) == 0) {
    SysList <- list()  # List which will stack and return all data
    #
    SysData <- toNumeric(na.exclude(dataSET[6:nrow(dataSET), ]), Order)
    #
    # Analyse data and return parameters
    model_pars <- summary(AQSys(SysData,  modelName))$parameters[, 1]
    # Select Model based on the user choice or standard value
    Fn <- ifelse(
      modelName %in% names(models_npars),
      AQSys.mathDesc(modelName),
      AQSys.err("0")
    )
    # define a straight line EQUATION
    Gn <- function (yMin, AngCoeff, xMAX, x) {
      yMin + AngCoeff * (x - xMAX)
    }
    # Add constant variable values to the equations
    modelFn <- function(x) Fn(model_pars, x)
    modelTl <- function(x) Gn(yMin, slope[1], xMAX, x)
    # Decide whether xmax will be calculated or use the provided value
    xMAX <- ifelse((!is.numeric(xmax) | is.null(xmax) | (xmax == "")),
                   max(SysData[, seq(1, ncol(SysData), 2)] * 1.1), xmax)
    #
    yMAX <- max(SysData[, 2]) # get ymax from the dataset
    SysL <- list()
    # go from xmin/1.5 to xmax in steps of 0.15
    x <- seq(min(SysData[, 1]) / 1.5, xMAX, ((xMAX - (min(SysData[, 1]) / 1.5)) / NP))
    BNDL <- setNames(data.frame(x, modelFn(x)), c("X", "Y"))
    BNDL["System"] <- "CriticalPointSeries"
    #
    TL <- 0
    dt <- 1
    reset_BNDL <- BNDL
    CrptFnd <- FALSE # Crital Point Found Logical variable
    DivFactor <- 25
    cat("\t - Calculating: [")
    while ((dt > tol) && !CrptFnd) {
      TL <- TL + 1
      yMAXTL <- yMAX + 1
      # The following lines tests successively different values of xmax until a valid ymax, 
      # which must be smaller than the experimental ymax.
      while (yMAXTL > yMAX) {
        yMin <- Fn(model_pars, xMAX) # EVALUATE REPLACE XMAX TO THE LIMIT OF SOLUBILITY?
        xRoots <- uniroot.all(function(x)(modelFn(x) - modelTl(x)), c(0, xMAX), tol = 0.1) 
        xMin <- min(xRoots) # As we started from the biggest root, we select the smallest
        yMAXTL <- modelFn(xMin) # calculate the y-value for the root found
        if (yMAXTL > yMAX) { # compare to the highest value allowed (it must be smaller than the maximum experimental Y)
          xMAX <- (xMAX - 0.001) # if bigger than physically possible, decrease xmax and recalculate
        }
      }
      # create and name a data.frame containing coordinates for the tieline calculated above
      SysL[[TL]] <- setNames(data.frame(c(xMAX, (xMAX + xMin) / 2, xMin), c(yMin, (yMAXTL + yMin) / 2, yMAXTL)), c("X", "Y"))
      SysL[[TL]]["System"] <- paste("CriticalPointSeries", "TL", TL, sep = "-")
      # check if compositions of both phases, as well as the global composition, are equal. If so, critical point was found.
      CrptFnd <- is.equal(SysL[[TL]], tol)
      # Bind the calculated tieline's data.frame to the output data.frame variable
      BNDL <- rbind(BNDL,  SysL[[TL]])
      # A monod-base equation to help convergence - 
      dt <- abs(SysL[[TL]][2, 1] - modelFn(SysL[[TL]][2, 1])) / ((2 * TL) / (xMAX / DivFactor))
      xMAX <- xMAX - dt
      #
      if (TL > 1000) {
        DivFactor <- DivFactor * 0.75
        SysL <- list()
        TL <- 0
        dt <- 1
        BNDL <- reset_BNDL
        xMAX <- ifelse((!is.numeric(xmax) | is.null(xmax) | (xmax == "")),
                       max(SysData[, seq(1, ncol(SysData), 2)] * 1.1), xmax)
        yMAX <- max(SysData[, 2])
        cat("#")
      }
    }
    cat("#]\n")
    # data.frame holding data regarding Critical Point convergence
    output_res <- setNames(data.frame(dt, TL), c("dt", "TL"))
    # Setting up data.frame to hold data from the global points
    GlobalPoints <- setNames(data.frame(matrix(ncol = 2, nrow = length(SysL))), c("XG", "YG"))
    # transfer data to specific globalpoint data.frame. It will be used to calculate other system compositions for a given tieline.
    for (TL in seq(1, length(SysL))) {
      GlobalPoints[TL, 1:2] <- SysL[[TL]][2, 1:2]
    }
    # Note that the criteria for the loops above means that at the end, all compositions are equal. Thus, the last tieline found
    XC <- SysL[[length(SysL) - 1]][["X"]][1]
    YC <- SysL[[length(SysL) - 1]][["Y"]][1]
    OUTPUT <- setNames(data.frame(XC, YC), c("XC", "YC"))
    #
    #
    #
    if (ext) {
      #OUTPUT_PLOT <- bndOrthPlot(subset(BNDL, BNDL$System == "CriticalPointSeries"), xlbl, ylbl)
      #return(subset(BNDL, BNDL$System == "CriticalPointSeries"))
      OUTPUT_PLOT <- AQSys.plot(SysData,
                                silent = TRUE,
                                modelName = modelName,
                                xmax = xmax,
                                xlbl = xlbl,
                                ylbl = ylbl)
      OUTPUT_PLOT <- OUTPUT_PLOT + geom_line(
        data = subset(BNDL, BNDL$System != "CriticalPointSeries"),
        aes_string(x = "X", y = "Y", group = "System"),
        colour = "red",
        alpha = 0.4
      ) + geom_point(
        data = subset(BNDL, BNDL$System != "CriticalPointSeries"),
        aes_string(x = "X", y = "Y"),
        colour = "black",
        size = 1,
        alpha = 1
      )
      OUTPUT_PLOT <- OUTPUT_PLOT + annotate(
        "point",
        x = XC,
        y = YC,
        colour = "black",
        bg = "gold",
        shape = 23,
        size = 2
      ) + theme(legend.position = "none")
      return(list(CriticalPoint=OUTPUT, Plot=OUTPUT_PLOT))
    } 
    return(OUTPUT)
    
  } else{
    # Return an error if an invalid dataset is provided.
    AQSys.err(9)
  }
}
####################################################################################################################

