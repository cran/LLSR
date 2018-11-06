####################################################################################################################
options(digits = 14)
####################################################################################################################
#' @import rootSolve
####################################################################################################################
#' @rdname AQSysEval
#' @name AQSysEval
#' @title AQSysEval
#' @description Import DB data from an Excel Worksheet.
#' @export AQSysEval
#' 
#' @param dataSET - Binodal Experimental data that will be used in the nonlinear fit
#' @param xmax Maximum value for the Horizontal axis' value (bottom-rich component). [type:double]
#' @param NP Number of points used to build the fitted curve. Default is 100. [type:Integer]
#' @param TLSlope The method assumes all tielines for a given ATPS are parallel, thus only one slope is required. [type:double]
#' @param modelName Character String specifying the nonlinear empirical equation to fit data.
#' The default method uses Merchuk's equation. Other mathematical descriptors can be listed using AQSysList(). [type:string]
#' @param convrgnceLines Magnify Plot's text to be compatible with High Resolution size [type:Logical]
#' @param nTL Number of tielines plotted for a given ATPS. Default is 3. [type:Integer]
#' @param nPoints Number of points chosen for a given tieline. Default is 3. [type:Integer]
#' @param tol limit of tolerance to reach to assume convergence. Default is 1e-5. [type:Integer]
#' @param xlbl Plot's Horizontal axis label.
#' @param ylbl Plot's Vertical axis label.
#' @param seriesNames Number of points used to build the fitted curve. Default is 100. [type:Integer]
#' @param save Save the generated plot in the disk using path and filename provided by the user. Default is FALSE. [type:Logical]
#' @param HR Magnify Plot's text to be compatible with High Resolution size [type:Logical]
#' @param autoname Number of points used to build the fitted curve. Default is FALSE. [type:Logical]
#' @param wdir The directory in which the plot file will be saved. [type:String]
#' @param silent save plot file without actually showing it to the user. Default is FALSE. [type:Logical]
#' 
#' @examples
#' \dontrun{
#' AQSysDB("C:/data.xls", ColDis = "xy", CAS = FALSE)
#'}
####################################################################################################################
# AQSysDB() is a simple approach that is ready to use any three-parameter equation
# and thus
#
# GIVEN A CURVE, A LINE SLOPE AND BOTH XMAX AND YMAX - FIND TWO POINTS (A LINE) WHICH LIES WITHIN THE
# CURVE LIMITS AND HAVE THE SAME SLOPE AS THE TIELINES

# PROCESS WILL BE INTERACTIVE BECAUSE A GIVEN XMAX MIGHT YIELD A YMAX ABOVE THE EXPERIMENTAL DATA OR A YMAX MIGHT YIELD A XMAX
# ALSO OUTSIDE THE EXPERIMENTAL RANGE
# SO THE POINT SELECTION WILL RUN UNTIL A SUITABLE PAIR IS FOUND
#dataSET <- dados[["db.data"]][, 93:104]
#
AQSysEval <- function(dataSET,
                      xmax = 30,
                      NP = 100,
                      TLSlope = -1.846924,
                      modelName = "merchuk",
                      convrgnceLines = FALSE,
                      nTL = 3,
                      nPoints = 3,
                      tol = 1e-5,
                      xlbl = "",
                      ylbl = "",
                      seriesNames = NULL,
                      save = FALSE,
                      HR = FALSE,
                      autoname = FALSE,
                      wdir = NULL,
                      silent = FALSE) {
  # Select which model will be used to generate the plot. Function return list of plots and respective number of parameters
  models_npars <- AQSysList(TRUE)
  # Divides the number of columns of the data set per two to calculate how many systems it hold. result must be even.
  nSys <- (ncol(dataSET) / 2)
  SysNames <- FALSE
  # Check data.frame validity and make an array of names for the systems if none is provided
  if ((ncol(dataSET) %% 2) == 0) {
    if (is.null(seriesNames) || !(length(seriesNames) == nSys)) {
      print(paste("The array seriesNames must have", nSys, "elements. Default names will be used instead."))
      seriesNames <- sapply(seq(1, nSys), function(x) paste("Series", x))
    } else {
      SysNames <- TRUE
    }
    SysList <- list()  # List which will stack and return all data
    PlotList <- list() # List which will stack and return all plots
    for (i in seq(1, nSys)) {
      #
      RawData <- dataSET[, (i * 2 - 1):(i * 2)]
      SysData <- LLSRxy(na.exclude(RawData[6:nrow(RawData), 1]),
                        na.exclude(RawData[6:nrow(RawData), 2]),
                        RawData[4, 2])
      # Analyse data and return parameters
      model_pars <- summary(AQSys(SysData,  modelName))$parameters[, 1]
      # Select Model based on the user choice or standard value
      Fn <- ifelse(
          modelName %in% names(models_npars),
          AQSys.mathDesc(modelName),
          AQSys.err("0")
        )
      # define a straight line EQ
      Gn <- function (yMin, TLSlope, xMAX, x) {
          yMin + TLSlope * (x - xMAX)
        }
      # Add constant variable values to the equations
      modelFn <- function(x) Fn(model_pars, x)
      modelTl <- function(x) Gn(yMin, TLSlope, xMAX, x)
      # Decide whether xmax will be calculated or use the provided value
      xMAX <- ifelse(is.null(xmax), max(SysData[, 1]), xmax)
      yMAX <- max(SysData[, 2]) # get ymax from the dataset
      SysList[[i]] <- list()
      SysL <- list()
      # go from xmin/1.5 to xmax in steps of 0.15
      x <- seq(min(SysData[, 1]) / 1.5, xMAX, ((xMAX - (min(SysData[, 1]) / 1.5)) / NP))
      BNDL <- setNames(data.frame(x, modelFn(x)), c("X", "Y"))
      BNDL["System"] <- seriesNames[i]
      #
      TL <- 0
      dt <- 1
      CrptFnd <- FALSE # Crital Point Found Logical variable
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
            xMAX <- (xMAX - 0.001) # if bigger than possible, decrease xmax and recalculate
          }
        }
        # create and name a data.frame containing coordinates for the tieline calculated above
        SysL[[TL]] <- setNames(data.frame(c(xMAX, (xMAX + xMin) / 2, xMin), c(yMin, (yMAXTL + yMin) / 2, yMAXTL)), c("X", "Y"))
        SysL[[TL]]["System"] <- paste(seriesNames[i], "TL", TL, sep = "-")
        # check if compositions of both phases, as well as the global composition, are equal. If so, critical point was found.
        CrptFnd <- is.equal(SysL[[TL]], tol)
        # Bind the calculated tieline's data.frame to the output data.frame variable
        BNDL <- rbind(BNDL,  SysL[[TL]])
        # A monod-base equation to help convergence
        dt <- abs(SysL[[TL]][2, 1] - modelFn(SysL[[TL]][2, 1])) / ((25 * TL) / (xMAX / 2))
        xMAX <- xMAX - dt
      }
      # data.frame holding data regarding Critical Point convergence
      output_res <- setNames(data.frame(dt, TL), c("dt", "TL"))
      # Setting up data.frame to hold data from the global points
      GlobalPoints <- setNames(data.frame(matrix(ncol = 2, nrow = length(SysL))), c("XG", "YG"))
      # transfer data to specific globalpoint data.frame. It will be used to calculate other system compositions for a given tieline.
      for (TL in seq(1, length(SysL))) {
        GlobalPoints[TL, 1:2] <- SysL[[TL]][2, 1:2]
      }
      SysL$GlobalPoints <- GlobalPoints
      # Note that the criteria for the loops above means that at the end, all compositions are equal. Thus, the last tieline found
      # essentially satisfy the definition of critical point. The lines below add such point in a specific dataframe for the system
      # under study/calculation
      XC <- SysL[[length(SysL)]][["X"]][1]
      YC <- SysL[[length(SysL)]][["Y"]][1]
      SysCvP <- setNames(data.frame(XC, YC), c("XC", "YC"))
      SysL$ConvergencePoint <- SysCvP
      # Print tieline's convergence data for the system
      # print(round(cbind(SysCvP, output_res), 4))
      # Max Tieline will be the first 'viable' tieline, i.e., which xmax yield a ymax within the experimental range
      maxTL <- SysL[[1]]
      # A procedure to find the minimum tieline was written as a standalone function
      minTL <- FindMinTL(SysCvP, GlobalPoints[1, ], max(maxTL["X"]), TLSlope, modelFn, tol)
      SysL$minTL <- minTL
      # when desired, a experimental design matrix with n tielines and m points is provided
      SysL$DOE <- seqTL(minTL, maxTL, TLSlope, modelFn, nTL, nPoints)
      # a hypothetical X=Y critical point is calculate -> probably deprecated and removed soon
      rootCritical <- uniroot(function(x)(modelFn(x) - x), c(0, xMAX))$root
      SysL$CriticalPoint <- rootCritical
      # prepare a plot with the system curve and all tielines. Convergence lines are included if requested.
      output_plot <- bndOrthPlot(subset(BNDL, BNDL$System == seriesNames[i]), xlbl, ylbl) 
      if (convrgnceLines){
        output_plot <- output_plot +
          geom_point(
            data = subset(BNDL, BNDL$System != seriesNames[i]),
            aes_string(x = "X", y = "Y"),
            colour = "red",
            size = 1,
            alpha = 0.4
          ) +
          geom_line(
            data = subset(BNDL, BNDL$System != seriesNames[i]),
            aes_string(x = "X", y = "Y", group = "System"),
            colour = "red",
            alpha = 0.4
          )
      } else 
        {
        output_plot <- output_plot  +
          geom_point(
            data = SysL$minTL,
            aes_string(x = "X", y = "Y"),
            colour = "red",
            size = 2,
            alpha = 1
          ) +
          geom_line(
            data = SysL$minTL,
            aes_string(x = "X", y = "Y"),
            colour = "red",
            size = 1, 
            alpha = 1
          ) +
          geom_point(
            data = maxTL,
            aes_string(x = "X", y = "Y"),
            colour = "red",
            size = 2,
            alpha = 1
          ) +
          geom_line(
            data = maxTL,
            aes_string(x = "X", y = "Y"),
            colour = "red",
            size = 1,
            alpha = 1
          ) + 
          geom_line(
            data = subset(SysL$DOE, SysL$DOE$Point %in% c("U", "G", "B")),
            aes_string(x = "X", y = "Y", group = "System"),
            colour = "#d11141",
            size = 3,
            alpha = 0.4
          ) +
          geom_point(
            data = subset(SysL$DOE, SysL$DOE$Point %in% "S"),
            aes_string(x = "X", y = "Y"),
            colour = "black",
            shape = 17,
            size = 3,
            alpha = 1
          )
        }
      output_plot <- output_plot +
        annotate(
          "point",
          x = XC,
          y = YC,
          colour = ifelse(convrgnceLines, "black", "red"),
          bg = ifelse(convrgnceLines, "gold", "red"),
          shape = 23,
          size = 2
        ) # + annotate("point", x = rootCritical,y = rootCritical, colour = "olivedrab2", shape = 18, size = 3)
      # Add plot to a list in order to be acessible to the user after method completes
      PlotList[[i]] <- output_plot
      # execute saving procedures
      if (autoname) {filename <- seriesNames[i]} else {wdir <- NULL}
      wdir <- saveConfig(output_plot, save, HR, filename, wdir, silent)
      if (silent == FALSE) {
        print(output_plot)
      }
      # add all data to a data.frame indexed by the system input order
      SysList[[i]] <- SysL
    }
    # return silently data.frame to the user
    invisible(list("data" = SysList, "plot" = output_plot))
  } else{
    # Return an error if an invalid dataset is provided.
    AQSys.err(9)
  }
}
