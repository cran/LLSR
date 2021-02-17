###############################################################################
options(digits = 14)
###############################################################################
#' @import rootSolve
###############################################################################
#' @rdname AQSysEval
#' @name AQSysEval
#' @title AQSysEval
#' @description The function perform a full ATPS characterization (parameters,
#'  tie-line boundaries and critical point), generating a brief report.
#' @export AQSysEval
#' @param dataSET - Binodal Experimental data that will be used in the nonlinear
#'  fit. [type:data.frame]
#' @param db A highly structure db containing data from previously analyzed 
#' data. LLSR database is used by default but user may input his own db if 
#' formatted properly.
#' @param xmax Maximum value for the Horizontal axis' value (bottom-rich
#'  component). [type:double]
#' @param ymax Maximum value for the vertical axis' value 
#' (bottom-rich component). [type:double]
#' @param NP Number of points used to build the fitted curve. 
#' Default is 100. [type:Integer]
#' @param slope The method assumes all tielines for a given ATPS are parallel,
#'  thus only one slope is required. [type:double]
#' @param modelName Character String specifying the nonlinear empirical
#'  equation to fit data.
#' The default method uses Merchuk's equation. Other mathematical descriptors 
#' can be listed using AQSysList(). [type:string]
#' @param convrgnceLines Magnify Plot's text to be compatible with High
#'  Resolution size [type:Logical]
#' @param nTL Number of tielines plotted for a given ATPS. 
#' Default is 3. [type:Integer]
#' @param nPoints Number of points chosen for a given tieline. 
#' Default is 3. [type:Integer]
#' @param tol limit of tolerance to reach to assume convergence.
#'  Default is 1e-5. [type:Integer]
#' @param xlbl Plot's Horizontal axis label. [type:String]
#' @param ylbl Plot's Vertical axis label. [type:String]
#' @param seriesNames Number of points used to build the fitted curve. 
#' Default is 100. [type:Integer]
#' @param save Save the generated plot in the disk using path and filename 
#' provided by the user. Default is FALSE. [type:Logical]
#' @param HR Magnify Plot's text to be compatible with High Resolution
#'  size [type:Logical]
#' @param autoname Number of points used to build the fitted curve. 
#' Default is FALSE. [type:Logical]
#' @param wdir The directory in which the plot file will be saved. [type:String]
#' @param silent save plot file without actually showing it to the user. 
#' Default is FALSE. [type:Logical]
# ' @param maxiter	- A positive integer specifying the maximum number of
#  iterations allowed.
#' 
#' @examples
#' \dontrun{
#' dataSET <- AQSearch.Binodal(db.uid='56b53a50f500c502fa4a65d197fc6d84')
#' xLabel <- "Ammonium Sulphate" 
#' yLabel <- "Poly(ethylene glycol) 2000" 
#' EvalData <- AQSysEval(dataSET2 , xlbl = xLabel, ylbl = yLabel)
#' }
#' @references 
#' KAUL, A. The Phase Diagram. In: HATTI-KAUL, R. (Ed.). 
#' Aqueous Two-Phase Systems: Methods and Protocols: Humana Press, v.11, 2000. 
#' cap. 2, p.11-21.  (Methods in Biotechnology). ISBN 978-0-89603-541-6.
#' (\href{https://link.springer.com/10.1385/1-59259-028-4:11}{SpringerLink})
#' 
AQSysEval <- function(dataSET,
                      db = LLSR::llsr_data,
                      xmax = NULL,
                      ymax = NULL,
                      NP = 100,
                      slope = NULL,
                      modelName = "merchuk",
                      convrgnceLines = FALSE,
                      nTL = 3,
                      nPoints = 3,
                      tol = 1e-4,
                      xlbl = "",
                      ylbl = "",
                      seriesNames = NULL,
                      save = FALSE,
                      HR = FALSE,
                      autoname = FALSE,
                      wdir = NULL,
                      silent = TRUE) {
  #
  if (is.null(slope)) {
    slope = as.numeric(findSlope(db, dataSET))
  } else if (!((ncol(dataSET) / 2) == length(slope))) {
    AQSys.err("11")
  }
  # Select which model will be used to generate the plot. Function return list
  #  of plots and respective number of parameters
  models_npars <- AQSysList(TRUE)
  # Divides the number of columns of the data set per two to calculate how many
  #  systems it hold. result must be even.
  nSys <- (ncol(dataSET) / 2)
  SysNames <- FALSE
  # Check data.frame validity and make an array of names for the systems if none
  #  is provided
  if ((ncol(dataSET) %% 2) == 0) {
    if (is.null(seriesNames) || !(length(seriesNames) == nSys)) {
      cat("[AQSysEval]\n")
      cat(
        paste(
          "\t - The array seriesNames must have",
          nSys,
          "element(s). \n \t - Default names will be used instead.\n"
        )
      )
      seriesNames <- sapply(seq(1, nSys), function(x) paste("Series", x))
    } else {
      SysNames <- TRUE
    }
    SysList <- list()  # List which will stack and return all data
    PlotList <- list() # List which will stack and return all plots
    for (i in seq(1, nSys)) {
      #
      RawData <- dataSET[, (i * 2 - 1):(i * 2)]
      SysOrder <- RawData[4, 2]
      SysData <- LLSRxy(na.exclude(RawData[6:nrow(RawData), 1]),
                        na.exclude(RawData[6:nrow(RawData), 2]),
                        SysOrder)
      # Analyse data and return parameters
      model_pars <- summary(AQSys(SysData,  modelName))$parameters[, 1]
      # Select Model based on the user choice or standard value
      Fn <- ifelse(
          modelName %in% names(models_npars),
          AQSys.mathDesc(modelName),
          AQSys.err("0")
        )
      # define a straight line EQUATION
      Gn <- function(yMin, AngCoeff, xMAX, x) {
          yMin + AngCoeff * (x - xMAX)
        }
      # Add constant variable values to the equations
      modelFn <- function(x) Fn(model_pars, x)
      modelTl <- function(x) Gn(yMin, slope[i], xMAX, x)
      # Decide whether xmax will be calculated or use the provided value
      yMAX <- ifelse(is.null(ymax), max(SysData[, 2])*0.9, ymax) 
      # get ymax from the dataset
      xMAX <- ifelse(is.null(xmax), 
                     max(SysData[, seq(1, ncol(SysData), 2)])*0.9, xmax)
      SysList[[i]] <- list()
      SysL <- list()
      #
      BNDL <- GenPlotSeries(SysData, xMAX, NP, modelFn, i, seriesNames)
      #
      xMAX <- ChckBndrs(modelFn, slope[i], BNDL, xMAX)
      reset_xMAX <- xMAX
      #
      BNDL <- GenPlotSeries(SysData, xMAX*1.15, NP, modelFn, i, seriesNames)
      #
      TL <- 0
      dt <- 1
      reset_BNDL <- BNDL
      CrptFnd <- FALSE # Critical Point Found Logical variable
      DivFactor <- 25
      #
      cat(paste0("\t - Calculating System #",  i," ["))
      while ((dt > tol) && !CrptFnd) {
        TL <- TL + 1
        yMAXTL <- yMAX + 1
        # The following lines tests successively different values of xmax until
        #  a valid ymax, 
        # which must be smaller than the experimental ymax.
        while (yMAXTL > yMAX) {
          yMin <- Fn(model_pars, xMAX) 
          # EVALUATE REPLACE XMAX TO THE LIMIT OF SOLUBILITY?
          xRoots <- uniroot.all(function(x)(modelFn(x) - modelTl(x)), 
                                c(0, xMAX), tol = 0.001) 
          xMin <- min(xRoots) 
          # As we started from the biggest root, we select the smallest
          yMAXTL <- modelFn(xMin) # calculate the y-value for the root found
          if ((yMAXTL > yMAX) | (length(xRoots) == 1)) { 
            # compare to the highest value allowed (it must be smaller than the
            #  maximum experimental Y)
            xMAX <- (xMAX - 0.001) 
            # if bigger than physically possible, decrease xmax and recalculate
          }
        }
        # create and name a data.frame containing coordinates for the tieline
        #  calculated above
        SysL[[TL]] <- setNames(data.frame(c(xMAX, (xMAX + xMin) / 2, xMin), 
                                          c(yMin, (yMAXTL + yMin) / 2, yMAXTL)),
                               c("X", "Y"))
        SysL[[TL]]["System"] <- paste(seriesNames[i], "TL", TL, sep = "-")
        # check if compositions of both phases, as well as the global 
        # composition, are equal. If so, critical point was found.
        CrptFnd <- is.equal(SysL[[TL]], 1e-3)
        # Bind the calculated tieline's data.frame to the output data.frame 
        # variable
        BNDL <- rbind(BNDL,  SysL[[TL]])
        # A monod-base equation to help convergence - 
        dt <- abs(SysL[[TL]][2, 1] - modelFn(SysL[[TL]][2, 1])) 
        xMAX <- xMAX - (dt / ((5 * TL) / (xMAX / DivFactor)))
        #
        if (TL == 1250) {
          DivFactor <- DivFactor * 0.9
          SysL <- list()
          TL <- 0
          dt <- 1
          BNDL <- reset_BNDL
          xMAX <- reset_xMAX
          yMAX <- ifelse(is.null(ymax), max(SysData[, 2])*1.1, ymax)
          cat("#")
        }
      }
      cat("#]\n")
      # data.frame holding data regarding Critical Point convergence
      output_res <- setNames(data.frame(dt, TL), c("dt", "TL"))
      # Setting up data.frame to hold data from the global points
      GlobalPoints <- setNames(
        data.frame(matrix(ncol = 2, nrow = length(SysL))), c("XG", "YG")
        )
      # transfer data to specific globalpoint data.frame. It will be used to
      #  calculate other system compositions for a given tieline.
      for (TL in seq(1, length(SysL))) {
        GlobalPoints[TL, 1:2] <- SysL[[TL]][2, 1:2]
      }
      SysL$GlobalPoints <- GlobalPoints
      # Note that the criteria for the loops above means that at the end, all
      #  compositions are equal. Thus, the last tieline found
      # essentially satisfy the definition of critical point. The lines below
      #  add such point in a specific dataframe for the system
      # under study/calculation
      XC <- SysL[[length(SysL) - 1]][["X"]][2]
      YC <- SysL[[length(SysL) - 1]][["Y"]][2]
      SysCvP <- setNames(data.frame(XC, YC), c("XC", "YC"))
      SysL$CriticalPoint <- SysCvP
      # Print tieline's convergence data for the system
      # print(round(cbind(SysCvP, output_res), 4))
      # Max Tieline will be the first 'viable' tieline, i.e., which xmax yield
      #  a ymax within the experimental range
      maxTL <- SysL[[1]]
      SysL$maxTL <- maxTL[, 1:2]
      # A procedure to find the minimum tieline was written as a standalone 
      # function
      minTL <- FindMinTL(SysCvP, GlobalPoints[1, ], max(maxTL["X"]), 
                         slope[i], modelFn, tol)
      SysL$minTL <- minTL
      # when desired, a experimental design matrix with n tielines and m points
      #  is provided
      SysL$DOE <- seqTL(minTL, maxTL, slope[i], modelFn, nTL, nPoints)
      # when desired, a experimental design matrix with n tielines and m points
      #  is provided
      SysL$TLL <- TLL(minTL, maxTL)
      # a hypothetical X=Y critical point is calculate -> probably deprecated
      #  and removed soon
      # rootCritical <- uniroot(function(x)(modelFn(x) - x), c(0, xMAX))$root
      # SysL$CriticalPoint <- rootCritical
      # add the Empirical Model parameters to the list, in order to allow 
      # subsequent use with no re-calculation
      SysL$PARs <- model_pars
      # add the Empirical Model's name for further use
      SysL$modelName <- modelName
      # prepare a plot with the system curve and all tielines. Convergence
      #  lines are included if requested.
      output_plot <- bndOrthPlot(dataSET = subset(
        BNDL, BNDL$System == seriesNames[i]), Order = SysOrder, xlbl = xlbl ,
        ylbl = ylbl) 
      if (convrgnceLines) {
        output_plot <- output_plot + 
          geom_line(
            data = subset(BNDL, BNDL$System != seriesNames[i]),
            aes_string(x = "X", y = "Y", group = "System"),
            colour = "red",
            alpha = 0.4
          ) +
          geom_point(
            data = subset(BNDL, BNDL$System != seriesNames[i]),
            aes_string(x = "X", y = "Y"),
            colour = "black",
            size = 1,
            alpha = 1
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
            data = subset(SysL$DOE, SysL$DOE$Point %in% c("T", "M", "B")),
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
        ) 
      # + annotate("point", x = rootCritical,y = rootCritical,
      #  colour = "olivedrab2", shape = 18, size = 3)
      # Add plot to a list in order to be acessible to the user after method 
      # completes
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
    invisible(list(
      "data" = if (length(SysList) > 1) {
        SysList
      } else{
        SysList[[1]]
      },
      "plot" = if (length(PlotList) > 1) {
        SysList
      } else{
        PlotList[[1]]
      }
    ))
  } else{
    # Return an error if an invalid dataset is provided.
    AQSys.err(9)
  }
}
