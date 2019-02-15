####################################################################################################################
options(digits = 14)
####################################################################################################################
#' @import rootSolve
####################################################################################################################
#' @rdname AQSysAnima
#' @name AQSysAnima
#' @title AQSysAnima
#' @description Import DB data from an Excel Worksheet.
#' @export AQSysAnima
#'
#' @param dataSET - Binodal Experimental data that will be used in the nonlinear fit
#' @param db A highly structure db containing data from previously analised data. LLSR database is used by default but user may input his own db if formatted properly.
#' @param xmax Maximum value for the Horizontal axis' value (bottom-rich component). [type:double]
#' @param NP Number of points used to build the fitted curve. Default is 100. [type:Integer]
#' @param slope The method assumes all tielines for a given ATPS are parallel, thus only one slope is required. [type:double]
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
# ' @param maxiter	- A positive integer specifying the maximum number of iterations allowed.
#'
#' @examples
#' \dontrun{
#' AQSysDB("C:/data.xls")
#'}
#'
#' @references
#' KAUL, A. The Phase Diagram. In: HATTI-KAUL, R. (Ed.). Aqueous Two-Phase Systems: Methods and Protocols: Humana Press, v.11, 2000. cap. 2, p.11-21.  (Methods in Biotechnology). ISBN 978-0-89603-541-6.
#' (\href{https://link.springer.com/10.1385/1-59259-028-4:11}{SpringerLink})
#'
AQSysAnima <- function(dataSET,
                      db = LLSR::llsr_data,
                      xmax = NULL,
                      NP = 100,
                      slope = NULL,
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
                      silent = TRUE) {
  #
  if (is.null(slope)) {
    slope = findSlope(db, dataSET)
  } else if (!((ncol(dataSET) / 2) == length(slope))) {
    AQSys.err("11")
  }
  # Select which model will be used to generate the plot. Function return list of plots and respective number of parameters
  models_npars <- AQSysList(TRUE)
  # Divides the number of columns of the data set per two to calculate how many systems it hold. result must be even.
  nSys <- (ncol(dataSET) / 2)
  SysNames <- FALSE
  image_points <- list()
  # Check data.frame validity and make an array of names for the systems if none is provided
  if (ncol(dataSET) == 2) {
    PIC <- 1
    seriesNames <- "PLOT"
    SysList <- list()  # List which will stack and return all data
    PlotList <- list() # List which will stack and return all plots
    
    #for (i in seq(1, nSys)) {
    #
    RawData <- dataSET[, 1:2]
    SysData <- LLSRxy(na.exclude(RawData[6:nrow(RawData), 1]),
                      na.exclude(RawData[6:nrow(RawData), 2]),
                      RawData[4, 2])
    # Analyse data and return parameters
    model_pars <-
      summary(AQSys(SysData,  modelName))$parameters[, 1]
    # Select Model based on the user choice or standard value
    Fn <- ifelse(modelName %in% names(models_npars),
                 AQSys.mathDesc(modelName),
                 AQSys.err("0"))
    # define a straight line EQUATION
    Gn <- function (yMin, AngCoeff, xMAX, x) {
      yMin + AngCoeff * (x - xMAX)
    }
    # Add constant variable values to the equations
    modelFn <- function(x) Fn(model_pars, x)
    modelTl <- function(x) Gn(yMin, slope[1], xMAX, x)
    # Decide whether xmax will be calculated or use the provided value
    xMAX <- ifelse(is.null(xmax), max(SysData[, seq(1, ncol(SysData), 2)] * 1.1), xmax)
    yMAX <- max(SysData[, 2]) # get ymax from the dataset
    SysList[[1]] <- list()
    SysL <- list()
    # go from xmin/1.5 to xmax in steps of 0.15
    x <- seq(min(SysData[, 1]) / 1.5, xMAX, ((xMAX - (min(SysData[, 1]) / 1.5)) / NP))
    BNDL <- setNames(data.frame(x, modelFn(x)), c("X", "Y"))
    BNDL["System"] <- seriesNames
    #
    TL <- 0
    dt <- 1
    reset_BNDL <- BNDL
    CrptFnd <- FALSE # Crital Point Found Logical variable
    DivFactor <- 5
    # prepare a plot with the system curve and all tielines. Convergence lines are included if requested.
    output_plot <- bndOrthPlot(subset(BNDL, BNDL$System == seriesNames), xlbl, ylbl)
    output_plot <- output_plot + annotate(
      "point",
      x = min(SysData[, 1]),
      y = yMAX,
      colour = "royalblue",
      shape = 17,
      size = 4
    )
    tl_x_vec <- NULL
    tl_y_vec <- NULL
    XY_limits <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("x_min", "y_max", "x_max", "y_min"))
    # yMAXTL <- yMAX + 1
    FOUND <- FALSE
    while ((dt > tol) && !CrptFnd) {
      TL <- TL + 1
      yMAXTL <- yMAX + 1
      # The following lines tests successively different values of xmax until a valid ymax,
      # which must be smaller than the experimental ymax.
      PT = 0
      while (yMAXTL > yMAX) {
        PT <- PT + 1
        yMin <- Fn(model_pars, xMAX) # EVALUATE REPLACE XMAX TO THE LIMIT OF SOLUBILITY?
        xRoots <- uniroot.all(function(x) (modelFn(x) - modelTl(x)), c(0, xMAX), tol = 0.1)
        xMin <- min(xRoots) # As we started from the biggest root, we select the smallest
        yMAXTL <- modelFn(xMin) # calculate the y-value for the root found
        if (yMAXTL > yMAX) {
          # compare to the highest value allowed (it must be smaller than the maximum experimental Y)
          xMAX <- xMAX - 0.005 # if bigger than physically possible, decrease xmax and recalculate
        } else {
          FOUND <- TRUE
          image_points$xMax <- ifelse(is.null(image_points$xMax), xMAX, image_points$xMax)
          image_points$yMin <- ifelse(is.null(image_points$yMin), yMin, image_points$yMin)
        }
        if (((PT %% 50) == 0) | FOUND){
          XY_limits <- rbind(XY_limits, setNames(c(xMin, yMAXTL, xMAX, yMin), c("x_min", "y_max", "x_max", "y_min")))
          if (is.null(names(image_points))){
            image_points <- list(XYLimits = XY_limits)
          } else {
            image_points$XYLimits <- XY_limits
          }
          #
          #
          #
          # PIC <- PIC + 1
          # saveImgSeq(output_plot, image_points, PIC, FOUND)
          #
          #
          #
        }
        #
      }
      # create and name a data.frame containing coordinates for the tieline calculated above
      SysL[[TL]] <- setNames(data.frame(c(xMAX, (xMAX + xMin) / 2, xMin), c(yMin, (yMAXTL + yMin) / 2, yMAXTL)), c("X", "Y"))
      SysL[[TL]]["System"] <- paste(seriesNames, "TL", TL, sep = "-")
      # check if compositions of both phases, as well as the global composition, are equal. If so, critical point was found.
      CrptFnd <- is.equal(SysL[[TL]], tol)
      # Bind the calculated tieline's data.frame to the output data.frame variable
      BNDL <- rbind(BNDL,  SysL[[TL]])
      # A monod-base equation to help convergence -
      dt <- abs(SysL[[TL]][2, 1] - modelFn(SysL[[TL]][2, 1])) / ((5 * TL) / (xMAX / DivFactor))
      xMAX <- xMAX - dt
      #
      if (TL > 1500) {
        tl_x_vec <- NULL
        tl_y_vec <- NULL
        DivFactor <- DivFactor * 0.5
        SysL <- list()
        TL <- 0
        dt <- 1
        BNDL <- reset_BNDL
        xMAX <- ifelse(is.null(xmax), max(SysData[, seq(1, ncol(SysData), 2)] * 1.1), xmax)
        yMAX <- max(SysData[, 2])
        cat("\t - Convergence not achieved. Reseting Search Factors...\n")
      }
      #
      image_points$tielines <- BNDL
      image_points$XYLimits <- XY_limits
      #
      PIC <- PIC + 1
      if (PIC >= 299){
        saveImgSeq(output_plot, image_points, PIC, FOUND)
      }
      
    }
    #
    TLLs <- SysL
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
    XC <- SysL[[length(SysL) - 1]][["X"]][1]
    YC <- SysL[[length(SysL) - 1]][["Y"]][1]
    SysCvP <- setNames(data.frame(XC, YC), c("XC", "YC"))
    # Max Tieline will be the first 'viable' tieline, i.e., which xmax yield a ymax within the experimental range
    maxTL <- SysL[[1]]
    minTL <- FindMinTL(SysCvP, GlobalPoints[1,], max(maxTL["X"]), slope[1], modelFn, tol)
    # execute saving procedures
    image_points$tielines <- BNDL
    image_points$XYLimits <- XY_limits
    image_points$boundaries <- TLL(minTL, maxTL)
    image_points$XYLimits <- XY_limits
    image_points$MAXTL <- maxTL
    image_points$MINTL <- minTL
    image_points$CriticalPoint <- SysCvP
    #
    PIC <- PIC + 1
    saveImgSeq(output_plot, image_points, PIC, FOUND)
  } else{
    # Return an error if an invalid dataset is provided.
    AQSys.err(9)
  }
}


saveImgSeq <- function(image_data, image_points, index, FOUND = FALSE) {
    # BKPlot <- image_data
    XYLimits <- image_points$XYLimits
    tielines <- image_points$tielines
    boundaries <- image_points$boundaries
    maxTL <- image_points$MAXTL
    minTL <- image_points$MINTL
    CriticalPoint <- image_points$CriticalPoint
    #
    # if (!("tielines" %in% names(image_points)) & !FOUND) {
      image_data <- image_data + geom_point(
        data = setNames(data.frame(
          c(XYLimits[, 1], XYLimits[, 3]),
          c(XYLimits[, 2], XYLimits[, 4])
        ), c("X", "Y")),
        aes_string(x = "X", y = "Y"),
        colour = "red",
        size = 2,
        alpha = 1
      )
    #} else
      if (FOUND) {
      # if (!("tielines" %in% names(image_points))) {
      image_data <- image_data + annotate(
        "point",
        x = as.numeric(image_points$xMax),
        y = as.numeric(image_points$yMin),
        colour = "royalblue",
        shape = 17,
        size = 4
      ) + geom_point(
        data = setNames(data.frame(
          c(XYLimits[, 1], XYLimits[, 3]),
          c(XYLimits[, 2], XYLimits[, 4])
        ), c("X", "Y")),
        aes_string(x = "X", y = "Y"),
        colour = "red",
        size = 2,
        alpha = 1
      )
    }
    #
    if ("tielines" %in% names(image_points)){
      image_data <- image_data + 
        geom_line(
          data = subset(tielines, tielines$System != "PLOT"),
          aes_string(x = "X", y = "Y", group = "System"),
          colour = "red",
          alpha = 0.4
        ) +
        geom_point(
          data = subset(tielines, tielines$System != "PLOT"),
          aes_string(x = "X", y = "Y"),
          colour = "red",
          size = 2,
          alpha = 1
        ) 
    }
    if ("CriticalPoint" %in% names(image_points)) {
      image_data <- image_data  +
        # geom_point(
        #   data = boundaries$minTL,
        #   aes_string(x = "X", y = "Y"),
        #   colour = "yellow",
        #   size = 2,
        #   alpha = 0.5
        # ) +
        # geom_line(
        #   data = boundaries$minTL,
        #   aes_string(x = "X", y = "Y"),
        #   colour = "yellow",
        #   size = 1, 
        #   alpha = 0.5
        # ) +
        # geom_point(
        #   data = maxTL,
        #   aes_string(x = "X", y = "Y"),
        #   colour = "yellow",
        #   size = 2,
        #   alpha = 0.5
        # ) +
        # geom_line(
        #   data = maxTL,
        #   aes_string(x = "X", y = "Y"),
        #   colour = "yellow",
        #   size = 1,
        #   alpha = 0.5
        # )  + 
        annotate(
          "point",
          x = CriticalPoint$XC,
          y = CriticalPoint$YC,
          colour = "black",
          bg = "gold",
          shape = 23,
          size = 2
        )
    }
    #
    wdir <- getwd()
    wdir <- file.path(wdir, "imgs")
    # do.call(file.remove, list(list.files(wdir, full.names = TRUE)))
    # file.remove(list.files(pattern=".png"))
    dir_and_file <- paste(wdir, paste("plot_(", index, ")",".png", sep = ""), sep = .Platform$file.sep)
    ggsave(
      filename = dir_and_file,
      plot = image_data,
      width = 21.14 / 2,
      height = 14.39 / 2
    )
    return(image_data)
  }

