options(digits = 14)
####################################################################################################################
#' @import rootSolve minpack.lm
#
merchuk <- function(XYdt, P1 = 10, P2 = 1, P3 = 0, ...) {
    # Output initialization
    FFn <- NULL
    # define header
    names(XYdt) <- c("XC", "YC")
    # define nonlinear system and solve it
    suppressWarnings(
      FFn <- tryCatch({
        nlsLM(
          YC ~ P1 * exp(P2 * (XC ^ (0.5)) - P3 * (XC ^ 3)),
          start = list(P1 = P1, P2 = P2, P3 = P3),
          data = XYdt,
          na.action = na.exclude
        )
      },
      error = function(e) {
        return(NULL)
      })
    )
    # return output from anaysis
    FFn
  }

murugesan <- function(XYdt, ...) {
  # Output initialization
  FFn <- NULL
  # define header
  names(XYdt) <- c("XC", "YC")
  ##  controll<-nls.control(maxiter=50,
  ##  tol=1e-10, minFactor = 1/1024,
  ##  printEval = FALSE, warnOnly = FALSE)
  # define nonlinear system and solve it
  suppressWarnings(
    FFn <- tryCatch({
      nls(
        YC ~ P1 + P2 * (XC) ^ 0.5 + P3 * XC,
        start = list(P1 = 50, P2 = 1, P3 = 0),
        data = XYdt,
        na.exclude
      )
    },
    error = function(e) {
      return(NULL)
    })
  )
  # return output from anaysis
  FFn
}

tello <- function(XYdt, ...) {
  # tello's method is highly dependent of guess values to obtain its parameters
  # so the method for calculating a more approximated guess value is described in
  # his article and implemented below
  #
  # Output initialization
  FFn <- NULL
  # calculate derivative of data
  tryCatch({
    df.sys <- diff(XYdt[, 1]) / diff(XYdt[, 2])
    nrow.sys <- nrow(XYdt)
    S <- smooth.spline(df.sys)$fit$coef[1:nrow.sys]
    XC <- XYdt[, 1]
    # merge data in a dataframe
    coef.est <- LLSRxy(XC, S)
    # define header
    names(coef.est) <- c("XC", "S")
    # define header
    names(XYdt) <- c("XC", "YC")
    # define nonlinear system and solve it using estimated parameters as guess
    suppressWarnings(
      FFn <- nlsLM(
        XC ~ exp((YC - P3) / P1) - P2,
        #YC ~ P1*log(XC + P2) + P3,
        start = list(P1 = -100, P2 = 10, P3 = 0),
        data = XYdt,
        na.action = na.exclude,
        control = nls.lm.control(maxiter = 25)
      )
    )
    # return output from anaysis
  },
  error = function(e) {
    return(NULL)
  })
  # return output from anaysis
  FFn
}

tang <- function(XYdt, P1 = 10, P2 = 0, P3 = 0, P4 = 0, ...) {
  # Output initialization
  FFn <- NULL
  # define header
  names(XYdt) <- c("XC", "YC")
  # define nonlinear system and solve it
  suppressWarnings(
    FFn <- tryCatch({
      nls(
        YC ~ exp(P1 + P2 * (XC ^ (0.5)) + P3 * XC + P4 * (XC ^ 2)),
        start = list(
          P1 = P1,
          P2 = P2,
          P3 = P3,
          P4 = P4
        ),
        data = XYdt,
        na.action = na.exclude
      )
    },
    error = function(e) {
      return(NULL)
    })
  )
  # return output from anaysis
  FFn
}
# 
chen <- function(XYdt, P1 = 10, P2 = 0, P3 = 0, P4 = 0,P5 = 0, ...) {
  # Output initialization
  FFn <- NULL
  # define header
  names(XYdt) <- c("XC", "YC")
  # define nonlinear system and solve it
  suppressWarnings(
    FFn <- tryCatch({
      nls(
        YC ~ exp(P1 + P2 * (XC ^ (0.5)) + P3 * XC + P4 * (XC ^ 2) + P5 * (XC ^ 3)),
        start = list(
          P1 = P1,
          P2 = P2,
          P3 = P3,
          P4 = P4,
          P5 = P5
        ),
        data = XYdt,
        na.action = na.exclude
      )
    },
    error = function(e) {
      return(NULL)
    })
  )
  # return output from anaysis
  FFn
}

# 
# xie <- function(XYdt, ...) {
#   # Output initialization
#   FFn <- NULL
#   # define header
#   names(XYdt) <- c("XC", "YC")
#   # define nonlinear system and solve it
#   FFn <- tryCatch({
#     nlsLM(
#       YC ~ exp(-(XC/P2) - (XC/P4))*(exp(XC/P4)*P1 + exp(XC/P2)*P3 + exp(XC/P2 + XC/P4)*P5), 
#       start = c(
#         P1 = 3,
#         P2 = 5,
#         P3 = 7,
#         P4 = 9,
#         P5 = 2
#       ),
#       data = XYdt,
#       na.action = na.exclude,
#       control = nls.lm.control(maxiter = 25)
#       )
#   },
#   error = function(e) {
#     return(NULL)
#   })
#   # return output from anaysis
#   FFn
# }
# 
# xueqiao <- function(XYdt, P1 = 0, P2 = 1, P3 = 1, P4 = 1, ...) {
#   # Output initialization
#   FFn <- NULL
#   # define header
#   names(XYdt) <- c("XC", "YC")
#   # define nonlinear system and solve it
#   FFn <- tryCatch({
#     nls(
#       YC ~ P1 * (-log(P2 * (XC / P3) + P4) / P2),
#       start = list(
#         P1 = P1,
#         P2 = P2,
#         P3 = P3,
#         P4 = P4,
#         P5 = P5
#       ),
#       data = XYdt,
#       na.action = na.exclude
#     )
#   },
#   error = function(e) {
#     return(NULL)
#   })
#   # return output from anaysis
#   FFn
# }
