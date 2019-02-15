####################################################################################################################
mathDescPair <- function(modelName) {
  #
  switch(modelName,
    "merchuk" = {
      FN <- "PARs[1] * exp(PARs[2] * x[$2]^0.5 - PARs[3]*x[$2]^3) - x[$1]"
    },
    "murugesan" = {
      FN <- "PARs[1] + PARs[2] * (x[$2]) ^ 0.5 + PARs[3] * x[$2] - x[$1]"
    },
    "tello" = {
      FN <- "PARs[1] * log(PARs[2] + x[$2]) + PARs[3] - x[$1]"
    },
    "tang" = {
      FN <- "exp(PARs[1] + PARs[2] * (x[$2] ^ (0.5)) + PARs[3] * x[$2] + PARs[4] * (x[$2] ^ 2)) - x[$1]"
    },
    "chen" = {
      FN <- "exp(PARs[1] + PARs[2] * (x[$2] ^ (0.5)) + PARs[3] * x[$2] + PARs[4] * (x[$2] ^ 2)+ PARs[5] * (x[$2] ^ 3)) - x[$1]"
    },
    AQSys.err("0")
  )
  # return chosen functions
  return(FN)
}
####################################################################################################################
AQSys.mathDesc <- function(mathDesc) {
  # Each switch option provides an equation that will be available to be used to
  # make plots and calculate compositions for a system with known parameters
  switch(
    mathDesc,
    "merchuk" = {
      Fn <- function(CoefSET, XC) {
        # equation's parameters
        P1 <- CoefSET[1]
        P2 <- CoefSET[2]
        P3 <- CoefSET[3]
        # merchuk's equation
        P1 * exp(P2 * (XC ^ (0.5)) - P3 * (XC ^ 3))
      }
    },
    "murugesan" = {
      Fn <- function(CoefSET, XC) {
        # equation's parameters
        P1 <- CoefSET[1]
        P2 <- CoefSET[2]
        P3 <- CoefSET[3]
        # murugesan's equation
        P1 + P2 * (XC) ^ 0.5 + P3 * XC
      }
    },
    "tello" = {
      Fn <- function(CoefSET, XC) {
        # equation's parameters
        P1 <- CoefSET[1]
        P2 <- CoefSET[2]
        P3 <- CoefSET[3]
        # tello's equation
        P1 * log(P2 + XC) + P3
      }
    },
    "tang" = {
      Fn <- function(CoefSET, XC) {
        # equation's parameters
        P1 <- CoefSET[1]
        P2 <- CoefSET[2]
        P3 <- CoefSET[3]
        P4 <- CoefSET[4]
        # tang's equation
        exp(P1 + P2 * (XC ^ (0.5)) + P3 * XC + P4 * (XC ^ 2))
      }
    },
    "chen" = {
      Fn <- function(CoefSET, XC) {
        # equation's parameters
        P1 <- CoefSET[1]
        P2 <- CoefSET[2]
        P3 <- CoefSET[3]
        P4 <- CoefSET[4]
        P5 <- CoefSET[5]
        # chen's equation
        exp(P1 + P2 * (XC ^ (0.5)) + P3 * XC + P4 * (XC ^ 2)+ P5 * (XC ^ 3))
      }
    },
    # "xie" = {
    #   Fn <- function(CoefSET, XC) {
    #     # equation's parameters
    #     P1 <- CoefSET[1]
    #     P2 <- CoefSET[2]
    #     P3 <- CoefSET[3]
    #     P4 <- CoefSET[4]
    #     P5 <- CoefSET[5]
    #     # xie's equation
    #     P1 * exp(- (Xc/P2)) + P3 * exp(- (Xc/P4)) + P5
    #   }
    # },
    # model needs more data than just the phase diagram
    # "xueqiao" = {
    #   Fn <- function(CoefSET, XC) {
    #     # equation's parameters
    #     P1 <- CoefSET[1]
    #     P2 <- CoefSET[2]
    #     P3 <- CoefSET[3]
    #     P4 <- CoefSET[4]
    #     # xueqiao's equation
    #     P1 * (-log(P2 * (XC / P3) + P4) / P2)
    #   }
    # },
    # if user selects an option not available, it triggers an error
    # (check AQSys.err.R for details)
    AQSys.err("0")
  )
  # return chosen function
  return(Fn)
}
####################################################################################################################
#' @rdname AQSysList
#' @export
#' @title Aqueous Systems Descriptors already implemented
#' @description The function returns a list of all mathematical descriptors available at the time.
#' @param npars Logic option to return a List variable containing the number of required parameters for each equation.
AQSysList <- function(npars=FALSE) {
  # a new entry in updte must be added for each new equation implemmented in AQSys.mathDesc()
  # updte entries' name must match AQSys.mathDesc switch options
  model_list <- list(
    "merchuk" = 3,
    "murugesan" = 3,
    "tello" = 3,
    "tang" = 4,
    "chen" = 5
    # ,"xie" = 5,
    # "xueqiao" = 4
  )
  updte <- if(npars) {
    model_list
  } else {
    names(model_list)
  }
  # return list
  return(updte)
}
####################################################################################################################
