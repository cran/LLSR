####################################################################################################################
AQSys.mathDesc <- function(mathDesc){
  switch(mathDesc,
         "merchuk"={
           Fn<-function(CoefSET,XC){
             P1<-CoefSET[1]
             P2<-CoefSET[2]
             P3<-CoefSET[3]
             P1*exp(P2*(XC^(0.5))-P3*(XC^3))
           }
         },
         "murugesan"={
           Fn <- function(CoefSET,XC){
             P1<-CoefSET[1]
             P2<-CoefSET[2]
             P3<-CoefSET[3]
             P1+P2*(XC)^0.5+P3*XC
           }
         },
         "tello"={
           Fn <- function(CoefSET,XC){
             P1<-CoefSET[1]
             P2<-CoefSET[2]
             P3<-CoefSET[3]
             P1*log(P2+XC)+P3
           }
         },
         AQSys.err("0")
  )
  return(Fn)
}
####################################################################################################################
#'@rdname AQSysList
#'@export AQSysList
#'@title Aqueous Systems Descriptors already implemented
#'@description The function returns a list of all mathematical descriptors available at the time.
AQSysList <- function(){
  AQSys.List()
}
####################################################################################################################
AQSys.List <- function(){
  updte <- c("merchuk", "murugesan", "tello")
  return(updte)
}
