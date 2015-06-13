####################################################################################################################
undigest <- function(db, str.md5){
  ans <- db$db.ref[grep(str.md5, db$db.ref[,2], fixed = TRUE), 1]
  #ans <- grep(str.md5, db$db.ref[,2], fixed = TRUE)
  ans
}
####################################################################################################################
db.check <- function(db){
  if (is.null(db)){
    AQSys.err("4")
  } else {
    if ( (is.list(db)) & (length(db) == 3) ){
      for (i in 1:length(db)){
        if (is.data.frame(db[[i]]) == FALSE) AQSys.err("3", k=names(db)[i])
      }
    } else { AQSys.err("7") }
  }
}
####################################################################################################################
