

#' prep function
#'
#' @param x is a data.table with the following columns:
#' - id: a unique identifier for participant
#' - age: age (in years, rounded down ) at baseline.
#' - time: days since baseline
#' - event: 0/1
#' @export

prep <- function(x){

  x <- as.data.table(x)

  if(uniqueN(x$id)!=length(x$id)){
    stop("duplicates detected in id. id must be unique")
  }

  any_decimal_age <- any(sapply(x$age%%1,function(z){isTRUE(all.equal(z,0))}))
  if(any_decimal_age){
    stop("decimals detected in age column. round these to integer years first")
    x[,age:=round(age)]
  }
  if(max(x$age,na.rm=TRUE)>120){
    stop("Ages above 120 detected. age column should be in units of years.")
  }
  if(length(setdiff(unique(x$event),c(0,1)))!=0){
    stop("values other than 0 or 1 detected in event column")
  }

  if(max(time,na.rm=TRUE)<50){
    warning("maximum time value is less than 50.
            Are you sure this column is in units of days since follow-up?
            Incorrect results will be returned if this column is not in units of days")
  }

  any_decimal_time <- any(sapply(x$time%%1,function(z){isTRUE(all.equal(z,0))}))
  if(any_decimal_time){
    stop("decimals detected in time column.
            Are you sure this column is in units of days since follow-up?
            Incorrect results will be returned if this column is not in units of days")
  }



  class()

  x[]

}

