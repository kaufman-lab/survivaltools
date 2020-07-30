#'Convert Dates/integers to the end of a specified interval
#'
#'\code{bin_time} codes values in x, according to which of the specified intervals they falls into,
#'to the end (or start or midpoint) of that interval
#'
#'Similar to cut but only for integers.
#' Unlike cut it returns the same class as the input (whereas cut returns factor representing the interval).
#' Easier to use than cut when dealing with a long sequence of intervals.
#' Note that if objects of class Date have decimal portions these will be silently ignored.
#'
#' @param x A vector of Date, IDates, or integers
#' @param y A two-column data.frame or data.table where the start column is the start Date, IDate, or integer
#'  of each period and the second column is the ending value (also Date, IDate, or integer) of each period.
#'  The intervals must be non-overlapping and have no gaps.
#' @param snapto a length-1 character vector equalling "start", "end", "mid", or "count". This designates
#'  whether the resulting vector of dates will be the start date, the end date, the midpoint (rounded down)
#'  of each period, or an integer ordered bin count.
#' @param nonmatcherror When TRUE, an error will occur if there are values of x not in intervals of y
#'  (ie outside the span of y). When FALSE, elements of the return vector corresponding to
#'  values of x not in intervals of y are returned as NA.
#' @return A vector of the same class as columns of y corresponding to elements of x
#' @examples
#'  bin_time(x=c(1L,6L,2L,14L),y=data.frame(seq(1L,21L,by=5L),seq(5L,25L,by=5L)))
#' @import data.table
#' @export
bin_time <-function(x,y,snapto="end",nonmatcherror=TRUE){

  stopifnot(class(x)[1]%in%c("integer","Date","IDate"))
  stopifnot(class(y[[1]])[1]%in%c("integer","Date","IDate"))
  stopifnot(class(y[[2]])[1]%in%c("integer","Date","IDate"))
  stopifnot(class(y[[1]])[1]==class(y[[1]])[1])

  stopifnot(length(y)==2)
  stopifnot(class(y)[1]%in%c("data.table","data.frame"))
  stopifnot(snapto %in% c("start","end","mid","count") )

  y <- as.data.table(y)


  if(nrow(y)!=nrow(stats::na.omit(y))){
    stop("intervals in y cannot be missing")
  }

  setnames(y,c("y1","y2"))
  setkey(y,y1) #order bins

  y[, iny:=TRUE]

 if(snapto=="start"){
   y[, out_value:=y1]
 }

  if(snapto=="end"){
    y[, out_value:=y2]
  }
  if(snapto=="mid"){
    y[, out_value:=floor((y2-y1)/2)]
  }

  if(snapto=="count"){
    y[, out_value:=1L:.N]
  }

  setkey(y,y1,y2)
  if(nrow(foverlaps(y,y))!=nrow(y)){
    stop("intervals in y must be non-overlapping")
  }

  if(!y[, list(V1=shift(y1,type="lead")-1,y2=y2)][, stats::na.omit(.SD)][,
      identical(as.numeric(V1),as.numeric(y2))]){
    stop("intervals in y must be comprehensive")
  }

  #check for overlaps in y

  xDT <- data.table(y1=x,y2=x)
  out <- foverlaps(xDT,y)


  if(any(is.na(out$iny)&nonmatcherror)){
    stop("there are values of x that do not fall in intervals of y")
  }
 out$out_value
}

