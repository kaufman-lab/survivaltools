

#' Melt multiple type columns to long format
#' Reshape to long a data.table with survival data for multiple event types.
#' @param x A data.table
#' @param id A length-1 character vector corresponding to a column of x that is a unique identifier.
#' @param measure.vars A list where each element is a length-n character vector corresponding to column names in x.
#'   The character vectors must be in the same order.  Each element of the list represents a set of
#'   columns that will be turned into a single column in the return.
#'  (such that measure.vars[[1]][1] corresponds to measure.vars[[2]][1],measure.vars[[1]][2]
#'  corresponds to measure.vars[[2]][2],etc).
#'  Typically, one element of measure.vars is a character vector corresponding to column names of follow-up dates/days and
#'  and another  element of measure.vars is a character vector corresponding to column names of binary event indicators.
#' @param value.name A character vector the length of measure.vars designating desired the column names
#'  in the return. Each element of value.name corresponds to the column that will be created from
#'  each element of measure.vars
#' @param type A vector the length of each element of measure.vars denoting the type of event.
#'  ie, type[1] denotes the type of event of measure.vars[[1]][1],measure.vars[[2]][1], etc and
#'  type[2] denotes the type of event of measure.vars[[1]][2],measure.vars[[2]][2], etc.
#'  type will become an indicator column in the return
#' @return
#' If there's only a single baseline date for all event types, ignore it for the purpose of this function
#' then merge it back in to the result. The same goes for any other variables in x that you want to be in y
#' that don't need to be melted.
#' @examples
#' x <- data.table(id=c(1,2,3),baseline_date=c(10,5,3),
#' mi_fu_date=c(15,8,20),mi_indicator=c(1,1,0),death_fu_date=c(26,8,20),death_indicator=c(0,1,0))
#' x_long <- melt_events(x,id="id",
#'  measure.vars=list(c("mi_fu_date","death_fu_date"),c("mi_indicator","death_indicator")),
#'  value.name=c("fu_date","event_indicator"),
#'  type=c("mi","death")
#' )
#' x_long <- x[,list(id,baseline_date)][x_long,on="id"] #merge baselin_date back in
#'
#' @export

melt_events <-function(x,id,measure.vars, value.name,type){
  if(length(unique(x[[id]]))!=nrow(x)){stop("duplicate values detected in id")}
  if(length(measure.vars)!=length(value.name)){
    stop("value.name must be the same length as measure.vars")
  }
  if(!all(lapply(measure.vars,length)==length(measure.vars[[1]]))){
    stop("Each element of measure.vars must be the same length")
  }
  if(length(type)!=length(measure.vars[[1]])){
    stop("type must be the same length as each element of measure.vars")
  }

  if(id=="...variable"){
  stop("id cannot be named ...variable. this is a reserved name internal to this funtion. rename this column")
    }

  if(id=="type"){
    stop("id cannot be named type this is a reserved name internal to this funtion. rename this column")
  }


  x_long <- melt(data=x,
                 id.vars=id,
                 measure.vars = measure.vars,
                 value.name = value.name,
                 variable.name = "...variable")

  key <- data.table(...variable=as.factor(as.character(1:length(measure.vars[[1]]))),
                    type=type)
  x_long[key,type:=type,on="...variable"]
  #x_long <- key[x_long,on="...variable"]
  x_long[, ...variable:=NULL]
  x_long[]
}
