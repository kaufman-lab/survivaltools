
#' Split all follow-up intervals at each occurence of an event
#'
#' Take a set of follow-up times represented as start/stop intervals and split each interval
#' at every occurance of an event so that time-varying variables can be merged in.
#'
#' The returned format where each participant has multiple rows each corresponding to
#' a non-overlapping interval of follow-up is sometimes known as the Andersen Gill formulation
#' of the Cox model. While the term "Andersen Gill" is more often used in the context of recurring events
#' this function was written not to deal with recurrent events but rather to facilitate merging
#' time-varying data into the analysis dataset more flexibly than allowed by the survival package.
#' As such, this funciton is untested for recurrent events so I've written stop() if it detects
#' recurrent events. It may be possible to modify the code to deal with recurrent events and time-varying
#' data, but this is not something I need so I've not considered that possiblity.
#'
#' Rather, the focus here is to deal with time-varying variables. The functions in the survival
#' package for dealing with time-varying variables are designed to deal with time-varying covariates
#' which change at a limited set of discrete time-points. This is the sort of data that's typical
#' in a prospective cohort study with multiple visits (ie, you collect blood pressure at the baseline visit
#' and at a small number of specific subsequent dates). It seems like the survival function was not designed
#' to deal with time-varying data which is changing constantly (as in the case of air pollution exposure).
#' Even when dealing with long-term averages, say a 2 year average prior to a given day, that running
#' 2 year average changes every day.  In my testing, attempting to provide the function in the survival package that
#' deals with time-varying data a dataset containing every two-year average (corresponding to every day of follow-up)
#' for every participant in a cohort leads to a crash. My guess is that the function is attempting
#' to split every individual's follow-up time at every day resulting in a large dataset that leads to
#' an out-of-memory error crash.
#'
#' However, even if exposure changes at every discrete time point in a survival dataset for every participant,
#' cutting everyone's follow-up time at every time point is unnecessary. For the purpose of accurate
#' cox estimation, it is only strictly necessary to cut follow-up intervals at the occurance of an event--
#' any event of the same type by any indivual in the cohort.
#'
#' (Actually, technically,
#'  if you're stratifying the cox model by a categorical covariate, it's probably only necessary to
#'  split the follow-up times of individuals only by the dates of events of other individuals in the same strata.
#'  but due to the fact that this would require a different analytic formulation of the dataset for different
#'  models stratified by different variables, this approach probably isn't worth it unless you encounter
#'  extreme computational problems.)
#'
#'  In any case, I've written this function to break up follow-up intervals at the date of every
#'  other event in the dataset.
#'
#'  Additionally, since when we're doing events analyses we're often considering multiple different types of
#'  events (each of which has their own set of event dates), this function has been written to simultaneously
#'  handle multiple event types in one function call via specification of a categorical variable
#'  (defining event types) in the group_vars argument. (Come to think of it, you could probably also put a
#'  categorical variable defining strata as discussed above in the group_vars variable.
#'  But you'd be careful in that you'd always need to stratify by that variable when fitting the cox model)
#'
#'
#'
#' @param x is a data.table
#' @param interval_vars length 2 character vector designating two column names in x:
#'  the first element of interval_vars designates the start time for follow-up (if using study time, create a column of zeros)
#'  the second element of interval_vars designates the start time for follow-up (ie censor or event date)
#' @param group_vars optionally,  groups_vars designates column name(s) in x in which the splitting algorithm will be repeated separately.
#'  (e.g. if you have different outcomes you can do this computation separately for each outcome with a single function
#'  call assuming that your dataset is formatted "long" with respect to outcome type (each particpant has one row per outcome type))
#' @param id_var a character vector designating a column name in x that identifies a row identifier that will be included in the final output
#' @param event_indicator character vector designating a column name in x that is of 0/1 indicator variable for events
#' @return returns an unkeyed data.table with columns of x plus two columns: start and end. the event_indicator
#' column will be equal to 1 only in a ppt's final interval (and only if that ppt had an event)
#' @export


split_by_event_times <- function(x,interval_vars,id_var,event_indicator,group_vars=NULL){

  if(!is.data.table(x)){stop("x must be a data.table. use setDT(x)")}

  if(any(c("actual_event_date","actual_event_date2","value")%in% c(group_vars,id_var,interval_vars) )){
    stop("actual_event_date/actual_event_date2/value are reserved names. rename your columns in x
         to avoid naming conflicts with temporary variables internal to this function")
  }

  if(any(c("start","end")%in% names(x) )){
    stop("start, end are reserved names. rename your columns in x to avoid naming conflicts
         with start and end which will be added to the return")
  }


  if(x[, sum(is.na(.SD)),.SDcols=c(interval_vars,id_var,event_indicator,group_vars)]!=0){
    stop("columns in x designated by arguments interval_vars,id_var,event_indicator,group_vars
         must all be nonmissing")
  }

  if(!x[,all(.SD[[interval_vars[2]]]>.SD[[interval_vars[1]]])]){
    stop("all interval ends must occur after interval starts")
  }

  if(!eval(parse(text=paste0("!x[,list(V1=sum(",event_indicator,")),by=c(id_var,group_vars)][,all(V1==1)]")))   ){
    stop("recurrent events within a ppt detected. This function was not tested for recurrent events")
  }

  if(!x[,.N,by=c(id_var,group_vars)][,all(N==1)]){
    stop("multiple rows per ppt detected. This function was not tested to deal with a dataset
         where particpants follow-up time has already been split into multiple rows.")
  }

  if(sum(x[[event_indicator]])==0){
    #there aren't any events in this dataset!
    out <- copy(x)
    out[, `:=`(start=.SD[[1]],end=.SD[[2]]),.SDcols=interval_vars]
    setcolorder(out, c(group_vars,id_var,"start","end",event_indicator))
    return(out)
  }


  #get unique event dates (excluding censor dates)
  #subset to events where event_indicator is equal to 1 and take unique values
  event_dates <- unique(x[J(1L), on=c(event_indicator),c(group_vars,interval_vars[2]),with=FALSE])
  setnames(event_dates,interval_vars[2],"actual_event_date")

  stopifnot(nrow(event_dates[is.na(actual_event_date)])==0)
  #calling this actual event date because these are all real events versus the event_date column in events_raw_long2 includes
  #dates which are not events and just the last day of follow-up
  event_dates[, actual_event_date2:=actual_event_date]


  setkeyv(event_dates,c(group_vars,"actual_event_date","actual_event_date2"))
  q <- foverlaps(x[,c(id_var, group_vars,interval_vars),with=FALSE],
                 event_dates,by.x=c(group_vars, interval_vars))

  #events[!q] are rows in events that are not in q. so this check ensures that every row in events ended up in q
  stopifnot(nrow(x[!q,on=c(group_vars,id_var,interval_vars)])==0)

  q[,actual_event_date2:=NULL]

  #error checks
  q[,stopifnot(all(!is.na(.SD[[interval_vars[1]]] )))] #there should be no missings in the baseline_date or the event_date. only actual_event can be missing
  q[,stopifnot(all(!is.na(.SD[[interval_vars[2]]] )))]
  q[,
    stopifnot(
      max(.SD[[interval_vars[1]]],.SD[[interval_vars[2]]],actual_event_date,na.rm=TRUE) == .SD[[interval_vars[2]]]
    ),by=c(group_vars,id_var)]

  q[,
    stopifnot(
      min(.SD[[interval_vars[1]]],.SD[[interval_vars[2]]],actual_event_date,na.rm=TRUE)==.SD[[interval_vars[1]]]
    ),by=c(group_vars,id_var)]


  #take a subject's baseline date, event_date, and any actual event dates that occur between those two,
  #stack them on top of each other, and take unique rows of dates, outcome_type, and idno (dropping/ignoring what type of date)
  q_long <- unique(
    melt(q, id.vars = c(group_vars,id_var),
         measure.vars=c("actual_event_date",interval_vars),value.name="value"
    )[,c(group_vars,id_var,"value"),with=FALSE]
  )

  #if there's a missing in ppt/outcome_type, it should be because an NA value was created
    #for actual_event_date during the foverlaps merge
  #this symbolizes that there aren't any other subject's events in that ppt's follow-up period and
   #it can therefore be safely dropped
  #make sure that if there's an NA for a idno/outcome_type group, it's exactly 1 NA and 2 nonmissing values
    #(that subject's baseline_date and that subject's "event_date"
  #(could be event or censor date) ) for that group.


  eval(parse(text=paste0(
    "q_long[q_long[is.na(value)],list(",paste0(c(id_var, group_vars),collapse=","),",value,isna_value=is.na(value)),",
    "on=c(id_var,group_vars)][,list(.N,V1=sum(isna_value)),by=c(id_var,group_vars)][,stopifnot(all(N==3)&V1==1)]"

  )))
  q_long <- q_long[!is.na(value)] #drop the rows with NA




  ##now that we have a list of all relevant dates for each idno/outcome_type combo,
  #use the value in the next row to define the end of the current row
  setkeyv(q_long,c(group_vars,id_var,"value")) #setkey sorts. You need to sort for shift to work correctly.
  q_long[, end:=shift(value,type="lead"),by=c(group_vars,id_var)]
  q_long <- q_long[!is.na(end)] #last row in each group is NA due to shift.
  #unnecessary to keep this row since the previous row already ends on the last date for that ppt
  setnames(q_long, "value","start")

  q_long[]



  ###q_long contains the intervals that will be used in the cox model merge it back into the events data to get the other columns back

  out <- q_long[x,on=c(group_vars,id_var),nomatch=NULL] #inner join
  #make sure there aren't any records lost to the inner join from either table. there shouldn't be.
  stopifnot(nrow(q_long[!x,on=c(group_vars,id_var)])==0)
  stopifnot(nrow(x[!q_long,on=c(group_vars,id_var)])==0)

  #if a ppt has an event, it happens only in the last interval for that outcome
  eval(parse(text=paste0("out[",interval_vars[2],"!=end,",event_indicator,":=0L] ")))

  setcolorder(out, c(group_vars,id_var,"start","end",event_indicator))


  #make sure the number of events hasn't changed in x compared to events_raw_long
  eval(parse(text=paste0(
  "stopifnot(
    all.equal(
      as.data.frame(x[,sum(",event_indicator,"),by=group_vars]),
      as.data.frame(out[,sum(",event_indicator,"),by=group_vars])
    )
  )"
  )))
  out[]
}
