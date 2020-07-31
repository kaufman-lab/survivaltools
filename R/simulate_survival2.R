## good overview of weibell, exponential, and cox
 #https://krex.k-state.edu/dspace/bitstream/handle/2097/8787/angelacrumer2011.pdf?sequence=3


#Comparison of Algorithms to Generate Event Times Conditional on Time-Dependent Covariates
#https://onlinelibrary.wiley.com/doi/pdf/10.1002/sim.3092

#https://www.cambridge.org/core/services/aop-cambridge-core/content/view/1945D7548766E76FB31C6C833976822E/S2049847018000195a.pdf/simulating_duration_data_for_the_cox_model.pdf

logit <- stats::binomial(link="logit")$linkfun
logistic <- stats::binomial(link="logit")$linkinv #logistic returns probability




##brief exploration of baseline hazards for my own reference
#this is not actually a function, just a code I want to keep nearby for reference
f <- function(){

  #note hazards can be above 1 but they must be positive
  #but a hazard above 1 are not meaningfully
    #different from each other when the time axis is discrete

  logistic(c(1,2,1,3)*c(.5,1.1,1,60))
  logistic(logit(logistic(c(1,2,1,3)))*c(.5,1.1,1,60))

#one example of a baseline hazard function: constant probablity of event in each day of 0.000015
baseline_hazard=function(x){0.000015}

#another example of a baseline hazard function: linear on the logit-scale
baseline_hazard=function(x){logistic(logit(0.000015)+x*3)}

#linear on the logit scale is nothing like a constant multiplier on the probabilty scale for large slopes:
baseline_hazard(1:10)/baseline_hazard(0:9)
log(baseline_hazard(1:10))/log(baseline_hazard(0:9))
log(baseline_hazard(1:10)- baseline_hazard(0:9))
log(baseline_hazard(1:10))- log(baseline_hazard(0:9))

#however, small slopes combined with small intercepts correspond to something like a constant
#approximately log-linear
baseline_hazard=function(x){logistic(logit(0.000015)+x*(1/75))}
baseline_hazard(1:10)/baseline_hazard(0:9)
log(baseline_hazard(1:10)/baseline_hazard(0:9))

#even so, this 1/75 is a rapid change in risk if x are units of days:
#after one year the risk increases ~130 fold
baseline_hazard(365)/baseline_hazard(0)

NULL

}








#' this is a title
#'
#' this is a summary
#'
#' note that the t_exit column in x represents the last day on which a ppt could be observed. typically
#' the end of follow-up, and may be the same day for everyone.
#'
#'
#' @param x a data.table containing at a mininum a unique identifier column named ppt_id, an integer time column
#' named time, a positive numeric column named baseline_hazard. Participants are assumed to be at risk if a time/participant combination exists as a row.
#' Within groups defined by participants, the values of time (or number of rows in each group) need not be the same.
#'  May also contain covariates.
#' @param covar_names a vector of column names in x corresponding to numeric covariates
#' @param covar_coefs A vector of coefficients (same length as covar_names)
#' @param return_hazard If TRUE, the return data.table will contain an attribute which itself is a data.table
#' containing ppt-specific hazard at each time-point.
#' @return a data.table
#' @export
simulate_survival <- function(x,
                              covar_names=NULL,
                              covar_coefs=NULL,
                              return_hazard=FALSE

){


  stopifnot(all(c("ppt_id","time","baseline_hazard")%in%names(x)))


  if(any(c("t_start","t_end","event")%in%names(x))){stop("t_start,t_end,event are reserved names and cannot be in x. check names of x")}
  if("hazard" %in% names(x)){stop("hazard is a reserved column name. don't include it in x")}
  stopifnot(length(covar_names)==length(covar_coefs))

  stopifnot(all(covar_names%in%names(x)))

  stopifnot(class(x[["time"]])=="integer")


  if(length(covar_names)){
    x[,hazard:=baseline_hazard*exp(as.matrix(.SD)%*%c(covar_coefs)),.SDcols=covar_names]
  }else{
    x[,hazard:=baseline_hazard]
  }
  on.exit(x[,hazard:=NULL][],add=TRUE)
  x[hazard >1L, hazard:=1L]
  x[,event:=rbinom(.N,1,prob = hazard)]
  on.exit(x[,event:=NULL][],add=TRUE)


  #which.max(event) returns the index of largest value of event (event should be 1 or 0).
  #if there are multiple events, it retuns the index of the first event

  #non-gforce optimized by-statements
  x[,`:=`(any_event=any(event)) ,by="ppt_id"]

  x[,`:=`(end=fifelse(any_event[1],time[which.max(event)],as.integer(NA))), by="ppt_id"]
  on.exit(x[,any_event:=NULL][],add=TRUE)
  on.exit(x[,end:=NULL][],add=TRUE)


  #gforce optimized by-statement
  out <- x[,list(t_start=min(time),
          eof=max(time),
          t_end=end[1],
          event=any_event[1]
  ),
  by="ppt_id"]


  out[is.na(t_end),t_end:=eof] #if end is non-missing, set it to the last observed time
  out[,eof:=NULL]
  out[,t_start:=t_start-1L] #participants are assumed to be at risk for every day provided in x
    #but interval formulation of cox model starts with an open interval on the left
    #so the cox formulation of the interval is actually 1 day before the first day at risk

  out[, event:=as.integer(event)]
  setkey(out, ppt_id, t_start)
  if(return_hazard){
    setattr(out,"hazard",x[,list(ppt_id,time,hazard)])
  }
  out
}

