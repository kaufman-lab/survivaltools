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
#' z must contain a row for every discrete time point in a ppt_id's range (as specified by t_start and t_exit in x).
#' keeping in mind that t_start, t_exit is open on the left and closed on the right as is typical for survival data.
#'
#'
#' @param x a data.table containing at a mininum a unique identifier column named ppt_id, an integer column
#'  representing time of study entry named t_start, and an integer vector representing time of study exit named t_exit.
#'  May also contain time-invariant covariates columns; these must be of class numeric (factors not implemented).
#' @param z a data.table containing a unique identifier column named id, an integer column representing time
#' named time, and a numeric column named baseline_hazard containing non-negative values. z may also contain
#' time-variant covariate columns; these must be of class numeric.
#' @param x_covar_names a vector of column names in x corresponding to numeric covariates
#' @param z_tv_covar_names a vector of column names in x corresponding to numeric time-varying covariates
#' @param covar_coefs A vector of coefficients (same length as x_covar_names)
#' @param tv_covar_coefs A vector of coefficients (same length as z_tv_covar_names)
#' @param return_hazard If TRUE, the return data.table will contain an attribute which itself is a data.table
#' containing ppt-specific hazard at each time-poin.
#' @export
simulate_survival <- function(x,
                              z,
                              x_covar_names=NULL,
                              z_tv_covar_names=NULL,
                              covar_coefs=NULL,
                              tv_covar_coefs=NULL,
                              return_hazard=FALSE

){

  original_names_in_x <- copy(names(x))
  stopifnot(all(c("ppt_id","time","baseline_hazard")%in%names(z)))
  stopifnot(all(c("ppt_id","t_start","t_exit")%in%names(x)))

  if(any(c("start","end","event")%in%names(x))){stop("start,end,event are reserved names and cannot be in x. check names of x")}

  stopifnot(length(x_covar_names)==length(covar_coefs))
  stopifnot(length(z_tv_covar_names)==length(tv_covar_coefs))

  stopifnot(all(x_covar_names%in%names(x)))
  stopifnot(all(z_tv_covar_names%in%names(z)))

  stopifnot(class(x[["t_start"]])=="integer")
  stopifnot(class(x[["t_exit"]])=="integer")
  stopifnot(class(z[["time"]])=="integer")

  stopifnot(is.numeric(z[["baseline_hazard"]]))
  stopifnot(all(z[["baseline_hazard"]]>=0))

  stopifnot(sum(duplicated(x[["ppt_id"]]))==0)

  ##state restore z
  stopifnot(all(! c("time2","survivaltools_zindex")%in%names(z) ))
  z[, time2:=time]
  z[,survivaltools_zindex:=1:.N]
  z_key <- key(z)
  on.exit(z[,time2:=NULL])
  on.exit(setorderv(z,"survivaltools_zindex"),add=TRUE)
  on.exit(setkeyv(z,z_key),add=TRUE)
  on.exit(z[,survivaltools_zindex:=NULL][],add=TRUE)


  ##state restore x
  stopifnot(all(! c("survivaltools_xindex")%in%names(x) ))
  x[,survivaltools_xindex:=1:.N]
  x_key <- key(x)
  on.exit(setorderv(x,"survivaltools_xindex"),add=TRUE)
  on.exit(setkeyv(x,x_key),add=TRUE)
  on.exit(x[,survivaltools_xindex:=NULL][],add=TRUE)




  ##check to make sure every day in x is represented in z
   #add one to the left interval in x since survival intervals are open on the left
   #and therefore we don't need data on that day

  stopifnot(all(! c("t_startplus1")%in%names(x) ))
  x[, t_startplus1:=t_start+1L]
  on.exit(x[,t_startplus1:=NULL][],add=TRUE)
  setkeyv(z,c("ppt_id","time","time2"))
  setkeyv(x,c("ppt_id","t_startplus1","t_exit"))
  j <- foverlaps(z,x,nomatch=NULL)
  all_days_in_x_represented_in_z <- j[,list(V1=.N==t_exit[1]-t_start[1]),by=c("ppt_id")][, all(V1)]
  if(!all_days_in_x_represented_in_z){
    stop("not every ppt-day combination indicated in the combination of ppt_id and t_start and t_exit values in x
                   are represented in z")
  }

  if("hazard" %in% names(j)){stop("hazard is a reserved column name. don't include it in x or z")}

  all_covars <- c(x_covar_names,z_tv_covar_names)
  if(length(all_covars)){
    j[,hazard:=baseline_hazard*exp(as.matrix(.SD)%*%c(covar_coefs,tv_covar_coefs)),.SDcols=all_covars]
  }else{
    j[,hazard:=baseline_hazard]
  }
  j[hazard >1L, hazard:=1L]
  if(any(j$hazard<0)){stop("hazard cannot be negative")}
  j[,event:=rbinom(.N,1,prob = hazard)]

  setkey(j,ppt_id,time)
  j[,any_event:=any(event) ,by="ppt_id"]
  out <- j[,list(start=t_start[1],
          end=fifelse(any_event[1],time[which.max(event)],t_exit[1]),
          event=as.integer(any_event[1])
  ),
  by="ppt_id"]

  stopifnot(nrow(out)==nrow(x))
  out2 <- out[x,on="ppt_id",nomatch=NULL][,c("ppt_id","start","end","event",
                                     setdiff(original_names_in_x,c("ppt_id","t_start","t_exit"))
                                     ),
                                  with=FALSE]
  if(return_hazard){
    setattr(out2,"hazard",j[,list(ppt_id,time,hazard)])
  }
  out2
}

