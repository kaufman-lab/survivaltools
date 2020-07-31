
#' expand a ppt table to long format
#'
#' take a table with one row per ppt and expand that to every combination of ppt and every possible time point
#'
#' @param ppt_tbl typically the return of generate_ppt_tbl, possibly with additional columns added.
#'  specifically: a data.table with columns ppt_id (assummed to be unique), t_start, t_exit, baseline_age,
#' and possibly other columns
#' @return a data.table where the number of rows will be the number of participants times the total number of unique integer
#'  time points in the set union of all the \[t_start,t_exit\] intervals in ppt_tbl.
#'  Ever ppt_id receives a row for every time point, whether or not the participant was at risk at that time.
#' it has the following columns:
#' - ppt_id
#' - time
#' - at_risk (indicator for whether participant was at risk during this time. this is just
#'   whether this time is in that participant's interval which for the purpose of the risk
#'  is considered open on the left and closed on the right.
#'   ie this is an indicator for whether time is in (t_start,t_exit\].
#' - age_tv time-varying age. calculated for convenience
#' - age_since65_tv time-varying age in units of days since age 65. assumes baseline_age was in units of days
#' - any other columns that were in ppt_tbl
#'@export
expand_ppt_table <- function(ppt_tbl){
  if(key(ppt_tbl)[1]=="ppt_id"){
    id <-  ppt_tbl$ppt_id
  }else{
    id <- sort(ppt_tbl$ppt_id)
  }
  z <- CJ(ppt_id=id, time=seq(min(ppt_tbl$t_start),max(ppt_tbl$t_exit)))


  stopifnot(!c("time","age_tv","age_since65_tv")%in%names(ppt_tbl))

  ppt_tbl_expanded <- ppt_tbl[z,on=c("ppt_id")]
  ppt_tbl_expanded[, at_risk:=between(time,t_start+1L,t_exit)]
  setattr(ppt_tbl_expanded, "sorted",c("ppt_id","time"))
  ppt_tbl_expanded[,age_since65_tv:=as.integer(floor(time-t_start+baseline_age-65*365.25))]

  ppt_tbl_expanded[,t_start:=NULL ]
  ppt_tbl_expanded[,t_exit:=NULL ]
  setcolorder(ppt_tbl_expanded,c("ppt_id","time","at_risk","baseline_age","age_since65_tv"))

  ppt_tbl_expanded[]
}




#' expand participant table, add hazard, and exposures
#'
#' calls expand_ppt_table to expand a participant table, then calls
#' generate_exposure_matrix and generate_baseline_hazard to add those columns.
#'
#'
#' If you need to use non-default options to generate_exposure_matrix or generate_baseline_hazard,
#' this function could easily be edited to pass arguments to those respective functions as a list.
#'
#' @param ppt_tbl typically the return of generate_ppt_tbl, possibly with additional columns added.
#'  specifically: a data.table with columns ppt_id (assummed to be unique), t_start, t_exit, baseline_age,
#' and possibly other columns
#' @param include_exposure if TRUE, adds an exposure column. Turn this off if you're not using the exposure column since
#' it takes time to generate this
#' @return a data.table where the number of rows will be the number of participants times the total number of unique integer
#'  time points in the set union of all the \[t_start,t_exit\] intervals in ppt_tbl.
#'  Ever ppt_id receives a row for every time point, whether or not the participant was at risk at that time.
#' it has the following columns:
#' - ppt_id
#' - time
#' - at_risk (indicator for whether participant was at risk during this time. this is just
#'   whether this time is in that participant's interval which for the purpose of the risk
#'  is considered open on the left and closed on the right.
#'   ie this is an indicator for whether time is in (t_start,t_exit\].
#' - age_tv time-varying age. calculated for convenience
#' - age_since65_tv time-varying age in units of days since age 65. assumes baseline_age was in units of days
#' - exposure
#' - baseline_hazard
#' - any other columns that were in ppt_tbl
#'@export
prep_tv_dataset <- function(ppt_tbl,include_exposure=TRUE){

  e <- expand_ppt_table(ppt_tbl)
  times <- e$time[e$ppt_id==1]
  ppts <- unique(e$ppt_id)
  exposures <- generate_exposure_matrix(times,ppts)
  if(include_exposure) e[, exposure:=as.vector(exposures)]
  e[, baseline_hazard:=rep(generate_baseline_hazard(times), times=uniqueN(ppt_id))]
  e[at_risk==TRUE][, at_risk:=NULL][]
}



