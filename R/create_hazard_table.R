#' Generate a baseline hazard for simulation
#'
#' A function to randomly generate baseline hazard as a function of time for each ppt
#'
#'
#'
#' The return table has a row for each combination of ppt and time,
#'  but the hazard is the same for all ppts.
#'
#' @param ppt_tbl A data.table; the return value of create_ppt_table
#' @param avg_daily_hazard A length-1 numeric vector signifying the desired average hazard
#' (ie probability of an event conditional on being at risk at that time)
#'   return a vector the of length n. Signifies the age at enrollment.
#' @param shape Controls the variabililty around avg_daily_hazard
#'   return a vector the of length n. Signifies calendar time at study entry (ie enrollment date).
#' @return a data.table with ppt_id, time, numeric column hazard,  integer columns age_tv (in units of days) age_since65_tv
#'  (both in units of days since 365.25*65)
#' @export
create_hazard_table <- function(ppt_tbl,avg_daily_hazard=0.000003,shape=7.5){
  x_all <- sort(ppt_tbl[,seq(min(t_start),max(t_exit))])
  n_points <- max(5, floor(length(x_all)*.01)) #select number of points to generate to build a spline

  #create a data.frame of n_points randomly selected times in the range plus the start and end times
  df <- data.frame(x=
                     c(range(x_all),
                       sample(setdiff(x_all,range(x_all)),n_points)
                     )
  )

  #for each time chosen, pick a value from the gamma distribution which will be the hazard at that point
  df$y <- pmin(rgamma(n=length(df$x),shape=shape,scale=avg_daily_hazard/shape),1)
  #x <- seq(0,.0001,by=0.000000001)
  #shape <- 7.5
  #plot(x,dgamma(x,shape=shape,scale=avg_daily_hazard/shape),type="l")
  #shape <- 4
  #lines(x,dgamma(x,shape=shape,scale=avg_daily_hazard/shape),type="l",col="red")
  #shape <- 10
  #lines(x,dgamma(x,shape=shape,scale=avg_daily_hazard/shape),type="l",col="green")

  #estimate a spline based on those points
  m <- lm(y~splines::bs(x,degree=5),data=df)
  #plot(y~x,data=df)
  #lines(x_all,predict(m, newdata=data.frame(x=x_all)))

  out <- data.table(time=x_all)
  out[, baseline_hazard:=pmin(1,pmax(0,predict(m,newdata=data.frame(x=x_all))))]

  stopifnot(!"original_index"%in%names(ppt_tbl))
  ppt_tbl[,original_index:=1:.N]
  ppt_tbl_key <- key(ppt_tbl)
  on.exit(setkeyv(ppt_tbl,ppt_tbl_key))
  on.exit(setorderv(ppt_tbl,"original_index"),add=TRUE)
  on.exit(ppt_tbl[,original_index:=NULL],add=TRUE)

  out[,time2:=time]

  setkey(out,time,time2)
  setkey(ppt_tbl,t_start,t_exit)
  out <- foverlaps(ppt_tbl, out,nomatch=NULL)[, list(ppt_id,time,baseline_hazard,baseline_age,t_start)]
  #create time-varying age from baseline age
  out[,age_tv:=time-t_start+baseline_age]
  out[,age_since65_tv:=as.integer(floor(age_tv-65*365.25))]


  out[, baseline_age:=NULL]
  out[, t_start:=NULL]

  setkey(out, ppt_id,time)
  out[]
  #it's possible, but unlikely, that the spline might predict below 0
  #even thought the data generating the splines is constrained to be above 0

}
