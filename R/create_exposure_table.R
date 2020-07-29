
#' A function to generate a time-varying variable for simulation
#'
#' A function to generate participant-specific random variables that vary smoothly over time.
#'
#'
#' First the function generates four temporal trends:
#' two temporal trends are random curves (splines fit to random points).
#' There is also a constant term and a linear term.
#'
#' Each participant has their own coefficient to each temporal trend.
#' the slope argument determines the slope of the linear trend
#'
#' @param ppt_tbl A data.table; the return value of create_ppt_table
#' @param slope the average slope for the linear temporal trend
#' @param slope_sd the standard deviation of the (normally distributed) distribution of participant-specific slopes
#' (which is centered around slope)
#' @param n_years Number of years over raw spline data to take average
#' @export


create_exposure_table <- function(ppt_tbl,slope=-0.0009,slope_sd=.0001, n_years=2L){
  #add an extra year to so you don't get weird divergence at the beginggin/end of the spline
  max_texit <-  max(ppt_tbl$t_exit)+365L
  min_tstart <- min(ppt_tbl$t_start)-365L*(n_years+1)
  #need n_years years of data prior to the minimum start date to calculate a n_years year avg

  x_all_original <- sort(ppt_tbl[,seq(min(t_start),max(t_exit))])
  x_all <- sort(ppt_tbl[,seq(min_tstart,max_texit)])
  n_points <- max(5, floor(length(x_all)*.05))
  #two splines
  df <- data.frame(x1=c(range(x_all),sample(setdiff(x_all,range(x_all)),n_points)),
                   x2=c(range(x_all),sample(setdiff(x_all,range(x_all)),n_points))
  )

  df$y1 <- rnorm(nrow(df),mean=0,sd=20)
  df$y2 <- rnorm(nrow(df),mean=0,sd=20)
  m1 <- lm(y1~splines::bs(x1,degree=30),data=df)
  m2 <- lm(y2~splines::bs(x2,degree=30),data=df)

  #a matrix of covariate data used to generate a single pm25 spline
  #covariates are: a constant, two splines, a linear term where
  #"zero" is the first observed day in ppt_tbl
  #having that be the zero makes picking a slope easier
  covar <- cbind(1,
                 predict(m1,newdata=data.frame(x1=x_all)),
                 predict(m2,newdata=data.frame(x2=x_all)),
                 x_all-min(x_all)
  )


  #plot(s1~x,data=df_timesplines, type="l",col="blue")
  #lines(s2~x,data=df_timesplines, col="green")

  #a matrix of beta-fields. ie: ppt-specific coefficients to splines
  #rows are ppt, columns are variables
  b <- cbind(rnorm(nrow(ppt_tbl),mean=18,sd=1),
             rnorm(nrow(ppt_tbl),sd=.4),
             rnorm(nrow(ppt_tbl),sd=.4),
             rnorm(nrow(ppt_tbl),slope,sd=slope_sd)
  )
  #linear decrease for everyone at slightly different rates


  q <- as.data.table(covar%*%t(b))
  setnames(q,as.character(sort(ppt_tbl$ppt_id)))
  q[, time:=x_all]
  out <- melt(data=q,
              id.vars="time",
              measure.vars=as.character(sort(ppt_tbl$ppt_id)),
              variable.name="ppt_id",
              value.name="pm25"
  )
  out[, ppt_id:=as.integer(ppt_id)]
  setcolorder(out,c("ppt_id","time","pm25"))
  setkey(out, ppt_id,time,pm25)
  #return(out[])

  mm <- frollmean(q[, as.character(sort(ppt_tbl$ppt_id)),with=FALSE],365*n_years)
  out_avg <- data.table(time=rep(x_all,length(ppt_tbl$ppt_id)),
                        ppt_id=rep(sort(ppt_tbl$ppt_id),each=length(x_all)),
                        pm25_avg=unlist(mm))

  out_final <- out[out_avg,on=c("ppt_id","time"),nomatch=NULL]
  out_final[]
}
