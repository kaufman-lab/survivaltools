
#' Generate a table containing time-varying exposure for use with simulate_survival
#'
#' A function to generate a table containing a participant-specific
#' time-dependent covariate that varies smoothly over time.
#' the time-dependent covariates differ between individuals but are correlated because
#' they are generated from a random linear combination of the same underlying splines.
#'
#'
#' @param ppt_tbl A data.table; the return value of create_ppt_table
#' @param mean The expected value of the curve at the midpoint of time.
#' @param slope the average slope of centered time.
#' @param slope_sd the standard deviation of the (normally distributed) distribution of participant-specific slopes
#' (which is centered around slope)
#' @param random_spline_points_range Passed to runif to generate points to create random splines.
#' not the actual range of the data around the mean because the spline can go beyond this range
#' and becasue these splines are multiplied by coefficients, but altering this will alter the
#' range of the variabilty around the linear trend determined by slope and mean.
#' @param degree  passed to generate_random_spline which is passed to bs
#' @param proportion_points_sampled determines number of points sampled to make the curve
#' @param dofoverlaps if TRUE, final dataset gets subset to ppt-specific follow-up periods,
#' otherwise the dataset is left at every combination of ppt and time whether or not the ppt
#' was being followed at that time (including times before and after the start of follow-up that
#' generated so that the actual used portion of the curve isn't on the edge which is variable.
#' @return
#' a data.table with the following columns: ppt_id, time, exposure, exposure_2yravg
#' @export




generate_exposure_table <- function(ppt_tbl,mean=15, slope=-0.0009,slope_sd=.0001,
                                    random_spline_points_range=25,
                                    degree=35,
                                    proportion_points_sampled=0.1,
                                    dofoverlaps=TRUE
                                    ){

  avg_duration <- 365*2

  #add an time beyond where you'll be using the curve
    #so that the edges of the curve which are variable are never utilized.
  x_range <- ppt_tbl[, c(min(t_start),max(t_exit))]
  lower <- (x_range[1]-avg_duration) - diff(x_range)*0.1
  upper <- x_range[2] + diff(x_range)*0.1
  x_range_widened <- c(lower,upper)
  n_points <- max(5, floor(diff(x_range_widened)*proportion_points_sampled)) #select number of points to generate to build a spline

  v <- random_spline_points_range/2

  ppts_all <- unique(ppt_tbl$ppt_id)

  z <- generate_curve_function(
     basis_functions=list(
       function(x){rep(1L,length(x))},
       function(x){x-mean(x_range)},
       generate_random_spline(x_range_widened, stats::runif(n_points,-v,v),degree=degree),
       generate_random_spline(x_range_widened, stats::runif(n_points,-v,v),degree=degree)
     ),
     coefs=list(
       stats::rnorm(length(ppts_all),mean=mean),
       stats::rnorm(length(ppts_all),slope,sd=slope_sd),
       stats::rnorm(length(ppts_all),sd=.3),
       stats::rnorm(length(ppts_all),sd=.3)
     ),
     id=ppts_all
   )

 x_all <- seq(x_range_widened[1],x_range_widened[2])


 exposures <- z(x_all,ppts_all)
 setnames(exposures, c("id","x","value"),c("ppt_id","time","exposure"))

 exposures[,exposure_2yravg:=frollmean(exposure,365L*2L),by=c("ppt_id")]

 if(dofoverlaps){

   exposures[, time2:=time]
   setkey(exposures,ppt_id, time,time2)

   out <- foverlaps(ppt_tbl, exposures,by.x=c("ppt_id","t_start","t_exit"),nomatch=NULL)

   out[, time2:=NULL]
   out[, t_exit:=NULL]
   out[, baseline_age:=NULL]
   out[, t_start:=NULL]

   out[]
 }else{
   setkey(exposures,ppt_id, time)
   exposures[]
 }



}


