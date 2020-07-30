#' Generate a table containing baseline hazard for use with simulate_survival
#'
#'
#' The return table has a row for each combination of ppt and time,
#'  but the hazard is the same for all ppts.
#'
#'  Degree affects the flexibility of the curve but the flexiblity of the curve
#'  is limited by the number of points that are sampled.
#'proportion_points_sampled
#'
#' @param ppt_tbl A data.table; the return value of create_ppt_table
#' @param daily_hazard_mean A length-1 numeric vector signifying the desired average hazard
#' (ie probability of an event conditional on being at risk at that time)
#'  Since the slope (specified in `daily_hazard_slope`) is on centered time, the expected
#'   value of the curve at the midpoint of the input times equals daily_hazard_mean.
#' @param daily_hazard_range parameters passed to runif to determine range of possible points used
#' to generate curve. not the actual possible range since the spline curve may exceed the range of the data.
#' @param daily_hazard_slope The slope of centered time on hazard curve (specify in units of hazard/day).
#' @param degree  passed to generate_random_spline which is passed to bs
#' @param proportion_points_sampled determines number of points sampled to make the random spline curve
#' @return a data.table with ppt_id, time (integer days since 1970-01-01),
#'  numeric column baseline_hazard,
#'  integer columns age_tv and age_since65_tv
#'  (in units of days since birthdate and days since age 65 respectively)
#'@examples
#' ppt_tbl <- generate_ppt_table(1000)
#' hazard_tbl <- generate_hazard_table(ppt_tbl)
#' hazard_tbl[ppt_id==1, plot(time,baseline_hazard, type="l")]
#' @export

generate_hazard_table <- function(ppt_tbl,
                                daily_hazard_mean=3e-06,
                                daily_hazard_range=1e-06,
                                daily_hazard_slope=-1e-10,
                                degree=20,
                                proportion_points_sampled=0.1
                                ){


  x_range <- ppt_tbl[, c(min(t_start),max(t_exit))]


  x_range_widened <- x_range + c(-1L, 1L)*diff(x_range)*0.1
  n_points <- max(5, floor(diff(x_range_widened)*proportion_points_sampled)) #select number of points to generate to build a spline

  z <- generate_curve_function(
    basis_functions=list(
      function(x){rep(1L,length(x))},
      function(x){x-mean(x_range)},
      generate_random_spline(x_range_widened,
                             stats::runif(n_points,-daily_hazard_range/2, daily_hazard_range/2),
                             degree=degree)
    ),
    coefs=list(daily_hazard_mean,daily_hazard_slope,1L)
  )

  x_all <- seq(x_range[1],x_range[2])

  haz <- z(x_all) #has columns x and value
  haz[, value:=pmin(1,pmax(0,value))]
  setnames(haz,c("x","value"),c("time","baseline_hazard"))

  ##expand hazard table to every participant by merging back to ppt_tbl on time
  haz[, time2:=time]
  setkey(haz,time,time2)

  out <- foverlaps(ppt_tbl, haz,by.x=c("t_start","t_exit"),nomatch=NULL)
  out[,age_tv:=time-t_start+baseline_age]
  out[,age_since65_tv:=as.integer(floor(age_tv-65*365.25))]

  out[, time2:=NULL]
  out[, t_exit:=NULL]
  out[, baseline_age:=NULL]
  out[, t_start:=NULL]

  out[]
}
