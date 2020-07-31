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
#' @param x A vector of time points
#' @param daily_hazard_mean A length-1 numeric vector signifying the desired average hazard
#' (ie probability of an event conditional on being at risk at that time)
#'  Since the slope (specified in `daily_hazard_slope`) is on centered time, the expected
#'   value of the curve at the midpoint of the input times equals daily_hazard_mean.
#' @param daily_hazard_range parameters passed to runif to determine range of possible points used
#' to generate curve. not the actual possible range since the spline curve may exceed the range of the data.
#' @param daily_hazard_slope The slope of centered time on hazard curve (specify in units of hazard/day).
#' @param degree  passed to generate_random_spline which is passed to bs
#' @param proportion_points_sampled determines number of points sampled to make the random spline curve
#' @return a vector of baseline hazards
#'@examples
#'
#' ppt_tbl <- generate_ppt_table(1000)
#' e <- expand_ppt_table(ppt_tbl)
#'
#' setkey(e, ppt_id, time)
#' times <- e[ppt_id==1,time]
#' e[, baseline_hazard:=rep(generate_baseline_hazard(times), times=uniqueN(ppt_id))]
#'
#'
#' @export

generate_baseline_hazard <- function(x,
                                daily_hazard_mean=3e-06,
                                daily_hazard_range=1e-06,
                                daily_hazard_slope=-1e-10,
                                degree=20,
                                proportion_points_sampled=0.1
                                ){


  x_range <- range(x)
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

  z(x)
}
