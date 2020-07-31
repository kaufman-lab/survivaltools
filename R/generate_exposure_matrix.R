
#' Generate a table containing time-varying exposure for use with simulate_survival
#'
#' A function to generate a table containing a participant-specific
#' time-dependent covariate that varies smoothly over time.
#' the time-dependent covariates differ between individuals but are correlated because
#' they are generated from a random linear combination of the same underlying splines.
#'
#' All combinations of x and ppt_id arguments are crossed to create a matrix.
#'
#' @param x A vector of time points
#' @param ppt_id A vector of ppt_ids
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
#' a matrix with length(x) rows and length(ppt_id) columns
#' @examples
#'
#' ppt_tbl <- generate_ppt_table(1000)
#' e <- expand_ppt_table(ppt_tbl)
#'
#' setkey(e, ppt_id, time)
#'
#' times <- e$time[e$ppt_id==1]
#' ppts <- unique(e$ppt_id)
#' exposures <- generate_exposure_matrix(times,ppts)
#' e[, exposure:=as.vector(exposures)]
#'
#' #plot a subset of participant time-series
#' matplot(times,exposures[,50:70] ,col=1:ncol(exposures ),type="l",lty=1)
#'
#' @export


generate_exposure_matrix <- function(x, ppt_id, mean=15, slope=-0.0009, slope_sd=.0001,
                                    random_spline_points_range=25,
                                    degree=10,
                                    proportion_points_sampled=0.1,
                                    dofoverlaps=TRUE
                                    ){

     #add an time beyond where you'll be using the curve
    #so that the edges of the curve which are variable are never utilized.
  x_range <- range(x)
  lower <- x_range[1] - diff(x_range)*0.15
  upper <- x_range[2] + diff(x_range)*0.15
  x_range_widened <- c(lower,upper)
  n_points <- max(5, floor(diff(x_range_widened)*proportion_points_sampled)) #select number of points to generate to build a spline

  v <- random_spline_points_range/2

  z <- generate_curve_function(
     basis_functions=list(
       function(x){rep(1L,length(x))},
       function(x){x-mean(x_range)},
       generate_random_spline(x_range_widened, stats::runif(n_points,-v,v),degree=degree),
       generate_random_spline(x_range_widened, stats::runif(n_points,-v,v),degree=degree)
     ),
     coefs=list(
       stats::rnorm(length(ppt_id),mean=mean),
       stats::rnorm(length(ppt_id),slope,sd=slope_sd),
       stats::rnorm(length(ppt_id),sd=.3),
       stats::rnorm(length(ppt_id),sd=.3)
     ),
     id=ppt_id
   )

 z(x,ppt_id)
}


