
#' Generate time-varying exposure data for simulation
#'
#' A function to generate participant-specific random variables that vary smoothly over time
#'
#' Defaults generate baseline ages as integers, randomly, with mean of ~64*365.25.
#' Default enrollment date is generated between 2001-01-01 and 2017-01-01.
#' Default t_exit (last possible day of enrollment) is 2017-01-01.
#'
#' @param n The desired number of participants (ie rows) of the output
#' @param baseline_age_expr An expression written as a function of n which will
#'   return a vector the of length n. Signifies the age at enrollment.
#' @param t_start An expression written as a function of n which will
#'   return a vector the of length n. Signifies calendar time at study entry (ie enrollment date).
#' @param t_exit An expression written as a function of n which will
#'   return a vector the of length n. Signifies calendar time for last possible date of follow-up.
#'
#' @return A data.table with columns ppt_id, time, pm25, and pm25_avg (where pm25_avg is just the
#' rolling average of pm25 where the size of the rolling window is defined by n_years)
#' @examples
#' generate_ppt_table(10,
#'                 baseline_age=as.integer(floor(55*365.25 + rbinom(n,size=30, prob=.3)*365.25)),
#'                 t_start=pmin(
#'                   as.integer(floor(as.numeric(
#'                    as.IDate("2000-01-01")+
#'                       abs(rnorm(n,300,sd=400))
#'                   )),
#'                   as.integer(as.IDate("2017-01-01"))
#'                   )),
#'                 t_exit=as.integer(as.IDate("2017-01-01"))
#')
#'
#' @export
generate_ppt_table <- function(n,
                             baseline_age_expr=as.integer(floor(
                               55*365.25 + rbinom(n,size=30, prob=.3)*365.25)),
                             t_start=pmin(
                               as.integer(floor(as.numeric(
                                 as.IDate("2000-01-01")+
                                   abs(rnorm(n,300,sd=400))
                               )),
                               as.integer(as.IDate("2017-01-01"))
                               )),
                             t_exit=as.integer(as.IDate("2017-01-01"))
                             ){

  out <-  data.table(
    ppt_id=1L:as.integer(n),
    baseline_age=eval(substitute(baseline_age_expr),envir=list(n=n), enclos = parent.frame()),
    t_start=eval(substitute(t_start),envir=list(n=n), enclos = parent.frame())
    ,
    t_exit=eval(substitute(t_exit),envir=list(n=n), enclos = parent.frame())
  )
  setattr(out, "sorted",c("ppt_id","t_start","t_exit"))
  out[]
}
