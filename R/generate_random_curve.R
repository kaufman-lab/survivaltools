#'
#' to increase the flexibility of the curve, you need to specify higher n_points and degee.
#'
#'
#'
#'
#' The random splines will have an expected value of 0 at all points. The returned curve will also
#' have an expected value of <intercept> at all points, unless you specify linear_slope to be non-zero,
#' in which case the curve will have an expected value of zero only at the midpoint of the domain (since
#' linear_slope is the slope of centered x)
#'
#'@param ppts NULL, or a unique vector identifying participant numbers.
#'@param x_domain A length-2 numeric vector specifying the inclusive interval over which the
#' curve should be defined.
#'@param random_trends Number of random temporal trends to generate
#'(in addition to the always-present constant and linear terms). Specify 0 for no random trends.
#'@param intercept The coefficient(s) for the constant trend. A vector either length of ppts or length-1.
#'@param linear_slope The coefficient(s) for the linear trend of centered x.
#' A vector either length of ppts or length-1.
#'@param random_spline_coefs Coefficient(s) for the random splines. A list with the length of the number of
#'random trends specified via random_trends.
#' Each element of this list can either be length-1 or the length of ppts.
#'@param n_points For the random splines, the number of points (not including end points) that will be
#'used to generate the spline.
#' Can be either length-1 (the same number of points will be used for all random_splines) or
#' the length of the number of random trends specified via, specifying that each random spline will be generated
#'  from a different number of points.
#'@random_points_expr For the random splines, an expression of n used to generate
#' the y value of the random points. See examples
#'@param degree The degee for the random splines.
#' Passed to bs function.
#' @examples
#'
#' z(ppts,x=1:(365*20))
#'
#' par(mfrow=c(3,2))
#' for(j in 1:3){
#'   z <- generate_random_curve_function(ppts,x_domain=c(1,365.25*20),random_trends=2,
#'                                       intercept=rnorm(length(ppts),mean=18),
#'                                       linear_slope=rnorm(length(ppts),-0.0009,sd=0.0001),
#'                                       random_spline_coefs=list(
#'                                         rnorm(length(ppts),sd=.4),
#'                                         rnorm(length(ppts),sd=.4)
#'                                       ),
#'                                       n_points=365,
#'                                       random_points_expr=rnorm(n,mean=0,sd=20),
#'                                       degree=30)
#'   x <- 1:(365*20)
#'   exposures <- z(ppts,x=x)
#'   for(i in 1:length(ppts)){
#'     if(i==1){
#'       exposures[J(ppts[i]), plot(x=x,y=pred,type="l",col=i,ylim=c(0,30)),on="ppt"]
#'     }
#'     else{
#'       exposures[J(ppts[i]), lines(x=x,y=pred,col=i,ylim=c(0,30)),on="ppt"]
#'     }
#'   }
#'
#' }
#' par(mfrow=c(1,1))
#'
#'
#'
#' max(5, floor(365*20)*.05)
#' shape <- 7.5
#' avg_daily_hazard=0.000003
#' curve(dgamma(x,shape=shape,scale=avg_daily_hazard/shape),from=0, to=0.0001)
#'
#' z2 <- generate_random_curve_function(ppts=NULL,x_domain=c(1,20),random_trends=1,
#'                                      intercept=0,
#'                                      linear_slope=-0.00000009,
#'                                      random_spline_coefs=1L,
#'                                      n_points=20,
#'                                      random_points_expr=rgamma(n,shape=shape,scale=avg_daily_hazard/shape),
#'                                      degree=10)
#'
#' z2(x=1:4)
#'
#'
#'@return if ppts is NULL, a function of x. if ppts is non-null, a function of ppts and x.
#' Either way the return function will return a data.table with columns: x, pred, and possibly ppts

generate_random_curve_function <- function(ppts, x_domain,
                                           random_trends,
                                           intercept=0L,
                                           linear_slope=0L,
                                           random_spline_coefs=NULL,
                                           n_points=NULL,
                                           random_points_expr=NULL,
                                           degree=NULL){
  sexpr <- substitute(random_points_expr)
  stopifnot(uniqueN(ppts)==length(ppts))
  stopifnot(is.numeric(random_trends))
  n_random_spline_trends <- random_trends

  if(n_random_spline_trends ==0){
    if(!is.null(n_points)) warning("n_points argument specifies the number of points
                                   used to generate the trends corresponding to random splines. Since
                                   random_trends==0, no random splines are generated and n_points is therefore ignored")
  }

  if(n_random_spline_trends > 0){
    if(is.null(n_points)) stop("n_points cannot be null when specifying random splines via random_trends ")
    if(is.null(degree)) stop("degree cannot be null when specifying radom splines via random_trends ")
    if(is.null(random_spline_coefs)) stop("random_spline_coefs cannot be null when specifying random splines via random_trends ")
  }

  stopifnot(length(n_points) %in% c(1L,n_random_spline_trends))

  if(!is.null(ppts)){
    stopifnot(length(intercept) %in% c(1L, length(ppts)))
    stopifnot(length(linear_slope) %in% c(1L, length(ppts)))
    if(!is.null(random_spline_coefs)){
      stopifnot(all(sapply(random_spline_coefs, length) %in% c(1L, length(ppts))))
    }
  }else{
    stopifnot(length(intercept) %in% c(1L))
    stopifnot(length(linear_slope) %in% c(1L))
    if(!is.null(random_spline_coefs)){

      stopifnot(all(sapply(random_spline_coefs, length) %in% c(1L)))
    }
  }



  n <- n_points+2L #add 2 for the end points

  if(n_random_spline_trends>0L){
    #for each random spline, generate a set of y values according to random_points_expr
    y_list <- replicate(n_random_spline_trends,
              expr=eval(sexpr, list(n=n)),
              simplify = FALSE)
    if(is.null(y_list[[1]])){stop("random_points_expr is returning NULL")}

    #for each random spline, generate a set of x values
    x_list <- replicate(n_random_spline_trends,
              expr=c(x_domain,runif(n_points,min=x_domain[1],max=x_domain[2])),
              simplify=FALSE
              )
    #there's a chance you could end up with duplicate x points, but this won't cause an error.

    bs_models <- mapply(x_list,y_list, FUN=function(x,y){
      lm(y~splines::bs(x,degree=degree))
    },SIMPLIFY = FALSE)

  }


  x_midpoint <- mean(x_domain)




  if(is.null(ppts)){
    out <- function(x){
      stopifnot(all(!is.na(x)))
      stopifnot(all(between(x,x_domain[1],x_domain[2])))
      if(n_random_spline_trends>0L){
        spline_preds <- lapply(bs_models, predict,newdata=data.frame(x=x))
      }else{
        spline_preds <- NULL
      }



      X <- do.call("cbind",c(list(1L,x-x_midpoint),spline_preds))
      b <- c(intercept, linear_slope, unlist(random_spline_coefs))
      data.table(x=x,pred=X%*%b)
    }
  }else{
    out <- function(ppt,x){
      stopifnot(all(!is.na(ppt)))
      stopifnot(all(!is.na(x)))
      if(identical(ppt,ppts)){
        all_ppts <- TRUE
      }else{
        all_ppts <- FALSE
        match_ppt <- match(ppt,ppts)
        stopifnot(all(!is.na(match_ppt)))
      }

      stopifnot(all(between(x,x_domain[1],x_domain[2])))
      if(n_random_spline_trends>0L){
        spline_preds <- lapply(bs_models, predict,newdata=data.frame(x=x))
      }else{
        spline_preds <- NULL
      }

      #X a matrix of covars that are a function of time
      #B is a matrix of coefficients that vary by ppt
      X <- do.call("cbind",c(list(1L,x-x_midpoint),spline_preds))
      B <- do.call("cbind", c(list(intercept), list(linear_slope), random_spline_coefs))
      if(!all_ppts) B <- B[match_ppt,]

      data.table(x=rep(x,times=length(ppt)),
                 ppt=rep(ppt,each=length(x)),
                 pred=as.vector(X%*%t(B))
                 )

    }
  }
  out

}


ppts <- paste0("pptid_",3:10)
z <- generate_random_curve_function(ppts,x_domain=c(1,365.25*20),random_trends=2,
                               intercept=rnorm(length(ppts),mean=18),
                               linear_slope=rnorm(length(ppts),-0.0009,sd=0.0001),
                               random_spline_coefs=list(
                                 rnorm(length(ppts),sd=.4),
                                 rnorm(length(ppts),sd=.4)
                               ),
                               n_points=365,
                               random_points_expr=rnorm(n,mean=0,sd=20),
                               degree=30)






