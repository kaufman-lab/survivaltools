#' Create a function from a set of basis functions and corresponding coefficients
#'
#' A function which returns a function which itself will return the linear combination of the provided
#' basis functions and coefficients. Coefficients may be vectors,
#'  allowing for id-specific linear combinations of those basis functions.
#'
#' Each set of coefficients may be a vector. The ordering of each coefficients
#' must correspond to each other. ie, if `coefs=list(c(1,2), c(3,4))` then 1 and 3 are taken
#' to be the coefficients for the first and second basis functions respectively for the first id and
#' 2 and 4 are taken to be the coefficients for the first and second basis functions respectively for the second id.
#' If any element of coefs has length greater than 1, then id must be specified because the return
#' function will have an id argument (which specifies which set of coefficients to use).
#'
#'@param basis_functions A list of functions of x that return a vector the same length as x
#'@param coefs A list of coefficients, the same length as basis_functions.
#' Each element of this list is either a vector of length-1 or the length of id.
#'@param id NULL, or a unique vector identifying participant numbers.
#'@return if id is NULL: a function of x which returns a vector indexed by x
#'if id is not NULL, a function of x and id (arguments in that order) which returns an m by n matrix
#' where m is length(x) and n is length(id).
#'@examples
#' set.seed(20)
#' f0 <- function(x){rep(1L,length(x))}
#' f1 <- function(x){x}
#' f2 <- generate_random_spline(c(1,100),stats::rnorm(20,sd=40),degree=10)
#' f3 <- generate_random_spline(c(1,100),stats::rnorm(20,sd=40),degree=10)
#'
#' #fixed coefficients
#' g <- generate_curve_function(list(f0,f1,f2,f3),coefs=list(30,5,1,1))
#' plot(1:100,g(1:100), type="l")
#'
#' #random coefficients
#' ids <- 1:10
#' g1 <- generate_curve_function(list(f0,f1,f2,f3),
#'                               coefs=list(stats::rnorm(length(ids), mean=30),
#'                                          stats::rnorm(length(ids),mean=5),
#'                                          stats::rnorm(length(ids)),
#'                                          stats::rnorm(length(ids))
#'                               ),
#'                               id=1:length(ids))
#' m1 <- g1(1:100,ids)
#' matplot(1:100,m1,col=1:ncol(m1),type="l",lty=1)
#'
#'
#'
#'  x <- 1:(365*20)
#'  x_range <- range(x)
#'  shape <- 7.5
#'  avg_daily_hazard=0.000003
#'  curve(dgamma(x,shape=shape,scale=avg_daily_hazard/shape),from=0, to=0.0001)
#'
#'  z <- generate_curve_function(
#'    basis_functions=list(
#'      function(x){x-mean(x_range)},
#'      generate_random_spline(x_range,
#'                             stats::rgamma(20,shape=shape,scale=avg_daily_hazard/shape),
#'                             degree=10)
#'    ),
#'    coefs=list(-5e-10,1)
#'  )
#'  plot(x,z(x),type="l")
#'
#'
#'
#'
#'  ppts <- 1:10
#'
#'  par(mfrow=c(3,2))
#'  for(j in 1:6){
#'    z2 <- generate_curve_function(
#'      basis_functions=list(
#'        function(x){rep(1L,length(x))},
#'        function(x){x-mean(x_range)},
#'        generate_random_spline(x_range, stats::rnorm(365,sd=20),degree=30),
#'        generate_random_spline(x_range,stats::rnorm(365,sd=20),degree=30)
#'      ),
#'      coefs=list(
#'        stats::rnorm(length(ppts),mean=18),
#'        stats::rnorm(length(ppts),-0.0009,sd=0.0001),
#'        stats::rnorm(length(ppts),sd=.4),
#'        stats::rnorm(length(ppts),sd=.4)
#'      ),
#'      id=ppts
#'    )
#'
#'  ppt_subset <-ppts[2:6]
#'  m2 <- z2(x=x,ppt_subset)
#'
#'  matplot(x,m2,col=1:ncol(m2),type="l",lty=1)
#' }
#'  par(mfrow=c(1,1))
#'
#'
#'
#'
#'
#'@export

generate_curve_function <- function(basis_functions, coefs,id=NULL){
  stopifnot(is.list(basis_functions))
  stopifnot(is.list(coefs))

  stopifnot(all(sapply(basis_functions, is.function)))
  stopifnot(identical(length(basis_functions),length(coefs)))


  coef_lengths <- sapply(coefs,length)

 if(any(coef_lengths)!=1L){
   if(!all(coef_lengths %in% c(1L, max(coef_lengths)))){
     stop("elements of coef which have length greater than 1 must all be the same length")
   }
 }

 n_id <- max(coef_lengths)

 if(n_id>1L & is.null(id)){ stop("id cannot be NULL when elements of coefs have length greater than 1")
 }

 if(!is.null(id)){
   if(length(id)!=n_id){
     stop("length(id) must be the same as max(sapply(coefs,length))")
   }
 }

 id_copy <- id


 if(n_id==1L){
   out <- function(x){
     stopifnot(all(!is.na(x)))

     #evaluate all the basis functions at x
     X <- sapply(basis_functions,function(f){f(x)})

     ##error check X
     if(is.list(X)){
       X_length <- sapply(X, length)
       if(!all(X_length[1]==X_length)){
         stop("basis functions are returning elements of different length")
       }else{
         stop("basis functions are being turned into a list by sapply for some reason.
              are basis functions all returning a vector?")
       }
     }

     coefs <- unlist(coefs) #coefs are all length-1L if n_id==1L
     as.vector(X%*% coefs)
   }

 }else{
   out <- function(x,id){

     stopifnot(all(!is.na(id)))
     stopifnot(all(!is.na(x)))

     match_id <- match(id,id_copy)

     #evaluate all the basis functions at x
     X <- sapply(basis_functions,function(f){f(x)})

     ##error check X:
     if(is.list(X)){
       X_length <- sapply(X, length)
       if(!all(X_length[1]==X_length)){
         stop("basis functions are returning elements of different length")
       }else{
         stop("basis functions are being turned into a list by sapply for some reason.
              are basis functions all returning a vector?")
       }
     }

     #X a matrix of covars that are a function of time
     #B is a matrix of coefficients that vary id

     B <- do.call("rbind", coefs)
     B <- B[,match_id]

    X%*%B

   }
 }
 out
}





#' Create a random curve using bs splines
#'
#' Generate a random set of x,y points then fit a spline through them.
#'
#'
#'
#' the length of random_points determines the number of points used
#' to fit the spline and therefore may limit the flexibility of the spline
#' even when degree is increased.
#'
#'Note that two of these points are allocated as end points (ie the values of x_range)
#'
#'bs will helpfully return a warning message when providing x values outside of the domain.
#'
#' The fitted/predicted values of the spline near the extent of
#' the input data can be quite variable.
#' This can be addressed by specifying a wider x_range than you plan on making predictions for.
#'
#'@param x_range A length-2 numeric vector specifying the inclusive interval over which the
#' curve should be defined.
#'@param random_points For the random splines, a vector of y values that will be used to generate the
#'spline. The length of random_points will determine the number of points used to generate the spline.
#'@param degree The degree for the random splines. Passed to bs function.
#'@return a function of x
#'@examples
#' set.seed(20)
#' f <- generate_random_spline(c(1,100),stats::rnorm(10),degree=4)
#' plot(1:100, f(1:100),type="l")
#'@export
generate_random_spline <- function(x_range, random_points, degree){
    stopifnot(is.numeric(x_range))
    stopifnot(length(x_range)==2L)

    n <- length(random_points)

    #generate a set of x values
    x <- c(x_range,stats::runif(n-2L,min=x_range[1],max=x_range[2]))
    #there's a chance you could end up with duplicate x points, but this won't cause an error.
    m <- stats::lm(random_points~splines::bs(x,degree=degree))

    rm(x) #not necessary due to scoping, but it makes the code clearer.
    function(x){
      stats::predict(m, newdata=list(x=x))
    }

}
