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
#'@return if id is NULL: a function of x which returns a data.table with columns "x" and "value".
#'if id is non-null, a function of both x and id which returns a data.table containing a row for every
#'combination of x and id,  with columns "x", "id",and "value".
#'@examples
#'set.seed(20)
#'f0 <- function(x){rep(1L,length(x))}
#'f1 <- function(x){x}
#'f2 <- generate_random_spline(c(1,100),20, stats::rnorm(n,sd=40),degree=10)
#'f3 <- generate_random_spline(c(1,100),20, stats::rnorm(n,sd=40),degree=10)
#'
#'#fixed coefficients
#'g <- generate_curve_function(list(f0,f1,f2,f3),coefs=list(30,5,1,1))
#'dt <- g(1:100)
#'plot(1:100,dt$value, type="l")
#'
#'#random coefficients
#'ids <- 1:10
#'g1 <- generate_curve_function(list(f0,f1,f2,f3),
#'                              coefs=list(stats::rnorm(length(ids), mean=30),
#'                                         stats::rnorm(length(ids),mean=5),
#'                                         stats::rnorm(length(ids)),
#'                                         stats::rnorm(length(ids))
#'                              ),
#'                              id=1:length(ids))
#'dt1 <- g1(1:100,ids)
#'plot(1,1,xlim=c(1,100),ylim=range(dt1$value),type="n")
#'sapply(ids,function(i){lines(dt1[id==i,x],dt1[id==i,value],col=i)})
#'
#'
#' x <- 1:(365*20)
#' x_range <- range(x)
#' shape <- 7.5
#' avg_daily_hazard=0.000003
#' curve(dgamma(x,shape=shape,scale=avg_daily_hazard/shape),from=0, to=0.0001)
#'
#' z <- generate_curve_function(
#'   basis_functions=list(
#'     function(x){x-mean(x_range)},
#'     generate_random_spline(x_range,n_points=20,
#'                            stats::rgamma(n,shape=shape,scale=avg_daily_hazard/shape),
#'                            degree=10)
#'   ),
#'   coefs=list(-5e-10,1)
#' )
#' plot(x,z(x)$value,type="l")
#'
#'
#'
#'
#' ppts <- 1:10
#'
#' par(mfrow=c(3,2))
#' for(j in 1:6){
#'   z2 <- generate_curve_function(
#'     basis_functions=list(
#'       function(x){rep(1L,length(x))},
#'       function(x){x-mean(x_range)},
#'       generate_random_spline(x_range,n_points=365, stats::rnorm(n,sd=20),degree=30),
#'       generate_random_spline(x_range,n_points=365, stats::rnorm(n,sd=20),degree=30)
#'     ),
#'     coefs=list(
#'       stats::rnorm(length(ppts),mean=18),
#'       stats::rnorm(length(ppts),-0.0009,sd=0.0001),
#'       stats::rnorm(length(ppts),sd=.4),
#'       stats::rnorm(length(ppts),sd=.4)
#'     ),
#'     id=ppts
#'   )
#'
#' exposures <- z2(ppts[2:6],x=x)
#' exposures_ppts <- unique(exposures$id)
#'
#' plot(1,1,type="n",xlim=x_range,ylim=range(exposures$value))
#' for(i in 1:length(exposures_ppts))
#'   exposures[J(exposures_ppts[i]), lines(x=x,y=value,col=i,ylim=c(0,30)),on="id"]
#' }
#' par(mfrow=c(1,1))
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
     data.table(x=x,value=as.vector(X%*% coefs))
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

     data.table(x=rep(x,times=length(id)),
                id=rep(id,each=length(x)),
                value=as.vector(X%*%B)
     )

   }
 }
 out
}





#' Create a random curve using bs splines
#'
#' Generate a random set of x,y points then fit a spline through them.
#'
#'@param x_range A length-2 numeric vector specifying the inclusive interval over which the
#' curve should be defined.
#'@param n_points For the random splines, the number of points (not including end points) that will be
#'used to generate the spline.
#'@param random_points_expr For the random splines, an expression of n used to generate
#' the y value of the random points. See examples
#'@param degree The degree for the random splines. Passed to bs function.
#'@return a function of x
#'@examples
#' set.seed(20)
#' f <- generate_random_spline(c(1,100),10, stats::rnorm(n),degree=4)
#' plot(1:100, f(1:100),type="l")
#'@export
generate_random_spline <- function(x_range, n_points, random_points_expr, degree){
    stopifnot(is.numeric(x_range))
    stopifnot(length(x_range)==2L)
    stopifnot(is.numeric(n_points))
    stopifnot(length(n_points)==1L)

    sexpr <- substitute(random_points_expr)
    n <- n_points+2L #add 2 for the end points

    #generate a set of y values according to random_points_expr
    y_points <- eval(sexpr, list(n=n))

    if(is.null(y_points)){stop("random_points_expr is returning NULL")}
    #generate a set of x values
    x <- c(x_range,stats::runif(n_points,min=x_range[1],max=x_range[2]))
    #there's a chance you could end up with duplicate x points, but this won't cause an error.
    m <- stats::lm(y_points~splines::bs(x,degree=degree))

    rm(x) #not necessary due to scoping, but it makes the code clearer.
    function(x){
      stats::predict(m, newdata=list(x=x))
    }

}
