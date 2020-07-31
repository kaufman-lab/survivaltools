test_that("generate_random_spline", {

  set.seed(20)
  f <- generate_random_spline(c(-3.2,25),rnorm(365,sd=20),degree=30)
  a <- f(seq(-3.2,25,by=.08))

  set.seed(20)
  y <- rnorm(365,sd=20)
  x <- c(c(-3.2,25),runif(363,min=-3.2,max=25))
  b <- predict(lm(y~splines::bs(x,degree=30)),newdata=list(x=seq(-3.2,25,by=.08)))

  expect_equal(a,b)
})



test_that("generate_curve_function", {


  x <- 1:(365*20)
  x_range <- range(x)
  set.seed(10)
  ppts <- 1:10

  f1 <- function(x){rep(1L,length(x))}
  f2 <- function(x){x-mean(x_range)}
  f3 <- generate_random_spline(x_range,rnorm(365,sd=20),degree=30)
  f4 <- generate_random_spline(x_range,rnorm(365,sd=20),degree=30)


  coef1 <- rnorm(length(ppts),mean=18)
  coef2 <- rnorm(length(ppts),-0.0009,sd=0.0001)
  coef3 <- rnorm(length(ppts),sd=.4)
  coef4 <- rnorm(length(ppts),sd=.4)



  z2 <- generate_curve_function(
    basis_functions=list(f1,f2,f3,f4),
    coefs=list(coef1,coef2,coef3,coef4),
    id=ppts
  )

  ppt_subset <- ppts[2:6]
  exposures <- z2(ppt_subset,x=x)


  v <- list()
  for(i in 1:length(ppt_subset)){
    current_ppt <- ppts %in% ppt_subset[i]
    v[[i]] <- coef1[current_ppt]*f1(x) + coef2[current_ppt]*f2(x) +
                           coef3[current_ppt]*f3(x) + coef4[current_ppt]*f4(x)
  }

  vv <- do.call("cbind",v)


  expect_true(all.equal(as.vector(vv), as.vector(exposures)))
})
