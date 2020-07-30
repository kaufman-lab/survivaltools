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

  exposures <- z2(ppts[2:6],x=x)

  exposures_ppts <- ppts[2:6]

  v <- list()
  for(i in 1:length(exposures_ppts)){
    current_ppt <- ppts %in% exposures_ppts[i]
    v[[i]] <- data.table(x=x,id=exposures_ppts[i],value=
                           coef1[current_ppt]*f1(x) + coef2[current_ppt]*f2(x) +
                           coef3[current_ppt]*f3(x) + coef4[current_ppt]*f4(x))
  }

  vv <- rbindlist(v)

  expect_true(all.equal(vv, exposures))
})
