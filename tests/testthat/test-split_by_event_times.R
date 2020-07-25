####manually construct a toy example z and solution z_expanded ####

##manually constructed toy example
z <-
  as.data.table(structure(
    list(
      ppt = 1:5,
      enroll_date = structure(c(12065, 12586,
                                11449, 11449, 11172), class = "Date"),
      fu_date = structure(c(14737,
                            13224, 15292, 15282, 11349), class = "Date"),
      event_indicator = c(0L,
                          1L, 1L, 1L, 0L),
      x = c(
        1.18991379345149,
        -0.448662309663428,-0.0366741417210273,
        1.23962758102562,
        0.861831762483491
      )
    ),
    row.names = c(NA,-5L),
    class = "data.frame"
  ))

#manually constructed toy solution
z_expanded <-
  structure(
    list(
      ppt = c(1L, 1L, 2L, 3L, 3L, 3L, 4L, 4L, 5L),
      start = structure(
        c(12065,
          13224, 12586, 11449, 13224, 15282, 11449, 13224, 11172),
        class = "Date"
      ),
      end = structure(
        c(13224, 14737, 13224, 13224, 15282, 15292,
          13224, 15282, 11349),
        class = "Date"
      ),
      event_indicator = c(0L,
                          0L, 1L, 0L, 0L, 1L, 0L, 1L, 0L),
      enroll_date = structure(
        c(12065,
          12065, 12586, 11449, 11449, 11449, 11449, 11449, 11172),
        class = "Date"
      ),
      fu_date = structure(
        c(14737, 14737, 13224, 15292, 15292,
          15292, 15282, 15282, 11349),
        class = "Date"
      ),
      x = c(
        1.18991379345149,
        1.18991379345149,
        -0.448662309663428,
        -0.0366741417210273,-0.0366741417210273,
        -0.0366741417210273,
        1.23962758102562,
        1.23962758102562,
        0.861831762483491
      )
    ),
    class = "data.frame",
    row.names = c(NA,-9L)
  )

test_that("basic splitting for toy example", {
  expect_equal(
    as.data.frame(split_by_event_times(z,interval_vars=c("enroll_date","fu_date"),
                         id_var="ppt",event_indicator="event_indicator")),
    z_expanded
  )
})


####compare cox result in binned model using bin count as time with using numeric date as time ####
###write a function that takes a raw unsplit survival dataset,
 #bins time,
 #splits the dataset (into andersen gill format) using the bintted time,
 #compares cox results when using binned time (ie dates) coerced to numeric as the time axis
   #with cox results when using integer count of bin froms from the origin bin
#the first uses days as the time scale and the second uses bins as the time scale
  #i expect these to be the same I just want to double check

#for the binned time that uses dates, this snaps to the end of the period (the default)

##RHS is formula of the form eg .~x1+x2
library(survival)
comparecox_bincount_bindates <- function(x,interval_vars,id_var,event_indicator,
                                         event_bins,RHS,group_vars=NULL){

  x <- copy(x)
  stopifnot(!"event_indicator" %in% setdiff(names(x),event_indicator))
  setnames(x, event_indicator,"event_indicator")

  stopifnot(!"ivar1" %in% setdiff(names(x),interval_vars))
  setnames(x, interval_vars[1],"ivar1")
  stopifnot(!"ivar2" %in% setdiff(names(x),interval_vars))
  setnames(x, interval_vars[2],"ivar2")

  interval_vars <- c("ivar1","ivar2")

  x[, ivar2:=bin_time(ivar2,event_bins)]
  x[, ivar1:=bin_time(ivar1,event_bins)]

  #remove any observations where the start date is also the end date
  x <- x[ivar2>ivar1]

  #### split the dataset into Anderson gill format based on dates####


  event_dataset_dates <- survivaltools::split_by_event_times(x,
                                                       interval_vars=interval_vars,
                                                       id_var=id_var,
                                                       event_indicator=event_indicator,
                                                       group_vars=group_vars)

  event_dataset_dates[, start:=as.numeric(start)]
  event_dataset_dates[, end:=as.numeric(end)]
  #even though it's time has already been binned/snapped to the end of the periods
    #it doesn't matter\



  x[, ivar1:=bin_time(ivar1,event_bins,snapto = "count")]
  x[, ivar2:=bin_time(ivar2,event_bins,snapto = "count")]


  event_dataset_counts <- survivaltools::split_by_event_times(x,
                                                             interval_vars=interval_vars,
                                                             id_var=id_var,
                                                             event_indicator=event_indicator,
                                                             group_vars=group_vars)



  coxph(update(Surv(start,end,event_indicator)~.,RHS),data=event_dataset_dates)
  list(
    coef(summary(coxph(update(Surv(start,end,event_indicator)~.,RHS),data=event_dataset_dates))),
    coef(summary(coxph(update(Surv(start,end,event_indicator)~.,RHS),data=event_dataset_counts)))
  )
}

event_bins <- data.table(start=seq(structure(7286L, class = c("IDate", "Date")),by=28,length=400),
                         end=seq(structure(7313L, class = c("IDate", "Date")),by=28,length=400)
)

q <- comparecox_bincount_bindates(x=z,interval_vars=c("enroll_date","fu_date"),
                             id_var="ppt",event_indicator="event_indicator",
                             event_bins=event_bins,RHS=.~x,group_vars=NULL)


test_that("bin counts vs numeric binned dates", {
  expect_equal(q[[1]],q[[2]]  )
})

###note that the cox results are very different between the binned model and the model that
 #uses the original unbinned dates:

z_expandedDT <- as.data.table(z_expanded)
z_expandedDT[, start:=as.numeric(start)]
z_expandedDT[, end:=as.numeric(end)]
coef(summary(coxph(Surv(start,end,event_indicator)~x,data=z_expandedDT)))
#this isn't surprising because the toy example was constructed to have a tie when binned
 #but not be tied on the raw date timescale
 #since there are only three events this tie makes a big difference
#the next step is to randomly generate several larger datasets and
  #compare binned vs unbinned coefficients to provide some certainty
 #that my binning approach isn't causing bias






stop()
#this doesn't work:
library(data.table)
library(survival)

#generate a dataset with time-invariant exposure and analyze it with a baseline timescale of age in days while adjusting
 #for calendar time

#start with some fake data I had already generated as an example for Jaime
x <- fread("
ppt enroll_date fu_date    event_indicator
1   2003-01-13  2010-05-08 0
2   2004-06-17  2006-03-17 1
3   2001-05-07  2011-11-14 1
4   2001-05-07  2011-11-04 1
5   2000-08-03  2001-01-27 0
6   2000-03-12  2008-01-27 0
7   2002-12-19  2004-08-29 1
")
x[, enroll_date:=as.IDate(enroll_date)]
x[, fu_date:=as.IDate(fu_date)]

#lots of shared ages at baseline to make this clear
x[, age_baseline:=c(25567,23947, 22507,23999,29367,22563,24389)]
x[, exposure:=rnorm(.N,mean=15)]
x[, time:=fu_date-enroll_date]
x[, age_fu:=age_baseline+time]

#valid:
summary(coxph(Surv(age_baseline,age_fu,event_indicator)~exposure,data=x))
summary(coxph(Surv(age_baseline,age_fu,event_indicator)~exposure+enroll_date,data=x))




#do some complicated Andersen Gill splitting using age as baseline
#x2 <- split_by_event_times(x,interval_vars=c("age_baseline","age_fu"),
#                    id_var="ppt",event_indicator="event_indicator")
#dput(x2)
x2 <- copy(structure(list(ppt = c(1L, 1L, 1L, 2L, 3L, 3L, 3L, 4L, 4L, 4L,4L, 5L, 6L, 6L, 6L, 7L, 7L), start = c(25567, 26350, 27832, 23947, 22507, 24585, 25008, 23999, 24585, 25008, 26350, 29367, 22563, 24585, 25008, 24389, 24585), end = c(26350, 27832, 28239, 24585, 24585, 25008, 26350, 24585, 25008, 26350, 27832, 29544, 24585, 25008, 25440, 24585, 25008), event_indicator = c(0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L), enroll_date = structure(c(12065L, 12065L, 12065L, 12586L, 11449L, 11449L, 11449L, 11449L, 11449L, 11449L, 11449L, 11172L, 11028L, 11028L, 11028L, 12040L, 12040L), class = c("IDate", "Date")), fu_date = structure(c(14737L, 14737L, 14737L, 13224L, 15292L, 15292L, 15292L, 15282L, 15282L, 15282L, 15282L, 11349L, 13905L, 13905L, 13905L, 12659L, 12659L), class = c("IDate", "Date")), age_baseline = c(25567, 25567, 25567, 23947, 22507, 22507, 22507, 23999, 23999, 23999, 23999, 29367, 22563, 22563, 22563, 24389, 24389), exposure = c(15.792449490405, 15.792449490405, 15.792449490405, 14.2628925303294, 14.7863559134546, 14.7863559134546, 14.7863559134546, 14.6779133344263, 14.6779133344263, 14.6779133344263, 14.6779133344263, 15.3640057296039, 15.7637972902047,15.7637972902047, 15.7637972902047, 15.2091881247596, 15.2091881247596), age_fu = c(28239, 28239, 28239, 24585, 26350, 26350, 26350,27832, 27832, 27832, 27832, 29544, 25440, 25440, 25440, 25008,
25008)), sorted = c("ppt", "start"), class = c("data.table", "data.frame"), row.names = c(NA, -17L)))

#in this new dataset, "start" and "end" are age in days
 # "start" and "end"  columns create an interval (start,end] which signifity the new tstart, tstop forumation (in units of age)
 #(start,end] are non-overlapping within ppt, no (start,end) interval never intersects any event date



x2[, days_since_enrollment:=(end-age_baseline)]
x2[, date_of_event:=enroll_date+days_since_enrollment]


summary(coxph(Surv(start,end,event_indicator)~exposure,data=x2)) #same as above
summary(coxph(Surv(start,end,event_indicator)~exposure+enroll_date,data=x2)) #same as above
summary(coxph(Surv(start,end,event_indicator)~exposure+date_of_event,data=x2))  #should be the same as the previous model but is slightly different

x2[, list(enroll_date,date_of_event,date_of_event)]
x2[, cor(as.numeric(enroll_date),as.numeric(date_of_event)),by=end]
