% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expand_ppt_table.R
\name{expand_ppt_table}
\alias{expand_ppt_table}
\title{expand a ppt table to long format}
\usage{
expand_ppt_table(ppt_tbl)
}
\arguments{
\item{ppt_tbl}{typically the return of generate_ppt_tbl, possibly with additional columns added.
specifically: a data.table with columns ppt_id (assummed to be unique), t_start, t_exit, baseline_age,
and possibly other columns}
}
\value{
a data.table where the number of rows will be the number of participants times the total number of unique integer
time points in the set union of all the [t_start,t_exit] intervals in ppt_tbl.
Ever ppt_id receives a row for every time point, whether or not the participant was at risk at that time.
it has the following columns:
\itemize{
\item ppt_id
\item time
\item at_risk (indicator for whether participant was at risk during this time. this is just
whether this time is in that participant's interval which for the purpose of the risk
is considered open on the left and closed on the right.
ie this is an indicator for whether time is in (t_start,t_exit].
\item age_tv time-varying age. calculated for convenience
\item age_since65_tv time-varying age in units of days since age 65. assumes baseline_age was in units of days
\item any other columns that were in ppt_tbl
}
}
\description{
take a table with one row per ppt and expand that to every combination of ppt and every possible time point
}
