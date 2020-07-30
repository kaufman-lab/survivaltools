#' extract coefficients from a cox model and turn them into a table
#'
#' take a single cox model and extract the coeffients from summary
#'
#'@param model A model returned by coxph
#'@return a data.table of coefficients from the models including a column "RHS" for the
#'RHS of the formula from the call
#'@export
extract_coef_table <- function(model){
  out <- as.data.table(coef(summary(model)),keep.rownames=TRUE)
  out[,RHS:=as.character(as.formula(model$call))[3]]
  setcolorder(out, "RHS")
  out[]
}



#' calculate average coefficients from a table
#'
#' Calculate average coefficients from a table containing results from several model runs
#'
#'
#'
#'@param results_table Tables of coefficients. Use rbindlist on a list of results from extract_coef_table
#'@return a data.table of summarized coefficients
#'@export
summarize_simulation_table <- function(results_table){
  out <- results_table[, list(mean_coef=mean(coef),
                              mean_se=mean(`se(coef)`),
                              sd_coef=sd(coef)
  ),
  keyby=c("rn","RHS")]
  out[]
}
