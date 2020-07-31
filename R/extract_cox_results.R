#' extract coefficients from a cox model and turn them into a table
#'
#' take a single cox model and extract the coeffients from summary
#'
#'@param model A model returned by coxph
#'@param RHS description of the model here. if NULL this will try to extrac the formulat from the call
#'but this may not work if you're fitting the model programmatically (ie via a loop)
#'@return a data.table of coefficients from the models including a column "RHS" for the
#'RHS of the formula from the call
#'@export
extract_coef_table <- function(model,RHS=NULL){
  out <- as.data.table(coef(summary(model)),keep.rownames=TRUE)
  if(is.null(RHS)){
    out[,RHS:=as.character(as.formula(model$call))[3]]
  }else{
    out[, RHS:=RHS]
  }
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
