#' Export model fit to compact format
#' 
#' Includes only the most important aspects (parameters). 
#' For example useful if exporting fits for prediction to JSON. 
#' 
#' @param fit a `disclapmixfit` object
#' 
#' @examples
#' data(danes)
#' db <- as.matrix(danes[rep(1:nrow(danes), danes$n), 1:(ncol(danes)-1)])
#' fit <- disclapmix(db, clusters = 4L)
#' str(export_compact_fit(fit), 1)
#' 
#' @export
export_compact_fit <- function(fit) {
  if (!is(fit, "disclapmixfit")) {
    stop("fit must be a disclapmixfit object")
  }
  
  x <- list(init_y_method = fit$init_y_method, 
            y = fit$y,
            tau = fit$tau,
            disclap_parameters = fit$disclap_parameters)
  
  x
}