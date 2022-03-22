#' Robust fitting
#' 
#' A wrapper around `disclapmix()` that tries to avoid errors. 
#' If `disclapmix()` fails, it tries one iteration with `glm_method` = 'glm.fit' 
#' and uses the resulting `v_matrix` as `init_v` in the next `disclapmix()` 
#' with `glm_method = 'internal_coef'`. 
#' Can sometimes avoid errors with SVD problems happening with 
#' `glm_method = 'internal_coef'` and `glm_method = 'internal_dev'`.
#' 
#' @inheritParams disclapmix
#' @param \dots Passed on to `disclapmix()`
#' 
#' @examples 
#' data(danes)
#' db <- as.matrix(danes[rep(1:nrow(danes), danes$n), 1:(ncol(danes)-1)])
#' fit <- disclapmix_robust(db, 3L)
#' fit
#' 
#' @export
disclapmix_robust <- function(x, clusters, ...) {
  verbose <- 0L
  dots <- list(...)
  if (!is.null(dots[["verbose"]])) {
    verbose <- dots[["verbose"]]
  }
  
  out <- tryCatch({
    disclapmix::disclapmix(x, clusters = clusters, ...)
  },
  
  error = function(cond) {
    if (verbose > 0L) {
      cat(as.character(Sys.time()), ": An error occured. ", 
        "Now trying glm.fit for first iteration to get sensible init_v.\n", sep = "")
    }
    
    # Ignore errors/warnings, not only 1 iteration
    fit0 <- suppressWarnings(disclapmix::disclapmix(x, clusters = clusters, iterations = 1L,
                                                    verbose = verbose,
                                                    eps = 1e-3,
                                                    init_y_method = "pam",
                                                    glm_method = "glm.fit",
                                                    glm_control_maxit = 100L,
                                                    glm_control_eps = 1e-4))
    
    new_init_v <- fit0$v_matrix

    fit <- disclapmix::disclapmix(x, clusters = clusters, init_v = fit0$v_matrix, ...)
    return(fit)
  })    
  
  return(out)
}

