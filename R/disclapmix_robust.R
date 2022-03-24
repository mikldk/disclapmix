#' Robust fitting
#' 
#' A wrapper around `disclapmix()` that tries to avoid errors. 
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
      cat(as.character(Sys.time()), ": An error occured. Trying different approaches.\n", sep = "")
    }
    
    ########################################################################
    # First try on only min(5*clusters, round(nrow(x)/2)) clusters
    # 10 times
    ########################################################################
    no_rnd_obs <- min(5L*clusters, round(nrow(x)/2))

    if (verbose > 0L) {
      cat(as.character(Sys.time()), ": Tries with ", no_rnd_obs, " random observations.\n", sep = "")
    }
    
    for (iter_try in seq_len(10L)) {
      if (verbose > 0L) {
        cat(as.character(Sys.time()), ": -------------------------------------------------------\n", sep = "")
        cat(as.character(Sys.time()), ": Try #", iter_try, " with ", no_rnd_obs, " random observations.\n", sep = "")
        cat(as.character(Sys.time()), ": -------------------------------------------------------\n", sep = "")
        cat(as.character(Sys.time()), ": > Start try #", iter_try, "\n", sep = "")
      }
      
      new_x <- x[sample(seq_len(nrow(x)), no_rnd_obs), ]
      
      fit_rnd <- tryCatch({
        disclapmix::disclapmix(new_x, clusters = clusters, ...)
      },
      error = function(cond) {
        return(NULL)
      })
      
      if (verbose > 0L) {
        cat(as.character(Sys.time()), ": < End try #", iter_try, "\n", sep = "")
      }
      
      if (is(fit_rnd, "disclapmixfit")) {
        fit_rnd_inity <- tryCatch({
          disclapmix::disclapmix(x, clusters = clusters, 
                                 init_y = fit_rnd$y, 
                                 init_y_method = NULL, ...)
        },
        error = function(cond) {
          return(NULL)
        })

        if (is(fit_rnd_inity, "disclapmixfit")) {
          return(fit_rnd_inity)
        }
      }
    }
    
    
    ########################################################################
    # Else, try glm.fit
    ########################################################################
    if (verbose > 0L) {
      cat(as.character(Sys.time()), ": Now trying glm.fit for first iteration to get sensible init_v.\n", sep = "")
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

