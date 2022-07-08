#' Adaptive fitting
#' 
#' A wrapper around `disclapmix_robust()` that instead of fitting one model 
#' for a given number of clusters, fits models until the best model (lowest marginal BIC) 
#' is in the interior (with margin `M`) of all number of clusters tried. 
#' 
#' E.g., the best model has 3 clusters and the margin `M = 5`, then 
#' this function ensures that models with 1, 2, ..., 3+5 = 8 clusters 
#' are fitted. If e.g. then 7 is better than 3, then it continues such that 
#' also models with up to 7+5 = 12 clusters are fitted. 
#' 
#' Note that models with 1-5 clusters are always fitted.
#' 
#' @inheritParams disclapmix
#' @param margin Fit models until there is at least this margin
#' @param criteria The slot to chose the best model from (small values indicate better model)
#' @param init_y_generator Function taking the number of clusters as input and returns `init_y` values
#' @param init_v_generator Function taking the number of clusters as input and returns `init_v` values
#' @param \dots Passed on to `disclapmix_robust()` (and further to `disclapmix()`)
#' 
#' @returns A list of all `disclapmix` fits
#' 
#' 
#' @examples 
#' data(danes)
#' db <- as.matrix(danes[rep(1:nrow(danes), danes$n), 1:(ncol(danes)-1)])
#' fits <- disclapmix_adaptive(db, margin = 5L)
#' fits
#' BICs <- sapply(fits, function(x) x$BIC_marginal)
#' BICs
#' ks <- sapply(fits, function(x) nrow(x$y)) # Always same as seq_along(fits)
#' ks
#' max_k <- max(ks)
#' best_k <- which.min(BICs)
#' max_k
#' best_k
#' max_k - best_k # = margin = 5
#' plot(ks, BICs, type = "b")
#' 
#' fits_clara <- disclapmix_adaptive(db, margin = 5L, init_y_method = "clara")
#' 
#' @export
disclapmix_adaptive <- function(x, 
                                label = "DL", 
                                margin = 5L, criteria = 'BIC_marginal', 
                                init_y_generator = NULL, 
                                init_v_generator = NULL, ...) {
  stopifnot(margin >= 1L)
  
  disclapmix_robust_call <- function(k, ...) {
    if (is.null(init_y_generator) && is.null(init_v_generator)) {
      return(disclapmix_robust(x = x, clusters = k, ...))
    } else if (!is.null(init_y_generator) && !is.null(init_v_generator)) {
      init_y <- init_y_generator(k)
      init_v <- init_v_generator(k)
      return(disclapmix_robust(x = x, clusters = k, init_y = init_y, init_v = init_v, init_y_method = NULL, ...))
    } else if (!is.null(init_y_generator)) {
      init_y <- init_y_generator(k)
      return(disclapmix_robust(x = x, clusters = k, init_y = init_y, init_y_method = NULL, ...))
    } else if (!is.null(init_v_generator)) {
      init_v <- init_v_generator(k)
      return(disclapmix_robust(x = x, clusters = k, init_v = init_v, ...))
    } else {
      stop("Unexpected")
    }
  }
  
  # Always fit 5
  old_ks <- seq_len(5L)
  
  dl_fits <- lapply(old_ks, function(k) {
    disclapmix_robust_call(k = k, ...)
  })
  dl_fits_BIC <- sapply(dl_fits, function(x) x[[criteria]])
  best_k <- which.min(dl_fits_BIC)
  
  max_k <- tail(old_ks, 1)
  
  while (best_k > (max_k - margin)) {
    max_k <- max_k + 1L
    dl_fits[[max_k]] <- disclapmix_robust_call(k = max_k, ...)
    dl_fits_BIC <- unlist(lapply(dl_fits, function(x) x[[criteria]]))
    best_k <- which.min(dl_fits_BIC)
  }
  
  # Check that the number of clusters are as expected:
  dl_fits_k <- unlist(lapply(dl_fits, function(x) nrow(x$y)))
  max_k <- max(dl_fits_k)
  stopifnot(best_k <= (max_k - margin))
  
  dl_fits <- lapply(dl_fits, function(x) {
    x$label <- label
    x
  })
  
  return(dl_fits)
}
