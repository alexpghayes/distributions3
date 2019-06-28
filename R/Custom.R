#' Create a custom discrete distribution
#'
#' @param ss A vector specifying the sample space.
#' @param probs A vector of same length as `ss` giving the probabilities of each outcome. Must sum to 1.
#' Defaults to `rep(1/length(ss), length(ss))`
#'
#' @return A `Custom` object.
#' @export
#'
#' @family discrete distributions
#'
#' @examples
#'
#' X <- Custom(ss = c(1,6,2), probs = c(0.4, 0.1, 0.5))
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 0.7)
#' log_pdf(X, 0.7)
#'
#' cdf(X, 0.7)
#' quantile(X, 0.7)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 0.7))
#'
Custom <- function(ss = NULL, probs = NULL){

  if(is.null(probs)){
    probs <- rep(1/length(ss), length.out = length(ss))
  }

  if(length(probs) != length(ss))
    stop("The number of probabilities does not match the size of the sample space!")

  if(sum(probs) != 1)
    stop("The specified probabilities do not sum to one!")

  d <- list(ss = ss, probs = probs)
  class(d) <- c("Custom", "distribution")

  return(d)
}

#' @export
print.Custom <- function(x, ...) {
  if(length(x$ss) > 5){
    ss_print <- paste(c(x$ss[1:3], '...', x$ss[length(x$ss) - c(1:0)]), collapse = ', ')
    probs_print <- paste(c(round(x$probs, 3)[1:3], '...', round(x$probs, 3)[length(x$probs) - c(1:0)]), collapse = ', ')
  } else {
    ss_print <- paste(x$ss, collapse = ', ')
    probs_print <- paste(round(x$probs, 3), collapse = ', ')
  }

  cat(glue::glue("Custom discrete distribution (Sample space = {[ss_print]}, Probabilities = {[probs_print]})",
                 .open = '[',
                 .close = ']'))
}

#' Draw a random sample from a Custom distribution
#'
#' @inherit Custom examples
#'
#' @param d A `Custom` object created by a call to [Custom()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector containing values from `outcomes` of length `n`.
#' @export
#'

random.Custom <- function(d, n = 1L, ...){
  sample(x = d$ss, size = n, prob = d$probs, replace = TRUE)
}

#' Evaluate the probability mass function of a Custom discrete distribution
#'
#' @inherit Custom examples
#' @inheritParams random.Custom
#'
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Custom <- function(d, x, ...) {
  d$probs[d$ss == x]
}

#' Evaluate the cumulative distribution function of a Custom discrete distribution
#'
#' @inherit Custom examples
#' @inheritParams random.Custom
#'
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Custom <- function(d, x, ...) {
  if(!is.numeric(d$ss))
    stop("The sample space is not numeric, and so the cdf is not well-defined.")

  if(!x %in% d$ss)
    stop("x is not in the sample space.")

  probs <- setNames(d$probs, d$ss)[order(d$ss)]
  setNames(cumsum(probs)[x], NULL)
}

#' Determine quantiles of a Custom discrete distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Custom examples
#' @inheritParams random.Custom
#'
#' @param p A vector of probabilites.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `p`.
#' @export
#'
quantile.Custom <- function(d, p, ...) {
  if(!is.numeric(d$ss))
    stop("The sample space is not numeric, and so quantiles are not well-defined.")

  if(any(p > 1) | any(p < 0))
    stop("All elements of p must be between 0 and 1.")

  probs <- setNames(d$probs, d$ss)[order(d$ss)]
  tmp_cdf <- cumsum(probs)

  qs <- as.numeric(names(sapply(p, function(x) tmp_cdf[tmp_cdf >= x][1])))

  return(qs)
}
