#' Methods for including distributions as vctrs in tibbles
#'
#' Methods for \code{\link[vctrs]{vec_proxy}} and \code{\link[vctrs]{vec_restore}}
#' from \pkg{vctrs} in order to include \code{distribution} objects in
#' \code{\link[tibble]{tibble}} objects.
#'
#' @details The methods for \code{\link[vctrs]{vec_proxy}} and
#'   \code{\link[vctrs]{vec_restore}} from \pkg{vctrs} are needed so that
#'   \code{distribution} objects can be included as a vector column in
#'   (and extracted from) \code{\link[tibble]{tibble}} data frames.
#'   \code{vec_proxy} simply adds the class \code{data.frame} which is the
#'   actual underlying data structure used by \code{distribution} objects.
#'   This way the number of rows etc. can be correctly determined. Conversely,
#'   \code{vec_restore} strips off the additional \code{data.frame} class and
#'   restores the original \code{distribution} classes. Users typically do not
#'   need to call \code{vec_proxy} and \code{vec_restore} directly.
#'
#' @param x,to Objects inheriting from \code{distribution}.
#' @param ... Currently not used.
#'
#' @return The `vec_proxy` method returns a `distribution` object which
#'   additionally inherits of `data.frame` while the `vec_restore` method
#'   restores the original `distribution` classes.
#'
#' @examples
#' \dontshow{ if(!requireNamespace("tibble")) {
#'   if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
#'     stop("not all packages required for the example are installed")
#'   } else q() }
#' }
#' ## Poisson GLM for FIFA 2018 goals data
#' data("FIFA2018", package = "distributions3")
#' m <- glm(goals ~ difference, data = FIFA2018, family = poisson)
#' 
#' ## Predict fitted Poisson distributions for teams with ability differences
#' ## of -1, 0, 1 (out-of-sample) using the new data as a data.frame
#' nd <- data.frame(difference = -1:1)
#' nd$dist <- prodist(m, newdata = nd)
#' nd
#'
#' ## Do the same using the new data as a tibble
#' library("tibble")
#' nt <- tibble(difference = -1:1)
#' nt$dist <- prodist(m, newdata = nt)
#' nt
#'
#' @keywords internal

#' @rdname vec_proxy.distribution
#' @exportS3Method vctrs::vec_proxy distribution
vec_proxy.distribution <- function(x, ...) {
  class(x) <- c(class(x), "data.frame")
  return(x)
}

#' @rdname vec_proxy.distribution
#' @exportS3Method vctrs::vec_restore distribution
vec_restore.distribution <- function(x, to, ...) {
  class(x) <- class(to)
  return(x)
}
