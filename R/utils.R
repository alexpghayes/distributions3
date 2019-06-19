#' Is an object a distribution?
#'
#' `is_distribution` tests if `x` inherits from `"distribution"`.
#'
#' @param x An object to test.
#'
#' @export
#'
#' @examples
#'
#' Z <- Normal()
#'
#' is_distribution(Z)
#' is_distribution(1L)
#'
is_distribution <- function(x) {
  inherits(x, "distribution")
}
