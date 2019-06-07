#' Is object a distribution?
#'
#' `is_distribution` tests if `x` is an distribution.
#' It simply checks that the object has a `"distribution"`
#' subclass.
#'
#' @param x An object to test.
#'
#' @export
#'
#' @examples
#' n <- normal()
#' is_distribution(n)

is_distribution <- function(x) {
  inherits(x, "distribution")
}
