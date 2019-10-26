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
is_distribution <- function(x) {
  inherits(x, "distribution")
}

#' Plot the CDF of a distribution
#'
#' A function to easily plot the CDF of a distribution using `ggplot2`. Requires `ggplot2` to be loaded.
#'
#' @param d A `distribution` object
#' @param limits either `NULL` (default) or a vector of length 2 that specifies the range of the x-axis
#' @param p If `limits` is `NULL`, the range of the x-axis will be the support of `d` if this is a bounded
#'   interval, or `quantile(d, p)` and `quantile(d, 1 - p)` if lower and/or upper limits of the support is
#'   `-Inf`/`Inf`.  Defaults to 0.001.
#' @param plot_theme specify theme of resulting plot using `ggplot2`. Default is `theme_minimal`
#'
#' @export
plot_cdf <- function(d, limits = NULL, p = 0.001,
                     plot_theme = ggplot2::theme_minimal, ...){

  if(!"ggplot2" %in% loadedNamespaces())
    stop("You must load ggplot2 for this function to work.")

  if(is.null(limits))
    limits <- support(d)

  if(limits[1] == -Inf){
    limits[1] <- quantile(d, p = p)
  }

  if(limits[2] == Inf){
    limits[2] <- quantile(d, p = 1-p)
  }

  if(class(d)[1] %in% c('Bernoulli', 'Binomial', 'Geometric', 'HyperGeometric',
                        'NegativeBinomial', 'Poisson')){
    plot_df <- data.frame(x = seq(limits[1], limits[2], by = 1))
    plot_df$y <- cdf(d, plot_df$x)

    out_plot <- ggplot2::ggplot(data = plot_df,
           aes(x = x, y = y)) +
      ggplot2::geom_bar(stat = 'identity', width = 1,
                        aes(color = I("black"),
                            fill = I("grey50"))) +
      plot_theme()
  }

  if(class(d)[1] %in% c('Beta', 'Cauchy', 'ChiSquare', 'Exponential',
                        'FisherF', 'Gamma', 'Logistic', 'LogNormal',
                        'Normal', 'StudentsT', 'Tukey', 'Uniform', 'Weibull')){
    plot_df <- data.frame(x = seq(limits[1], limits[2], by = 0.001))
    plot_df$y <- cdf(d, plot_df$x)

    out_plot <- ggplot2::ggplot(data = plot_df,
                    aes(x = x, y = y)) +
      ggplot2::geom_line() +
      plot_theme()
  }

  return(out_plot)

}


#' Plot the CDF of a distribution
#'
#'
#'
#' @export
plot_pdf <- function(d, limits = NULL, p = 0.001,
                     plot_theme = ggplot2::theme_minimal, ...){

  if(!"ggplot2" %in% loadedNamespaces())
    stop("You must load ggplot2 for this function to work.")


  if(is.null(limits))
    limits <- support(d)

  if(limits[1] == -Inf){
    limits[1] <- quantile(d, p = p)
  }

  if(limits[2] == Inf){
    limits[2] <- quantile(d, p = 1-p)
  }

  if(class(d)[1] %in% c('Bernoulli', 'Binomial', 'Geometric', 'HyperGeometric',
                        'NegativeBinomial', 'Poisson')){
    plot_df <- data.frame(x = seq(limits[1], limits[2], by = 1))
    plot_df$y <- pdf(d, plot_df$x)

    out_plot <- ggplot(data = plot_df,
                    aes(x = x, y = y)) +
      geom_bar(stat = 'identity', width = 1,
                        aes(color = I("black"),
                            fill = I("grey50"))) +
      plot_theme()
  }

  if(class(d)[1] %in% c('Beta', 'Cauchy', 'ChiSquare', 'Exponential',
                        'FisherF', 'Gamma', 'Logistic', 'LogNormal',
                        'Normal', 'StudentsT', 'Tukey', 'Uniform', 'Weibull')){
    plot_df <- data.frame(x = seq(limits[1], limits[2], by = 0.001))
    plot_df$y <- pdf(d, plot_df$x)

    out_plot <- ggplot(data = plot_df,
                    aes(x = x, y = y)) +
      geom_line() +
      plot_theme()
  }

  return(out_plot)

}

