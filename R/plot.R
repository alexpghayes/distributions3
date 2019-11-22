#' Plot the p.m.f, p.d.f or c.d.f. of a univariate distribution
#'
#' Plot method for an object inheriting from class \code{"distribution"}.
#' By default the probability density function (p.d.f.), for a continuous
#' variable, or the probability mass function (p.m.f.), for a discrete
#' variable, is plotted.  The cumulative distribution function (c.d.f.)
#' will be plotted if \code{cdf = TRUE}.  Multiple functions are included
#' in the plot if any of the parameter vectors in \code{x} has length greater
#' than 1.  See the argument \code{all}.
#'
#' @param x an object of class \code{c("name", "distribution")}, where
#'   \code{"name"} is the name of the distribution.
#' @param cdf A logical scalar.  If \code{cdf = TRUE} then the cumulative
#'   distribution function (c.d.f.) is plotted.  Otherwise, the probability
#'   density function (p.d.f.), for a continuous variable, or the probability
#'   mass function (p.m.f.), for a discrete variable, is plotted.
#' @param p A numeric vector.  If \code{xlim} is not passed in \code{...}
#'   then \code{p} is the fallback option for setting the range of values
#'   over which the p.m.f, p.d.f. or c.d.f is plotted.  See **Details**.
#' @param len An integer scalar.  If \code{x} is a continuous distribution
#'   object then \code{len} is the number of values at which the p.d.f or
#'   c.d.f. is evaluated to produce the plot.  The larger \code{len} is the
#'   smoother is the curve.
#' @param all A logical scalar.  If \code{all = TRUE} then a separate
#'   distribution is plotted for all the combinations of parameter
#'   values present in the parameter vectors present in \code{x}.  These
#'   combinations are generated using \code{\link{expand.grid}}.  If
#'   \code{all = FALSE} then the number of distributions plotted is equal to
#'   the maximum of the lengths of these parameter vectors, with shorter
#'   vectors recycled to this length if necessary using \code{\link{rep_len}}.
#' @param legend_args A list of arguments to be passed to
#'   \code{\link[graphics]{legend}}.  In particular, the argument \code{x}
#'   (perhaps in conjunction with \code{legend_args$y}) can be used to set the
#'   position of the legend.  If \code{legend_args$x} is not supplied then
#'   \code{"bottomright"} is used if \code{cdf = TRUE} and \code{"topright"} if
#'   \code{cdf = FALSE}.
#' @param ...  Further arguments to be passed to \code{\link[graphics]{plot}},
#'   \code{\link[stats:ecdf]{plot.ecdf}} and \code{\link[graphics]{lines}},
#'   such as \code{xlim, ylim, xlab, ylab, main, lwd, lty, col, pch}.
#' @details If \code{xlim} is passed in \code{...} then this determines the
#'   range of values of the variable to be plotted on the horizontal axis.
#'   If \code{x} is a discrete distribution object then the values for which
#'   the p.m.f. or c.d.f. is plotted is the smallest set of consecutive
#'   integers that contains both components of \code{xlim}.  Otherwise,
#'   \code{xlim} is used directly.
#'
#'   If \code{xlim} is not passed in \code{...} then the range of values spans
#'   the support of the distribution, with the following proviso: if the
#'   lower (upper) endpoint of the distribution is \code{-Inf} (\code{Inf})
#'   then the lower (upper) limit of the plotting range is set to the
#'   \code{p[1]}\% (\code{p[2]}\%) quantile of the distribution.
#'
#'   If the name of \code{x} is a single upper case letter then that name is
#'   used to labels the axes of the plot.  Otherwise, \code{x} and
#'   \code{P(X = x)} or \code{f(x)} are used.
#'
#'   A legend is included only if at least one of the parameter vectors in
#'   \code{x} has length greater than 1.
#'
#'   Plots of c.d.f.s are produced using calls to
#'   \code{\link[stats]{approxfun}} and \code{\link[stats:ecdf]{plot.ecdf}}.
#' @return An object with the same class as \code{x}, in which the parameter
#'   vectors have been expanded to contain a parameter combination for each
#'   function plotted.
#' @examples
#' B <- Binomial(20, 0.7)
#' plot(B)
#' plot(B, cdf = TRUE)
#'
#' B2 <- Binomial(20, c(0.1, 0.5, 0.9))
#' plot(B2, legend_args = list(x = "top"))
#' x <- plot(B2, cdf = TRUE)
#' x$size
#' x$p
#'
#' X <- Poisson(2)
#' plot(X)
#' plot(X, cdf = TRUE)
#'
#' G <- Gamma(c(1,3), 1:2)
#' plot(G)
#' plot(G, all = TRUE)
#' plot(G, cdf = TRUE)
#'
#' C <- Cauchy()
#' plot(C, p = c(1, 99), len = 10000)
#' plot(C, cdf = TRUE, p = c(1, 99))
#' @export
plot.distribution <- function(x, cdf = FALSE, p = c(0.1, 99.9), len = 1000,
                              all = FALSE, legend_args = list(), ...) {
  if (!is_distribution(x)) {
    stop("use only with \"distribution\" objects")
  }
  # Extract the name of the distribution
  distn <- class(x)[1]
  # Supported distributions
  discrete <- c("Bernoulli", "Binomial", "Categorical", "Geometric",
                "HyperGeometric", "Multinomial", "NegativeBinomial", "Poisson")
  continuous <-  c("Beta", "Cauchy", "ChiSquare", "Erlang", "Exponential",
                   "Frechet", "FisherF", "GEV", "Gamma", "GP", "Gumbel",
                   "Laplace", "LogNormal", "Logistic", "Normal", "StudentsT",
                   "Tukey", "Uniform", "Weibull")
  # Check that the distribution is recognised
  if (!(distn %in% c(discrete, continuous))) {
    stop("This distribution is not supported")
  }
  # Indicator of whether or not the distribution is discrete
  x_is_discrete <- distn %in% discrete
  # The number of parameters and their names
  np <- length(x)
  par_names <- names(x)
  # Expands the vector(s) of parameters to create a parameter combination for
  # each function to be plotted.
  par_lengths <- sapply(x, length)
  if (all) {
    n_distns <- prod(par_lengths)
    xx <- as.list(expand.grid(x, KEEP.OUT.ATTRS = FALSE))
  } else {
    n_distns <- max(par_lengths)
    xx <- lapply(x, rep_len, n_distns)
  }
  class(xx) <- class(x)
  # Create a title for the plot
  # If n_distns = 1 then place the parameter values in the title
  # If n_distns > 1 then place the parameter values in the legend
  if (n_distns > 1) {
    my_main <- paste0(distn, " (")
    if (np > 1) {
      for (i in 1:(np - 1)) {
        my_main <- paste0(my_main, par_names[i], ", ")
      }
    }
    my_main <- paste0(my_main, par_names[np], ")")
  } else {
    my_main <- paste0(distn, " (")
    if (np > 1) {
      for (i in 1:(np - 1)) {
        my_main <- paste0(my_main, x[i], ", ")
      }
    }
    my_main <- paste0(my_main, x[np], ")")
  }
  if (cdf) {
    my_main <- paste(my_main, "c.d.f.")
  } else {
    if (x_is_discrete) {
      my_main <- paste(my_main, "p.m.f.")
    } else {
      my_main <- paste(my_main, "p.d.f.")
    }
  }
  # Extract user-supplied arguments for graphics::plot()
  user_args <- list(...)
  # If xlim is supplied then use it.
  # Otherwise, use default values (but no -Inf or Inf)
  if (is.null(user_args[['xlim']])) {
    my_xlim <- quantile(x, rep(0:1, each = n_distns))
    my_xlim <- ifelse(is.finite(my_xlim), my_xlim,
                      quantile(x, rep(p / 100, each = n_distns)))
    my_xlim <- range(my_xlim)
  } else {
    my_xlim <- user_args$xlim
  }
  # Set x and y axis labels. If the variable name is a single upper case letter
  # then use that name in the labels.  Otherwise, use the generic x and
  # P(X = x) or f(x).
  variable_name <- substitute(x)
  if (nchar(variable_name) == 1 & toupper(variable_name) == variable_name) {
    my_xlab <- tolower(substitute(x))
    if (cdf) {
      my_ylab <- paste0("F(", my_xlab, ")")
    } else {
      if (x_is_discrete) {
        my_ylab <- paste0("P(", substitute(x), " = ", my_xlab, ")")
      } else {
        my_ylab <- paste0("f(", my_xlab, ")")
      }
    }
  } else {
    my_xlab <- "x"
    my_ylab <- "P(X = x)"
  }
  # Function to create the legend text
  create_legend_text <- function(x, n_distns) {
    leg_text <- numeric(n_distns)
    for (i in 1:n_distns) {
      text_i <- lapply(x, "[[", i)
      leg_text[i] <- paste0(text_i, collapse = ", ")
    }
    return(leg_text)
  }
  # Plot function for discrete cdf with default arguments
  discrete_cdf_plot <- function(x, xvals, ..., ylim = c(0, 1), xlab = my_xlab,
                                ylab = my_ylab, lwd = 2, main = my_main,
                                pch = 16, col = 1:n_distns) {
    col <- rep_len(col, n_distns)
    probs <- matrix(NA, ncol = n_distns, nrow = length(xvals), byrow = TRUE)
    for (i in 1:n_distns) {
      new_xx <- lapply(xx, "[[", i)
      class(new_xx) <- class(x)
      probs[, i] <- cdf(new_xx, xvals)
    }
    rval <- stats::approxfun(xvals, probs[, 1], method = "constant",
                             yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    plot(rval, ylim = ylim, xlab = xlab, ylab = ylab, axes = FALSE, lwd = lwd,
         main = main, col = col[1], pch = pch, ...)
    if (n_distns > 1) {
      for (i in 2:n_distns) {
        rval <- stats::approxfun(xvals, probs[, i], method = "constant",
                                 yleft = 0, yright = 1, f = 0,
                                 ties = "ordered")
        class(rval) <- c("ecdf", "stepfun", class(rval))
        graphics::lines(rval, lwd = lwd, col = col[i], pch = pch, ...)
      }
      # Add a legend
      if (is.null(legend_args[['legend']])) {
        legend_args$legend <- create_legend_text(xx, n_distns)
      }
      if (is.null(legend_args[['title']])) {
        legend_args$title <- paste0(par_names, collapse = ", ")
      }
      if (is.null(legend_args[['col']])) {
        legend_args$col <- col
      }
      if (is.null(legend_args[['lwd']])) {
        legend_args$lwd <- lwd
      }
      if (is.null(legend_args[['lty']])) {
        legend_args$lty <- 1
      }
      if (is.null(legend_args[['pch']])) {
        legend_args$pch <- pch
      }
      do.call(graphics::legend, legend_args)
    }
    graphics::axis(1, at = xvals)
    graphics::axis(2)
  }
  # Plot function for discrete pmf with default arguments
  discrete_pmf_plot <- function(x, xvals, ..., xlab = my_xlab, ylab = my_ylab,
                                lwd = 2, main = my_main, pch = 16,
                                col = 1:n_distns) {
    yvals <- matrix(NA, ncol = n_distns, nrow = length(xvals), byrow = TRUE)
    for (i in 1:n_distns) {
      new_xx <- lapply(xx, "[[", i)
      class(new_xx) <- class(x)
      yvals[, i] <- pmf(new_xx, xvals)
    }
    graphics::matplot(xvals, yvals, type = "p", xlab = xlab, ylab = ylab,
                      axes = FALSE, lwd = lwd, main = main, pch = pch,
                      col = col, ...)
    graphics::axis(1, at = xvals)
    graphics::axis(2)
    # If n_distns > 1 then add a legend
    if (n_distns > 1) {
      if (is.null(legend_args[['legend']])) {
        legend_args$legend <- create_legend_text(xx, n_distns)
      }
      if (is.null(legend_args[['title']])) {
        legend_args$title <- paste0(par_names, collapse = ", ")
      }
      if (is.null(legend_args[['col']])) {
        legend_args$col <- col
      }
      if (is.null(legend_args[['pch']])) {
        legend_args$pch <- pch
      }
      do.call(graphics::legend, legend_args)
    }
  }
  # Plot function for discrete distributions with defaults
  continuous_plot <- function(x, xvals, ..., xlab = my_xlab, ylab = my_ylab,
                              lwd = 2, lty = 1, main = my_main,
                              col = 1:n_distns) {
    yvals <- matrix(NA, ncol = n_distns, nrow = length(xvals), byrow = TRUE)
    for (i in 1:n_distns) {
      new_xx <- lapply(xx, "[[", i)
      class(new_xx) <- class(x)
      if (cdf) {
        yvals[, i] <- cdf(new_xx, xvals)
      } else {
        yvals[, i] <- pdf(new_xx, xvals)
      }
    }
    graphics::matplot(xvals, yvals, type = "l", xlab = xlab, ylab = ylab,
                      axes = FALSE, lwd = lwd, lty = lty, main = main, ...)
    graphics::axis(1)
    graphics::axis(2)
    graphics::box(bty = "l")
    if (cdf) {
      graphics::abline(h = 0:1, lty = 2)
    }
    # If n_distns > 1 then add a legend
    if (n_distns > 1) {
      if (is.null(legend_args[['legend']])) {
        legend_args$legend <- create_legend_text(xx, n_distns)
      }
      if (is.null(legend_args[['title']])) {
        legend_args$title <- paste0(par_names, collapse = ", ")
      }
      if (is.null(legend_args[['col']])) {
        legend_args$col <- col
      }
      if (is.null(legend_args[['lwd']])) {
        legend_args$lwd <- lwd
      }
      if (is.null(legend_args[['lty']])) {
        legend_args$lty <- lty
      }
      do.call(graphics::legend, legend_args)
    }
  }
  # Start the legend. If legend$x hasn't been supplied then set defaults.
  if (is.null(legend_args[['x']])) {
    if (cdf) {
      legend_args[['x']] <- "bottomright"
    } else {
      legend_args[['x']] <- "topright"
    }
  }
  if (x_is_discrete) {
    # Assume integer support
    xvals <- floor(my_xlim[1]):ceiling(my_xlim[2])
    if (cdf) {
      discrete_cdf_plot(x, xvals, ...)
    } else {
      discrete_pmf_plot(x, xvals, ...)
    }
  } else {
    xvals <- seq(my_xlim[1], my_xlim[2], length.out = len)
    continuous_plot(x, xvals, ...)
  }
  return(invisible(xx))
}
