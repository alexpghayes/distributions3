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
                     plot_theme = NULL){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("the ggplot2 package is needed. Please install it.", call. = FALSE)
  }
  if (is.null(plot_theme)) {
    plot_theme <- ggplot2::theme_minimal
  }
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
           ggplot2::aes_string(x = "x", y = "y")) +
      ggplot2::geom_bar(stat = 'identity', width = 1,
                        ggplot2::aes(color = I("black"),
                            fill = I("grey50"))) +
      plot_theme()
  }

  if(class(d)[1] %in% c('Beta', 'Cauchy', 'ChiSquare', 'Exponential',
                        'FisherF', 'Gamma', 'Logistic', 'LogNormal',
                        'Normal', 'StudentsT', 'Tukey', 'Uniform', 'Weibull')){
    plot_df <- data.frame(x = seq(limits[1], limits[2], length.out = 5000))
    plot_df$y <- cdf(d, plot_df$x)

    out_plot <- ggplot2::ggplot(data = plot_df,
                    ggplot2::aes_string(x = "x", y = "y")) +
      ggplot2::geom_line() +
      plot_theme()
  }

  return(out_plot)

}


#' Plot the PDF of a distribution
#'
#' A function to easily plot the PDF of a distribution using `ggplot2`. Requires `ggplot2` to be loaded.
#'
#' @param d A `distribution` object
#' @param limits either `NULL` (default) or a vector of length 2 that specifies the range of the x-axis
#' @param p If `limits` is `NULL`, the range of the x-axis will be the support of `d` if this is a bounded
#'   interval, or `quantile(d, p)` and `quantile(d, 1 - p)` if lower and/or upper limits of the support is
#'   `-Inf`/`Inf`.  Defaults to 0.001.
#' @param plot_theme specify theme of resulting plot using `ggplot2`. Default is `theme_minimal`
#'
#' @export
plot_pdf <- function(d, limits = NULL, p = 0.001,
                     plot_theme = NULL){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("the ggplot2 package is needed. Please install it.", call. = FALSE)
  }
  if (is.null(plot_theme)) {
    plot_theme <- ggplot2::theme_bw
  }
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

    out_plot <- ggplot2::ggplot(data = plot_df,
                       ggplot2::aes_string(x = "x", y = "y")) +
     ggplot2::geom_bar(stat = 'identity', width = 1,
               ggplot2::aes(color = I("black"),
                   fill = I("grey50"))) +
      #xlab("x") +
      plot_theme()
  }

  if(class(d)[1] %in% c('Beta', 'Cauchy', 'ChiSquare', 'Exponential',
                        'FisherF', 'Gamma', 'Logistic', 'LogNormal',
                        'Normal', 'StudentsT', 'Tukey', 'Uniform', 'Weibull')){
    plot_df <- data.frame(x = seq(limits[1], limits[2], length.out = 5000))
    plot_df$y <- pdf(d, plot_df$x)

    out_plot <- ggplot2::ggplot(data = plot_df,
                    ggplot2::aes_string(x = "x", y = "y")) +
      ggplot2::geom_line() +
      plot_theme()
  }

  out_plot$mapping$d <- class(d)[1]

  for(i in seq_along(d))
    out_plot$mapping[[paste0("param", i)]] <- d[[i]]

  return(out_plot)
}


#' @rdname geom_auc
#' @export
stat_auc <- function(mapping = NULL,
                     data = NULL,
                     geom = "auc",
                     position = "identity",
                     na.rm = FALSE,
                     show.legend = NA,
                     inherit.aes = TRUE,
                     from = -Inf,
                     to = Inf,
                     annotate = FALSE,
                     digits = 3,
                     ...) {

  ggplot2::layer(
    geom = geom,
    stat = StatAuc,
    data = data,
    mapping = mapping,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      from = from,
      to = to,
      annotate = annotate,
      digits = digits,
      ...
    )
  )
}


#' @rdname geom_auc
#' @format NULL
#' @usage NULL
#' @export
StatAuc <- ggplot2::ggproto("StatAuc", ggplot2::Stat,

  compute_group = function(data,
                           scales,
                           from = -Inf,
                           to = Inf,
                           annotate = FALSE,
                           digits = 3) {


    ## set data outside interval to zero
    data[data$x < from | data$x > to, 'y'] <- 0

    if (!isFALSE(annotate)) {
      ## compute auc for label based on original implementation
      n_params <- sum(grepl("param", names(data)))
      d <- do.call(eval(parse(text = paste0("function(...) ", data$d[1], "(...)"))),
                   args = lapply(paste0("param", 1:n_params), function(x) data[[x]][1]))

      data$label <- paste0(
        "P(", from, "< X < ", to, ") = ",
        round(cdf(d, to) - cdf(d, from), digits = digits)
      )
    } else {
      data$label <- NA
    }

    return(data)
  },

  required_aes = c("x", "y")
)


#' Fill out area under the curve for a plotted PDF
#'
#' @param from Left end-point of interval
#' @param to Right end-point of interval
#' @param annotate Should P() be added in the upper left corner as an annotation?
#' Works also with a colour character, e.g., "red".
#' @param digits Number of digits shown in annotation
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_area
#'
#' @export
#'
#' @examples
#'
#' X <- Normal()
#'
#' plot_pdf(X) + geom_auc(to = -0.645)
#' plot_pdf(X) + geom_auc(from = -0.645, to = 0.1, annotate = TRUE)
geom_auc <- function(mapping = NULL,
                     data = NULL,
                     stat = "auc",
                     position = "identity",
                     na.rm = FALSE,
                     show.legend = NA,
                     inherit.aes = TRUE,
                     from = -Inf,
                     to = Inf,
                     annotate = FALSE,
                     digits = 3,
                     ...) {

  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomAuc,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      from = from,
      to = to,
      annotate  = annotate,
      digits = digits,
      ...
    )
  )
}

#' @rdname geom_auc
#' @format NULL
#' @usage NULL
#' @export
GeomAuc <- ggplot2::ggproto("GeomAuc", ggplot2::GeomArea,

  required_aes = c("x", "y", "label"),

  default_aes = ggplot2::aes(colour = NA, fill = "grey40", size = 0.5, linetype = 1,
    alpha = NA),

  draw_panel = function(data,
                        panel_params,
                        coord,
                        na.rm = FALSE,
                        flipped_aes = FALSE,
                        outline.type = c("both", "upper", "lower", "full"),
                        parse = FALSE,
                        check_overlab = FALSE,
                        annotate = FALSE) {

    outline.type <- match.arg(outline.type)

    if (!isFALSE(annotate)) {

      annotate <- if (isTRUE(annotate)) "black" else annotate[1]

      data_annotate <- data.frame(
        x = -Inf,
        y = Inf,
        label = data$label[1],
        colour = annotate,
        size = 3.88 * data$size[1] * 2,
        angle = 0,
        hjust = -0.1,
        vjust = 2,
        alpha = NA,
        family = "",
        fontface = 1,
        lineheight = 1.2
      )

      grid::grobTree(
        ggplot2::GeomArea$draw_group(data,
                                 panel_params,
                                 coord,
                                 na.rm = na.rm,
                                 flipped_aes = flipped_aes,
                                 outline.type = outline.type),
        ggplot2::GeomText$draw_panel(data_annotate,
                                     panel_params,
                                     coord,
                                     parse = parse,
                                     na.rm = na.rm,
                                     check_overlap = check_overlab)
      )
    } else {
        ggplot2::GeomArea$draw_group(data,
                                 panel_params,
                                 coord,
                                 na.rm = na.rm,
                                 flipped_aes = flipped_aes,
                                 outline.type = outline.type)
    }

  }

)


