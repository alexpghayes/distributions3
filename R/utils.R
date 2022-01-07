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
#'
#' @examples
#'
#' N1 <- Normal()
#' plot_cdf(N1) 
#' 
#' N2 <- Normal(0, c(1, 2))
#' plot_cdf(N2)
#'
#' B1 <- Binomial(10, 0.2)
#' plot_cdf(B1) 
#' 
#' B2 <- Binomial(10, c(0.2, 0.5))
#' plot_cdf(B2)
plot_cdf <- function(d, limits = NULL, p = 0.001,
                     plot_theme = NULL){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("the ggplot2 package is needed. Please install it.", call. = FALSE)
  }
  if (is.null(plot_theme)) {
    plot_theme <- ggplot2::theme_minimal
  }

  ## get limits
  if(is.null(limits))
    limits[1] <- min(support(d, drop = FALSE)[, 1])
    limits[2] <- max(support(d, drop = FALSE)[, 2])

  if(limits[1] == -Inf){
    limits[1] <- min(quantile(d, p = p))
  }

  if(limits[2] == Inf){
    limits[2] <- max(quantile(d, p = 1 - p))
  }

  if(class(d)[1] %in% c('Bernoulli', 'Binomial', 'Geometric', 'HyperGeometric',
                        'NegativeBinomial', 'Poisson')){

    ## span x and compute cdf
    x <- seq(limits[1], limits[2], by = 1)
    y <- data.frame(t(cdf(d, x, drop = FALSE)))
    names(y) <- NULL

    plot_df <- list()
    for (i in seq_along(y)) {
      plot_df[[i]] <- data.frame("x" = x, "y" = y[i], "group" = i)
    }
    plot_df <- do.call("rbind", plot_df)

    ## set names
    if (!is.null(names(d))) {
      plot_df$group <- factor(plot_df$group, levels = 1L:length(d), labels = names(d))
    }

    ## actual plot
    out_plot <- ggplot2::ggplot(
        data = plot_df,
        ggplot2::aes_string(x = "x", y = "y")
      ) +
      ggplot2::geom_bar(
        stat = 'identity',
        width = 1,
        ggplot2::aes(
          color = I("black"),
          fill = I("grey50"))
        ) +
      ggplot2::facet_grid(group ~ .) +
      plot_theme()
  }

  if(class(d)[1] %in% c('Beta', 'Cauchy', 'ChiSquare', 'Exponential',
                        'FisherF', 'Gamma', 'Logistic', 'LogNormal',
                        'Normal', 'StudentsT', 'Tukey', 'Uniform', 'Weibull')){

    ## span x and compute cdf
    x <- seq(limits[1], limits[2], length.out = 5000)
    y <- data.frame(t(cdf(d, x, drop = FALSE)))
    names(y) <- NULL

    plot_df <- list()
    for (i in seq_along(y)) {
      plot_df[[i]] <- data.frame("x" = x, "y" = y[i], "group" = i)
    }
    plot_df <- do.call("rbind", plot_df)

    ## set names
    if (!is.null(names(d))) {
      plot_df$group <- factor(plot_df$group, levels = 1L:length(d), labels = names(d))
    }

    ## actual plot
    out_plot <- ggplot2::ggplot(
        data = plot_df,
        ggplot2::aes_string(x = "x", y = "y")
      ) +
      ggplot2::geom_line() +
      ggplot2::facet_grid(group ~ .) +
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
#'
#' @examples
#'
#' N1 <- Normal()
#' plot_pdf(N1) 
#' 
#' N2 <- Normal(0, c(1, 2))
#' plot_pdf(N2)
#'
#' B1 <- Binomial(10, 0.2)
#' plot_pdf(B1) 
#' 
#' B2 <- Binomial(10, c(0.2, 0.5))
#' plot_pdf(B2)
plot_pdf <- function(d, limits = NULL, p = 0.001,
                     plot_theme = NULL){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("the ggplot2 package is needed. Please install it.", call. = FALSE)
  }
  if (is.null(plot_theme)) {
    plot_theme <- ggplot2::theme_bw
  }

  ## get limits
  if(is.null(limits))
    limits[1] <- min(support(d, drop = FALSE)[, 1])
    limits[2] <- max(support(d, drop = FALSE)[, 2])

  if(limits[1] == -Inf){
    limits[1] <- min(quantile(d, p = p))
  }

  if(limits[2] == Inf){
    limits[2] <- max(quantile(d, p = 1 - p))
  }

  if(class(d)[1] %in% c('Bernoulli', 'Binomial', 'Geometric', 'HyperGeometric',
                        'NegativeBinomial', 'Poisson')){

    ## span x and compute pdf
    x <- seq(limits[1], limits[2], by = 1)
    y <- data.frame(t(pdf(d, x, drop = FALSE)))
    names(y) <- NULL

    plot_df <- list()
    for (i in seq_along(y)) {
      plot_df[[i]] <- data.frame("x" = x, "y" = y[i], "group" = i)
    }
    plot_df <- do.call("rbind", plot_df)

    ## set names
    if (!is.null(names(d))) {
      plot_df$group <- factor(plot_df$group, levels = 1L:length(d), labels = names(d))
    }

    ## actual plot
    out_plot <- ggplot2::ggplot(
        data = plot_df,
        ggplot2::aes_string(x = "x", y = "y")
      ) +
      ggplot2::geom_bar(
        stat = 'identity', 
        width = 1,
        ggplot2::aes(
          color = I("black"),
          fill = I("grey50"))
        ) +
      ggplot2::facet_grid(group ~ .) +
      plot_theme()
  }

  if(class(d)[1] %in% c('Beta', 'Cauchy', 'ChiSquare', 'Exponential',
                        'FisherF', 'Gamma', 'Logistic', 'LogNormal',
                        'Normal', 'StudentsT', 'Tukey', 'Uniform', 'Weibull')){

    ## span x and compute pdf
    x <- seq(limits[1], limits[2], length.out = 5000)
    y <- data.frame(t(pdf(d, x, drop = FALSE)))
    names(y) <- NULL

    plot_df <- list()
    for (i in seq_along(y)) {
      plot_df[[i]] <- data.frame("x" = x, "y" = y[i], "group" = i)
    }
    plot_df <- do.call("rbind", plot_df)

    ## set names
    if (!is.null(names(d))) {
      plot_df$group <- factor(plot_df$group, levels = 1L:length(d), labels = names(d))
    }

    ## actual plot
    out_plot <- ggplot2::ggplot(
        data = plot_df,
        ggplot2::aes_string(x = "x", y = "y")
      ) +
      ggplot2::geom_line() +
      ggplot2::facet_grid(group ~ .) +
      plot_theme()
  }

  ## save parameters for "stat_auc"
  out_plot$mapping$d <- class(d)[1]
  for(i in seq_along(colnames(d)))
    out_plot$mapping[[paste0("param", i)]] <- rep(as.matrix(d)[, i], each = length(x))

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
#' N1 <- Normal()
#' plot_pdf(N1) + geom_auc(to = -0.645)
#' plot_pdf(N1) + geom_auc(from = -0.645, to = 0.1, annotate = TRUE)
#' 
#' N2 <- Normal(0, c(1, 2))
#' plot_pdf(N2) + geom_auc(to = 0)
#' plot_pdf(N2) + geom_auc(from = -2, to = 2, annotate = TRUE)
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


apply_dpqr <- function(d,
                       FUN,
                       at = NULL,
                       drop = TRUE,
                       type_prefix = "x",
                       ...) {

  # -------------------------------------------------------------------
  # SANITY CHECKS
  # -------------------------------------------------------------------
  stopifnot(is_distribution(d))
  stopifnot(is.function(FUN))
  stopifnot(is.null(at) || is.numeric(at))
  stopifnot(is.logical(drop))
  stopifnot(is.character(type_prefix))

  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## * Is 'at' some kind of 'data'
  ## * or missing altogether ('none')
  attype <- if (is.null(at) || names(formals(FUN))[1L] == "d") {
    "none"
  } else {
    "data"
  }

  # -------------------------------------------------------------------
  # PREPARE OUTPUT CONDITIONAL ON `attype`
  # -------------------------------------------------------------------
  ## If 'at' is missing: prediction is just a transformation of the parameters
  if (attype == "none") {

    rval <- FUN(d, ...)
    if (is.null(dim(rval))) names(rval) <- rownames(d)

  ## Otherwise 'at' is 'data':
  ## set up a function that suitably expands 'at' (if necessary)
  ## and then evaluates it at the predicted parameters ('data')
  } else {
    FUN4 <- function(at, d, ...) {
      n <- NROW(d)
      if (!is.data.frame(at)) {
        if (length(at) == 1L) at <- rep.int(as.vector(at), n)  ## as vector (case 1)
        if (length(at) != n) at <- rbind(at)  ## as matrix (case 2)
      }
      if (is.matrix(at) && NROW(at) == 1L) {  ## case 2
        at <- matrix(rep(at, each = n), nrow = n)
        rv <- FUN(as.vector(at), d = d[rep(1L:n, ncol(at))], ...)
        rv <- matrix(rv, nrow = n)

        if (length(rv != 0L)) {
          rownames(rv) <- rownames(d)
            if (all(at[1L, ] == 1L)) {
              colnames(rv) <- paste(
                type_prefix,
                seq_along(at[1L, ]), 
                sep = "_"
              )
            } else {
              colnames(rv) <- paste(
                type_prefix,
                make_suffix(at[1L, ], digits = pmax(3L, getOption("digits") - 3L)),
                sep = "_"
              )
          }
        }
      } else {  ## case 1
        rv <- FUN(at, d = d, ...)
        names(rv) <- rownames(d)
      }
      return(rv)
    }

    rval <- FUN4(at, d = d, ...)
  }

  # -------------------------------------------------------------------
  # RETURN
  # -------------------------------------------------------------------
  ## return a data.frame (drop=FALSE) or should it be dropped
  ## to a vector if possible (drop=TRUE): default = TRUE
  if (drop) {
    if (!is.null(dim(rval)) && NROW(rval) == 1L) {
      rval <- unname(drop(rval))
    }
  } else {
    if (is.null(dim(rval))) {
      rval <- as.matrix(rval)
      if (ncol(rval) == 1L) {
        colnames(rval) <- paste(
          type_prefix,
          make_suffix(unique(at), digits = pmax(3L, getOption("digits") - 3L)),
          sep = "_"
        )
      }
    }
    if (!inherits(rval, "data.frame")) rval <- as.data.frame(rval)
  }

  return(rval)
}


## Methods ---------------------------------------------------------------------

#' @export
dim.distribution <- function(x) NULL

#' @export
length.distribution <- function(x) length(unclass(x)[[1L]])

#' @export
`[.distribution` <- function(x, i) {
  cl <- class(x)
  class(x) <- "data.frame"
  x <- x[i, , drop = FALSE]
  class(x) <- cl
  return(x)
} 

#' @export
format.distribution <- function(x, digits = getOption("digits") - 3L, ...) {
  cl <- class(x)[1L]
  n <- names(x)
  if(is.null(attr(x, "row.names"))) attr(x, "row.names") <- 1L:length(x)
  class(x) <- "data.frame"
  f <- sprintf("%s distribution (%s)", cl, apply(rbind(apply(as.matrix(x), 2L, format, digits = digits, ...)), 1L, function(p) paste(names(x), "=", as.vector(p), collapse = ", ")))
  setNames(f, n)
}

#' @export
print.distribution <- function(x, digits = getOption("digits") - 3L, ...) {
  print(format(x, digits = digits), ...)
  invisible(x)
}

#' @export
names.distribution <- function(x) {
  n <- attr(x, "row.names")
  if(identical(n, seq_along(x))) NULL else n
}

#' @export
`names<-.distribution` <- function(x, value) {
  cl <- class(x)
  class(x) <- "data.frame"
  rownames(x) <- value
  class(x) <- cl
  return(x)
}

#' @export
dimnames.distribution <- function(x) {
  list(
    attr(x, "rownames"),
    names(unclass(x))
  )
}

## (a) Data frame of parameters
as_data_frame_parameters <- function(x, ...) {
  class(x) <- "data.frame"
  return(x)
}

## (b) Data frame with distribution column
as_data_frame_column <- function(x, ...) {
  d <- data.frame(x = seq_along(x))
  rownames(d) <- names(x)
  d$x <- x
  names(d) <- deparse(substitute(x))
  return(d)
}

## Convention: "as.data.frame" uses version (b) and "as.matrix" uses version (a)

#' @export
as.data.frame.distribution <- as_data_frame_column

#' @export
as.matrix.distribution <- function(x, ...) {
  x <- as_data_frame_parameters(x, ...)
  as.matrix(x)
}

#' @export
as.list.distribution <- function(x, ...) {
  x <- as_data_frame_parameters(x, ...)
  as.list(x)
}

#' @export
c.distribution <- function(...) {
  x <- list(...)
  cl <- class(x[[1L]])
  x <- lapply(x, function(d) {
    class(d) <- "data.frame"
    d
  })
  x <- do.call("rbind", x)
  class(x) <- cl
  return(x)
}

#' @export
summary.distribution <- function(object, ...) {
  cat(sprintf("%s distribution:", class(object)[1L]), "\n")
  class(object) <- "data.frame"
  summary(object, ...)
}

make_suffix <- function(x, digits = 3) {
  rval <- sapply(x, format, digits = digits)
  nok <- duplicated(rval)
  while(any(nok) && digits < 10) {
    digits <- digits + 1
    rval[nok] <- sapply(x[nok], format, digits = digits)
    nok <- duplicated(rval)
  }
  nok <- duplicated(rval) | duplicated(rval, fromLast = TRUE)
  if(any(nok)) rval[nok] <- make.unique(rval[nok], sep = "_") 
  return(rval)
}

make_support <- function(min, max, drop = TRUE) {
  rval <- cbind(min = min, max = max)
  if (drop && NROW(rval) == 1L) unname(rval[1L, ]) else rval 
}

