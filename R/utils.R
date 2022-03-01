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


# -------------------------------------------------------------------
# HELPER FUNCTION FOR VECTORIZATION OF DISTRIBUTION OBJECTS
# -------------------------------------------------------------------

apply_dpqr <- function(d,
                       FUN,
                       at,
                       drop = TRUE,
                       type = NULL,
                       ...) {

  ## sanity checks
  stopifnot(
    is_distribution(d),
    is.function(FUN),
    is.numeric(at),
    is.logical(drop),
    is.character(type)
  )

  ## basic properties:
  ## rows n = number of distributions
  ## columns k = number of arguments at || number of random replications
  rnam <- names(d)
  n <- length(d)
  k <- if (type == "random") as.numeric(at) else length(at)

  ## make names (not needed for random numbers)
  anam <- if (type == "random") NULL else make_suffix(at, digits = pmax(3L, getOption("digits") - 3L))

  ## handle different types of "at"
  if (type != "random") {
    if (k == 0L) {
      return(matrix(numeric(0L), nrow = n, ncol = 0L, dimnames = list(rnam, NULL)))
    } else if (k == 1L) {
      at <- rep.int(as.vector(at), n)
    } else if (k == n && is.null(dim(at))) {
      k <- 1L
    } else {
      at <- as.vector(at)
      k <- length(at)
    }
  }

  ## "at" labels
  cnam <- paste(substr(type, 1L, 1L), if (type == "random") seq_len(k) else anam, sep = "_")

  ## handle zero-length distribution vector
  if (n == 0L) return(matrix(numeric(0L), nrow = 0L, ncol = k, dimnames = list(NULL, cnam)))

  ## call FUN
  if(type == "random") {
    rval <- if (n == 1L) {
      FUN(at, d = d, ...)
    } else {
      replicate(at, FUN(n, d = d))
    }
  } else {
    rval <- if (k == 1L) {
      FUN(at, d = d, ...)
    } else {
      vapply(at, FUN, numeric(n), d = d, ...)
    }
  }

  ## handle dimensions
  if (k == 1L && drop) {
    rval <- as.vector(rval)
    names(rval) <- rnam
  } else if (length(anam) > k) {
    cnam <- type
    rval <- matrix(rval, nrow = n, ncol = k, dimnames = list(rnam, cnam))
  } else if (drop) {
    rval <- drop(matrix(rval, nrow = n, ncol = k, dimnames = list(rnam, cnam)))
  } else {
    rval <- matrix(rval, nrow = n, ncol = k, dimnames = list(rnam, cnam))
  }

  return(rval)
}


# -------------------------------------------------------------------
# METHODS FOR DISTRIBUTION OBJECTS
# -------------------------------------------------------------------

#' @export
dim.distribution <- function(x) NULL

#' @export
length.distribution <- function(x) length(unclass(x)[[1L]])

#' @export
`[.distribution` <- function(x, i) {
  cl <- class(x)
  nm <- names(x)
  class(x) <- "data.frame"
  x <- x[i, , drop = FALSE]
  class(x) <- cl
  if (is.null(nm)) attr(x, "row.names") <- seq_along(x)
  return(x)
}

#' @export
format.distribution <- function(x, digits = getOption("digits") - 3L, ...) {
  cl <- class(x)[1L]
  if (length(x) < 1L) {
    return(character(0))
  }
  n <- names(x)
  if (is.null(attr(x, "row.names"))) attr(x, "row.names") <- 1L:length(x)
  class(x) <- "data.frame"
  f <- sprintf("%s distribution (%s)", cl, apply(rbind(apply(as.matrix(x), 2L, format, digits = digits, ...)), 1L, function(p) paste(names(x), "=", as.vector(p), collapse = ", ")))
  setNames(f, n)
}

#' @export
print.distribution <- function(x, digits = getOption("digits") - 3L, ...) {
  if (length(x) < 1L) {
    cat(sprintf("%s distribution of length zero\n", class(x)[1L]))
  } else {
    print(format(x, digits = digits), ...)
  }
  invisible(x)
}

#' @export
names.distribution <- function(x) {
  n <- attr(x, "row.names")
  if (identical(n, seq_along(x))) NULL else n
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
  rval <- format(x, digits = digits, trim = TRUE, drop0trailing = TRUE)
  nok <- duplicated(rval)
  while (any(nok) && digits < 10) {
    digits <- digits + 1
    rval[nok] <- format(x[nok], digits = digits, trim = TRUE, drop0trailing = TRUE)
    nok <- duplicated(rval)
  }
  nok <- duplicated(rval) | duplicated(rval, fromLast = TRUE)
  if (any(nok)) rval[nok] <- make.unique(rval[nok], sep = "_")
  return(rval)
}

make_support <- function(min, max, d, drop = TRUE) {
  rval <- matrix(c(min, max), ncol = 2, dimnames = list(names(d), c("min", "max")))
  if (drop && NROW(rval) == 1L) rval[1L, , drop = TRUE] else rval
}

make_positive_integer <- function(n) {
  n <- if (length(n) > 1L) length(n) else suppressWarnings(try(as.integer(n), silent = TRUE))
  if (inherits(n, "try-error") || is.na(n) || n < 0L) {
    stop("Invalid arguments")
  }
  n
}
