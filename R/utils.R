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
        if (length(at) == 1L) at <- rep.int(as.vector(at), n) ## as vector (case 1)
        if (length(at) != n) at <- rbind(at) ## as matrix (case 2)
      }
      if (is.matrix(at) && NROW(at) == 1L) { ## case 2
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
      } else { ## case 1
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
  class(x) <- "data.frame"
  x <- x[i, , drop = FALSE]
  class(x) <- cl
  return(x)
}

#' @export
format.distribution <- function(x, digits = getOption("digits") - 3L, ...) {
  cl <- class(x)[1L]
  n <- names(x)
  if (is.null(attr(x, "row.names"))) attr(x, "row.names") <- 1L:length(x)
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
  rval <- sapply(x, format, digits = digits)
  nok <- duplicated(rval)
  while (any(nok) && digits < 10) {
    digits <- digits + 1
    rval[nok] <- sapply(x[nok], format, digits = digits)
    nok <- duplicated(rval)
  }
  nok <- duplicated(rval) | duplicated(rval, fromLast = TRUE)
  if (any(nok)) rval[nok] <- make.unique(rval[nok], sep = "_")
  return(rval)
}

make_support <- function(min, max, drop = TRUE) {
  rval <- cbind(min = min, max = max)
  if (drop && NROW(rval) == 1L) unname(rval[1L, ]) else rval
}
