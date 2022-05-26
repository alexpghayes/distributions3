prodist <- function(object, ...) {
  UseMethod("prodist")
}

prodist.lm <- function(object, ...) {
  mu <- predict(object, ...)
  sigma <- sqrt(sum(residuals(object)^2)/df.residual(object))
  Normal(mu = mu, sigma = sigma)
}

prodist.glm <- function(object, ..., dispersion = NULL) {
  mu <- predict(object, type = "response", ...)
  phi <- dispersion
  if(is.null(phi)) { 
    phi <- if(object$family$family %in% c("poisson", "binomial")) {
      1
    } else {
      sum((object$weights * object$residuals^2)[object$weights > 0])/ df.residual(object)
    }
  }
  switch(object$family$family,
    "gaussian" = Normal(mu = mu, sigma = sqrt(phi)),
    "poisson" = Poisson(lambda = mu),
    "binomial" = Binomial(size = 1, p = mu), ## FIXME predict size?
    "Gamma" = distributions3::Gamma(shape = 1/phi, rate = 1/(phi * mu)),
    "inverse.gaussian" = stop("inverse Gaussian distributions3 object not implemented yet"),
    "quasi" = stop("quasi family is not associated with a full probability distribution"),
    "quasibinomial" = stop("quasi-Poisson family is not associated with a full probability distribution"),
    "quasipoisson" = stop("quasi-Binomial family is not associated with a full probability distribution"),
    stop(sprintf("%s family not supported yet", object$family$family))
  )
}

prodist.negbin <- function(object, ...) {
  mu <- predict(object, type = "response", ...)
  NegativeBinomial(mu = mu, size = object$theta)
}

prodist.polr <- function(object, ...) {
  p <- predict(object, type = "prob", ...)
  Multinomial(size = 1, p = p) ## FIXME predict size?
}

prodist.multinom <- function(object, ...) {
  p <- predict(object, type = "prob", ...)
  Multinomial(size = 1, p = p) ## FIXME predict size?
}

prodist.Arima <- function(object, ...) {
  p <- predict(object, ...)
  n <- stats::.preformat.ts(p$pred)
  n <- if(is.null(dim(n))) {
    format(stats::time(p$pred))
  } else {
    t(outer(rownames(n), colnames(n), paste))[t(n) != ""]
  }
  Normal(mu = setNames(as.numeric(p$pred), n), sigma = as.numeric(p$se))
}
