#' Extracting fitted or predicted probability distributions from models
#' 
#' Generic function with methods for various model classes for extracting
#' fitted (in-sample) or predicted (out-of-sample) probability `distributions3`
#' objects.
#' 
#' To facilitate making probabilistic forecasts based on regression and time
#' series model objects, the function \code{prodist} extracts fitted or
#' predicted probability \code{distribution} objects. Currently, methods are
#' provided for objects fitted by \code{\link[stats]{lm}},
#' \code{\link[stats]{glm}}, and \code{\link[stats]{arima}} in base R as
#' well as \code{glm.nb} from the \pkg{MASS} package. All methods essentially
#' proceed in two steps: First, the standard \code{\link[stats]{predict}}
#' method for these model objects is used to compute fitted (in-sample, default)
#' or predicted (out-of-sample) distribution parameters. Typically, this includes
#' the mean plus further parameters describing scale, dispersion, shape, etc.).
#' Second, the \code{distributions} objects are set up using the generator
#' functions from \pkg{distributions3}.
#' 
#' @aliases prodist.lm prodist.glm prodist.negbin prodist.Arima
#' 
#' @param object A model object.
#' @param ... Arguments passed on to methods, typically for calling the
#' underlying \code{\link[stats]{predict}} methods, e.g., \code{newdata} for
#' \code{\link[stats]{lm}} or \code{\link[stats]{glm}} objects or \code{n.ahead}
#' for \code{\link[stats]{arima}} objects.
#' 
#' @return An object inheriting from \code{distribution}.
#' 
#' @seealso \code{\link[stats]{predict}}, \code{\link[stats]{lm}},
#' \code{\link[stats]{glm}}, \code{\link[stats]{arima}}
#' 
#' @keywords distribution
#' 
#' @examples
#' 
#' ## Model: Linear regression
#' ## Fit: lm
#' ## Data: 1920s cars data
#' data("cars", package = "datasets")
#' 
#' ## Stopping distance (ft) explained by speed (mph)
#' reg <- lm(dist ~ speed, data = cars)
#' 
#' ## Extract fitted normal distributions (in-sample, with constant variance)
#' pd <- prodist(reg)
#' head(pd)
#' 
#' ## Compute corresponding medians and 90% interval
#' qd <- quantile(pd, c(0.05, 0.5, 0.95))
#' head(qd)
#' 
#' ## Visualize observations with predicted quantiles
#' plot(dist ~ speed, data = cars)
#' matplot(cars$speed, qd, add = TRUE, type = "l", col = 2, lty = 1)
#' 
#' 
#' ## Model: Poisson generalized linear model
#' ## Fit: glm
#' ## Data: FIFA 2018 World Cup data
#' data("FIFA2018", package = "distributions3")
#' 
#' ## Number of goals per team explained by ability differences
#' poisreg <- glm(goals ~ difference, data = FIFA2018, family = poisson)
#' summary(poisreg)
#' ## Interpretation: When the ratio of abilities increases by 1 percent,
#' ## the expected number of goals increases by around 0.4 percent
#' 
#' ## Predict fitted Poisson distributions for teams with equal ability (out-of-sample)
#' nd <- data.frame(difference = 0)
#' prodist(poisreg, newdata = nd)
#' 
#' ## Extract fitted Poisson distributions (in-sample)
#' pd <- prodist(poisreg)
#' head(pd)
#' 
#' ## Extract log-likelihood from model object
#' logLik(poisreg)
#' 
#' ## Replicate log-likelihood via distributions object
#' sum(log_pdf(pd, FIFA2018$goals))
#' log_likelihood(pd, FIFA2018$goals)
#' 
#' 
#' ## Model: Autoregressive integrated moving average model
#' ## Fit: arima
#' ## Data: Quarterly approval ratings of U.S. presidents (1945-1974)
#' data("presidents", package = "datasets")
#' 
#' ## ARMA(2,1) model
#' arma21 <- arima(presidents, order = c(2, 0, 1))
#' 
#' ## Extract predicted normal distributions for next two years
#' p <- prodist(arma21, n.ahead = 8)
#' p
#' 
#' ## Compute median (= mean) forecast along with 80% and 95% interval
#' quantile(p, c(0.5, 0.1, 0.9, 0.025, 0.975))
#' 
#' @export
prodist <- function(object, ...) {
  UseMethod("prodist")
}

#' @export
prodist.lm <- function(object, ...) {
  mu <- predict(object, ...)
  sigma <- sqrt(sum(residuals(object)^2)/df.residual(object))
  Normal(mu = mu, sigma = sigma)
}

#' @export
prodist.glm <- function(object, ..., dispersion = NULL) {
  ## fitted means
  mu <- predict(object, type = "response", ...)

  ## dispersion parameter phi
  phi <- dispersion
  if(is.null(phi)) { 
    phi <- if(object$family$family %in% c("poisson", "binomial")) {
      1
    } else {
      sum((object$weights * object$residuals^2)[object$weights > 0])/ df.residual(object)
    }
  }
  
  ## size for binomial distributions is available only in-sample:
  if(object$family$family == "binomial") {
    size <- if("newdata" %in% names(list(...))) 1L else object$prior.weights
  }
  
  ## set up distributions object (if possible)
  switch(object$family$family,
    "gaussian" = Normal(mu = mu, sigma = sqrt(phi)),
    "poisson" = Poisson(lambda = mu),
    "binomial" = Binomial(size = size, p = mu),
    "Gamma" = distributions3::Gamma(shape = 1/phi, rate = 1/(phi * mu)),
    "inverse.gaussian" = stop("inverse Gaussian distributions3 object not implemented yet"), ## FIXME: could use SuppDists for this
    "quasi" = stop("quasi family is not associated with a full probability distribution"),
    "quasibinomial" = stop("quasi-Poisson family is not associated with a full probability distribution"),
    "quasipoisson" = stop("quasi-Binomial family is not associated with a full probability distribution"),
    stop(sprintf("%s family not supported yet", object$family$family))
  )
}

#' @export
prodist.negbin <- function(object, ...) {
  mu <- predict(object, type = "response", ...)
  NegativeBinomial(mu = mu, size = object$theta)
}

#' @export
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



## Further examples requiring other packages ---------------
## 
## library("MASS")
## m2 <- glm.nb(goals ~ difference, data = FIFA2018)
## pd2 <- prodist(m2)
## 
## data("cns", package = "faraway")
## cns$CNS <- cns$An + cns$Sp + cns$Other
## cns_logit <- glm(cbind(CNS, NoCNS) ~ Water + Work, data = cns, family = binomial)
## logLik(cns_logit)
## log_likelihood(prodist(cns_logit), cns$CNS)
## 
## data("SwissLabor", package = "AER")
## swiss_logit <- glm(participation ~ . + I(age^2), data = SwissLabor, family = binomial)
## 
## 
## Multinomial and ordinal models --------------------------
## FIXME: Multinomial() distribution not yet vectorized
## 
## prodist.polr <- function(object, ...) {
##   p <- predict(object, type = "prob", ...)
##   Multinomial(size = 1, p = p)
## }
## 
## prodist.multinom <- function(object, ...) {
##   p <- predict(object, type = "prob", ...)
##   Multinomial(size = 1, p = p) ## FIXME predict size?
## }
##
## library("nnet")
## data("BankWages", package = "AER")
## bwmale <- subset(BankWages, gender == "male")
## bw_mnl <- multinom(job ~ education + minority, data = bwmale, trace = FALSE)
## 
## library("MASS")
## bw_olm <- polr(job ~ education + minority, data = bwmale, Hess = TRUE)
## summary(bw_olm)
