context("test-prodist")

## simulate data
set.seed(0)
d <- data.frame(x = runif(10, -1, 1))
d$y <- rpois(10, exp(-0.5 + 1.5 * d$x))
nd <- data.frame(x = 0)

test_that("prodist() works for lm", {
  m_lm <- lm(y ~ x, data = d)
  expect_equal(prodist(m_lm), Normal(mu = fitted(m_lm), sigma = summary(m_lm)$sigma))
  expect_equal(prodist(m_lm, newdata = nd), Normal(mu = predict(m_lm, newdata = nd), sigma = summary(m_lm)$sigma))
})

test_that("prodist() works for glm", {
  m_glm <- glm(y ~ x, data = d, family = poisson)
  expect_equal(prodist(m_glm), Poisson(lambda = fitted(m_glm)))
  expect_equal(prodist(m_glm, newdata = nd), Poisson(lambda = predict(m_glm, newdata = nd, type = "response")))
})

test_that("prodist() works for arima", {
  m_ar <- arima(d$y, order = c(1, 0, 0))
  p_ar <- predict(m_ar, n.ahead = 3)
  expect_equal(prodist(m_ar, n.ahead = 3),
    Normal(mu = setNames(as.numeric(p_ar$pred), as.character(time(p_ar$pred))), sigma = as.numeric(p_ar$se)))
})
