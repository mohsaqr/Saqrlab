# ===========================================================================
# Survival simulation -- Cox proportional-hazards generative model
# Returns saqr_sim with $data and $params
# ===========================================================================


#' Simulate Survival Data (Cox Proportional Hazards)
#'
#' @description Generate right-censored survival data from a Cox
#'   proportional-hazards generative model. The hazard for subject \eqn{i} is
#'   \deqn{h(t \mid X_i) = h_0(t) \exp(X_i \beta),}
#'   where \eqn{h_0(t)} is a parametric baseline hazard (Weibull, exponential,
#'   or Gompertz). Event times are drawn by inverse-CDF sampling given each
#'   subject's linear predictor, then administrative/exponential censoring is
#'   applied so that roughly \code{censoring_rate} of observations are
#'   censored. Designed so that \code{coxph(Surv(time, status) ~ ., data)}
#'   recovers \code{betas} at large \code{n}.
#'
#' @param n Integer. Number of subjects. Default: 300.
#' @param n_covariates Integer. Number of covariates \eqn{p}. Default: 2.
#'   Ignored when \code{betas} is supplied (its length wins).
#' @param betas Numeric vector or NULL. True log hazard ratios, one per
#'   covariate. When NULL, a moderate alternating sequence is auto-generated
#'   (e.g. \code{c(0.5, -0.3, 0.5, -0.3, ...)}).
#' @param baseline Character. Baseline hazard family: \code{"weibull"}
#'   (default), \code{"exponential"}, or \code{"gompertz"}.
#' @param lambda Positive numeric. Baseline scale parameter. Default: 0.1.
#' @param shape Positive numeric. Baseline shape parameter. For Weibull this is
#'   the Weibull shape; for Gompertz the level parameter; ignored for
#'   exponential. Default: 1.
#' @param censoring_rate Numeric in \[0, 1). Target proportion of censored
#'   observations. Default: 0.3.
#' @param covariate_type Character. Covariate distribution: \code{"normal"}
#'   (default, standard normal) or \code{"binary"} (Bernoulli(0.5), coded 0/1).
#' @param seed Integer or NULL. Random seed.
#'
#' @return A \code{\link{saqr_sim}} object with:
#'   \describe{
#'     \item{\code{$data}}{data.frame with columns \code{time} (numeric event
#'       or censoring time), \code{status} (1 = event, 0 = censored), and
#'       \code{X1}...\code{Xp} (covariates).}
#'     \item{\code{$params}}{list with \code{betas} (true log hazard ratios),
#'       \code{baseline}, \code{lambda}, \code{shape}, \code{covariate_type},
#'       and \code{realized_censoring_rate}.}
#'   }
#'
#' @details Inverse-CDF event times use the closed-form survival inverse for
#'   each baseline. With \eqn{U \sim \mathrm{Unif}(0,1)} and linear predictor
#'   \eqn{\eta = X\beta}:
#'   Weibull \eqn{t = (-\log U / (\lambda e^{\eta}))^{1/\mathrm{shape}}};
#'   exponential \eqn{t = -\log U / (\lambda e^{\eta})};
#'   Gompertz \eqn{t = (1/\mathrm{shape}) \log(1 - \mathrm{shape}\log U /
#'   (\lambda e^{\eta}))}.
#'
#' @examples
#' r <- simulate_survival(n = 400, n_covariates = 2, seed = 1)
#' print(r)
#' r$params$betas
#' r$params$realized_censoring_rate
#' head(as.data.frame(r))
#'
#' # Binary covariates, exponential baseline, explicit effects
#' r2 <- simulate_survival(n = 500, betas = c(0.8, -0.5),
#'                         baseline = "exponential",
#'                         covariate_type = "binary", seed = 42)
#' summary(r2)
#'
#' @export
simulate_survival <- function(n = 300, n_covariates = 2, betas = NULL,
                              baseline = c("weibull", "exponential",
                                           "gompertz"),
                              lambda = 0.1, shape = 1,
                              censoring_rate = 0.3,
                              covariate_type = c("normal", "binary"),
                              seed = NULL) {
  baseline <- match.arg(baseline)
  covariate_type <- match.arg(covariate_type)

  stopifnot(
    is.numeric(n), length(n) == 1L, n >= 5L,
    is.numeric(n_covariates), length(n_covariates) == 1L, n_covariates >= 1L,
    is.numeric(lambda), length(lambda) == 1L, lambda > 0,
    is.numeric(shape), length(shape) == 1L, shape > 0,
    is.numeric(censoring_rate), length(censoring_rate) == 1L,
    censoring_rate >= 0, censoring_rate < 1
  )

  n <- as.integer(n)

  # Resolve betas (its length governs the number of covariates).
  if (is.null(betas)) {
    p <- as.integer(n_covariates)
    base_pattern <- c(0.5, -0.3)
    betas <- rep(base_pattern, length.out = p)
  } else {
    stopifnot(is.numeric(betas), length(betas) >= 1L)
    p <- length(betas)
  }

  if (!is.null(seed)) set.seed(as.integer(seed))

  # Covariate matrix X (n x p), vectorised across all subjects at once.
  x_values <- if (covariate_type == "binary") {
    stats::rbinom(n * p, size = 1L, prob = 0.5)
  } else {
    stats::rnorm(n * p)
  }
  x_mat <- matrix(x_values, nrow = n, ncol = p)
  colnames(x_mat) <- paste0("X", seq_len(p))

  # Linear predictor and per-subject hazard multiplier.
  eta <- drop(x_mat %*% betas)
  hazard_mult <- exp(eta)

  # Inverse-CDF event times from the chosen baseline (fully vectorised).
  u <- stats::runif(n)
  neg_log_u <- -log(u)
  event_time <- switch(
    baseline,
    weibull     = (neg_log_u / (lambda * hazard_mult))^(1 / shape),
    exponential = neg_log_u / (lambda * hazard_mult),
    gompertz    = (1 / shape) *
      log1p(shape * neg_log_u / (lambda * hazard_mult))
  )

  # Calibrate exponential censoring so ~censoring_rate are censored.
  # Censoring time C ~ Exp(rate = cens_rate_param); status = 1 if T <= C.
  if (censoring_rate <= 0) {
    cens_time <- rep(Inf, n)
  } else {
    cens_rate_param <- calibrate_censoring_rate(
      event_time = event_time, target = censoring_rate
    )
    cens_time <- stats::rexp(n, rate = cens_rate_param)
  }

  status <- as.integer(event_time <= cens_time)
  obs_time <- pmin(event_time, cens_time)
  realized_cens <- mean(status == 0L)

  df <- data.frame(time = obs_time, status = status)
  df <- cbind(df, as.data.frame(x_mat))

  saqr_sim(
    data   = df,
    params = list(
      betas                   = stats::setNames(betas, paste0("X", seq_len(p))),
      baseline                = baseline,
      lambda                  = lambda,
      shape                   = shape,
      covariate_type          = covariate_type,
      realized_censoring_rate = realized_cens
    ),
    type = "survival",
    seed = seed,
    call = match.call()
  )
}


# ---------------------------------------------------------------------------
# Calibrate the exponential censoring rate to hit a target censoring fraction.
# P(censored) = P(C < T) where C ~ Exp(rate). A 1-D root-find on `rate`
# against the empirical event times gives a well-calibrated value.
# ---------------------------------------------------------------------------
#' @keywords internal
#' @noRd
calibrate_censoring_rate <- function(event_time, target) {
  # Fraction censored as a function of the censoring rate, evaluated against
  # the realized event times: E[ P(C < t_i) ] = mean(1 - exp(-rate * t_i)).
  cens_fraction <- function(rate) mean(1 - exp(-rate * event_time))

  # Bracket the root: rate near 0 -> ~0 censored; large rate -> ~1 censored.
  lo <- 1e-8
  hi <- 1
  # Expand the upper bound until it over-shoots the target.
  while (cens_fraction(hi) < target && hi < 1e8) hi <- hi * 2

  if (cens_fraction(hi) < target) {
    return(hi)  # cannot reach target; return the largest tried.
  }
  stats::uniroot(
    function(rate) cens_fraction(rate) - target,
    lower = lo, upper = hi
  )$root
}
