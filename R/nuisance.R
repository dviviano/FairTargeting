# Nuisance estimation ----------------------------------------------------------

.family_object <- function(family) {
  switch(
    family,
    gaussian = stats::gaussian(),
    binomial = stats::binomial(),
    .abort("Unsupported family. Use 'gaussian' or 'binomial'.")
  )
}

.infer_outcome_family <- function(y) {
  uy <- sort(unique(stats::na.omit(y)))
  if (length(uy) <= 2 && all(uy %in% c(0, 1))) "binomial" else "gaussian"
}

.fit_glm_matrix <- function(x, y, family = "gaussian") {
  x <- as.matrix(x)
  y <- as.numeric(y)
  y_mean <- mean(y)

  if (length(y) < 2 || anyNA(y_mean)) {
    return(list(type = "mean", value = 0, family = family))
  }

  if (family == "binomial" && length(unique(y)) < 2) {
    return(list(type = "mean", value = .clip_prob(y_mean), family = family))
  }

  fit <- try(
    stats::glm.fit(
      x = .add_intercept(x),
      y = y,
      family = .family_object(family)
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error") || is.null(fit$coefficients)) {
    return(list(type = "mean", value = y_mean, family = family))
  }

  coefs <- as.numeric(fit$coefficients)
  coefs[!is.finite(coefs)] <- 0

  list(
    type = "glm_matrix",
    coefficients = coefs,
    family = family,
    fitted = fit,
    value = y_mean
  )
}

.predict_glm_matrix <- function(object, x) {
  x <- as.matrix(x)

  if (identical(object$type, "mean")) {
    return(rep(object$value, nrow(x)))
  }

  eta <- as.numeric(.add_intercept(x) %*% object$coefficients)

  if (object$family == "binomial") {
    return(stats::plogis(eta))
  }

  eta
}

.estimate_nuisance_glm <- function(dat,
                                   outcome,
                                   treatment,
                                   sensitive,
                                   covariates,
                                   cross_fit,
                                   n_folds,
                                   seed,
                                   outcome_family,
                                   trim,
                                   min_obs) {
  y <- as.numeric(dat[[outcome]])
  d <- dat[[treatment]]
  s <- dat[[sensitive]]
  n <- nrow(dat)

  if (!isTRUE(cross_fit)) {
    n_folds <- 1
  }

  n_folds <- min(as.integer(n_folds), n)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  fold_id <- if (n_folds == 1) {
    rep(1L, n)
  } else {
    sample(rep(seq_len(n_folds), length.out = n))
  }

  cov_spec <- .make_design_spec(
    dat,
    sensitive = sensitive,
    covariates = covariates,
    include_sensitive_first = FALSE
  )

  x_cov <- .design_matrix(dat, cov_spec)

  prop_spec <- .make_design_spec(
    dat,
    sensitive = sensitive,
    covariates = covariates,
    include_sensitive_first = TRUE
  )

  x_prop <- .design_matrix(dat, prop_spec)

  e_hat <- rep(NA_real_, n)
  mu_hat11 <- mu_hat10 <- mu_hat01 <- mu_hat00 <- rep(NA_real_, n)
  fold_models <- vector("list", n_folds)

  fit_cell <- function(train, dd, ss) {
    idx <- train & d == dd & s == ss

    if (sum(idx) < min_obs) {
      val <- if (sum(train & d == dd) >= min_obs) {
        mean(y[train & d == dd])
      } else {
        mean(y[train])
      }

      list(type = "mean", value = val, family = outcome_family)
    } else {
      .fit_glm_matrix(
        x_cov[idx, , drop = FALSE],
        y[idx],
        family = outcome_family
      )
    }
  }

  for (k in seq_len(n_folds)) {
    hold <- fold_id == k
    train <- !hold

    if (n_folds == 1) {
      train <- rep(TRUE, n)
      hold <- rep(TRUE, n)
    }

    prop_model <- .fit_glm_matrix(
      x_prop[train, , drop = FALSE],
      d[train],
      family = "binomial"
    )

    e_hat[hold] <- .predict_glm_matrix(
      prop_model,
      x_prop[hold, , drop = FALSE]
    )

    models <- list(
      mu_hat11 = fit_cell(train, 1, 1),
      mu_hat10 = fit_cell(train, 0, 1),
      mu_hat01 = fit_cell(train, 1, 0),
      mu_hat00 = fit_cell(train, 0, 0)
    )

    mu_hat11[hold] <- .predict_glm_matrix(
      models$mu_hat11,
      x_cov[hold, , drop = FALSE]
    )

    mu_hat10[hold] <- .predict_glm_matrix(
      models$mu_hat10,
      x_cov[hold, , drop = FALSE]
    )

    mu_hat01[hold] <- .predict_glm_matrix(
      models$mu_hat01,
      x_cov[hold, , drop = FALSE]
    )

    mu_hat00[hold] <- .predict_glm_matrix(
      models$mu_hat00,
      x_cov[hold, , drop = FALSE]
    )

    fold_models[[k]] <- list(
      propensity = prop_model,
      outcome = models
    )
  }

  list(
    e = .clip_prob(e_hat, eps = trim),
    mu_hat11 = mu_hat11,
    mu_hat10 = mu_hat10,
    mu_hat01 = mu_hat01,
    mu_hat00 = mu_hat00,
    fold_id = fold_id,
    models = list(
      folds = fold_models,
      covariate_design = cov_spec,
      propensity_design = prop_spec
    )
  )
}

.estimate_nuisance_glmnet <- function(dat,
                                      outcome,
                                      treatment,
                                      sensitive,
                                      covariates,
                                      cross_fit,
                                      n_folds,
                                      seed,
                                      outcome_family,
                                      trim) {
  .require_glmnet()

  y <- as.numeric(dat[[outcome]])
  d <- dat[[treatment]]
  s <- dat[[sensitive]]

  x_spec <- .make_design_spec(
    dat,
    sensitive = sensitive,
    covariates = covariates,
    include_sensitive_first = FALSE
  )

  X <- .design_matrix(dat, x_spec)
  X_int <- X

  X_reg <- cbind(
    X,
    d,
    d * X_int * s,
    d * X_int * (1 - s),
    s * X_int,
    d * X_int
  )

  prop_spec <- .make_design_spec(
    dat,
    sensitive = sensitive,
    covariates = covariates,
    include_sensitive_first = TRUE
  )

  X_prop <- .design_matrix(dat, prop_spec)

  n_folds <- min(as.integer(n_folds), nrow(dat))

  if (isTRUE(cross_fit) && n_folds < 2) {
    .abort("`n_folds` must be at least 2 when `cross_fit = TRUE`.")
  }

  if (isTRUE(cross_fit)) {
    mu <- cross_fitting_mean(
      Y = y,
      X_reg = X_reg,
      X = X,
      X_int = X_int,
      seeds = seed %||% 123,
      K = n_folds,
      family = outcome_family
    )

    e_hat <- cross_fitting_propensity(
      D = d,
      X = X_prop,
      seeds = seed %||% 123,
      K = n_folds
    )

    fold_id <- rep(NA_integer_, nrow(dat))
  } else {
    if (!is.null(seed)) {
      set.seed(seed)
    }

    y_fit <- if (outcome_family == "binomial") as.factor(y) else y

    mean_fit <- glmnet::glmnet(
      y = y_fit,
      x = as.matrix(X_reg),
      family = outcome_family,
      lambda = exp(-12)
    )

    pred_mean <- function(XX) {
      as.numeric(stats::predict(
        mean_fit,
        newx = as.matrix(XX),
        type = if (outcome_family == "binomial") "response" else "link"
      ))
    }

    mu <- cbind(
      m_hat11 = pred_mean(cbind(X, 1, X_int, X_int * 0, X_int, X_int)),
      m_hat10 = pred_mean(cbind(X, 0, 0 * X_int, X_int * 0, X_int, 0 * X_int)),
      m_hat01 = pred_mean(cbind(X, 1, X_int * 0, X_int, 0 * X_int, X_int)),
      m_hat00 = pred_mean(cbind(X, 0, 0 * X_int, X_int * 0, 0 * X_int, 0 * X_int))
    )

    prop_fit <- glmnet::cv.glmnet(
      y = as.factor(d),
      x = as.matrix(X_prop),
      family = "binomial"
    )

    e_hat <- as.numeric(stats::predict(
      prop_fit,
      newx = as.matrix(X_prop),
      type = "response"
    ))

    fold_id <- rep(1L, nrow(dat))
  }

  list(
    e = .clip_prob(as.numeric(e_hat), eps = trim),
    mu_hat11 = as.numeric(mu[, "m_hat11"]),
    mu_hat10 = as.numeric(mu[, "m_hat10"]),
    mu_hat01 = as.numeric(mu[, "m_hat01"]),
    mu_hat00 = as.numeric(mu[, "m_hat00"]),
    fold_id = fold_id,
    models = list(
      covariate_design = x_spec,
      propensity_design = prop_spec,
      mean_design_columns = colnames(X_reg)
    )
  )
}

#' Estimate nuisance functions for Fair Policy Targeting.
#'
#' Estimates propensity scores and conditional outcome means used to construct
#' doubly robust scores for Fair Policy Targeting. The default `method =
#' "glmnet"` uses penalized regression with cross-fitting. A base-`glm` method is
#' also available for smaller problems and debugging workflows.
#'
#' @param data Data frame.
#' @param outcome Outcome column name.
#' @param treatment Binary treatment column name.
#' @param sensitive Binary sensitive-attribute column name.
#' @param covariates Character vector of covariate column names.
#' @param cross_fit Logical. If `TRUE`, nuisance functions are estimated
#'   out-of-fold.
#' @param n_folds Number of folds used when `cross_fit = TRUE`.
#' @param seed Optional fold-generation seed.
#' @param outcome_family Outcome model family. One of `"auto"`, `"gaussian"`,
#'   or `"binomial"`.
#' @param trim Propensity-score trimming bound.
#' @param method Nuisance estimation method. One of `"glmnet"` or `"glm"`.
#' @param min_obs Minimum observations in a treatment-by-sensitive cell for the
#'   `glm` method.
#'
#' @return An object of class `"fairpolicy_nuisance"`.
#' @export
estimate_nuisance <- function(data,
                              outcome,
                              treatment,
                              sensitive,
                              covariates,
                              cross_fit = TRUE,
                              n_folds = 5,
                              seed = NULL,
                              outcome_family = c("auto", "gaussian", "binomial"),
                              trim = 1e-3,
                              method = c("glmnet", "glm"),
                              min_obs = 5) {
  .validate_required_data(data, outcome, treatment, sensitive, covariates)

  outcome_family <- match.arg(outcome_family)
  method <- match.arg(method)

  dat <- as.data.frame(data)
  y <- as.numeric(dat[[outcome]])

  dat[[treatment]] <- .as_binary(dat[[treatment]], treatment)
  dat[[sensitive]] <- .as_binary(dat[[sensitive]], sensitive)

  .validate_original_vectors(
    y,
    dat[[treatment]],
    dat[[sensitive]]
  )

  if (outcome_family == "auto") {
    outcome_family <- .infer_outcome_family(y)
  }

  if (!is.numeric(trim) || length(trim) != 1 || trim <= 0 || trim >= 0.5) {
    .abort("`trim` must be a single number in (0, 0.5).")
  }

  n_folds <- as.integer(n_folds)

  if (n_folds < 1) {
    .abort("`n_folds` must be at least 1.")
  }

  fit <- if (method == "glmnet") {
    .estimate_nuisance_glmnet(
      dat,
      outcome,
      treatment,
      sensitive,
      covariates,
      cross_fit,
      n_folds,
      seed,
      outcome_family,
      trim
    )
  } else {
    .estimate_nuisance_glm(
      dat,
      outcome,
      treatment,
      sensitive,
      covariates,
      cross_fit,
      n_folds,
      seed,
      outcome_family,
      trim,
      min_obs
    )
  }

  s <- dat[[sensitive]]
  p_s1 <- .clip_prob(mean(s), eps = trim)

  m1 <- ifelse(s == 1, fit$mu_hat11, fit$mu_hat01)
  m0 <- ifelse(s == 1, fit$mu_hat10, fit$mu_hat00)

  if (outcome_family == "binomial") {
    for (nm in c("mu_hat11", "mu_hat10", "mu_hat01", "mu_hat00")) {
      fit[[nm]] <- .clip_prob(fit[[nm]], eps = trim)
    }

    m1 <- .clip_prob(m1, eps = trim)
    m0 <- .clip_prob(m0, eps = trim)
  }

  structure(
    list(
      e = .clip_prob(fit$e, eps = trim),
      propensity1 = .clip_prob(fit$e, eps = trim),
      propensity2 = p_s1,
      p_s1 = p_s1,
      p_s0 = 1 - p_s1,
      mu_hat11 = as.numeric(fit$mu_hat11),
      mu_hat10 = as.numeric(fit$mu_hat10),
      mu_hat01 = as.numeric(fit$mu_hat01),
      mu_hat00 = as.numeric(fit$mu_hat00),
      m1 = as.numeric(m1),
      m0 = as.numeric(m0),
      m = data.frame(
        mu_hat11 = as.numeric(fit$mu_hat11),
        mu_hat10 = as.numeric(fit$mu_hat10),
        mu_hat01 = as.numeric(fit$mu_hat01),
        mu_hat00 = as.numeric(fit$mu_hat00),
        m1 = as.numeric(m1),
        m0 = as.numeric(m0),
        check.names = FALSE
      ),
      fold_id = fit$fold_id,
      cross_fit = isTRUE(cross_fit),
      outcome_family = outcome_family,
      trim = trim,
      method = method,
      models = fit$models,
      columns = list(
        outcome = outcome,
        treatment = treatment,
        sensitive = sensitive,
        covariates = covariates
      ),
      call = match.call()
    ),
    class = "fairpolicy_nuisance"
  )
}
