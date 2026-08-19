# Doubly robust scores ---------------------------------------------------------

.validate_nuisance_components <- function(nuisance, n) {
  required <- c(
    "e",
    "propensity2",
    "mu_hat11",
    "mu_hat10",
    "mu_hat01",
    "mu_hat00",
    "m1",
    "m0"
  )

  missing <- setdiff(required, names(nuisance))

  if (length(missing) > 0) {
    .abort(sprintf(
      "`nuisance` is missing required components: %s.",
      paste(missing, collapse = ", ")
    ))
  }

  for (nm in setdiff(required, "propensity2")) {
    if (length(nuisance[[nm]]) != n) {
      .abort(sprintf(
        "`nuisance$%s` must have one value per observation.",
        nm
      ))
    }
  }

  invisible(TRUE)
}

#' Compute doubly robust welfare scores.
#'
#' Constructs the score components used for Fair Policy Targeting, including
#' group-specific baseline terms, treatment-effect scores, and counterfactual
#' conditional means for envy-based fairness criteria.
#'
#' @param outcome Numeric outcome vector.
#' @param treatment Binary treatment vector.
#' @param sensitive Binary sensitive-attribute vector.
#' @param nuisance Object returned by `estimate_nuisance()`.
#' @param e Optional treatment propensity scores, used only if `nuisance` is
#'   `NULL`.
#' @param m Optional nuisance data frame or list with `m1`, `m0`, `mu_hat11`,
#'   `mu_hat10`, `mu_hat01`, and `mu_hat00`.
#' @param p_s1 Optional probability that the sensitive attribute equals 1.
#' @param cost_treatment Treatment cost subtracted from treated outcomes.
#' @param trim Propensity-score trimming bound.
#'
#' @return An object of class `"fairtargeting_scores"`.
#' @export
compute_dr_scores <- function(outcome,
                              treatment,
                              sensitive,
                              nuisance = NULL,
                              e = NULL,
                              m = NULL,
                              p_s1 = NULL,
                              cost_treatment = 0,
                              trim = 1e-3) {
  Y <- as.numeric(outcome)
  D <- .as_binary(treatment, "treatment")
  S <- .as_binary(sensitive, "sensitive")
  n <- length(Y)

  .validate_original_vectors(Y, D, S)

  if (!is.numeric(cost_treatment) ||
      length(cost_treatment) != 1 ||
      is.na(cost_treatment)) {
    .abort("`cost_treatment` must be a single numeric value.")
  }

  if (!is.null(nuisance)) {
    if (!inherits(nuisance, "fairtargeting_nuisance")) {
      .abort("`nuisance` must be an object returned by `estimate_nuisance()`.")
    }

    .validate_nuisance_components(nuisance, n)

    e <- nuisance$e
    p_s1 <- nuisance$propensity2
    m <- nuisance$m
    trim <- nuisance$trim %||% trim
  }

  if (is.null(e) || is.null(m)) {
    .abort("Supply either `nuisance` or both `e` and `m`.")
  }

  e <- .clip_prob(e, eps = trim)

  if (length(e) != n) {
    .abort("`e` must have one value per observation.")
  }

  if (is.null(p_s1)) {
    p_s1 <- mean(S)
  }

  propensity2 <- .clip_prob(p_s1, eps = trim)
  p_s0 <- 1 - propensity2

  m <- as.data.frame(m, check.names = FALSE)

  required <- c(
    "m1",
    "m0",
    "mu_hat11",
    "mu_hat10",
    "mu_hat01",
    "mu_hat00"
  )

  missing <- setdiff(required, names(m))

  if (length(missing) > 0) {
    .abort(sprintf(
      "`m` is missing required components: %s.",
      paste(missing, collapse = ", ")
    ))
  }

  for (nm in required) {
    if (length(m[[nm]]) != n) {
      .abort(sprintf("`m$%s` must have one value per observation.", nm))
    }

    m[[nm]] <- as.numeric(m[[nm]])
  }

  m1 <- m$m1
  m0 <- m$m0

  G_i1 <- (Y - m0) * (1 - D) * S / ((1 - e) * propensity2) +
    m0 * S / propensity2

  G_i2 <- (Y - cost_treatment - m1) * D * S / (e * propensity2) +
    m1 * S / propensity2

  g_i_S <- G_i2 - G_i1

  G_i12 <- (Y - m0) * (1 - D) * (1 - S) /
    ((1 - e) * (1 - propensity2)) +
    m0 * (1 - S) / (1 - propensity2)

  G_i22 <- (Y - cost_treatment - m1) * D * (1 - S) /
    (e * (1 - propensity2)) +
    m1 * (1 - S) / (1 - propensity2)

  g_i_S2 <- G_i22 - G_i12
  all_g_i <- g_i_S + g_i_S2

  structure(
    list(
      n = n,
      Y = Y,
      D = D,
      S = S,
      outcome = Y,
      treatment = D,
      sensitive = S,
      propensity1 = e,
      e = e,
      propensity2 = propensity2,
      p_s1 = propensity2,
      p_s0 = p_s0,
      m1 = m1,
      m0 = m0,
      mu_hat11 = m$mu_hat11,
      mu_hat10 = m$mu_hat10,
      mu_hat01 = m$mu_hat01,
      mu_hat00 = m$mu_hat00,
      G_i1 = as.numeric(G_i1),
      G_i2 = as.numeric(G_i2),
      G_i12 = as.numeric(G_i12),
      G_i22 = as.numeric(G_i22),
      g_i_S = as.numeric(g_i_S),
      g_i_S2 = as.numeric(g_i_S2),
      all_g_i = as.numeric(all_g_i),
      tau_s1 = as.numeric(g_i_S),
      tau_s0 = as.numeric(g_i_S2),
      cost_treatment = cost_treatment,
      trim = trim
    ),
    class = "fairtargeting_scores"
  )
}

#' Estimate group-specific welfare for a policy.
#'
#' @param scores Object returned by `compute_dr_scores()`.
#' @param policy_s1 Numeric vector of policy values evaluated with the sensitive
#'   attribute set to 1.
#' @param policy_s0 Numeric vector of policy values evaluated with the sensitive
#'   attribute set to 0. If `NULL`, `policy_s1` is reused.
#' @param include_baseline Logical. If `TRUE`, include the no-treatment baseline
#'   terms in the welfare calculation.
#' @param scale Scale for scalar output. `"sum"` returns sample-sum welfare and
#'   `"mean"` divides by the number of observations.
#'
#' @return Named numeric vector `c(W1, W0)`.
#' @export
estimate_group_welfare <- function(scores,
                                   policy_s1,
                                   policy_s0 = NULL,
                                   include_baseline = FALSE,
                                   scale = c("sum", "mean")) {
  if (!inherits(scores, "fairtargeting_scores")) {
    .abort("`scores` must be produced by `compute_dr_scores()`.")
  }

  scale <- match.arg(scale)

  policy_s1 <- as.numeric(policy_s1)
  policy_s0 <- if (is.null(policy_s0)) policy_s1 else as.numeric(policy_s0)

  if (length(policy_s1) != scores$n || length(policy_s0) != scores$n) {
    .abort(
      "Policy vectors must have length equal to the number of observations in `scores`."
    )
  }

  if (any(!is.finite(policy_s1)) || any(!is.finite(policy_s0))) {
    .abort("Policy vectors must be finite.")
  }

  w1 <- sum(scores$g_i_S * policy_s1)
  w0 <- sum(scores$g_i_S2 * policy_s0)

  if (isTRUE(include_baseline)) {
    w1 <- w1 + sum(scores$G_i1)
    w0 <- w0 + sum(scores$G_i12)
  }

  if (scale == "mean") {
    w1 <- w1 / scores$n
    w0 <- w0 / scores$n
  }

  c(W1 = w1, W0 = w0)
}

# Fairness criteria ------------------------------------------------------------

#' Compute Fair Policy Targeting unfairness criteria.
#'
#' Computes unfairness criteria for fitted or candidate policy values using the
#' score objects produced by `compute_dr_scores()`.
#'
#' @param scores Object returned by `compute_dr_scores()`.
#' @param policy_s1 Numeric vector of policy values evaluated with the sensitive
#'   attribute set to 1.
#' @param policy_s0 Numeric vector of policy values evaluated with the sensitive
#'   attribute set to 0. If `NULL`, `policy_s1` is reused.
#' @param distance Fairness criterion. One of `"welfare"`,
#'   `"relative_welfare"`, `"parity"`, or `"envy"`.
#' @param fairness Optional alternative name for `distance`.
#' @param two_directions Logical. If `TRUE`, welfare, relative-welfare, and
#'   parity unfairness are computed as absolute two-sided disparities. If
#'   `FALSE`, signed disparities are used.
#' @param scale Scale for scalar output. `"sum"` returns sample-sum criteria;
#'   `"mean"` divides scalar criteria by the number of observations.
#'
#' @return A list of class `"fairtargeting_fairness"` containing the selected
#'   unfairness criterion and related group-specific components.
#' @export
compute_fairness <- function(scores,
                             policy_s1,
                             policy_s0 = NULL,
                             distance = NULL,
                             fairness = NULL,
                             two_directions = TRUE,
                             scale = c("sum", "mean")) {
  if (!inherits(scores, "fairtargeting_scores")) {
    .abort("`scores` must be produced by `compute_dr_scores()`.")
  }
  distance <- .normalize_distance(distance, fairness)
  scale <- match.arg(scale)
  policy_s1 <- as.numeric(policy_s1)
  policy_s0 <- if (is.null(policy_s0)) policy_s1 else as.numeric(policy_s0)
  if (length(policy_s1) != scores$n || length(policy_s0) != scores$n) {
    .abort("Policy vectors must have length equal to the number of observations in `scores`.")
  }
  if (any(!is.finite(policy_s1)) || any(!is.finite(policy_s0))) {
    .abort("Policy vectors must be finite.")
  }

  S <- scores$S
  p <- scores$propensity2

  welfare1 <- sum(scores$g_i_S * policy_s1) + sum(scores$G_i1)
  welfare0 <- sum(scores$g_i_S2 * policy_s0) + sum(scores$G_i12)
  welfare_signed <- welfare0 - welfare1
  welfare_unfairness <- if (isTRUE(two_directions)) abs(welfare1 - welfare0) else welfare_signed

  rel1 <- sum(scores$g_i_S * policy_s1)
  rel0 <- sum(scores$g_i_S2 * policy_s0)
  relative_signed <- rel0 - rel1
  relative_unfairness <- if (isTRUE(two_directions)) abs(rel1 - rel0) else relative_signed

  parity1 <- sum(S * policy_s1 / p)
  parity0 <- sum((1 - S) * policy_s0 / (1 - p))
  parity_signed <- parity0 - parity1
  parity_unfairness <- if (isTRUE(two_directions)) abs(parity1 - parity0) else parity_signed

  envy1 <- sum(scores$mu_hat01 * policy_s1 * S / p) +
    sum(scores$mu_hat00 * (1 - policy_s1) * S / p) - rel0
  envy0 <- sum(scores$mu_hat11 * policy_s0 * (1 - S) / (1 - p)) +
    sum(scores$mu_hat10 * (1 - policy_s0) * (1 - S) / (1 - p)) - rel1
  envy_unfairness <- envy1 + envy0

  divisor <- if (scale == "mean") scores$n else 1
  objective <- switch(
    distance,
    welfare = welfare_unfairness,
    relative_welfare = relative_unfairness,
    parity = parity_unfairness,
    envy = envy_unfairness
  ) / divisor

  structure(
    list(
      distance = distance,
      fairness = distance,
      objective = as.numeric(objective),
      two_directions = isTRUE(two_directions),
      scale = scale,
      welfare1 = welfare1 / divisor,
      welfare0 = welfare0 / divisor,
      welfare_signed = welfare_signed / divisor,
      welfare_unfairness = welfare_unfairness / divisor,
      relative_welfare1 = rel1 / divisor,
      relative_welfare0 = rel0 / divisor,
      relative_welfare_signed = relative_signed / divisor,
      relative_welfare_unfairness = relative_unfairness / divisor,
      parity1 = parity1 / divisor,
      parity0 = parity0 / divisor,
      parity_signed = parity_signed / divisor,
      parity_unfairness = parity_unfairness / divisor,
      envy_s1 = envy1 / divisor,
      envy_s0 = envy0 / divisor,
      envy_unfairness = envy_unfairness / divisor,
      prediction_parity = (parity1 - parity0) / scores$n,
      abs_prediction_parity = abs(parity1 - parity0) / scores$n,
      welfare_relative = c(W1 = rel1 / divisor, W0 = rel0 / divisor),
      welfare_level = c(W1 = welfare1 / divisor, W0 = welfare0 / divisor)
    ),
    class = "fairtargeting_fairness"
  )
}

.original_fairness_components <- function(scores, distance, two_directions = TRUE) {
  distance <- .normalize_distance(distance)
  S <- scores$S
  p <- scores$propensity2
  constant1 <- constant2 <- 0

  if (distance == "welfare") {
    component1 <- scores$g_i_S
    component2 <- scores$g_i_S2
    constant1 <- sum(scores$G_i1)
    constant2 <- sum(scores$G_i12)
  } else if (distance == "relative_welfare") {
    component1 <- scores$g_i_S
    component2 <- scores$g_i_S2
    distance <- "welfare"
  } else if (distance == "envy") {
    component1 <- scores$mu_hat01 * S / p - scores$mu_hat00 * S / p - scores$g_i_S
    component2 <- scores$mu_hat11 * (1 - S) / (1 - p) - scores$mu_hat10 * (1 - S) / (1 - p) - scores$g_i_S2
  } else if (distance == "parity") {
    component1 <- S / p
    component2 <- (1 - S) / (1 - p)
  }

  list(
    distance_for_model = distance,
    component1 = component1,
    component2 = component2,
    constant1 = constant1,
    constant2 = constant2,
    two_directions = isTRUE(two_directions)
  )
}
