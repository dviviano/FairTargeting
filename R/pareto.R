# Pareto frontier and second-stage optimization --------------------------------

.policy_counterfactuals <- function(beta,
                                    X,
                                    S,
                                    no_parity_constraint = FALSE,
                                    probabilistic = FALSE,
                                    threshold_probabilistic = FALSE,
                                    probs = NULL) {
  X <- as.matrix(X)
  beta <- as.matrix(beta)
  n <- nrow(X)

  if (isTRUE(no_parity_constraint)) {
    XX1 <- as.matrix(cbind(1, 1, X[, -1, drop = FALSE]))
    XX0 <- as.matrix(cbind(1, 0, X[, -1, drop = FALSE]))
  } else {
    XX1 <- as.matrix(cbind(1, X))
    XX0 <- as.matrix(cbind(1, X))
  }

  pi1 <- pi0 <- matrix(NA_real_, nrow = nrow(beta), ncol = n)
  xi1 <- xi0 <- matrix(NA_real_, nrow = nrow(beta), ncol = n)

  for (j in seq_len(nrow(beta))) {
    sc1 <- as.numeric(XX1 %*% beta[j, ])
    sc0 <- as.numeric(XX0 %*% beta[j, ])

    if (isTRUE(threshold_probabilistic)) {
      pp <- as.numeric(probs[j, ])

      xi1[j, ] <- as.numeric(sc1 > 0)
      xi0[j, ] <- as.numeric(sc0 > 0)
      pi1[j, ] <- ifelse(sc1 > 0, pp[2] + pp[1], pp[2])
      pi0[j, ] <- ifelse(sc0 > 0, pp[2] + pp[1], pp[2])
    } else if (isTRUE(probabilistic)) {
      pi1[j, ] <- sc1
      pi0[j, ] <- sc0
      xi1[j, ] <- sc1
      xi0[j, ] <- sc0
    } else {
      xi1[j, ] <- as.numeric(sc1 > 0)
      xi0[j, ] <- as.numeric(sc0 > 0)
      pi1[j, ] <- xi1[j, ]
      pi0[j, ] <- xi0[j, ]
    }
  }

  list(
    policy_s1 = pi1,
    policy_s0 = pi0,
    xi_s1 = xi1,
    xi_s0 = xi0
  )
}

.frontier_candidate_table <- function(frontier,
                                      scores,
                                      alpha_seq,
                                      X,
                                      no_parity_constraint = FALSE,
                                      probabilistic = FALSE,
                                      threshold_probabilistic = FALSE) {
  cf <- .policy_counterfactuals(
    beta = frontier$beta,
    X = X,
    S = scores$S,
    no_parity_constraint = no_parity_constraint,
    probabilistic = probabilistic,
    threshold_probabilistic = threshold_probabilistic,
    probs = frontier$probs %||% NULL
  )

  W1 <- as.numeric(cf$policy_s1 %*% scores$g_i_S)
  W0 <- as.numeric(cf$policy_s0 %*% scores$g_i_S2)
  W1_level <- W1 + sum(scores$G_i1)
  W0_level <- W0 + sum(scores$G_i12)

  data.frame(
    grid_index = seq_along(alpha_seq),
    alpha = alpha_seq,
    W1 = W1,
    W0 = W0,
    W1_level = W1_level,
    W0_level = W0_level,
    Wbar = as.numeric(alpha_seq * W1 + (1 - alpha_seq) * W0),
    objective = as.numeric(frontier$objective),
    check.names = FALSE
  )
}

#' Filter Pareto-efficient welfare points.
#'
#' @param x A two-column numeric matrix or data frame.
#' @param w1 Name or index of the first welfare column.
#' @param w0 Name or index of the second welfare column.
#' @param tolerance Numeric tolerance for dominance checks.
#'
#' @return A data frame containing only nondominated rows.
#' @export
filter_pareto_frontier <- function(x,
                                   w1 = 1,
                                   w0 = 2,
                                   tolerance = 0) {
  x <- as.data.frame(x)

  if (nrow(x) == 0L) {
    return(x)
  }

  W1 <- as.numeric(x[[w1]])
  W0 <- as.numeric(x[[w0]])

  if (any(!is.finite(W1)) || any(!is.finite(W0))) {
    .abort("Pareto welfare columns must contain only finite numeric values.")
  }

  dominated <- rep(FALSE, nrow(x))

  for (i in seq_len(nrow(x))) {
    for (j in seq_len(nrow(x))) {
      if (i == j) next

      weakly_better <- W1[j] >= W1[i] - tolerance &&
        W0[j] >= W0[i] - tolerance

      strictly_better <- W1[j] > W1[i] + tolerance ||
        W0[j] > W0[i] + tolerance

      if (weakly_better && strictly_better) {
        dominated[i] <- TRUE
        break
      }
    }
  }

  x[!dominated, , drop = FALSE]
}

#' Estimate a discretized Pareto frontier.
#'
#' Estimates the alpha-grid Pareto frontier for Fair Policy Targeting using the
#' selected policy class and Gurobi optimization backend. The resulting frontier
#' values are used by the second-stage fairness optimization in
#' `fit_fair_policy()`.
#'
#' @param data Training data.
#' @param scores Object returned by `compute_dr_scores()`.
#' @param sensitive Sensitive-attribute column name.
#' @param covariates Policy covariate column names.
#' @param policy_class One of `"maxscore"`, `"probabilistic"`, or
#'   `"threshold_probabilistic"`.
#' @param capacity Optional treatment capacity, supplied as a fraction or count.
#' @param alpha_grid Optional explicit alpha grid.
#' @param discretization Number of alpha-grid points if `alpha_grid` is `NULL`.
#' @param solver Optimization backend. The default and supported exact backend
#'   is `"gurobi"`.
#' @param gurobi_params Optional Gurobi parameter list.
#' @param no_parity_constraint Logical. Whether to exclude the sensitive
#'   attribute from the policy rule while retaining it for group-specific
#'   evaluation.
#' @param additional_fairness_constraint Logical. Whether to include an
#'   additional fairness constraint.
#' @param parity_constraint Constraint direction for the additional parity
#'   constraint.
#' @param max_time Gurobi time limit.
#' @param numcores Number of Gurobi threads.
#' @param tolerance_frontier Optimization tolerance for frontier estimation.
#' @param open_alpha_grid Logical. If `TRUE`, use an open alpha grid that
#'   excludes 0 and 1.
#' @param ... Additional arguments reserved for future extensions.
#'
#' @return An object of class `"fairpolicy_frontier"`.
#' @export
estimate_pareto_frontier <- function(data,
                                     scores,
                                     sensitive,
                                     covariates,
                                     policy_class = c(
                                       "maxscore",
                                       "probabilistic",
                                       "threshold_probabilistic"
                                     ),
                                     capacity = NULL,
                                     alpha_grid = NULL,
                                     discretization = NULL,
                                     solver = "gurobi",
                                     gurobi_params = NA,
                                     no_parity_constraint = FALSE,
                                     additional_fairness_constraint = FALSE,
                                     parity_constraint = ">=",
                                     max_time = 300,
                                     numcores = 1,
                                     tolerance_frontier = 10^(-3),
                                     open_alpha_grid = FALSE,
                                     ...) {
  if (!inherits(scores, "fairpolicy_scores")) {
    .abort("`scores` must be produced by `compute_dr_scores()`.")
  }

  policy_class <- match.arg(policy_class)
  solver <- match.arg(solver, choices = "gurobi")

  if (solver != "gurobi") {
    .abort(
      "`estimate_pareto_frontier()` uses `solver = \"gurobi\"` for exact frontier optimization."
    )
  }

  .require_gurobi()

  dat <- as.data.frame(data)
  dat[[sensitive]] <- .as_binary(dat[[sensitive]], sensitive)

  design_spec <- .make_design_spec(
    dat,
    sensitive = sensitive,
    covariates = covariates,
    include_sensitive_first = TRUE
  )

  X <- .design_matrix(dat, design_spec)

  .validate_original_vectors(scores$Y, scores$D, scores$S, X)

  if (!identical(as.integer(X[, 1]), as.integer(scores$S))) {
    .abort(
      "The policy design matrix must contain the sensitive attribute as its first column."
    )
  }

  if (is.null(alpha_grid)) {
    alpha_grid <- make_alpha_grid(
      n_grid = discretization,
      n = scores$n,
      open = open_alpha_grid
    )
  }

  alpha_grid <- as.numeric(alpha_grid)

  if (any(!is.finite(alpha_grid)) || any(alpha_grid < 0 | alpha_grid > 1)) {
    .abort("`alpha_grid` must contain finite values between 0 and 1.")
  }

  discretization <- length(alpha_grid)
  max_treated_units <- .resolve_capacity_original(capacity, scores$n)

  probabilistic <- identical(policy_class, "probabilistic")
  threshold_probabilistic <- identical(policy_class, "threshold_probabilistic")

  frontier <- .estimate_Pareto_frontier_core(
    Y = scores$Y,
    X = X,
    D = scores$D,
    S = scores$S,
    propensity1 = scores$propensity1,
    propensity2 = scores$propensity2,
    scale_Y = TRUE,
    discretization = discretization,
    cost_treatment = scores$cost_treatment,
    params = gurobi_params,
    max_treated_units = max_treated_units,
    numcores = numcores,
    maxtime = max_time,
    alpha_seq = alpha_grid,
    m1 = scores$m1,
    m0 = scores$m0,
    additional_fairness_constraint = additional_fairness_constraint,
    parity_constraint = parity_constraint,
    probabilistic = probabilistic,
    threshold_probabilistic = threshold_probabilistic,
    parallel = FALSE,
    tolerance = tolerance_frontier
  )

  candidates <- .frontier_candidate_table(
    frontier,
    scores,
    alpha_grid,
    X,
    no_parity_constraint = no_parity_constraint,
    probabilistic = probabilistic,
    threshold_probabilistic = threshold_probabilistic
  )

  efficient_candidates <- filter_pareto_frontier(
    candidates,
    w1 = "W1",
    w0 = "W0",
    tolerance = tolerance_frontier
  )

  structure(
    list(
      candidates = candidates,
      efficient_candidates = efficient_candidates,
      frontier = frontier,
      alpha_grid = alpha_grid,
      frontier_values = as.numeric(frontier$objective),
      n = scores$n,
      policy_class = policy_class,
      capacity = capacity,
      max_treated_units = max_treated_units,
      solver = solver,
      gurobi_params = gurobi_params,
      no_parity_constraint = isTRUE(no_parity_constraint),
      additional_fairness_constraint = isTRUE(additional_fairness_constraint),
      parity_constraint = parity_constraint,
      design_spec = design_spec,
      X = X,
      scores = scores
    ),
    class = "fairpolicy_frontier"
  )
}

# Pareto utilities -------------------------------------------------------------

.threshold_probability_from_score <- function(score, parameters) {
  parameters <- as.numeric(parameters)
  if (length(parameters) < 2L) {
    .abort("`parameters` must contain at least two probability values.")
  }

  ifelse(
    score > 0,
    sum(parameters[(length(parameters) - 1L):length(parameters)]),
    parameters[length(parameters)]
  )
}

#' Compute nondominated welfare points from two candidate sets.
#'
#' Stacks two two-column welfare matrices and removes rows that are strictly
#' dominated in both welfare dimensions. For weak-Pareto filtering with a
#' tolerance, use [filter_pareto_frontier()].
#'
#' @param welf1 Numeric matrix/data frame with two welfare columns.
#' @param welf0 Numeric matrix/data frame with two welfare columns.
#'
#' @return A numeric matrix containing the nondominated rows.
#' @export
compute_pareto_frontier <- function(welf1, welf0) {
  total_welfare <- rbind(as.matrix(welf1), as.matrix(welf0))

  if (ncol(total_welfare) != 2L) {
    .abort("`welf1` and `welf0` must have exactly two columns after row-binding.")
  }

  is_dominated <- rep(FALSE, nrow(total_welfare))

  for (i in seq_len(nrow(total_welfare))) {
    for (j in seq_len(nrow(total_welfare))) {
      if (i == j) next

      if (total_welfare[i, 1] < total_welfare[j, 1] &&
          total_welfare[i, 2] < total_welfare[j, 2]) {
        is_dominated[i] <- TRUE
        break
      }
    }
  }

  total_welfare[!is_dominated, , drop = FALSE]
}
