# Main estimator ---------------------------------------------------------------

.policy_from_solver_result <- function(result, design_spec, policy_class,
                                       no_parity_constraint = FALSE,
                                       threshold_probabilistic = FALSE,
                                       X_train = NULL) {
  beta <- as.numeric(result$beta)
  probs <- result$probs %||% NULL

  structure(
    list(
      policy_class = policy_class,
      method = "gurobi",
      design_spec = design_spec,
      beta = beta,
      probs = probs,
      no_parity_constraint = isTRUE(no_parity_constraint),
      probabilistic = identical(policy_class, "probabilistic"),
      threshold_probabilistic = isTRUE(threshold_probabilistic),
      fitted_policy = as.numeric(result$policies),
      selected_alpha = result$alpha,
      solver_result = result$results,
      X_train = X_train
    ),
    class = "fairtargeting_policy"
  )
}

.tree_policy_from_row <- function(row, design_spec, length_tree,
                                  no_parity_constraint, X_train = NULL) {
  num_variables <- 2^length_tree - 1
  leaves <- 2^length_tree
  tree_vector_length <- 2 * num_variables + 2 + leaves
  tree_vector <- as.numeric(unlist(row[seq_len(tree_vector_length)]))

  structure(
    list(
      policy_class = "tree",
      method = "exhaustive_tree",
      design_spec = design_spec,
      tree_vector = tree_vector,
      length_tree = length_tree,
      no_parity_constraint = isTRUE(no_parity_constraint),
      probabilistic = TRUE,
      threshold_probabilistic = FALSE,
      X_train = X_train,
      selected_alpha = tree_vector[2 * num_variables + 2]
    ),
    class = "fairtargeting_policy"
  )
}

.predict_training_counterfactuals <- function(policy, X) {
  if (identical(policy$policy_class, "tree")) {
    X <- as.matrix(X)

    if (isTRUE(policy$no_parity_constraint)) {
      X1 <- as.matrix(cbind(1, X[, -1, drop = FALSE]))
      X0 <- as.matrix(cbind(0, X[, -1, drop = FALSE]))
    } else {
      X1 <- X0 <- X
    }

    return(list(
      policy_s1 = as.numeric(
        compute_predictions_tree(
          policy$tree_vector,
          X1,
          length_tree = policy$length_tree
        )
      ),
      policy_s0 = as.numeric(
        compute_predictions_tree(
          policy$tree_vector,
          X0,
          length_tree = policy$length_tree
        )
      )
    ))
  }

  beta <- matrix(policy$beta, nrow = 1)

  cf <- .policy_counterfactuals(
    beta = beta,
    X = X,
    S = X[, 1],
    no_parity_constraint = policy$no_parity_constraint,
    probabilistic = policy$probabilistic,
    threshold_probabilistic = policy$threshold_probabilistic,
    probs = if (is.null(policy$probs)) NULL else matrix(policy$probs, nrow = 1)
  )

  list(
    policy_s1 = as.numeric(cf$policy_s1[1, ]),
    policy_s0 = as.numeric(cf$policy_s0[1, ])
  )
}

#' Fit a Fair Policy Targeting rule.
#'
#' Estimates a Fair Policy Targeting rule using a two-stage optimization
#' procedure. First, the method estimates a discretized empirical Pareto
#' frontier over an alpha grid. Second, it solves a fairness optimization problem
#' subject to approximate Pareto-frontier constraints.
#'
#' The selected policy is not restricted to be one of the alpha-grid maximizers.
#' It may be any feasible policy in the model class that satisfies the
#' approximate Pareto-frontier constraints.
#'
#' @param data Data frame.
#' @param outcome Outcome column name.
#' @param treatment Binary treatment column name.
#' @param sensitive Binary sensitive-attribute column name. Value 1 denotes the
#'   sensitive group.
#' @param covariates Character vector of policy covariates. The policy design
#'   matrix is constructed with the sensitive attribute as the first column,
#'   followed by the listed covariates.
#' @param policy_class One of `"maxscore"`, `"probabilistic"`,
#'   `"threshold_probabilistic"`, or `"tree"`.
#' @param distance Fairness criterion. One of `"welfare"`,
#'   `"relative_welfare"`, `"parity"`, or `"envy"`.
#' @param fairness Optional alternative name for `distance`.
#' @param two_directions Logical. If `TRUE`, welfare, relative-welfare, and
#'   parity unfairness are computed as two-sided absolute disparities. If
#'   `FALSE`, signed disparities are used.
#' @param capacity Optional treatment capacity. Values from 0 to 1 are treated
#'   as fractions; values greater than 1 are treated as counts.
#' @param discretization Number of alpha-grid points. Defaults to
#'   `floor(sqrt(n))`.
#' @param n_grid Optional alternative name for `discretization`.
#' @param alpha_seq Optional explicit alpha grid. By default, the grid is
#'   `seq(0, 1, length.out = discretization)`.
#' @param lambda Approximate Pareto-frontier slack parameter.
#' @param solver Optimization backend. The default and supported exact backend
#'   is `"gurobi"`.
#' @param cross_fit Logical. Whether to cross-fit nuisance functions.
#' @param n_folds Number of folds for nuisance cross-fitting.
#' @param seed Optional random seed.
#' @param outcome_family Outcome model family. One of `"auto"`, `"gaussian"`,
#'   or `"binomial"`.
#' @param nuisance Optional nuisance object returned by `estimate_nuisance()`.
#' @param nuisance_method Nuisance estimation method. One of `"glmnet"` or
#'   `"glm"`.
#' @param trim Propensity trimming bound.
#' @param cost_treatment Optional treatment cost.
#' @param probabilistic Optional logical flag for probabilistic policy classes.
#'   If `NULL`, it is inferred from `policy_class`.
#' @param threshold_probabilistic Optional logical flag for
#'   threshold-probabilistic policy classes. If `NULL`, it is inferred from
#'   `policy_class`.
#' @param no_parity_constraint Logical. Whether to exclude the sensitive
#'   attribute from the policy rule while retaining it for group-specific
#'   evaluation.
#' @param noparity_constraint Optional alternative name for `no_parity_constraint`.
#' @param additional_fairness_constraint Logical. Whether to include an
#'   additional fairness constraint in the optimization problem.
#' @param parity_constraint Constraint direction for the additional parity
#'   constraint.
#' @param UnFairness_bound Optional bound for constrained fairness extensions.
#' @param max_time Gurobi time limit for the second-stage fairness optimization.
#' @param max_time_frontier Gurobi time limit for each alpha-frontier problem.
#' @param gurobi_params Optional Gurobi parameter list.
#' @param numcores Number of Gurobi threads.
#' @param tolerance_frontier Optimization tolerance for frontier problems.
#' @param tolerance_optimization Optimization tolerance for the second-stage
#'   problem.
#' @param open_alpha_grid Logical. If `TRUE`, use an open alpha grid that
#'   excludes 0 and 1.
#' @param tree_num_splits Number of split candidates for tree policies.
#' @param tree_length Tree depth parameter.
#' @param tree_epsilon_optimality Epsilon-optimality pruning value for tree
#'   search.
#' @param ... Reserved for future extensions.
#'
#' @return An S3 object of class `"fairtargeting"`.
#' @export
fit_fair_policy <- function(data,
                            outcome,
                            treatment,
                            sensitive,
                            covariates,
                            policy_class = c(
                              "maxscore",
                              "probabilistic",
                              "threshold_probabilistic",
                              "tree"
                            ),
                            distance = NULL,
                            fairness = NULL,
                            two_directions = TRUE,
                            capacity = NULL,
                            discretization = NULL,
                            n_grid = NULL,
                            alpha_seq = NULL,
                            lambda = 0,
                            solver = "gurobi",
                            cross_fit = TRUE,
                            n_folds = 5,
                            seed = NULL,
                            outcome_family = c("auto", "gaussian", "binomial"),
                            nuisance = NULL,
                            nuisance_method = c("glmnet", "glm"),
                            trim = 1e-3,
                            cost_treatment = 0,
                            probabilistic = NULL,
                            threshold_probabilistic = NULL,
                            no_parity_constraint = NULL,
                            noparity_constraint = NULL,
                            additional_fairness_constraint = FALSE,
                            parity_constraint = ">=",
                            UnFairness_bound = Inf,
                            max_time = 300,
                            max_time_frontier = max_time,
                            gurobi_params = NA,
                            numcores = 1,
                            tolerance_frontier = 10^(-3),
                            tolerance_optimization = 10^(-6),
                            open_alpha_grid = FALSE,
                            tree_num_splits = 4,
                            tree_length = 2,
                            tree_epsilon_optimality = 0.001,
                            ...) {
  policy_class <- match.arg(policy_class)
  solver <- match.arg(solver, choices = "gurobi")
  outcome_family <- match.arg(outcome_family)
  nuisance_method <- match.arg(nuisance_method)
  distance <- .normalize_distance(distance, fairness)
  no_parity_constraint <- .resolve_no_parity(
    no_parity_constraint,
    noparity_constraint
  )

  if (!is.numeric(lambda) || length(lambda) != 1 || is.na(lambda) || lambda < 0) {
    .abort("`lambda` must be a single non-negative number.")
  }
  if (!tree_length %in% 1:3) {
    .abort("`tree_length` must be 1, 2, or 3.")
  }

  .validate_required_data(data, outcome, treatment, sensitive, covariates)

  dat <- as.data.frame(data)
  dat[[treatment]] <- .as_binary(dat[[treatment]], treatment)
  dat[[sensitive]] <- .as_binary(dat[[sensitive]], sensitive)

  .validate_original_vectors(
    dat[[outcome]],
    dat[[treatment]],
    dat[[sensitive]]
  )

  if (solver != "gurobi") {
    .abort("`fit_fair_policy()` uses `solver = \"gurobi\"` for exact Fair Policy Targeting optimization.")
  }

  .require_gurobi()

  if (is.null(nuisance)) {
    nuisance <- estimate_nuisance(
      data = dat,
      outcome = outcome,
      treatment = treatment,
      sensitive = sensitive,
      covariates = covariates,
      cross_fit = cross_fit,
      n_folds = n_folds,
      seed = seed,
      outcome_family = outcome_family,
      trim = trim,
      method = nuisance_method
    )
  } else if (!inherits(nuisance, "fairtargeting_nuisance")) {
    .abort(
      "`nuisance` must be NULL or an object returned by `estimate_nuisance()`."
    )
  }

  scores <- compute_dr_scores(
    outcome = dat[[outcome]],
    treatment = dat[[treatment]],
    sensitive = dat[[sensitive]],
    nuisance = nuisance,
    cost_treatment = cost_treatment,
    trim = trim
  )

  if (is.null(discretization)) {
    discretization <- n_grid %||% floor(sqrt(scores$n))
  }

  if (is.null(alpha_seq)) {
    alpha_seq <- make_alpha_grid(
      n_grid = discretization,
      n = scores$n,
      open = open_alpha_grid
    )
  }

  discretization <- length(alpha_seq)
  max_treated_units <- .resolve_capacity_original(capacity, scores$n)

  design_spec <- .make_design_spec(
    dat,
    sensitive = sensitive,
    covariates = covariates,
    include_sensitive_first = TRUE
  )

  X <- .design_matrix(dat, design_spec)

  if (!identical(as.integer(X[, 1]), as.integer(scores$S))) {
    .abort(
      "The policy design matrix must contain the sensitive attribute as its first column."
    )
  }

  if (is.null(probabilistic)) {
    probabilistic <- identical(policy_class, "probabilistic")
  }

  if (is.null(threshold_probabilistic)) {
    threshold_probabilistic <- identical(policy_class, "threshold_probabilistic")
  }

  if (isTRUE(threshold_probabilistic)) {
    policy_class <- "threshold_probabilistic"
    probabilistic <- TRUE
  } else if (isTRUE(probabilistic)) {
    policy_class <- "probabilistic"
  }

  if (identical(policy_class, "tree")) {
    tree_results <- estimate_fairness_optimal_tree(
      Y = scores$Y,
      X = X,
      D = scores$D,
      S = scores$S,
      propensity1 = scores$propensity1,
      propensity2 = scores$propensity2,
      scale_Y = TRUE,
      discretization = discretization,
      cost_treatment = cost_treatment,
      params = gurobi_params,
      mu_hat11 = scores$mu_hat11,
      mu_hat01 = scores$mu_hat01,
      mu_hat00 = scores$mu_hat00,
      mu_hat10 = scores$mu_hat10,
      max_treated_units = max_treated_units,
      maxtime = max_time_frontier,
      alpha_seq = alpha_seq,
      noparity_constraint = no_parity_constraint,
      additional_fairness_constraint = additional_fairness_constraint,
      m0 = scores$m0,
      m1 = scores$m1,
      numcores = numcores,
      two_directions = two_directions,
      tolerance = tolerance_optimization,
      num_splits = tree_num_splits,
      length_tree = tree_length,
      epsilon_optimality = tree_epsilon_optimality,
      parallel = FALSE,
      parity_constraint = parity_constraint
    )

    selected_row <- switch(
      distance,
      parity = tree_results$best_parity,
      envy = tree_results$best_envy,
      relative_welfare = tree_results$best_relative_welfare,
      welfare = tree_results$best_welfare
    )

    policy <- .tree_policy_from_row(
      selected_row,
      design_spec,
      tree_length,
      no_parity_constraint,
      X_train = X
    )

    cf <- .predict_training_counterfactuals(policy, X)

    unfairness <- compute_fairness(
      scores,
      cf$policy_s1,
      cf$policy_s0,
      distance = distance,
      two_directions = two_directions
    )

    frontier <- structure(
      list(
        tree_results = tree_results,
        alpha_grid = alpha_seq,
        policy_class = "tree",
        design_spec = design_spec,
        X = X,
        scores = scores
      ),
      class = "fairtargeting_frontier"
    )

    selected_alpha <- policy$selected_alpha
    selected_grid_index <- which.min(abs(alpha_seq - selected_alpha))
  } else {
    frontier <- estimate_pareto_frontier(
      data = dat,
      scores = scores,
      sensitive = sensitive,
      covariates = covariates,
      policy_class = policy_class,
      capacity = capacity,
      alpha_grid = alpha_seq,
      discretization = discretization,
      solver = "gurobi",
      gurobi_params = gurobi_params,
      no_parity_constraint = no_parity_constraint,
      additional_fairness_constraint = additional_fairness_constraint,
      parity_constraint = parity_constraint,
      max_time = max_time_frontier,
      numcores = numcores,
      tolerance_frontier = tolerance_frontier,
      open_alpha_grid = open_alpha_grid
    )

    result <- Est_objective_estimandMaxscore(
      Y = scores$Y,
      X = X,
      D = scores$D,
      S = scores$S,
      propensity1 = scores$propensity1,
      propensity2 = scores$propensity2,
      scale_Y = TRUE,
      discretization = discretization,
      cost_treatment = cost_treatment,
      params = gurobi_params,
      mu_hat11 = scores$mu_hat11,
      mu_hat01 = scores$mu_hat01,
      mu_hat00 = scores$mu_hat00,
      mu_hat10 = scores$mu_hat10,
      max_treated_units = max_treated_units,
      maxtime1 = max_time,
      maxtime2 = max_time_frontier,
      alpha_seq = alpha_seq,
      m1 = scores$m1,
      m0 = scores$m0,
      quick_run = FALSE,
      no_parity_constraint = no_parity_constraint,
      additional_fairness_constraint = additional_fairness_constraint,
      parity_constraint = parity_constraint,
      frontier = frontier$frontier,
      unique_values = 1 - no_parity_constraint,
      distance = distance,
      probabilistic = isTRUE(probabilistic),
      numcores = numcores,
      threshold_probabilistic = isTRUE(threshold_probabilistic),
      two_directions = two_directions,
      parallel = FALSE,
      tolerance_frontier = tolerance_frontier,
      tolerance_optimization = tolerance_optimization,
      frontier_slack = lambda * sqrt(scores$n)
    )

    policy <- .policy_from_solver_result(
      result = result$result,
      design_spec = design_spec,
      policy_class = policy_class,
      no_parity_constraint = no_parity_constraint,
      threshold_probabilistic = isTRUE(threshold_probabilistic),
      X_train = X
    )

    cf <- .predict_training_counterfactuals(policy, X)

    unfairness <- compute_fairness(
      scores,
      cf$policy_s1,
      cf$policy_s0,
      distance = distance,
      two_directions = two_directions
    )

    selected_alpha <- result$result$alpha
    selected_grid_index <- match(
      selected_alpha[1],
      alpha_seq,
      nomatch = NA_integer_
    )
  }

  out <- list(
    call = match.call(),
    policy = policy,
    policy_parameters = policy,
    policy_class = policy_class,
    distance = distance,
    fairness = distance,
    two_directions = isTRUE(two_directions),
    capacity = capacity,
    max_treated_units = max_treated_units,
    lambda = lambda,
    frontier_slack = lambda * sqrt(scores$n),
    alpha_grid = alpha_seq,
    selected_alpha = selected_alpha,
    selected_grid_index = selected_grid_index,
    frontier = frontier,
    welfare_by_group = list(
      relative = unfairness$welfare_relative,
      level = unfairness$welfare_level
    ),
    unfairness = unfairness,
    scores = scores,
    nuisance = nuisance,
    training_data = dat[, unique(c(sensitive, covariates)), drop = FALSE],
    columns = list(
      outcome = outcome,
      treatment = treatment,
      sensitive = sensitive,
      covariates = covariates
    ),
    design_spec = design_spec,
    X = X,
    solver = solver,
    backend = "gurobi",
    no_parity_constraint = no_parity_constraint,
    additional_fairness_constraint = isTRUE(additional_fairness_constraint),
    parity_constraint = parity_constraint,
    UnFairness_bound = UnFairness_bound
  )

  structure(out, class = "fairtargeting")
}

# Prediction methods -----------------------------------------------------------

.predict_threshold_probabilistic <- function(score, probs) {
  probs <- as.numeric(probs)

  if (length(probs) != 2L) {
    .abort("`probs` must have length 2 for threshold-probabilistic policies.")
  }

  as.numeric(ifelse(score > 0, probs[2] + probs[1], probs[2]))
}

.predict_policy_values <- function(policy,
                                   newdata,
                                   type = c("response", "score"),
                                   sensitive_value = NULL) {
  type <- match.arg(type)

  if (!inherits(policy, "fairtargeting_policy")) {
    .abort("`policy` must be a fairtargeting policy object.")
  }

  x <- .design_matrix(
    newdata,
    policy$design_spec,
    sensitive_value = sensitive_value
  )

  if (identical(policy$method, "exhaustive_tree")) {
    score <- as.numeric(
      compute_predictions_tree(
        policy$tree_vector,
        x,
        length_tree = policy$length_tree
      )
    )

    return(score)
  }

  if (!identical(policy$method, "gurobi")) {
    .abort("Prediction is available for policies fitted by `fit_fair_policy()`.")
  }

  score <- as.numeric(cbind(1, x) %*% policy$beta)

  if (type == "score") {
    return(score)
  }

  if (isTRUE(policy$threshold_probabilistic)) {
    return(.predict_threshold_probabilistic(score, policy$probs))
  }

  if (isTRUE(policy$probabilistic)) {
    return(score)
  }

  as.numeric(score > 0)
}

#' Predict treatment policies from a fitted fairtargeting object.
#'
#' @param object A fitted object returned by `fit_fair_policy()`.
#' @param newdata Optional new data. If `NULL`, predictions are returned for the
#'   training data stored in the object.
#' @param type Prediction type. `"response"` returns policy values and `"score"`
#'   returns the underlying score when available.
#' @param sensitive_value Optional counterfactual sensitive value, either 0 or 1.
#' @param ... Unused.
#'
#' @return Numeric vector of policy values or scores.
#' @method predict fairtargeting
#' @export
predict.fairtargeting <- function(object,
                               newdata = NULL,
                               type = c("response", "score"),
                               sensitive_value = NULL,
                               ...) {
  type <- match.arg(type)

  if (!inherits(object, "fairtargeting")) {
    .abort("`object` must have class 'fairtargeting'.")
  }

  if (is.null(newdata)) {
    newdata <- object$training_data
  }

  .predict_policy_values(
    object$policy,
    newdata = newdata,
    type = type,
    sensitive_value = sensitive_value
  )
}

#' Predict treatment policies from a fairtargeting policy object.
#'
#' @param object A policy object stored in a fitted `fairtargeting` object.
#' @param newdata New data.
#' @param type Prediction type. `"response"` returns policy values and `"score"`
#'   returns the underlying score when available.
#' @param sensitive_value Optional counterfactual sensitive value, either 0 or 1.
#' @param ... Unused.
#'
#' @return Numeric vector of policy values or scores.
#' @method predict fairtargeting_policy
#' @export
predict.fairtargeting_policy <- function(object,
                                      newdata,
                                      type = c("response", "score"),
                                      sensitive_value = NULL,
                                      ...) {
  type <- match.arg(type)

  .predict_policy_values(
    object,
    newdata = newdata,
    type = type,
    sensitive_value = sensitive_value
  )
}

# Print methods ----------------------------------------------------------------

#' Print a fitted fairtargeting object.
#'
#' @param x A fitted `fairtargeting` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @method print fairtargeting
#' @export
print.fairtargeting <- function(x, ...) {
  cat("Fair Policy Targeting fit\n")
  cat("  policy class: ", x$policy_class, "\n", sep = "")
  cat("  backend:      ", x$backend %||% x$solver, "\n", sep = "")
  cat("  distance:     ", x$distance %||% x$fairness, "\n", sep = "")
  cat("  two-sided:    ", x$two_directions, "\n", sep = "")
  cat("  capacity:     ", if (is.null(x$capacity)) "none" else x$capacity, "\n", sep = "")
  cat("  alpha grid:   ", length(x$alpha_grid), " values", "\n", sep = "")

  if (!is.null(x$selected_alpha)) {
    cat("  selected alpha:", paste(x$selected_alpha, collapse = ", "), "\n")
  }

  cat("  unfairness:   ", signif(x$unfairness$objective, 5), "\n", sep = "")

  invisible(x)
}

#' Print a fairtargeting Pareto frontier object.
#'
#' @param x A `fairtargeting_frontier` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @method print fairtargeting_frontier
#' @export
print.fairtargeting_frontier <- function(x, ...) {
  cat("Fair Policy Targeting Pareto frontier\n")
  cat("  policy class: ", x$policy_class, "\n", sep = "")
  cat("  grid size:    ", length(x$alpha_grid), "\n", sep = "")
  cat("  backend:      ", x$solver %||% "gurobi", "\n", sep = "")

  if (!is.null(x$efficient_candidates)) {
    cat("  nondominated candidates:", nrow(x$efficient_candidates), "\n", sep = " ")
    print(utils::head(x$efficient_candidates), row.names = FALSE)
  } else if (!is.null(x$candidates)) {
    print(utils::head(x$candidates), row.names = FALSE)
  }

  invisible(x)
}

# Summary methods --------------------------------------------------------------

#' Summarize a fitted fairtargeting object.
#'
#' @param object A fitted `fairtargeting` object.
#' @param ... Unused.
#'
#' @return An object of class `summary.fairtargeting`.
#' @method summary fairtargeting
#' @export
summary.fairtargeting <- function(object, ...) {
  if (!inherits(object, "fairtargeting")) {
    .abort("`object` must have class 'fairtargeting'.")
  }

  structure(
    list(
      policy_class = object$policy_class,
      backend = object$backend %||% object$solver,
      distance = object$distance,
      two_directions = object$two_directions,
      capacity = object$capacity,
      alpha_grid = object$alpha_grid,
      selected_alpha = object$selected_alpha,
      unfairness = object$unfairness,
      welfare_by_group = object$welfare_by_group,
      no_parity_constraint = object$no_parity_constraint,
      additional_fairness_constraint = object$additional_fairness_constraint
    ),
    class = "summary.fairtargeting"
  )
}

#' Print a fairtargeting summary object.
#'
#' @param x A `summary.fairtargeting` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @method print summary.fairtargeting
#' @export
print.summary.fairtargeting <- function(x, ...) {
  cat("Summary of Fair Policy Targeting fit\n")
  cat("  policy class: ", x$policy_class, "\n", sep = "")
  cat("  backend:      ", x$backend %||% "unknown", "\n", sep = "")
  cat("  distance:     ", x$distance, "\n", sep = "")
  cat("  two-sided:    ", x$two_directions, "\n", sep = "")

  if (!is.null(x$capacity)) {
    cat("  capacity:     ", x$capacity, "\n", sep = "")
  }

  if (!is.null(x$alpha_grid)) {
    cat("  alpha grid:   ", length(x$alpha_grid), " values\n", sep = "")
  }

  if (!is.null(x$selected_alpha)) {
    cat("  selected alpha:", paste(x$selected_alpha, collapse = ", "), "\n")
  }

  cat("\nUnfairness objective:\n")
  print(x$unfairness$objective)

  cat("\nWelfare by group (relative):\n")
  print(x$welfare_by_group$relative)

  cat("\nWelfare by group (level):\n")
  print(x$welfare_by_group$level)

  invisible(x)
}
