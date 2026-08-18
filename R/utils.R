# Utilities and package helpers ------------------------------------------------

utils::globalVariables(c("i", "k", "g_i"))


`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.abort <- function(message, call. = FALSE) {
  stop(message, call. = call.)
}

.require_gurobi <- function() {
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    .abort(paste(
      "The Gurobi optimization backend requires Gurobi and the `gurobi` R",
      "package. Install Gurobi and its R bindings before running",
      "the exact Fair Policy Targeting optimizer."
    ))
  }

  invisible(TRUE)
}

.call_gurobi <- function(model, params = list()) {
  .require_gurobi()
  gurobi_ns <- asNamespace("gurobi")
  solver <- get("gurobi", envir = gurobi_ns)
  solver(model, params = params)
}

.require_glmnet <- function() {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    .abort(paste(
      "The `glmnet` nuisance-estimation backend requires the `glmnet` R",
      "package. Install `glmnet`, or call `estimate_nuisance(..., method = 'glm')`."
    ))
  }

  invisible(TRUE)
}

.check_columns <- function(data, columns, what = "columns") {
  missing <- setdiff(columns, names(data))

  if (length(missing) > 0) {
    .abort(sprintf(
      "Missing %s in `data`: %s.",
      what,
      paste(missing, collapse = ", ")
    ))
  }

  invisible(TRUE)
}

.check_no_missing <- function(data, columns) {
  bad <- columns[vapply(data[columns], function(x) anyNA(x), logical(1))]

  if (length(bad) > 0) {
    .abort(sprintf(
      "Missing values are not allowed in: %s. Remove or impute them before fitting a fair policy.",
      paste(bad, collapse = ", ")
    ))
  }

  invisible(TRUE)
}

.as_binary <- function(x, name) {
  if (is.logical(x)) {
    return(as.integer(x))
  }

  if (is.factor(x)) {
    if (nlevels(x) != 2) {
      .abort(sprintf("`%s` must have exactly two levels.", name))
    }

    return(as.integer(x == levels(x)[2]))
  }

  if (is.character(x)) {
    ux <- sort(unique(x))

    if (length(ux) != 2) {
      .abort(sprintf("`%s` must have exactly two distinct values.", name))
    }

    return(as.integer(x == ux[2]))
  }

  if (is.numeric(x) || is.integer(x)) {
    ux <- sort(unique(stats::na.omit(x)))

    if (!all(ux %in% c(0, 1))) {
      .abort(sprintf(
        "`%s` must be binary and coded as 0/1, logical, character, or a two-level factor.",
        name
      ))
    }

    return(as.integer(x))
  }

  .abort(sprintf(
    "`%s` must be binary and coded as 0/1, logical, character, or a two-level factor.",
    name
  ))
}

.validate_required_data <- function(data,
                                    outcome,
                                    treatment,
                                    sensitive,
                                    covariates) {
  if (!is.data.frame(data)) {
    .abort("`data` must be a data.frame or tibble-like object.")
  }

  if (!is.character(outcome) || length(outcome) != 1) {
    .abort("`outcome` must be a single column name.")
  }

  if (!is.character(treatment) || length(treatment) != 1) {
    .abort("`treatment` must be a single column name.")
  }

  if (!is.character(sensitive) || length(sensitive) != 1) {
    .abort("`sensitive` must be a single column name.")
  }

  if (!is.character(covariates) || length(covariates) < 1) {
    .abort("`covariates` must be a non-empty character vector of column names.")
  }

  needed <- unique(c(outcome, treatment, sensitive, covariates))

  .check_columns(data, needed)
  .check_no_missing(data, needed)

  y <- data[[outcome]]

  if (!is.numeric(y) && !is.integer(y) && !is.logical(y)) {
    .abort("`outcome` must be numeric, integer, or logical.")
  }

  invisible(TRUE)
}

.validate_vectors <- function(Y, D, S, X = NULL) {
  n <- length(Y)

  if (length(D) != n || length(S) != n) {
    .abort("`Y`, `D`, and `S` must have the same length.")
  }

  D <- .as_binary(D, "D")
  S <- .as_binary(S, "S")

  if (!all(c(0, 1) %in% unique(S))) {
    .abort("`S` must contain both groups 0 and 1.")
  }

  if (!all(c(0, 1) %in% unique(D))) {
    .abort("`D` must contain both treatment states 0 and 1.")
  }

  if (!is.null(X)) {
    X <- as.matrix(X)

    if (nrow(X) != n) {
      .abort("`X` must have one row per observation in `Y`, `D`, and `S`.")
    }

    if (ncol(X) < 1) {
      .abort("`X` must have at least one column.")
    }

    if (anyNA(X) || any(!is.finite(X))) {
      .abort(
        "`X` must contain only finite, non-missing numeric values after model-matrix construction."
      )
    }
  }

  invisible(TRUE)
}

.validate_original_vectors <- .validate_vectors

.clip_prob <- function(x, eps = 1e-3) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

.resolve_capacity <- function(capacity, n) {
  if (is.null(capacity)) {
    return(list(max_count = n, fraction = 1, supplied = FALSE))
  }

  if (!is.numeric(capacity) || length(capacity) != 1 || is.na(capacity)) {
    .abort("`capacity` must be NULL or a single numeric value.")
  }

  if (capacity < 0) {
    .abort("`capacity` must be non-negative.")
  }

  if (capacity <= 1) {
    max_count <- floor(capacity * n)

    return(list(
      max_count = max_count,
      fraction = capacity,
      supplied = TRUE
    ))
  }

  max_count <- floor(capacity)

  if (max_count > n) {
    .abort(
      "When `capacity` is larger than 1 it is interpreted as a count and cannot exceed `nrow(data)`."
    )
  }

  list(
    max_count = max_count,
    fraction = max_count / n,
    supplied = TRUE
  )
}

.resolve_capacity_count <- function(capacity, n) {
  .resolve_capacity(capacity, n)$max_count
}

.resolve_capacity_original <- .resolve_capacity_count

.normalize_distance <- function(distance = NULL, fairness = NULL) {
  if (!is.null(distance) && !is.null(fairness)) {
    if (!identical(distance, fairness)) {
      .abort("Use either `distance` or `fairness`, or set both to the same value.")
    }
  }

  x <- distance %||% fairness %||% "envy"

  aliases <- c(
    parity_abs = "parity",
    welfare_disparity = "welfare",
    welfare_abs = "welfare",
    relative_welfare_abs = "relative_welfare"
  )

  x <- if (x %in% names(aliases)) unname(aliases[[x]]) else x

  choices <- c("welfare", "relative_welfare", "parity", "envy")

  if (!x %in% choices) {
    .abort(sprintf(
      "`distance` must be one of: %s.",
      paste(choices, collapse = ", ")
    ))
  }

  x
}

.resolve_no_parity <- function(no_parity_constraint = NULL,
                               noparity_constraint = NULL) {
  if (!is.null(no_parity_constraint) &&
      !is.null(noparity_constraint) &&
      !identical(isTRUE(no_parity_constraint), isTRUE(noparity_constraint))) {
    .abort(
      "`no_parity_constraint` and `noparity_constraint` cannot disagree."
    )
  }

  isTRUE(no_parity_constraint %||% noparity_constraint %||% FALSE)
}

.make_design_spec <- function(data,
                              sensitive,
                              covariates,
                              include_sensitive_first = TRUE) {
  variables <- unique(c(
    if (isTRUE(include_sensitive_first)) sensitive else character(0),
    covariates
  ))

  .check_columns(data, variables)

  form <- stats::reformulate(
    termlabels = variables,
    response = NULL,
    intercept = FALSE
  )

  trm <- stats::terms(form, data = data)

  mf <- stats::model.frame(
    trm,
    data = data[, variables, drop = FALSE],
    na.action = stats::na.pass
  )

  mm <- stats::model.matrix(trm, data = mf)

  if (ncol(mm) < 1) {
    .abort("The policy design matrix has zero columns. Check `covariates`.")
  }

  if (isTRUE(include_sensitive_first)) {
    s_col <- sensitive

    if (!s_col %in% colnames(mm)) {
      .abort(
        "The policy design matrix must contain the sensitive attribute as its first column."
      )
    }

    mm <- cbind(
      mm[, s_col, drop = FALSE],
      mm[, setdiff(colnames(mm), s_col), drop = FALSE]
    )
  }

  factor_cols <- names(mf)[vapply(mf, is.factor, logical(1))]
  xlevels <- lapply(mf[factor_cols], levels)

  list(
    variables = variables,
    sensitive = sensitive,
    include_sensitive_first = isTRUE(include_sensitive_first),
    terms = trm,
    columns = colnames(mm),
    contrasts = attr(mm, "contrasts"),
    xlevels = xlevels
  )
}

.design_matrix <- function(newdata, spec, sensitive_value = NULL) {
  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  .check_columns(newdata, spec$variables)

  dat <- newdata[, spec$variables, drop = FALSE]

  if (!is.null(sensitive_value) && spec$sensitive %in% names(dat)) {
    if (!sensitive_value %in% c(0, 1)) {
      .abort("`sensitive_value` must be 0, 1, or NULL.")
    }

    dat[[spec$sensitive]] <- sensitive_value
  }

  mf <- stats::model.frame(
    spec$terms,
    data = dat,
    na.action = stats::na.pass,
    xlev = spec$xlevels
  )

  mm <- stats::model.matrix(
    spec$terms,
    data = mf,
    contrasts.arg = spec$contrasts
  )

  missing <- setdiff(spec$columns, colnames(mm))

  if (length(missing) > 0) {
    zeros <- matrix(0, nrow = nrow(mm), ncol = length(missing))
    colnames(zeros) <- missing
    mm <- cbind(mm, zeros)
  }

  mm[, spec$columns, drop = FALSE]
}

.add_intercept <- function(x) {
  cbind(`(Intercept)` = 1, as.matrix(x))
}

#' Create an alpha grid for the Pareto frontier.
#'
#' Creates a grid of alpha values used to estimate the discretized Pareto
#' frontier. By default, the grid includes both endpoints, 0 and 1. Set
#' `open = TRUE` to use an open grid that excludes the endpoints.
#'
#' @param n_grid Integer grid size. If `NULL`, `n` must be supplied.
#' @param n Optional sample size used when `n_grid` is `NULL`.
#' @param open Logical. If `TRUE`, use an open grid in `(0, 1)`.
#'
#' @return Numeric vector of alpha values.
#' @export
make_alpha_grid <- function(n_grid = NULL, n = NULL, open = FALSE) {
  if (is.null(n_grid)) {
    if (is.null(n)) {
      .abort("Supply either `n_grid` or `n`.")
    }

    n_grid <- floor(sqrt(n))
  }

  if (!is.numeric(n_grid) || length(n_grid) != 1 || is.na(n_grid)) {
    .abort("`n_grid` must be a single positive integer.")
  }

  n_grid <- as.integer(n_grid)

  if (n_grid < 1) {
    .abort("`n_grid` must be at least 1.")
  }

  if (isTRUE(open)) {
    seq(
      from = 1 / (n_grid + 1),
      to = n_grid / (n_grid + 1),
      length.out = n_grid
    )
  } else {
    seq(from = 0, to = 1, length.out = n_grid)
  }
}

#' Simulate data for a minimal Fair Policy Targeting example.
#'
#' @param n Number of observations.
#' @param seed Optional random seed.
#'
#' @return A data.frame with outcome `Y`, treatment `D`, sensitive attribute
#'   `S`, and covariates `X1`, `X2`, and `X3`.
#' @export
simulate_fairpolicy_data <- function(n = 300, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.numeric(n) || length(n) != 1 || n < 10) {
    .abort("`n` must be a single integer of at least 10.")
  }

  x1 <- stats::rnorm(n)
  x2 <- stats::rbinom(n, 1, 0.5)
  x3 <- stats::runif(n, -1, 1)

  s <- stats::rbinom(
    n,
    1,
    stats::plogis(-0.2 + 0.5 * x1 - 0.3 * x2)
  )

  e <- stats::plogis(-0.1 + 0.4 * x1 - 0.2 * x2 + 0.3 * s)

  d <- stats::rbinom(n, 1, e)

  tau <- 0.6 + 0.7 * x1 - 0.4 * x2 - 0.3 * s + 0.2 * x3

  y0 <- 0.5 + 0.5 * x1 + 0.2 * x2 - 0.4 * s +
    stats::rnorm(n, sd = 0.7)

  y <- y0 + d * tau

  data.frame(
    Y = y,
    D = d,
    S = s,
    X1 = x1,
    X2 = x2,
    X3 = x3
  )
}
