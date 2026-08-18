# fairpolicy

Fair Policy Targeting for treatment allocation problems.

`fairpolicy` implements tools for estimating welfare-improving policy rules under fairness and Pareto-efficiency considerations. The package follows the Fair Policy Targeting framework of Viviano and Bradic (2022) and provides routines for nuisance estimation, doubly robust welfare scoring, Pareto-frontier approximation, and Gurobi-based fair policy optimization.

The package is designed for settings with

- an observed outcome;
- a binary treatment or allocation indicator;
- a binary sensitive attribute;
- policy covariates used to construct treatment rules.

The typical workflow is:

1. prepare a data set with outcome, treatment, sensitive group, and covariates;
2. estimate nuisance functions for treatment assignment and conditional outcomes;
3. construct doubly robust welfare scores for the two groups;
4. approximate the empirical Pareto frontier over an alpha grid;
5. solve the fair policy targeting problem;
6. inspect the selected policy, selected alpha value, welfare by group, unfairness criterion, and predictions.

## Reference

This package implements methods from:

Viviano, D., & Bradic, J. (2022). *Fair Policy Targeting*. arXiv:2005.12395v3.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("ostasovskyi/fairPolicy")
```

Load the package:

```r
library(fairpolicy)
```

The package uses base R and `stats` for lightweight utilities. Nuisance estimation with penalized regression uses `glmnet`; parallel utilities use `foreach` and `doParallel`.

```r
install.packages(c("glmnet", "foreach", "doParallel"))
```

Policy optimization is solved with Gurobi. To fit Fair Policy Targeting rules, install the Gurobi Optimizer and the Gurobi R package. Gurobi is not installed from CRAN. After installing and licensing Gurobi, install its R package from the Gurobi installation directory. On Linux or WSL, the command usually has the following form, with the path adjusted to your installed version:

```sh
R CMD INSTALL /opt/gurobi/linux64/R/gurobi_*.tar.gz
```

Check that the R interface is available:

```r
library(gurobi)
```

## Quick start

The following example creates a small synthetic data set and fits a fair policy using a maximum-score policy rule. The example uses the base-`glm` nuisance backend for readability. For larger or high-dimensional applications, use the default `nuisance_method = "glmnet"`.

```r
library(fairpolicy)

n <- 120
df <- data.frame(
  S = rbinom(n, 1, 0.45),
  X1 = rnorm(n),
  X2 = rnorm(n)
)

df$D <- rbinom(
  n,
  1,
  plogis(-0.2 + 0.5 * df$X1 - 0.3 * df$S)
)

tau <- 0.4 + 0.6 * df$X1 - 0.5 * df$S
df$Y <- 1 + 0.5 * df$X1 + 0.2 * df$X2 + tau * df$D + rnorm(n)

fit <- fit_fair_policy(
  data = df,
  outcome = "Y",
  treatment = "D",
  sensitive = "S",
  covariates = c("X1", "X2"),
  policy_class = "maxscore",
  distance = "parity",
  two_directions = TRUE,
  capacity = 0.50,
  discretization = 10,
  lambda = 0,
  solver = "gurobi",
  nuisance_method = "glm",
  n_folds = 5,
  seed = 123,
  max_time = 300
)

fit
summary(fit)

head(predict(fit, newdata = df))
head(predict(fit, newdata = df, type = "score"))
```

The returned object stores the fitted policy, the alpha grid used to approximate the Pareto frontier, the selected alpha value, welfare summaries by group, the selected unfairness criterion, doubly robust scores, and nuisance estimates.

## Main function

### `fit_fair_policy()`

`fit_fair_policy()` is the main user-facing function. It estimates nuisance functions when needed, constructs doubly robust welfare scores, approximates a Pareto frontier, and solves the second-stage fairness optimization problem.

```r
fit <- fit_fair_policy(
  data = df,
  outcome = "Y",
  treatment = "D",
  sensitive = "S",
  covariates = c("X1", "X2"),
  policy_class = "maxscore",
  distance = "parity",
  capacity = 0.50,
  solver = "gurobi"
)
```

Important arguments are grouped below.

| Argument group | Inputs |
|---|---|
| Data | `data`, `outcome`, `treatment`, `sensitive`, `covariates` |
| Policy class | `policy_class = "maxscore"`, `"probabilistic"`, `"threshold_probabilistic"`, or `"tree"` |
| Fairness criterion | `distance = "welfare"`, `"relative_welfare"`, `"parity"`, or `"envy"` |
| Pareto approximation | `discretization`, `n_grid`, `alpha_seq`, `lambda`, `open_alpha_grid` |
| Capacity | `capacity`, supplied as a treatment fraction or treatment count |
| Optimization | `solver = "gurobi"`, `max_time`, `max_time_frontier`, `gurobi_params`, `numcores` |
| Nuisance estimation | `nuisance`, `nuisance_method`, `cross_fit`, `n_folds`, `seed`, `outcome_family`, `trim` |
| Fairness extensions | `two_directions`, `no_parity_constraint`, `additional_fairness_constraint`, `parity_constraint` |

Main outputs include:

| Output | Meaning |
|---|---|
| `policy` | Fitted policy object |
| `policy_class` | Policy class used in optimization |
| `alpha_grid` | Alpha grid used for Pareto-frontier approximation |
| `selected_alpha` | Alpha value selected by the second-stage problem |
| `frontier` | Pareto-frontier approximation object |
| `unfairness` | Selected unfairness criterion and related components |
| `welfare_by_group` | Relative and level welfare for the two groups |
| `scores` | Doubly robust score object used by the optimizer |
| `nuisance` | Nuisance estimates used to construct scores |

## Step-by-step workflow

The complete workflow can also be run manually.

### 1. Estimate nuisance functions

```r
nuisance <- estimate_nuisance(
  data = df,
  outcome = "Y",
  treatment = "D",
  sensitive = "S",
  covariates = c("X1", "X2"),
  method = "glm",
  cross_fit = TRUE,
  n_folds = 5,
  seed = 123
)
```

`estimate_nuisance()` estimates treatment propensity scores, the sensitive-group probability, and conditional outcome functions. It returns an object of class `fairpolicy_nuisance`.

### 2. Construct doubly robust scores

```r
scores <- compute_dr_scores(
  outcome = df$Y,
  treatment = df$D,
  sensitive = df$S,
  nuisance = nuisance
)

str(scores)
```

The score object contains the group-specific doubly robust welfare components used for Pareto-frontier estimation and fairness calculations.

### 3. Approximate the Pareto frontier

```r
frontier <- estimate_pareto_frontier(
  data = df,
  scores = scores,
  sensitive = "S",
  covariates = c("X1", "X2"),
  policy_class = "maxscore",
  capacity = 0.50,
  alpha_grid = seq(0, 1, length.out = 10),
  solver = "gurobi",
  max_time = 300
)

frontier
head(frontier$candidates)
head(frontier$efficient_candidates)
```

Each alpha value controls the relative weight placed on the two group-specific welfare objectives. The second-stage fair policy optimization searches over feasible policies satisfying approximate Pareto-frontier constraints and minimizes the selected unfairness criterion.

## Policy classes

The package supports four policy classes:

```r
policy_class = "maxscore"
policy_class = "probabilistic"
policy_class = "threshold_probabilistic"
policy_class = "tree"
```

- `"maxscore"` fits deterministic linear threshold rules.
- `"probabilistic"` fits linear probabilistic treatment rules.
- `"threshold_probabilistic"` fits threshold rules with treatment probabilities on either side of the threshold.
- `"tree"` searches over simple tree-based policy rules.

## Fairness criteria

The main fairness criteria are:

```r
distance = "welfare"
distance = "relative_welfare"
distance = "parity"
distance = "envy"
```

For welfare, relative-welfare, and parity criteria, `two_directions = TRUE` computes absolute two-sided disparities. If `two_directions = FALSE`, the signed disparity is used. The argument `fairness` is accepted as an alternative name for `distance`.

## Alpha grid and Pareto slack

The alpha grid defines the discretized approximation to the Pareto frontier:

```r
alpha_seq <- seq(0, 1, length.out = discretization)
```

If `alpha_seq` is not supplied, the package constructs a grid using `discretization`; when `discretization` is not supplied, the default is based on `floor(sqrt(n))`. The slack parameter `lambda` controls approximate Pareto feasibility in the second-stage optimization problem.

## Capacity constraints

Treatment capacity can be supplied as either a fraction or a count:

```r
capacity = 0.40  # treat at most 40 percent of the sample
capacity = 150   # treat at most 150 observations
```

If `capacity` is `NULL`, no binding capacity limit is imposed beyond the sample size.

## Prediction and summaries

After fitting, use standard S3 methods:

```r
fit
summary(fit)

policy_values <- predict(fit, newdata = df)
policy_scores <- predict(fit, newdata = df, type = "score")

head(policy_values)
head(policy_scores)
```

Counterfactual sensitive-attribute predictions can be requested with `sensitive_value`:

```r
head(predict(fit, newdata = df, sensitive_value = 1))
head(predict(fit, newdata = df, sensitive_value = 0))
```

## Input requirements

The data supplied to `fit_fair_policy()` should contain:

- a numeric, integer, or logical outcome variable;
- a binary treatment indicator;
- a binary sensitive attribute;
- at least one policy covariate;
- no missing values in the variables used for fitting.

Binary variables may be coded as `0`/`1`, logical values, character values with two distinct categories, or two-level factors. By convention, `S = 1` denotes the sensitive group.