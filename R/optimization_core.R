# Optimization core ------------------------------------------------------------
#
# Core optimization routines for Fair Policy Targeting. This file contains
# model-construction helpers, Pareto-frontier estimation, maximum-score policy
# optimization, probabilistic policy optimization, and exhaustive tree-search
# routines used by the package-facing wrappers.


# ---- from helpers.R ----
## This script contains helpers for the main library

## Inputs:
## Y: outcome of interest
## X_reg: X used in the regression function
## X_reg as a structure cbind(X, D, D * X_int * S, D * X_int * (1 - S), S * X_int, D * X_int)
## X covariates (without interaction terms)
## X_int: interaction terms
## seeds: seeds for splitting folds
## K: number of folds
## family: either binomial or gaussian
## Output: estimated functions on different subgroups
cross_fitting_mean <- function(Y, X_reg, X, X_int, seeds = 123, K = 5, family = 'binomial'){

  set.seed(seeds)
  ii <- sample(rep(1:K, length= dim(X)[1]))
  m_hat11 <- m_hat10 <- m_hat01 <- m_hat00 <- rep(NA, dim(X)[1])
  for(i in 1:K){

    train <- ii != i
    hold <- ii == i
    if(family == 'binomial'){
    my_reg <- glmnet::glmnet(y = as.factor(Y)[train], x = as.matrix(X_reg)[train, ], family = family, lambda = exp(-12))
    } else {
      my_reg <- glmnet::glmnet(y = Y[train], x = as.matrix(X_reg)[train, ], family = family, lambda = exp(-12))
    }
    XX <- cbind(X, 1,  X_int ,  X_int * 0, X_int , X_int)
    if(family == 'binomial'){
    m_hat11[hold] <- stats::predict(my_reg, newx = as.matrix(XX)[hold,], type = 'response')
    } else {
      m_hat11[hold] <- stats::predict(my_reg, newx = as.matrix(XX)[hold,])
    }
    XX <- cbind(X, 0,  0 * X_int ,  X_int * 0, X_int, 0 * X_int)
    if(family == 'binomial'){
    m_hat10[hold] <- stats::predict(my_reg, newx = as.matrix(XX)[hold,], type = 'response')
    } else {
      m_hat10[hold] <- stats::predict(my_reg, newx = as.matrix(XX)[hold,])
    }

    XX <- cbind(X, 1,  X_int *0 ,  X_int, 0  * X_int, X_int )
    if(family == 'binomial'){
    m_hat01[hold] <- stats::predict(my_reg, newx = as.matrix(XX)[hold,], type = 'response')
    } else {
      m_hat01[hold] <- stats::predict(my_reg, newx = as.matrix(XX)[hold,])
    }
    XX <- cbind(X, 0,  0 * X_int ,  X_int * 0, 0 * X_int, 0*X_int)
    if(family == 'binomial'){
    m_hat00[hold] <- stats::predict(my_reg, newx = as.matrix(XX)[hold,], type = 'response')
    } else {
      m_hat00[hold] <- stats::predict(my_reg, newx = as.matrix(XX)[hold,])
    }


  }
  return(cbind(m_hat11, m_hat10, m_hat01, m_hat00))
}

## Estimate the propensity score via cross fitting
## Inputs: D (treatment assignment)
##         X (covariates)
cross_fitting_propensity <- function(D,X, seeds = 123, K = 5){
  set.seed(seeds)
  ii <- sample(rep(1:K, length= dim(X)[1]))
  propensity1 <- rep(NA, dim(X)[1])
  for(i in 1:K){
    train <- ii != i
    hold <- ii == i
    propensity_score <- glmnet::cv.glmnet(y = as.factor(D)[train], x = as.matrix(X)[train,], family = 'binomial')
    propensity1[hold] <- stats::predict(propensity_score, newx = as.matrix(X)[hold,], type = 'response')
  }
  return(propensity1)
}

## The function initialize the quadratic program

## Inputs: Y outcome
##         X : covariates for targeting, the first entry is assumed to be the sensity attribute
##         D: treatment assignment
##         S: sensitive attribute
##         propensity1: probability of treatment
##         propensity2 : probability of sensitive attribute (note: this are vector with the predict probs with n entries)
##         B: coefficients upper bounds
##         params: parameters of the program
##         tolerance_constraint: tolerance coefficient for the MILP program
##         scale_Y: whethere Y is rescaled(default F)
##         cost_treatment (default 0)
##         alpha: importance weights
##         g_i1: objective for the sensitive group
##         g_i2: objective for the opposite group
##         max_treated units: maximum number of treated individuals
##         maxtime: maximum time
##         m1, m0: predicted conditional means
##         no_parity_constraint: whether an additional constraint on the welfare is included or not (use or not use S for prediction?)
##         distance: unfairness measure
##         constant1, constant2: offset constants which depend on the unfairness measure
##         two_directions: only active if distance != envy. If T it takes fairness as absolute value, otherwise it takes the
##                         difference between the priviledge and sensitive group

## return: gurobi model
create_model_quadratic_program <- function(Y, X, D, S, propensity1, propensity2,B=1, params = NA,
                                           tolerance_constraint = 10**(-7),  scale_Y = F,
                                           cost_treatment  = 0, alpha = 1/2, g_i1, g_i2, max_treated_units,
                                           maxtime = 300, m1 = 0, m0 = 0,
                                           no_parity_constraint = F,
                                           distance = 'envy', constant1 = NA,constant2 = NA,
                                           two_directions = T) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  ## Recall: S is the first column of X
  XX1 <- cbind(1, 1, X[,-1])
  XX0 <- cbind(1, 0, X[,-1])
  if(no_parity_constraint == F){
    XX1 <- cbind(1, X)
    XX0 <- cbind(1, X)
  }
  colnames(XX1) <- colnames(XX0) <- paste0('V', c(1:(p+1)))
  XX <- rbind(XX1, XX0)
  # Prepare input
  C <- B*max(apply(XX, 1, function(x) sum(abs(x))))
  XX <- as.matrix(XX/C)
  model  <- list()
  ## Specify the constraints
  ## Consider the vector of (z_1, ..., z_n,beta_0 , beta_1, ..., beta_p)
  if(distance == 'envy'){
  model$obj <- c(g_i1, g_i2, rep(0, p + 1))
  model$Q <- diag(c(rep(0, n), rep(0, n), rep(0, p + 1)))
  model$A<- rbind(cbind(diag(1, nrow = 2*n), -XX), cbind(diag(1, nrow = 2*n), -XX),
                  c(rep(1,n)*S, rep(1, n)*(1 - S), rep(0, p +1)))
  model$modelsense<-'min'
  model$rhs<- c(rep(1 - tolerance_constraint, dim(model$A)[1]/2),
                rep(tolerance_constraint, dim(model$A)[1]/2), max_treated_units)
  model$sense<- c(rep('<=', dim(model$A)[1]/2), rep('>', dim(model$A)[1]/2), '<=')
  model$vtype<- c(rep('B', 2*n), rep('C', p+1))
  model$ub<- c(rep(1,2*n), rep(B,1+p))
  model$lb<- c(rep(0,2*n), rep(-B,p + 1))
  }
  if(distance == 'welfare' | distance == 'relative_welfare'){
  # Welfare program has two more variables to guarantee that the absolute value in the objective
  # function holds
  # the last two entries are two binary variables one indicating W_1 - W_0 >= 0, the other W_0 - W_1 >= 0
  ## note: if distance was specified as relative welfare than automatically these constants terms are zero

  model$obj <- c(rep(0, 2*n + p + 1), constant1 - constant2, constant2 - constant1)
  Q_matrix <- matrix(0, nrow= 2*n + p + 3, ncol = 2*n + p + 3)
  Q_matrix[1:n, 2 * n + p + 2] <- g_i1
  Q_matrix[(n + 1):(2*n), 2 * n + p + 2] <- -g_i2

  Q_matrix[1:n, 2 * n + p + 3] <- -g_i1
  Q_matrix[(n + 1):(2*n), 2 * n + p + 3] <- g_i2

  model$A<- rbind(cbind(diag(1, nrow = 2*n), -XX, 0, 0),
                  cbind(diag(1, nrow = 2*n), -XX, 0, 0),
                  c(rep(1,n)*S, rep(1, n)*(1 - S), rep(0, p +3)))
  model$Q <- Q_matrix

  model$modelsense<-'min'
  ## tolerance_constraint enters here
  model$rhs<- c(rep(1 - tolerance_constraint, dim(model$A)[1]/2),
                rep(tolerance_constraint, dim(model$A)[1]/2),
                max_treated_units)
  model$sense<- c(rep('<=', dim(model$A)[1]/2), rep('>', dim(model$A)[1]/2), '<=')
  model$vtype<- c(rep('B', 2*n), rep('C', p+1), rep('B', 2))
  # Put bounds on the parameter space (If commented, parameter space = real line)
  model$ub<- c(rep(1,2*n), rep(B,1+p), rep(1, 2))
  model$lb<- c(rep(0,2*n), rep(-B,p + 1), rep(0, 2))

  ## Objective if do not use absolute value for unfairness
  saved_objective_one_direction <- c(-g_i1, g_i2, rep(0, 3 + p))
  }

  if(distance == 'parity'){

    ## constant1, constant2 here should be 0s
    model$obj <- c(rep(0, 2*n + p + 1), 0, 0)
    Q_matrix <- matrix(0, nrow= 2*n + p + 3, ncol = 2*n + p + 3)
    Q_matrix[1:n, 2 * n + p + 2] <- S/mean(S)
    Q_matrix[(n + 1):(2*n), 2 * n + p + 2] <- -(1 - S)/(1 - mean(S))

    Q_matrix[1:n, 2 * n + p + 3] <- -S/mean(S)
    Q_matrix[(n + 1):(2*n), 2 * n + p + 3] <- (1 - S)/(1 - mean(S))

    model$A<- rbind(cbind(diag(1, nrow = 2*n), -XX, 0, 0),
                    cbind(diag(1, nrow = 2*n), -XX, 0, 0),
                    c(rep(1,n)*S, rep(1, n)*(1 - S), rep(0, p +3)))
    model$Q <- Q_matrix

    model$modelsense<-'min'
    ## tolerance_constraint enters here
    model$rhs<- c(rep(1 - tolerance_constraint, dim(model$A)[1]/2),
                  rep(tolerance_constraint, dim(model$A)[1]/2), max_treated_units)
    model$sense<- c(rep('<=', dim(model$A)[1]/2), rep('>', dim(model$A)[1]/2), '<=')
    model$vtype<- c(rep('B', 2*n), rep('C', p+1), rep('B', 2))
    # Put bounds on the parameter space (If commented, parameter space = real line)
    model$ub<- c(rep(1,2*n), rep(B,1+p), rep(1, 2))
    model$lb<- c(rep(0,2*n), rep(-B,p + 1), rep(0, 2))
    saved_objective_one_direction <- c(-S/mean(S),(1 - S)/(1 - mean(S)), rep(0, 3 + p))
  }
  ## Case where we take differences without absolute values in the objective (unfairness)
  if(two_directions == F & distance != 'envy') {
    model$obj <- saved_objective_one_direction
  }
  return(list(model, g_i1, g_i2))
}
### Same function but with probabilistic assignment (see Appendix B.3)

create_model_quadratic_program_probabilistic <- function(Y, X, D, S, propensity1, propensity2,B=1,
                                                         params = NA,  tolerance_constraint = 10**(-7),  scale_Y = F,
                                                         cost_treatment  = 0, alpha = 1/2, g_i1, g_i2,
                                                         max_treated_units, maxtime = 300, m1 = 0, m0 = 0,
                                                         no_parity_constraint = F, distance = 'envy',
                                                         constant1 = NA,constant2 = NA, two_directions = T) {
  n <- dim(X)[1]
  p <- dim(X)[2]

  XX1 <- cbind(1, 1, X[,-1])
  XX0 <- cbind(1, 0, X[,-1])
  if(no_parity_constraint == F){
    XX1 <- cbind(1, X)
    XX0 <- cbind(1, X)
  }
  colnames(XX1) <- colnames(XX0) <- paste0('V', c(1:(p+1)))
  XX <- rbind(XX1, XX0)
  # Prepare input
  XX <- as.matrix(XX)
  model  <- list()
  if(distance == 'envy'){
    model$obj <- c(g_i1, g_i2, rep(0, p + 1))
    model$Q <- diag(c(rep(0, n), rep(0, n), rep(0, p + 1)))
    model$A<- rbind(cbind(diag(1, nrow = 2*n), -XX), c(rep(1,n)*S, rep(1, n)*(1 - S), rep(0, p +1)))
    model$modelsense<-'min'
    model$rhs<- c(rep(0, dim(model$A)[1] - 1), max_treated_units)
    model$sense<- c(rep('=', dim(model$A)[1] - 1), '<=')
    model$vtype<- c(rep('C', 2*n), rep('C', p+1))
    model$ub<- c(rep(1,2*n), rep(B,1+p))
    model$lb<- c(rep(0,2*n), rep(-B,p + 1))
  }
  if(distance == 'welfare' | distance == 'relative_welfare'){
    # Welfare program has two more variables to guarantee that the absolute value in the objective
    # function holds
    # the last two entries are two binary variables one indicating W_1 - W_0 >= 0, the other W_0 - W_1 >= 0
    model$obj <- c(rep(0, 2*n + p + 1), constant1 - constant2, constant2 - constant1)
    Q_matrix <- matrix(0, nrow= 2*n + p + 3, ncol = 2*n + p + 3)
    Q_matrix[1:n, 2 * n + p + 2] <- g_i1
    Q_matrix[(n + 1):(2*n), 2 * n + p + 2] <- -g_i2

    Q_matrix[1:n, 2  * n + p + 3] <- -g_i1
    Q_matrix[(n + 1):(2*n), 2 * n + p + 3] <- g_i2

    model$A<- rbind(cbind(diag(1, nrow = 2*n), -XX, 0, 0),
                    c(rep(1,n)*S, rep(1, n)*(1 - S), rep(0, p +3)))
    model$Q <- Q_matrix

    model$modelsense<-'min'
    ## tolerance_constraint enters here
    model$rhs<- c(rep(0, dim(model$A)[1] - 1), max_treated_units)
    model$sense<- c(rep('=', dim(model$A)[1] - 1), '<=')
    model$vtype<- c(rep('C', 2*n), rep('C', p+1), rep('B', 2))
    # Put bounds on the parameter space (If commented, parameter space = real line)
    model$ub<- c(rep(1,2*n), rep(B,1+p), rep(1, 2))
    model$lb<- c(rep(0,2*n), rep(-B,p + 1), rep(0, 2))
    ## Objective if do not use absolute value
    saved_objective_one_direction <- c(-g_i1, g_i2, rep(0, 3 + p))
  }

  if(distance == 'parity'){
    model$obj <- c(rep(0, 2*n + p + 1), constant1 - constant2, constant2 - constant1)
    Q_matrix <- matrix(0, nrow= 2*n + p + 3, ncol = 2*n + p + 3)
    Q_matrix[1:n, 2 * n + p + 2] <- S/mean(S)
    Q_matrix[(n + 1):(2*n), 2 * n + p + 2] <- -(1 - S)/(1 - mean(S))

    Q_matrix[1:n, 2 * n + p + 3] <- -S/mean(S)
    Q_matrix[(n + 1):(2*n), 2 * n + p + 3] <- (1 - S)/(1 - mean(S))

    model$A<- rbind(cbind(diag(1, nrow = 2*n), -XX, 0, 0),
                    c(rep(1,n)*S, rep(1, n)*(1 - S), rep(0, p +3)))
    model$Q <- Q_matrix

    model$modelsense<-'min'
    ## tolerance_constraint enters here
    model$rhs<- c(rep(0, dim(model$A)[1] - 1), max_treated_units)
    model$sense<- c(rep('=', dim(model$A)[1] - 1), '<=')
    model$vtype<- c(rep('C', 2*n), rep('C', p+1), rep('B', 2))
    # Put bounds on the parameter space (If commented, parameter space = real line)
    model$ub<- c(rep(1,2*n), rep(B,1+p), rep(1, 2))
    model$lb<- c(rep(0,2*n), rep(-B,p + 1), rep(0, 2))
    saved_objective_one_direction <- c(-S/mean(S),(1 - S)/(1 - mean(S)), rep(0, 3 + p))
  }
  if(two_directions == F & distance != 'envy') {
    model$obj <- saved_objective_one_direction
  }
  return(list(model, g_i1, g_i2))
}


## Same function as before with probabilistic assignment and different probabilities for different thresholds
## you have xi1, ...xin, beta_1, ..., beta_p, indicator1, indicator2, prob1, prob2, prob1, ..., prob2
## each n is two times, first time for S= 1, second time for S = 0
create_model_quadratic_program_threshold_probabilistic <- function(Y, X, D, S, propensity1, propensity2,B=1, params = NA,
                                           tolerance_constraint = 10**(-7),  scale_Y = F,
                                           cost_treatment  = 0, alpha = 1/2, g_i1, g_i2, max_treated_units,
                                           maxtime = 300, m1 = 0, m0 = 0,
                                           no_parity_constraint = F,
                                           distance = 'envy', constant1 = NA,constant2 = NA,
                                           two_directions = T) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  ## Recall: S is the first column of X
  XX1 <- cbind(1, 1, X[,-1])
  XX0 <- cbind(1, 0, X[,-1])
  if(no_parity_constraint == F){
    XX1 <- cbind(1, X)
    XX0 <- cbind(1, X)
  }
  colnames(XX1) <- colnames(XX0) <- paste0('V', c(1:(p+1)))

  XX <- rbind(XX1, XX0)
  # Prepare input
  C <- B*max(apply(XX, 1, function(x) sum(abs(x))))
  XX <- as.matrix(XX/C)
  model  <- list()
  if(distance == 'envy'){
    model$obj <- c(rep(0, 2 * n + p + 3), g_i1, g_i2)
    model$Q <- diag(c(rep(0, n), rep(0, n), rep(0, p + 3), rep(0, 2 * n)))
    model$A<- rbind(cbind(diag(1, nrow = 2*n), -XX, 0, 0, diag(0, nrow = 2 * n)),
                    cbind(diag(1, nrow = 2*n), -XX, 0, 0, diag(0, nrow = 2 * n)),
                    c(rep(0, p + 1 + 2 * n), rep(1,n)*S, rep(1, n)*(1 - S)),
                    c(rep(0,2 * n), rep(0, p + 1), 1, 1 , rep(0, 2 * n))
                    )
    model$modelsense<-'min'
    model$rhs<- c(rep(1 - tolerance_constraint, (dim(model$A)[1] - 2)/2),
                  rep(tolerance_constraint, (dim(model$A)[1] - 2)/2), max_treated_units, 1)
    model$sense<- c(rep('<=', (dim(model$A)[1] -2)/2), rep('>', (dim(model$A)[1] -2)/2), '<=', '<=')
    model$vtype<- c(rep('B', 2*n), rep('C', p+ 3 + 2 * n))
    model$ub<- c(rep(1,2*n), rep(B,1+p), 1, 1, rep(1, 2 * n))
    model$lb<- c(rep(0,2*n), rep(-B,p + 1), 0, 0, rep(0, 2 * n))
  }
  if(distance == 'welfare' | distance == 'relative_welfare'){
    # Welfare program has two more variables to guarantee that the absolute value in the objective
    # function holds
    # the last two entries are two binary variables one indicating W_1 - W_0 >= 0, the other W_0 - W_1 >= 0
    ## note: if distance was specified as relative welfare than automatically these constants terms are zero

    model$obj <- c(rep(0, 2*n + p + 1), constant1 - constant2, constant2 - constant1, rep(0, 2 * n + 2))
    Q_matrix <- matrix(0, nrow= 2*n + p + 3 + 2 + 2 * n, ncol = 2*n + p + 3 + 2 + 2 * n)
    Q_matrix[(2 * n + p + 3 + 2 + 1):(2 * n + p + 3 + 2 + n), 2 * n + p + 2] <- g_i1
    Q_matrix[(2 * n + p + 3 + 2 + 1 + n):(2 * n + p + 3 + 2 + 2*n), 2 * n + p + 2] <- -g_i2

    Q_matrix[(2 * n + p + 3 + 2 + 1):(2 * n + p + 3 + 2 + n), 2 * n + p + 3] <- -g_i1
    Q_matrix[(2 * n + p + 3 + 2 + 1 + n):(2 * n + p + 3 + 2 + 2*n), 2 * n + p + 3] <- g_i2

    model$A<- rbind(cbind(diag(1, nrow = 2*n), -XX, 0, 0, 0 ,0 , diag(0, nrow = 2*n)),
                    cbind(diag(1, nrow = 2*n), -XX, 0, 0, 0, 0, diag(0, nrow = 2*n)),
                    c(rep(0, 2 * n), rep(0, p +5), rep(1,n)*S, rep(1, n)*(1 - S) ),
                    c(rep(0, 2 * n), rep(0, p +3),1, 1, rep(0, 2 * n ) )
                    )
    model$Q <- Q_matrix

    model$modelsense<-'min'
    ## tolerance_constraint enters here
    model$rhs<- c(rep(1 - tolerance_constraint, (dim(model$A)[1] - 2)/2),
                  rep(tolerance_constraint, (dim(model$A)[1] - 2)[1]/2),
                  max_treated_units, 1)
    model$sense<- c(rep('<=', (dim(model$A)[1] - 2)/2), rep('>', (dim(model$A)[1] - 2)/2), '<=', '<=')
    model$vtype<- c(rep('B', 2*n), rep('C', p+1), rep('B', 2), rep('C', 2 * n + 2))
    # Put bounds on the parameter space (If commented, parameter space = real line)
    model$ub<- c(rep(1,2*n), rep(B,1+p), rep(1, 2), rep(1, 2 + 2 * n))
    model$lb<- c(rep(0,2*n), rep(-B,p + 1), rep(0, 2), rep(0, 2 + 2 * n))
    saved_objective_one_direction <- c(rep(0, 2 * n), rep(0, 5 + p), -g_i1,  g_i2)
  }

  if(distance == 'parity'){

    ## constant1, constant2 here should be 0s
    model$obj <- c(rep(0, 2*n + p + 1), 0, 0, 0, 0, rep(0, 2 * n))
    Q_matrix <- matrix(0, nrow= 4*n + p + 5, ncol = 4*n + p + 5)
    Q_matrix[(2 * n + p + 3 + 2 + 1):(2 * n + p + 3 + 2 + n), 2 * n + p + 2] <- S/mean(S)
    Q_matrix[(2 * n + p + 3 + 2 + 1 + n):(2 * n + p + 3 + 2 + 2*n), 2 * n + p + 2] <- -(1 - S)/(1 - mean(S))

    Q_matrix[(2 * n + p + 3 + 2 + 1):(2 * n + p + 3 + 2 + n), 2 * n + p + 3] <- -S/mean(S)
    Q_matrix[(2 * n + p + 3 + 2 + 1 + n):(2 * n + p + 3 + 2 + 2*n), 2  * n + p + 3] <- (1 - S)/(1 - mean(S))

    model$A<- rbind(cbind(diag(1, nrow = 2*n), -XX, 0, 0, 0 ,0 , diag(0, nrow = 2*n)),
                    cbind(diag(1, nrow = 2*n), -XX, 0, 0, 0, 0, diag(0, nrow = 2*n)),
                    c(rep(0, 2 * n), rep(0, p +5), rep(1,n)*S, rep(1, n)*(1 - S) ),
                    c(rep(0, 2 * n), rep(0, p +3),1, 1, rep(0, 2 * n) )
    )

    model$Q <- Q_matrix

    model$modelsense<-'min'
    ## tolerance_constraint enters here
    model$rhs<- c(rep(1 - tolerance_constraint, (dim(model$A)[1] - 2)/2),
                  rep(tolerance_constraint, (dim(model$A)[1] - 2)/2), max_treated_units, 1)
    model$sense<- c(rep('<=',(dim(model$A)[1] - 2)/2), rep('>', (dim(model$A)[1] - 2)/2), '<=', '<=')
    model$vtype<- c(rep('B', 2*n), rep('C', p+1), rep('B', 2), rep('C', 2 * n + 2))
    # Put bounds on the parameter space (If commented, parameter space = real line)
    model$ub<- c(rep(1,2*n), rep(B,1+p), rep(1, 2), rep(1, 2 + 2* n))
    model$lb<- c(rep(0,2*n), rep(-B,p + 1), rep(0, 2), rep(0, 2 + 2 * n))
    saved_objective_one_direction <- c(rep(0, 2 * n), rep(0, 5 + p), -S/mean(S), (1 - S)/(1 - mean(S)))
  }
  ## Return objective when we do not consider the absolute value for unfairness
  if(two_directions == F & distance != 'envy') {
    model$obj <- saved_objective_one_direction
  }
  return(list(model, g_i1, g_i2))
}

## Inputs: as before
## output: estimate the maximum score with deterministic assignment
Est_max_score <- function(Y, X, D, S, propensity1, propensity2,
                          B=1, params = NA, model_only = F, tolerance_constraint = 10**(-7),
                          scale_Y = F,
                          cost_treatment  = 0, alpha = 1/2, g_i = NA,
                          max_treated_units, maxtime = 300, m1 = 0, m0 = 0,
                          additional_fairness_constraint = F,
                          parity_constraint = '>=', cores = 1, warm_start = NA) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  XX <- cbind(1, X)
  g_i_S <- rep(0, n)
  g_i_S2 <- rep(0,n)
  if(is.na(g_i)[1]){
    if(is.na(S)[1]){
      G_i1 <- (Y - m0) * (1-D) / ((1-propensity1)) + m0
      G_i2 <- (Y - cost_treatment - m1) * D / (propensity1) + m1
      g_i = G_i2 - G_i1
    }
    else{
      ## Propensity for the sensitive attribute
      G_i1 <- (Y - m0) * (1-D) * S / ((1-propensity1)*propensity2) + m0 * S/propensity2
      G_i2 <- (Y - cost_treatment - m1) * D * S/ (propensity1*propensity2) + m1 * S/propensity2
      g_i_S = G_i2 - G_i1

      G_i1 <- (Y - m0) * (1-D) * (1 - S) / ((1-propensity1)*(1 - propensity2)) + m0 * (1 - S)/(1 - propensity2)
      G_i2 <- (Y - cost_treatment - m1) * D * (1 - S)/ (propensity1*(1 - propensity2)) + m1 * (1 - S)/(1 - propensity2)
      g_i_S2 = G_i2 - G_i1
      g_i = alpha * g_i_S + (1 - alpha) * g_i_S2
    }
  }

  maximum_X <- max(apply(XX, 1, function(h) max(abs(h))))
  XX <- XX/maximum_X
  unique_values <- unique(XX)
  index_unique <- apply(XX,1, function(y) which(apply(unique_values, 1, function(x) all(y == x))))
  g_i_before <- g_i
  g_i <- sapply(c(1:dim(unique_values)[1]), function(x) sum(g_i[which(x == index_unique)]))
  g_i_S <- sapply(c(1:dim(unique_values)[1]), function(x) sum(g_i_S[which(x == index_unique)]))
  g_i_S2 <- sapply(c(1:dim(unique_values)[1]), function(x) sum(g_i_S2[which(x == index_unique)]))

  XX <- unique_values
  n_original <- n
  n <- length(g_i)
  n_indexes <- sapply(c(1:n), function(x) sum(index_unique == x))
  C <- B*max(apply(XX, 1, function(x) sum(abs(x))))
  XX <- as.matrix(XX/C)
  model  <- list()
  AA <- rbind(cbind(diag(1, nrow = n), -XX), cbind(diag(1, nrow = n), -XX), c(n_indexes, rep(0, p +1)))

  model$obj<- c(g_i, rep(0, p + 1))
  model$modelsense<-'max'
  rhs <- c(rep(1 - tolerance_constraint, dim(AA)[1]/2), rep(tolerance_constraint, dim(AA)[1]/2), max_treated_units)
  sense <- c(rep('<=', dim(AA)[1]/2), rep('>', dim(AA)[1]/2), '<=')
  model$vtype<- c(rep('B', n), rep('C', p+1))
  # Put bounds on the parameter space (If commented, parameter space = real line)
  model$ub<- c(rep(1,n), rep(B,1+p))
  model$lb<- c(rep(0,n), rep(-B,p + 1))
  if(additional_fairness_constraint){
    AA <-  rbind(AA, c(g_i_S - g_i_S2, rep(0, p + 1)))
    rhs <- c(rhs, 0)
    sense <- c(sense, parity_constraint)
  }
  model$rhs <- rhs
  model$sense <- sense
  model$A <- AA
  if(is.na(params)[1]){
    params<- list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9, TimeLimit = maxtime, BarConvTol = exp(-2),
                  Threads = cores, Disconnected=0, Heuristics = 0, NodefileStart = 0.5)}
  if(is.na(warm_start)[1] == F) model$start = warm_start
  if(model_only) return(model)
  .require_gurobi()
  result<- .call_gurobi(model, params = params)
  beta_hat <- result$x[(n+1):(n+p+1)]
  pi_est <- apply(cbind(1,X), 1, function(x) ifelse(x%*%beta_hat > 0, 1, 0))
  return(list(obj_est = sum(g_i_before*pi_est), g_i = g_i, result = result$x, pi = pi_est, beta = beta_hat))
}

## Estimate a linear probability model
## with linear programming (note: this leads to fast computations due to lack of integer variables)

Est_max_score_probabilistic <- function(Y, X, D, S, propensity1, propensity2, B=1, params = NA, model_only = F, tolerance_constraint = 10**(-7),
                                        scale_Y = F,
                                        cost_treatment  = 0, alpha = 1/2, g_i = NA,
                                        max_treated_units, maxtime = 300, m1 = 0, m0 = 0,
                                        additional_fairness_constraint = F,
                                        parity_constraint = '>=', cores = 1, intercept = T) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  XX = X
  if(intercept) XX <- cbind(1, X)
  p <- dim(XX)[2]
  g_i_S <- rep(0, n)
  g_i_S2 <- rep(0,n)

  if(is.na(g_i)[1]){
    if(is.na(S)[1]){
      G_i1 <- (Y - m0) * (1-D) / ((1-propensity1)) + m0
      G_i2 <- (Y - cost_treatment - m1) * D / (propensity1) + m1
      g_i = G_i2 - G_i1
    }
    else{
      ## Propensity for the sensitive attribute
      G_i1 <- (Y - m0) * (1-D) * S / ((1-propensity1)*propensity2) + m0 * S/propensity2
      G_i2 <- (Y - cost_treatment - m1) * D * S/ (propensity1*propensity2) + m1 * S/propensity2
      g_i_S = G_i2 - G_i1

      G_i1 <- (Y - m0) * (1-D) * (1 - S) / ((1-propensity1)*(1 - propensity2)) + m0 * (1 - S)/(1 - propensity2)
      G_i2 <- (Y - cost_treatment - m1) * D * (1 - S)/ (propensity1*(1 - propensity2)) + m1 * (1 - S)/(1 - propensity2)
      g_i_S2 = G_i2 - G_i1
      g_i = alpha * g_i_S + (1 - alpha) * g_i_S2
    }
  }


  unique_values <- unique(XX)
  index_unique <- apply(XX,1, function(y) which(apply(unique_values, 1, function(x) all(y == x))))
  g_i_before <- g_i
  g_i_S_before = g_i_S
  g_i_S2_before = g_i_S2

  g_i <- sapply(c(1:dim(unique_values)[1]), function(x) sum(g_i[which(x == index_unique)]))
  g_i_S <- sapply(c(1:dim(unique_values)[1]), function(x) sum(g_i_S[which(x == index_unique)]))
  g_i_S2 <- sapply(c(1:dim(unique_values)[1]), function(x) sum(g_i_S2[which(x == index_unique)]))

  XX <- unique_values
  n_original <- n
  n <- length(g_i)
  n_indexes <- sapply(c(1:n), function(x) sum(index_unique == x))

  XX <- as.matrix(XX)
  model  <- list()
  AA <- rbind(cbind(diag(1, nrow = n), -XX), c(n_indexes, rep(0, p )))

  model$obj<- c(g_i, rep(0, p ))
  model$modelsense<-'max'
  rhs <- c(rep(0, dim(AA)[1] - 1), max_treated_units)
  sense <- c(rep('=', dim(AA)[1] - 1 ), '<=')
  ## State here that variables are continuous
  model$vtype<- c(rep('C', n), rep('C', p))
  ## Use here bounds on [0,1] for each variable z_i
  model$ub<- c(rep(1,n), rep(B,p))
  model$lb<- c(rep(0,n), rep(-B,p ))
  if(additional_fairness_constraint){
    AA <-  rbind(AA, c(g_i_S - g_i_S2, rep(0, p )))
    rhs <- c(rhs, 0)
    sense <- c(sense, parity_constraint)
  }
  model$rhs <- rhs
  model$sense <- sense
  model$A <- AA
  if(is.na(params)[1]){
    params<- list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9, TimeLimit = maxtime, BarConvTol = exp(-2),
                  Threads = cores, Disconnected=0, Heuristics = 0, NodefileStart = 0.5, OutputFlag=0)}
  if(model_only) return(model)
  .require_gurobi()
  result<- .call_gurobi(model, params = params)
  beta_hat <- result$x[(n+1):(n+p)]
  if(intercept) {
    pi_est <- apply(cbind(1,X), 1, function(x) x%*%beta_hat)
  } else {
    pi_est <- apply(X, 1, function(x) x%*%beta_hat)
  }
  welfare_1 = sum(g_i_S_before * pi_est)
  welfare_0 = sum(g_i_S2_before * pi_est)
  return(list(obj_est = sum(g_i_before*pi_est), g_i = g_i, result = result$x, pi = pi_est, beta = beta_hat,
              welfare_1 = welfare_1, welfare_0 = welfare_0))
}


## A decision rule that assign probability beta1 to individuals below a certain threshold and beta0 to individuals above a certain threshold
## last three variables *(beta_1, beta_2, gamma_1, ... gamma_n) are the probabilities that we assign
## with gamma_i= xi beta_2 + beta_3
##      beta_2 + beta_3 <= 1
##              xi \in {0,1}
##              xi = 1\{X\beta + beta_0 S > 0} (hyperplane cut)


Est_threshold_probabilistic <- function(Y, X, D, S, propensity1, propensity2,
                          B=1, params = NA, model_only = F, tolerance_constraint = 10**(-7),
                          scale_Y = F,
                          cost_treatment  = 0, alpha = 1/2, g_i = NA,
                          max_treated_units, maxtime = 300, m1 = 0, m0 = 0,
                          additional_fairness_constraint = F,
                          parity_constraint = '>=', cores = 1) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  XX <- cbind(1, X)
  g_i_S <- rep(0, n)
  g_i_S2 <- rep(0,n)
  if(is.na(g_i)[1]){
    if(is.na(S)[1]){
      G_i1 <- (Y - m0) * (1-D) / ((1-propensity1)) + m0
      G_i2 <- (Y - cost_treatment - m1) * D / (propensity1) + m1
      g_i = G_i2 - G_i1
    }
    else{
      ## Propensity for the sensitive attribute
      G_i1 <- (Y - m0) * (1-D) * S / ((1-propensity1)*propensity2) + m0 * S/propensity2
      G_i2 <- (Y - cost_treatment - m1) * D * S/ (propensity1*propensity2) + m1 * S/propensity2
      g_i_S = G_i2 - G_i1

      G_i1 <- (Y - m0) * (1-D) * (1 - S) / ((1-propensity1)*(1 - propensity2)) + m0 * (1 - S)/(1 - propensity2)
      G_i2 <- (Y - cost_treatment - m1) * D * (1 - S)/ (propensity1*(1 - propensity2)) + m1 * (1 - S)/(1 - propensity2)
      g_i_S2 = G_i2 - G_i1
      g_i = alpha * g_i_S + (1 - alpha) * g_i_S2 ## Assume that the policy function does not depend on the sensitive attribute
    }
  }

  maximum_X <- max(apply(XX, 1, function(h) max(abs(h))))
  XX <- XX/maximum_X
  unique_values <- unique(XX)
  index_unique <- apply(XX,1, function(y) which(apply(unique_values, 1, function(x) all(y == x))))
  g_i_before <- g_i
  g_i <- sapply(c(1:dim(unique_values)[1]), function(x) sum(g_i[which(x == index_unique)]))
  g_i_S <- sapply(c(1:dim(unique_values)[1]), function(x) sum(g_i_S[which(x == index_unique)]))
  g_i_S2 <- sapply(c(1:dim(unique_values)[1]), function(x) sum(g_i_S2[which(x == index_unique)]))

  XX <- unique_values
  n_original <- n
  n <- length(g_i)
  n_indexes <- sapply(c(1:n), function(x) sum(index_unique == x))
  C <- B*max(apply(XX, 1, function(x) sum(abs(x))))
  XX <- as.matrix(XX/C)
  model  <- list()
  AA <- rbind(cbind(diag(1, nrow = n), -XX, 0, 0, diag(0, nrow = n)), cbind(diag(1, nrow = n), -XX, 0, 0, diag(0, nrow = n)),
              c(rep(0, n), rep(0, p +1), 0, 0, n_indexes),
              c(rep(0, p +1 + n), 1, 1, rep(0, n))
              )
  model$obj<- c(rep(0, p + n + 3), g_i)
  model$modelsense<-'max'
  rhs <- c(rep(1 - tolerance_constraint, dim(AA)[1]/2 - 1), rep(tolerance_constraint, dim(AA)[1]/2 - 1), max_treated_units, 1)
  sense <- c(rep('<=', dim(AA)[1]/2 - 1), rep('>', dim(AA)[1]/2 -1 ), '<=', '<=')
  model$vtype<- c(rep('B', n), rep('C', p + 3 + n))
  # Put bounds on the parameter space (If commented, parameter space = real line)
  model$ub<- c(rep(1,n), rep(B,1+p), rep(1, n + 2))
  model$lb<- c(rep(0,n), rep(-B,p + 1), rep(0, n + 2))
  if(additional_fairness_constraint){
    AA <-  rbind(AA, c(rep(0, n), rep(0, p + 3), g_i_S - g_i_S2))
    rhs <- c(rhs, 0)
    sense <- c(sense, parity_constraint)
  }
  model$rhs <- rhs
  model$sense <- sense
  model$A <- AA
  quadratic_const <- list()
  for(i in 1:n){
    quadcon <- list()
    quadcon$Qc <- matrix(0, nrow = 2 * n + p + 3, ncol = 2 * n + p + 3)
    quadcon$Qc[i, n + p + 2] <- -1
    quadcon$q <- rep(0, 2 * n + p + 3)
    quadcon$q[n + p + 3 + i] <- 1
    quadcon$q[n + p + 3] <- -1
    quadcon$rhs = 0
    quadcon$sense = '='
    quadratic_const[[i]] <- quadcon
  }
  model$quadcon <- quadratic_const
  if(is.na(params)[1]){
    params<- list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9, TimeLimit = maxtime, BarConvTol = exp(-2),
                  Threads = cores, Disconnected=0, Heuristics = 0, NodefileStart = 0.5)
    }
  if(model_only) return(model)
  .require_gurobi()
  result<- .call_gurobi(model, params = params)
  beta_hat <- result$x[(n+1):(n+p+1)]
  probs <- result$x[(n+p+2):(n+p+3)]
  pi_est <- apply(cbind(1,X), 1, function(x) ifelse(x%*%beta_hat > 0, probs[2] + probs[1], probs[2]))
  return(list(obj_est = sum(g_i_before*pi_est), g_i = g_i, result = result$x, pi = pi_est, beta = beta_hat, probs = probs))
}


# ---- from library.R ----
## Preliminaries:

## Functions support 4 notions of unfairness
## welfare: distance between E[Y(pi(X)) | S = 1] - E[Y(pi(X)) | S = 0]
## relative_welfare: E[(Y(1) - Y(0))\pi(X) | S = 1] - E[(Y(1) - Y(0))\pi(X) | S = 0]
## parity: E[\pi(X) | S = 1] - E[\pi(X) | S = 0]
## envy: see main text
## It also support notions in **absolute values**: in this case you need to specify two_directions = T


## Estimate fairness using the maximum score function
## Inputs : Y,X,D,S,propensity1 as above
##          p_s: probability P(S = 1)
##          scale_Y: whether the outcome is rescaled
##          discretization: number of elements in the grid
##          params: parameters for gurobi
##          frontier_objective: estimated frontier for initialization
##          mu_hatSD : conditional mean function on group (S,D)
##          warm_start: starting value
##          alpha_seq: sequence of alphas for the grid
##          no_parity_constrant = F: the policy is constant in the sensitive attribute
##          additional_fairness_constraint: if an additional constraint on the welfare of female >= welfare males is imposed
##          unique_values : whether we can collapse unique values to make computation faster (do not change the default)
##          distance: unfairness criterion
##          m0, m1: conditional mean functions
##          whether a probabilistic rule is chosen
##          numcores: number of cores for parallelization
##          two_directions: use the absolute value for the construction of the unfairness?
##                          only active if distance != envy
##                          note: you may set two_directions = F when the quadratic program is not PSD
##                          note: two_directions = F: the structure of the program is the same as two_directions = T
##                                                    the program still contains two binary variables
##                                                    indicating whether the unfairness is positive or negative
##                                                    but these variables do not enter in the objective function
## Return: results from the model

## structure of the variables
## (xi_11, ..., xi_n1, xi_10, ..., xi_n0, betas, constraint1, constraint2, u1, ..., uN)
## constraint1,2 not used for envy: indicate absolute values for unfairness and are binary (see below)

Est_fairnessMaxScore <- function(Y, X, D, S, propensity1, p_s,  scale_Y = T,
                                 discretization = floor(sqrt(length(Y))),
                                 cost_treatment = 0, params=NA,  frontier_objective,
                                 mu_hat11 = NA, mu_hat01 = NA, mu_hat00 = NA, mu_hat10 = NA,
                                 all_g_i, max_treated_units, maxtime = 300, warm_start = NA, alpha_seq,
                                 noparity_constraint = F,
                                 additional_fairness_constraint = F, unique_values = 1 - noparity_constraint,
                                 distance = 'envy', m0 = 0, m1 = 0,
                                 probabilistic = F, numcores = 10, two_directions = T, tolerance = 10**(-6),
                                 return_model = F){

  propensity2 <- p_s
  constant1 <- constant2 <- 0
  ## Construct the relative welfares
  G_i1 <- (Y - m0) * (1-D) * S / ((1-propensity1)*propensity2) + m0 * S/propensity2
  G_i2 <- (Y - cost_treatment - m1) * D * S/ (propensity1*propensity2) + m1 * S/propensity2
  g_i_S = G_i2 - G_i1
  G_i12 <- (Y - m0) * (1-D) * (1 - S) / ((1-propensity1)*(1 - propensity2)) + m0 * (1 - S)/(1 - propensity2)
  G_i22 <- (Y - cost_treatment - m1) * D * (1 - S)/ (propensity1*(1 - propensity2)) + m1 * (1 - S)/(1 - propensity2)
  g_i_S2 = G_i22 - G_i12

  ## Different effects for distance being either welfare or envy
  if(distance == 'welfare'){
    component1 <- g_i_S
    component2 <- g_i_S2
    constant1 <- sum(G_i1)
    constant2 <- sum(G_i12)
  }

  if(distance == 'relative_welfare'){
    component1 <- g_i_S
    component2 <- g_i_S2
    ## Consider **relative** improvement in welfare
    constant1 <- 0
    constant2 <- 0
    distance = 'welfare'
  }

  if(distance == 'envy'){
    component1 <- mu_hat01*S/propensity2 -  mu_hat00*S/propensity2 - g_i_S ## policy 1
    component2 <- mu_hat11*(1 - S)/(1 - propensity2) -  mu_hat10*(1 - S)/(1 - propensity2) - g_i_S2 ## policy
  }

  if(distance == 'parity'){
    component1 <- S/propensity2
    component2 <- (1 - S)/(1 - propensity2)
  }

  p <- dim(X)[2]
  model  <- list()
  ## Compute constraints and objective function
  if(probabilistic == F){
    initialize <- create_model_quadratic_program(Y = NA, X, D = NA, S, propensity1 = NA,
                                                 propensity2 = NA,B=1, params = NA,
                                                 tolerance_constraint = tolerance,  scale_Y = NA,
                                                 cost_treatment  = NA, alpha = NA, g_i1 = component1, g_i2 = component2,
                                                 max_treated_units = max_treated_units,
                                                 no_parity_constraint =noparity_constraint,
                                                 constant1 = constant1, constant2 = constant2, distance = distance,
                                                 two_directions = two_directions)
  } else {

    initialize <- create_model_quadratic_program_probabilistic(Y = NA, X, D = NA, S, propensity1 = NA,
                                                               propensity2 = NA,B=1, params = NA,
                                                               tolerance_constraint = tolerance,  scale_Y = NA,
                                                               cost_treatment  = NA, alpha = NA, g_i1 = component1, g_i2 = component2,
                                                               max_treated_units = max_treated_units,
                                                               no_parity_constraint =noparity_constraint,
                                                               constant1 = constant1, constant2 = constant2, distance = distance,
                                                               two_directions = two_directions)

  }
  init_Q <- initialize[[1]]$Q
  init_Q <- cbind(init_Q, matrix(0, nrow = nrow(init_Q), ncol = discretization))
  init_Q <- rbind(init_Q, matrix(0, nrow = discretization, ncol = ncol(init_Q)))
  if(two_directions) model$Q <- init_Q

  ## Add additional constraints on the function being pareto optimal
  model_init <- initialize[[1]]
  g_i <- initialize[[2]]
  g_i2 <- initialize[[3]]

  model$quadcon <- list()
  ## acc: it counts the number of binary variables we use
  ## acc = 2 for three notions: these two variables are binary and indicate
  ##         whether the unfairness is positive or negative
  ##         are positioned after betas (see created_model_quadratic_program for description)
  acc <- 0
  if(distance == 'welfare' | distance == 'parity' | distance == 'relative_welfare') acc <- 2
  ## Qudratic constraints
  ## indicating Pareto efficiency
  for(j in 1:discretization){
    constraint_q <- list()
    constraint_q$Qc <- matrix(0, nrow = length(model_init$obj) + discretization, ncol =length(model_init$obj) + discretization )
    constraint_q$Qc[1:(length(g_i)),j + 2*length(g_i) + dim(X)[2] + 1 + acc ] <- S*alpha_seq[j]*all_g_i ## consider the g_is for S == 1, assume all_g_i contains the true g_i
    constraint_q$Qc[(length(g_i) + 1):(2*length(g_i)),j + 2*length(g_i) + dim(X)[2] + 1 + acc ] <- (1 - S)*(1 - alpha_seq[j])*all_g_i ## consider the g_s for S==0
    constraint_q$q <- rep(0, length(model_init$obj) + discretization)

    constraint_q$q[length(model_init$obj) + j] <- -frontier_objective[j]
    constraint_q$rhs <- -10**(-8)
    constraint_q$sense <- '>'
    model$quadcon[[j]] <- constraint_q
  }

  ## Add additional quadratic constraint
  ## to indicate whether the welfare is positive or negative
  ## this is used when absolute value is considered
  if(distance == 'welfare'){
    constraint_q <- list()
    constraint_q$Qc <- matrix(0, nrow = length(model_init$obj) + discretization, ncol =length(model_init$obj) + discretization)
    constraint_q$Qc[1:(length(g_i)),2*length(g_i) + dim(X)[2] + 2 ] <- g_i
    constraint_q$Qc[(length(g_i) + 1):(2*length(g_i)),2*length(g_i) + dim(X)[2] + 2 ] <- -g_i2
    constraint_q$q <- c(rep(0, length(model_init$obj) - 2), constant1 - constant2, 0, rep(0, discretization))
    constraint_q$rhs <- -10**(-8)
    constraint_q$sense <- '>'
    model$quadcon[[j + 1]] <- constraint_q


    constraint_q <- list()
    constraint_q$Qc <- matrix(0, nrow = length(model_init$obj) + discretization, ncol =length(model_init$obj) + discretization)
    constraint_q$Qc[1:(length(g_i)),2*length(g_i) + dim(X)[2] + 3 ] <- -g_i
    constraint_q$Qc[(length(g_i) + 1):(2*length(g_i)),2*length(g_i) + dim(X)[2] + 3 ] <- g_i2
    constraint_q$q <- c(rep(0, length(model_init$obj) - 2), 0, -constant1 + constant2, rep(0, discretization))
    constraint_q$rhs <- -10**(-8)
    constraint_q$sense <- '>'
    model$quadcon[[j + 2]] <- constraint_q
  }

  ## Same with parity as above
  if(distance == 'parity'){
    constraint_q <- list()
    constraint_q$Qc <- matrix(0, nrow = length(model_init$obj) + discretization, ncol =length(model_init$obj) + discretization)
    constraint_q$Qc[1:(length(g_i)),2*length(g_i) + dim(X)[2] + 2 ] <- S/propensity2
    constraint_q$Qc[(length(g_i) + 1):(2*length(g_i)),2*length(g_i) + dim(X)[2] + 2 ] <- -(1- S)/(1 - propensity2)
    constraint_q$q <- c(rep(0, length(model_init$obj) - 2), constant1 - constant2, 0, rep(0, discretization))
    constraint_q$rhs <- -10**(-8)
    constraint_q$sense <- '>'
    model$quadcon[[j + 1]] <- constraint_q


    constraint_q <- list()
    constraint_q$Qc <- matrix(0, nrow = length(model_init$obj) + discretization, ncol =length(model_init$obj) + discretization)
    constraint_q$Qc[1:(length(g_i)),2*length(g_i) + dim(X)[2] + 3 ] <- -S/propensity2
    constraint_q$Qc[(length(g_i) + 1):(2*length(g_i)),2*length(g_i) + dim(X)[2] + 3 ] <- (1 - S)/(1 - propensity2)
    constraint_q$q <- c(rep(0, length(model_init$obj) - 2), 0, -constant1 + constant2, rep(0, discretization))
    constraint_q$rhs <- -10**(-8)
    constraint_q$sense <- '>'
    model$quadcon[[j + 2]] <- constraint_q
  }

  model$obj<- c(model_init$obj, rep(0, discretization))
  model$modelsense<- 'min'
  model$vtype<- c(model_init$vtype, rep('B', discretization))
  model$ub<- c(model_init$ub, rep(1, discretization))
  model$lb<- c(model_init$lb, rep(0, discretization))

  ## Construct the matrix A
  init_A <- model_init$A
  init_A <- cbind(init_A, matrix(0, nrow = nrow(init_A), ncol = discretization))
  init_A <- rbind(init_A, c(rep(0,length(model_init$obj) ), rep(1, discretization)))

  ## Use the warm start:
  ## this is the best allocation obtained from the estimated frontier

  if(is.na(warm_start)[1] == F) {
    add_vec <- rep(0, discretization)
    add_vec[warm_start[length(warm_start)]] <- 1
    model$start = c(warm_start[-length(warm_start)], add_vec)
  }
  rhs<- c(model_init$rhs, 0.9)
  sense<- c(model_init$sense, '>')

  if(additional_fairness_constraint){
    ## Additional constraint if true
    add_constraint <- c(S*all_g_i, -(1 - S)*all_g_i, rep(0, ncol(init_A) - 2*length(all_g_i)))
    init_A <- rbind(init_A, add_constraint)
    rhs<- c(rhs, 0)
    sense<- c(sense, '>=')

  }
  model$A<- init_A
  model$rhs <- rhs
  model$sense <- sense


  if(is.na(params)[1]){
    params<- list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9, TimeLimit =maxtime, BarConvTol = exp(-2), Threads = numcores, Disconnected=0,Heuristics=0, NodefileStart = 0.5)
  }
  if(return_model) return(list(model = model, initialize = initialize, frontier_objective = frontier_objective))
  ## return the results
  .require_gurobi()
  result<- .call_gurobi(model, params = params)
  n <- 2*length(g_i)
  beta_hat <- result$x[(n+1):(n+p+1)]
  if(distance == 'envy') alpha <- result$x[((n+p+2)):length(result$x)]
  if(distance %in% c('welfare', 'parity', 'relative_welfare') ) alpha <- result$x[((n+p+4)):length(result$x)]
  to_report <- alpha_seq[alpha == 1]
  ## Since X was rescaled by the same constant we can use X (not XX) here
  pi_est <- apply(cbind(1,X), 1, function(x) ifelse(x%*%beta_hat > 0, 1, 0))
  if(probabilistic) pi_est <- apply(cbind(1,X), 1, function(x) x%*%beta_hat)
  return(list(obj_est = result$objval, policies = pi_est, beta = beta_hat, alpha = to_report, results = result))
}

## Estimate maximum score with probabilistic assignments
##          note: here we assign two different probability if above or below a certain cutoff
## Inputs : Y,X,D,S,propensity1 as above
##          p_s: probability P(S = 1)
##          scale_Y: whether the outcome is rescaled
##          discretization: number of elements in the grid
##          params: parameters for gurobi
##          frontier_objective: estimated frontier for initialization
##          mu_hatSD : conditional mean function on group (S,D)
##          warm_start: starting value
##          alpha_seq: sequence of alphas for the grid
##          no_parity_constrant = F: the policy is constant in the sensitive attribute
##          additional_fairness_constraint: if an additional constraint on the welfare of female >= welfare males is imposed
##          unique_values : whether we can collapse unique values to make computation faster (do not change the default)
##          distance: unfairness criterion
##          m0, m1: conditional mean functions
##          whether a probabilistic rule is chosen
##          numcores: number of cores for parallelization
##          two directions: consider the fairness in absolute value? (note: this require PSD programming if T)
##          return_model: if true return the gurobi model
## Return: results from the model


## structure of the variables
## (xi_11, ..., xi_n1, xi_10, ..., xi_n0, betas, constraint1, constraint2, prob1, prob2, gamma11, ..., gamma1n,
##                                                                         gamma10, ..., gamma0n, u1, ..., uN)
## gamma_is: it is the probability of treatment for individual i if she had S = s
## gamma_is = xi_is prob1 + prob2
## prob1 + prob2 <= 1
## xi_is: as in the maximum score function

## constraint1,2 not used for envy: indicate absolute values for unfairness and are binary (see below)
Est_fairnessMaxScore_threshold_probabilistic <- function(Y, X, D, S, propensity1, p_s,  scale_Y = T,
                                 discretization = floor(sqrt(length(Y))),
                                 cost_treatment = 0, params=NA,  frontier_objective,
                                 mu_hat11, mu_hat01, mu_hat00, mu_hat10,
                                 all_g_i, max_treated_units, maxtime = 300, warm_start = NA, alpha_seq,
                                 noparity_constraint = F,
                                 additional_fairness_constraint = F, unique_values = 1 - noparity_constraint,
                                 distance = 'envy', m0 = NA, m1 = NA,
                                 probabilistic = T, numcores = 10,
                                 threshold_probabilistic = T,
                                 two_directions = T, return_model = F, tolerance = 10**(-6)){

  propensity2 <- p_s
  constant1 <- constant2 <- 0
  ## Construct the relative welfares
  G_i1 <- (Y - m0) * (1-D) * S / ((1-propensity1)*propensity2) + m0 * S/propensity2
  G_i2 <- (Y - cost_treatment - m1) * D * S/ (propensity1*propensity2) + m1 * S/propensity2
  g_i_S = G_i2 - G_i1
  G_i12 <- (Y - m0) * (1-D) * (1 - S) / ((1-propensity1)*(1 - propensity2)) + m0 * (1 - S)/(1 - propensity2)
  G_i22 <- (Y - cost_treatment - m1) * D * (1 - S)/ (propensity1*(1 - propensity2)) + m1 * (1 - S)/(1 - propensity2)
  g_i_S2 = G_i22 - G_i12

  ## Different effects for distance being either welfare or envy
  if(distance == 'welfare'){
   component1 <- g_i_S
   component2 <- g_i_S2
   constant1 <- sum(G_i1)
   constant2 <- sum(G_i12)
  }

  if(distance == 'relative_welfare'){
    component1 <- g_i_S
    component2 <- g_i_S2
    ## Consider **relative** improvement in welfare
    constant1 <- 0
    constant2 <- 0
    distance = 'welfare'
    }

  if(distance == 'envy'){
    component1 <- mu_hat01*S/propensity2 -  mu_hat00*S/propensity2 - g_i_S ## policy 1
    component2 <- mu_hat11*(1 - S)/(1 - propensity2) -  mu_hat10*(1 - S)/(1 - propensity2) - g_i_S2 ## policy
  }

  if(distance == 'parity'){
    component1 <- S/propensity2
    component2 <- (1 - S)/(1 - propensity2)
  }

  p <- dim(X)[2]
  model  <- list()
  ## Compute constraints and objective function

    initialize <- create_model_quadratic_program_threshold_probabilistic(Y = NA, X, D = NA, S, propensity1 = NA,
                                                               propensity2 = NA,B=1, params = NA,
                                                               tolerance_constraint = tolerance,  scale_Y = NA,
                                                               cost_treatment  = NA, alpha = NA, g_i1 = component1, g_i2 = component2,
                                                               max_treated_units = max_treated_units,
                                                               no_parity_constraint =noparity_constraint,
                                                               constant1 = constant1, constant2 = constant2, distance = distance,
                                                               two_directions = two_directions)


  init_Q <- initialize[[1]]$Q
  init_Q <- cbind(init_Q, matrix(0, nrow = nrow(init_Q), ncol = discretization))
  init_Q <- rbind(init_Q, matrix(0, nrow = discretization, ncol = ncol(init_Q)))
  ## if two directions we have a quadratic program
  ## otherwise we have a linear program
  if(distance != 'envy' & two_directions) model$Q <- init_Q

  ## Add additional constraints on the function being pareto optimal
  model_init <- initialize[[1]]
  g_i <- initialize[[2]]
  g_i2 <- initialize[[3]]

  model$quadcon <- list()
  acc <- 2
  if(distance == 'welfare' | distance == 'parity' | distance == 'relative_welfare') acc <- 4
  for(j in 1:discretization){
    constraint_q <- list()
    constraint_q$Qc <- matrix(0, nrow = length(model_init$obj) + discretization,
                              ncol =length(model_init$obj) + discretization )

    starting <- 2 * length(g_i) + dim(X)[2] + acc + 2
    constraint_q$Qc[starting:(starting + length(g_i) - 1),j + 4*length(g_i) + dim(X)[2] + 1 + acc ] <- S*alpha_seq[j]*all_g_i ## consider the g_is for S == 1, assume all_g_i contains the true g_i
    constraint_q$Qc[(starting + length(g_i)):(starting + 2*length(g_i) - 1),j + 4*length(g_i) + dim(X)[2] + 1 + acc ] <- (1 - S)*(1 - alpha_seq[j])*all_g_i ## consider the g_s for S==0
    constraint_q$q <- rep(0, length(model_init$obj) + discretization)


    constraint_q$q[length(model_init$obj) + j] <- -frontier_objective[j]
    constraint_q$rhs <- -10**(-8)
    constraint_q$sense <- '>'
    model$quadcon[[j]] <- constraint_q
  }

  ## Add additional quadratic constraint for welfare distance
  if(distance == 'welfare'){
  constraint_q <- list()
  constraint_q$Qc <- matrix(0, nrow = length(model_init$obj) + discretization, ncol =length(model_init$obj) + discretization)
  constraint_q$Qc[starting:(length(g_i) + starting - 1),2*length(g_i) + dim(X)[2] + 2 ] <- g_i
  constraint_q$Qc[(length(g_i) + starting):(2*length(g_i)  + starting - 1),2*length(g_i) + dim(X)[2] + 2 ] <- -g_i2
  constraint_q$q <- c(rep(0, 2 * length(g_i) + dim(X)[2] + 1 ), constant1 - constant2, 0,
                      rep(0, 2 + 2 * length(g_i)),
                      rep(0, discretization))
  constraint_q$rhs <- -10**(-8)
  constraint_q$sense <- '>'
  model$quadcon[[j + 1]] <- constraint_q


  constraint_q <- list()
  constraint_q$Qc <- matrix(0, nrow = length(model_init$obj) + discretization, ncol =length(model_init$obj) + discretization)
  constraint_q$Qc[starting:(length(g_i) + starting - 1),2*length(g_i) + dim(X)[2] + 3 ] <- -g_i
  constraint_q$Qc[(length(g_i) + starting):(2*length(g_i)  + starting - 1),2*length(g_i) + dim(X)[2] + 3 ] <- g_i2
  constraint_q$q <- c(rep(0, 2 * length(g_i) + dim(X)[2] + 1 ), 0, -constant1 + constant2, rep(0, 2 + 2 * length(g_i)), rep(0, discretization))
  constraint_q$rhs <- -10**(-8)
  constraint_q$sense <- '>'
  model$quadcon[[j + 2]] <- constraint_q
  }

  if(distance == 'parity'){
    constraint_q <- list()
    constraint_q$Qc <- matrix(0, nrow = length(model_init$obj) + discretization, ncol =length(model_init$obj) + discretization)
    constraint_q$Qc[starting:(length(g_i) + starting - 1),2*length(g_i) + dim(X)[2] + 2 ] <- S/propensity2
    constraint_q$Qc[(length(g_i) + starting):(2*length(g_i)  + starting - 1),2*length(g_i) + dim(X)[2] + 2 ] <- -(1- S)/(1 - propensity2)
    constraint_q$q <- c(rep(0, 2 * length(g_i) + dim(X)[2] + 1 ), constant1 - constant2, 0,  rep(0, 2 + 2 * length(g_i)), rep(0, discretization))
    constraint_q$rhs <- -10**(-8)
    constraint_q$sense <- '>'
    model$quadcon[[j + 1]] <- constraint_q


    constraint_q <- list()
    constraint_q$Qc <- matrix(0, nrow = length(model_init$obj) + discretization, ncol =length(model_init$obj) + discretization)
    constraint_q$Qc[starting:(length(g_i) + starting - 1),2*length(g_i) + dim(X)[2] + 3 ] <- -S/propensity2
    constraint_q$Qc[(length(g_i) + starting):(2*length(g_i)  + starting - 1),2*length(g_i) + dim(X)[2] + 3 ] <- (1 - S)/(1 - propensity2)
    constraint_q$q <- c(rep(0, 2 * length(g_i) + dim(X)[2] + 1 ), 0, -constant1 + constant2, rep(0, 2 + 2 * length(g_i)),  rep(0, discretization))
    constraint_q$rhs <- -10**(-8)
    constraint_q$sense <- '>'
    model$quadcon[[j + 2]] <- constraint_q
  }


  ## Add quadratic constraints on the xi_i and beta_i
  quadratic_const <- list()
  k = 1
  n = length(g_i)
  for(i in 1:n){

    ## Impose the constraints only on the S = 1 since the first n components generate zero effect on the objective for those with S = 0
    if(S[i] == 1){
    quadcon <- list()
    quadcon$Qc <- matrix(0, nrow = 4 * n + p + 5 + discretization, ncol = 4 * n + p + 5 + discretization)
    quadcon$Qc[i, 2 * n + p + 4] <- -1
    quadcon$q <- rep(0, 4 * n + p + 5 +  discretization)
    quadcon$q[2 * n + p + 5 + i] <- 1
    quadcon$q[2 * n + p + 5] <- -1
    quadcon$rhs = 0
    quadcon$sense = '='
    quadratic_const[[k]] <- quadcon
    k = k + 1
    }
  }

  for(i in (n+1):(2 * n)){

    ## Impose the constraints only on the S = 0 since the second n components generates no effect on the objective for S = 1
    if(S[i - n] == 0){
      quadcon <- list()
      quadcon$Qc <- matrix(0, nrow = 4 * n + p + 5+ discretization, ncol = 4 * n + p + 5+ discretization)
      quadcon$Qc[i, 2 * n + p + 4] <- -1
      quadcon$q <- rep(0, 4 * n + p + 5  + discretization)
      quadcon$q[2 * n + p + 5 + i] <- 1
      quadcon$q[2 * n + p + 5] <- -1
      quadcon$rhs = 0
      quadcon$sense = '='
      quadratic_const[[k]] <- quadcon
      k = k + 1
    }
  }

  model$quadcon <- append(model$quadcon, quadratic_const)


  model$obj<- c(model_init$obj, rep(0, discretization))
  model$modelsense<- 'min'
  model$vtype<- c(model_init$vtype, rep('B', discretization))
  model$ub<- c(model_init$ub, rep(1, discretization))
  model$lb<- c(model_init$lb, rep(0, discretization))

  ## Construct the matrix A
  init_A <- model_init$A
  init_A <- cbind(init_A, matrix(0, nrow = nrow(init_A), ncol = discretization))
  init_A <- rbind(init_A, c(rep(0,length(model_init$obj) ), rep(1, discretization)))


  if(is.na(warm_start)[1] == F) {
    add_vec <- rep(0, discretization)
    add_vec[warm_start[length(warm_start)]] <- 1
    model$start = c(warm_start[-length(warm_start)], add_vec)
  }
  rhs<- c(model_init$rhs, 0.9)
  sense<- c(model_init$sense, '>')

  if(additional_fairness_constraint){
   ## Additional constraint if true
   add_constraint <- rep(0, ncol(init_A))
   add_constraint[(2 * length(all_g_i) + dim(X)[1] + 1 + acc + 1):(4 * length(all_g_i) + dim(X)[1] + 1 + acc + 1) ]  <- c(S*all_g_i, -(1 - S)*all_g_i)
   init_A <- rbind(init_A, add_constraint)
   rhs<- c(rhs, 0)
   sense<- c(sense, '>=')

  }
  model$A<- init_A
  model$rhs <- rhs
  model$sense <- sense


  if(is.na(params)[1]){
    params<- list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9, TimeLimit =maxtime, BarConvTol = exp(-2), Threads = numcores, Disconnected=0,Heuristics=0, NodefileStart = 0.5)
    }
  if(return_model) return(model)
  .require_gurobi()
  result<- .call_gurobi(model, params = params)
  n <- 2*length(g_i)
  beta_hat <- result$x[(n+1):(n+p+1)]
  ## Probability of treatment
  probs <- result$x[(n+4 + p):(n+p+5)]
  if(distance == 'envy') alpha <- result$x[((n+p+4 + n)):length(result$x)]
  if(distance %in% c('welfare', 'parity', 'relative_welfare') ) alpha <- result$x[((n+p+6 + n)):length(result$x)]
  ## alpha active on the frontier
  to_report <- alpha_seq[alpha > 0.9]
  ## Since X was rescaled by the same constant we can use X (not XX) here
  pi_est <- apply(cbind(1,X), 1, function(x) ifelse(x%*%beta_hat > 0, probs[2] + probs[1], probs[2]))

  return(list(obj_est = result$objval, policies = pi_est, beta = beta_hat, alpha = to_report, probs = probs, results = result))
}



## Estimate the Pareto frontier with MILP
## argument: as before
## output: list

.estimate_Pareto_frontier_core <- function(Y, X, D, S, propensity1, propensity2, scale_Y = TRUE,                                             discretization = floor(sqrt(length(Y))),  cost_treatment = 0, params=NA,
                                             max_treated_units = length(Y),
                                             numcores = 10, maxtime = 300, alpha_seq = seq(from = 0, to = 1, length = discretization),
                                             m1 = 0, m0 = 0, additional_fairness_constraint = F,
                                             parity_constraint = '>=', probabilistic = F,
                                             threshold_probabilistic = F, parallel = T, tolerance = 10**(-3)){
  if (!requireNamespace('foreach', quietly = TRUE)) .abort('Parallel frontier code requires the `foreach` package.')
  foreach <- foreach::foreach
  `%do%` <- foreach::`%do%`
  `%dopar%` <- foreach::`%dopar%`
  if(parallel){
    if (!requireNamespace('doParallel', quietly = TRUE)) .abort('Parallel frontier code requires the `doParallel` package.')
    doParallel::registerDoParallel(numcores)
    res <- foreach(i = alpha_seq, .combine = rbind, .export = c('Y', 'D', 'X', 'S', 'propensity1',
                                                              'propensity2', 'params', 'scale_Y',
                                                              'cost_treatment', 'max_treated_units',
                                                              'maxtime', 'm1', 'm0', 'additional_fairness_constraint',
                                                              'parity_constraint'))%dopar%{
                  if(probabilistic == F & threshold_probabilistic == F){
                  result <- Est_max_score(Y, X, D, S, propensity1, propensity2,B=1, params = params, tolerance_constraint = tolerance,
                                          scale_Y = scale_Y, additional_fairness_constraint = additional_fairness_constraint,
                                          parity_constraint = parity_constraint,
                  cost_treatment  = cost_treatment, alpha = i, max_treated_units = max_treated_units, maxtime = maxtime, m1 = m1, m0 = m0)

                  } else if (probabilistic == T & threshold_probabilistic == F) {
                    result <- Est_max_score_probabilistic(Y, X, D, S, propensity1, propensity2,B=1, params = params, tolerance_constraint = tolerance,
                                            scale_Y = scale_Y, additional_fairness_constraint = additional_fairness_constraint,
                                            parity_constraint = parity_constraint,
                                            cost_treatment  = cost_treatment, alpha = i, max_treated_units = max_treated_units, maxtime = maxtime, m1 = m1, m0 = m0)
                  } else {
                    result <- Est_threshold_probabilistic(Y, X, D, S, propensity1, propensity2,B=1, params = params, tolerance_constraint = tolerance,
                                                          scale_Y = scale_Y, additional_fairness_constraint = additional_fairness_constraint,
                                                          parity_constraint = parity_constraint,
                                                          cost_treatment  = cost_treatment, alpha = i, max_treated_units = max_treated_units, maxtime = maxtime, m1 = m1, m0 = m0)

                  }
                  c(result[[2]], result[[1]], result[[3]], result[[4]], length(result[[2]]), length(result[[3]]))
                                                              }

  } else {

    res <- foreach(i = alpha_seq, .combine = rbind, .export = c('Y', 'D', 'X', 'S', 'propensity1',
                                                                'propensity2', 'params', 'scale_Y',
                                                                'cost_treatment', 'max_treated_units',
                                                                'maxtime', 'm1', 'm0', 'additional_fairness_constraint',
                                                                'parity_constraint'))%do%{
                                                                  if(probabilistic == F & threshold_probabilistic == F){
                             result <- Est_max_score(Y, X, D, S, propensity1, propensity2,B=1, params = params, tolerance_constraint = tolerance,
                                 scale_Y = scale_Y, additional_fairness_constraint = additional_fairness_constraint,
                                     parity_constraint = parity_constraint, cost_treatment  = cost_treatment, alpha = i,
                                 max_treated_units = max_treated_units, maxtime = maxtime, m1 = m1, m0 = m0, cores = numcores) } else if (probabilistic == T & threshold_probabilistic == F) {
                       result <- Est_max_score_probabilistic(Y, X, D, S, propensity1, propensity2,B=1, params = params, tolerance_constraint =tolerance,
                                                                                                                                    scale_Y = scale_Y, additional_fairness_constraint = additional_fairness_constraint,
                                                                                                                                    parity_constraint = parity_constraint,
                                                                                                                                    cost_treatment  = cost_treatment, alpha = i, max_treated_units = max_treated_units, maxtime = maxtime, m1 = m1, m0 = m0,
                                                             cores = numcores)
                                 } else {   result <- Est_threshold_probabilistic(Y, X, D, S, propensity1, propensity2,B=1, params = params, tolerance_constraint =tolerance,
                                                                                  scale_Y = scale_Y, additional_fairness_constraint = additional_fairness_constraint,
                                                                                  parity_constraint = parity_constraint,
                                                                                  cost_treatment  = cost_treatment, alpha = i, max_treated_units = max_treated_units, maxtime = maxtime, m1 = m1, m0 = m0, cores = numcores)

                                                                                            }
                                                                  c(result[[2]], result[[1]],
                                                                    result[[3]], result[[4]], length(result[[2]]),
                                                                    length(result[[3]]))



                                                                }

  }

  if(threshold_probabilistic == F){
  nn <- res[1, dim(res)[2] - 1]
  nn2 <- res[1, dim(res)[2]]
  res <- res[, -dim(res)[2]]
  aa1 <- res[, 1:nn]
  aa2 <- res[, nn + 1]
  return(list(g_i = res[, 1:nn], objective = res[, nn + 1], results = res[, c((nn + 2):(nn + 1 + nn2))],
              policies = res[, c((nn + 2 + nn2):(dim(res)[2]-1))], beta = res[, c((nn + 2 + nn):(nn + 1 + nn2))]))
  } else {
    nn <- res[1, dim(res)[2] - 1]
    nn2 <- res[1, dim(res)[2]]
    res <- res[, -dim(res)[2]]
    aa1 <- res[, 1:nn]
    aa2 <- res[, nn + 1]
    return(list(g_i = res[, 1:nn], objective = res[, nn + 1], results = res[, c((nn + 2):(nn + 1 + nn2))],
                policies = res[, c((nn + 2 + nn2):(dim(res)[2]-1))], beta = res[, c((nn + 1 + nn + 1):(2 * nn + 1 + dim(X)[2] + 1 ))],
                probs = res[, c((2 * nn + 2 + dim(X)[2] + 1 ):(2 * nn + 2 + dim(X)[2] + 2 ) )]
                ))
    }
}


## Wrapper function for MILP
## additional arguments
##                      quick_run: return approximate value based on Pareto frontier calculation?
##                      two_directions: if T uses unfairness in abs value else use unfairness without absolute value (dff between advantage and dis group)
##                      two_directions = F can be used if the quadratic program is not PSD
##                      parallel : compute the frontier in parallel (more RAM requirement)
##                      tolerance_frontier: tolerance used to compute the frontier (dafault is 10**(-3) which may be conservative). You may increase it till 10**(-6) but check
##                                          solutions if you do.
##                      tolerance_optimization: tolerance for the optimization, must be smaller than the tolerance for the frontier
Est_objective_estimandMaxscore <- function(Y, X, D, S, propensity1, propensity2, scale_Y = T,
                                           discretization = floor(sqrt(length(Y))),cost_treatment = 0, params=NA,
                                           mu_hat11, mu_hat01, mu_hat00, mu_hat10, model_only = F,
                                           max_treated_units = length(Y),
                                           maxtime1 = 300, maxtime2 = 100,
                                           alpha_seq = seq(from = 0, to = 1, length = discretization),
                                           m1 = 0, m0 = 0, quick_run = F, no_parity_constraint = F,
                                           additional_fairness_constraint = F,
                                           parity_constraint = '>=', frontier = NA, unique_values = 1 - no_parity_constraint,
                                           distance = 'envy', probabilistic = F, numcores = 10, threshold_probabilistic = F,
                                           two_directions = T, parallel = T, tolerance_frontier = 10**(-3), tolerance_optimization = 10**(-6),
                                           return_frontier = F, frontier_slack = 0){
  if(is.na(frontier)[1]){

  frontier <- .estimate_Pareto_frontier_core(Y = Y, X = X, D = D, S = S, propensity1= propensity1, propensity2 = propensity2,
                                      scale_Y = scale_Y,  discretization = discretization,
                                       cost_treatment = cost_treatment,
                                       params = params, max_treated_units = max_treated_units,
                                       maxtime = maxtime2, alpha_seq = alpha_seq, m1 = m1, m0 = m0,
                                       additional_fairness_constraint = additional_fairness_constraint,
                                       parity_constraint = parity_constraint, probabilistic = probabilistic,
                                      numcores = numcores, threshold_probabilistic = threshold_probabilistic,
                                      parallel = parallel, tolerance = tolerance_frontier)

  }
  if(return_frontier) return(frontier)
  frontier_objective <- frontier[[2]] - frontier_slack
  results_frontier <- frontier[[4]] ## Store the policies
  results_frontier_collapsed <- frontier[[3]]

  G_i1 <- (Y - m0) * (1-D) * S / ((1-propensity1)*propensity2) + m0 * S/propensity2
  G_i2 <- (Y - cost_treatment - m1) * D * S/ (propensity1*propensity2) + m1 * S/propensity2
  g_i_S = G_i2 - G_i1

  G_i12 <- (Y - m0) * (1-D) * (1 - S) / ((1-propensity1)*(1 - propensity2)) + m0 * (1 - S)/(1 - propensity2)
  G_i22 <- (Y - cost_treatment - m1) * D * (1 - S)/ (propensity1*(1 - propensity2)) + m1 * (1 - S)/(1 - propensity2)
  g_i_S2 = G_i22 - G_i12


  all_g_i = g_i_S + g_i_S2 ## Save the welfare criterion

  ## Do not consider a probabilistic threshold here
  if(threshold_probabilistic == F){
  beta <- frontier[[5]]
  if(unique_values == F){
  XX1 <- as.matrix(cbind(1, 1, X[,-1]))
  XX0 <- as.matrix(cbind(1, 0, X[,-1]))
  } else {
  XX1 <- as.matrix(cbind(1, X))
  XX0 <- as.matrix(cbind(1, X))
  }

  ## Compute a warm-start for the MILP
  policy1 <- t(apply(beta, 1, function(x) sapply(XX1%*%x, function(y) ifelse(y > 0, 1, 0))))
  policy0 <- t(apply(beta, 1, function(x) sapply(XX0%*%x, function(y) ifelse(y > 0, 1, 0))))


  ## Compute the distance for the welfare-based fairness

  if(distance == 'welfare'){
    welfare1 <- apply(policy1, 1, function(x) sum(g_i_S*x) + sum(G_i1))
    welfare0 <- apply(policy0, 1, function(x) sum(g_i_S2*x) + sum(G_i12))
    objective_warm_starts_welfare <-  welfare0 - welfare1
    if(two_directions) objective_warm_starts_welfare <- abs(welfare1 - welfare0)
    least_unfair = which(objective_warm_starts_welfare == min(objective_warm_starts_welfare))
    least_unfair = least_unfair[which.min(abs(least_unfair - discretization/2))]
    indicator <- ifelse(welfare1[least_unfair] -  welfare0[least_unfair] > 0, 1, 0)
    warm_start = c(policy1[least_unfair,], policy0[least_unfair,], beta[least_unfair,], indicator, 1 - indicator, least_unfair)
  }

  if(distance == 'relative_welfare'){
    welfare1 <- apply(policy1, 1, function(x) sum(g_i_S*x))
    welfare0 <- apply(policy0, 1, function(x) sum(g_i_S2*x))
    objective_warm_starts_welfare <-  welfare0 - welfare1
    if(two_directions) objective_warm_starts_welfare <- abs(welfare1 - welfare0)
    least_unfair = which(objective_warm_starts_welfare == min(objective_warm_starts_welfare))
    least_unfair = least_unfair[which.min(abs(least_unfair - discretization/2))]
    indicator <- ifelse(welfare1[least_unfair] -  welfare0[least_unfair] > 0, 1, 0)
    warm_start = c(policy1[least_unfair,], policy0[least_unfair,], beta[least_unfair,], indicator, 1 - indicator, least_unfair)
  }

  if(distance == 'parity'){
    w1 <- apply(policy1, 1, function(x) sum(S*x/mean(S)))
    w0 <- apply(policy0, 1, function(x) sum((1 - S)*x/(1 - mean(S))))
    objective_warm_starts_welfare <- w0 - w1
    if(two_directions)  objective_warm_starts_welfare <- abs(w1 - w0)
    least_unfair = which(objective_warm_starts_welfare == min(objective_warm_starts_welfare))
    least_unfair = least_unfair[which.min(abs(least_unfair - discretization/2))]
    indicator <- ifelse(w1[least_unfair] -  w0[least_unfair] > 0, 1, 0)
    warm_start = c(policy1[least_unfair,], policy0[least_unfair,], beta[least_unfair,], indicator, 1 - indicator, least_unfair)

  }

  ## Compute the distance for the envy-based fairness

  if(distance == 'envy'){
  welfare1 <- apply(policy1, 1, function(x) sum(g_i_S*x))
  welfare0 <- apply(policy0, 1, function(x) sum(g_i_S2*x))
  objective_warm_starts1 <- apply(policy1, 1, function(x) sum(mu_hat01*x*S)/propensity2 +  sum(mu_hat00*(1 - x)*S)/propensity2) -  welfare0
  objective_warm_starts2 <- apply(policy0, 1, function(x) sum(mu_hat11*x*(1 - S))/(1 - propensity2) +  sum(mu_hat10*(1 - x)*(1 - S))/(1 - propensity2)) - welfare1
  objective_warm_starts_envy <- objective_warm_starts1 + objective_warm_starts2
  least_unfair = which(objective_warm_starts_envy == min(objective_warm_starts_envy))
  least_unfair = least_unfair[which.min(abs(least_unfair - discretization/2))]
  warm_start = c(policy1[least_unfair,], policy0[least_unfair,], beta[least_unfair,], least_unfair)
  }
  } else {
    ## Compute warm start for the probabilistic with threshold (note: warm start only imporves computational time but not performance)

    beta <- frontier[[5]]
    probs <- frontier[[6]]
    if(unique_values == F){
      XX1 <- as.matrix(cbind(1, 1, X[,-1]))
      XX0 <- as.matrix(cbind(1, 0, X[,-1]))
    } else {
      XX1 <- as.matrix(cbind(1, X))
      XX0 <- as.matrix(cbind(1, X))
    }

    ## Compute a warm-start for the binary indicator and the probabilistic assignments
    xi1 <- t(apply(cbind(beta,probs),
                   1, function(x) sapply(XX1%*%(x[1:dim(beta)[2]]), function(y) ifelse(y > 0, 1, 0))))
    policy1 <- t(apply(cbind(beta,probs),
                       1, function(x) sapply(XX1%*%(x[1:dim(beta)[2]]), function(y) ifelse(y > 0, x[dim(beta)[2] + 1] +
                                                                                           x[dim(beta)[2] + 2], x[dim(beta)[2] + 2]
                                                                                         ))))
    xi0 <- t(apply(cbind(beta,probs),
                   1, function(x) sapply(XX0%*%(x[1:dim(beta)[2]]), function(y) ifelse(y > 0, 1, 0))))
    policy0 <- t(apply(cbind(beta,probs),
                       1, function(x) sapply(XX0%*%(x[1:dim(beta)[2]]), function(y) ifelse(y > 0, x[dim(beta)[2] + 1] +
                                                                                           x[dim(beta)[2] + 2], x[dim(beta)[2] + 2]))))


    ## Compute the distance for the welfare-based fairness

    if(distance == 'welfare'){
      welfare1 <- apply(policy1, 1, function(x) sum(g_i_S*x) + sum(G_i1))
      welfare0 <- apply(policy0, 1, function(x) sum(g_i_S2*x) + sum(G_i12))
      objective_warm_starts_welfare <- welfare0 - welfare1
      if(two_directions) objective_warm_starts_welfare <- abs(welfare1 - welfare0)
      least_unfair = which(objective_warm_starts_welfare == min(objective_warm_starts_welfare))
      least_unfair = least_unfair[which.min(abs(least_unfair - discretization/2))]
      indicator <- ifelse(welfare1[least_unfair] -  welfare0[least_unfair] > 0, 1, 0)
      warm_start = c(xi1[least_unfair,], xi0[least_unfair,], beta[least_unfair,], indicator, 1 - indicator, probs[least_unfair,],
                     policy1[least_unfair,], policy0[least_unfair,], least_unfair)
    }

    if(distance == 'relative_welfare'){
      welfare1 <- apply(policy1, 1, function(x) sum(g_i_S*x))
      welfare0 <- apply(policy0, 1, function(x) sum(g_i_S2*x))
      objective_warm_starts_welfare <- welfare0 - welfare1
      if(two_directions) objective_warm_starts_welfare <- abs(welfare1 - welfare0)
      least_unfair = which(objective_warm_starts_welfare == min(objective_warm_starts_welfare))
      least_unfair = least_unfair[which.min(abs(least_unfair - discretization/2))]
      indicator <- ifelse(welfare1[least_unfair] -  welfare0[least_unfair] > 0, 1, 0)
      warm_start = c(xi1[least_unfair,], xi0[least_unfair,], beta[least_unfair,], indicator, 1 - indicator, probs[least_unfair,],
                     policy1[least_unfair,], policy0[least_unfair,], least_unfair)
    }

    if(distance == 'parity'){
      w1 <- apply(policy1, 1, function(x) sum(S*x/mean(S)))
      w0 <- apply(policy0, 1, function(x) sum((1 - S)*x/(1 - mean(S))))
      objective_warm_starts_welfare <- w0 - w1
      if(two_directions) objective_warm_starts_welfare <- abs(w1 - w0)
      least_unfair = which(objective_warm_starts_welfare == min(objective_warm_starts_welfare))
      least_unfair = least_unfair[which.min(abs(least_unfair - discretization/2))]
      indicator <- ifelse(w1[least_unfair] -  w0[least_unfair] > 0, 1, 0)
      warm_start = c(xi1[least_unfair,], xi0[least_unfair,], beta[least_unfair,], indicator, 1 - indicator, probs[least_unfair,],
                     policy1[least_unfair,], policy0[least_unfair,], least_unfair)

    }

    ## Compute the distance for the envy-based fairness

    if(distance == 'envy'){
      welfare1 <- apply(policy1, 1, function(x) sum(g_i_S*x))
      welfare0 <- apply(policy0, 1, function(x) sum(g_i_S2*x))
      objective_warm_starts1 <- apply(policy1, 1, function(x) sum(mu_hat01*x*S)/propensity2 +  sum(mu_hat00*(1 - x)*S)/propensity2) -  welfare0
      objective_warm_starts2 <- apply(policy0, 1, function(x) sum(mu_hat11*x*(1 - S))/(1 - propensity2) +  sum(mu_hat10*(1 - x)*(1 - S))/(1 - propensity2)) - welfare1
      objective_warm_starts_envy <- objective_warm_starts1 + objective_warm_starts2
      least_unfair = which(objective_warm_starts_envy == min(objective_warm_starts_envy))
      least_unfair = least_unfair[which.min(abs(least_unfair - discretization/2))]
      warm_start = c(xi1[least_unfair,], xi0[least_unfair,], beta[least_unfair,], probs[least_unfair,],
                     policy1[least_unfair,], policy0[least_unfair,], least_unfair)
    }

  }

  if(quick_run) return(list(result = warm_start, frontier = frontier))
  if(threshold_probabilistic == F) {
  result <- Est_fairnessMaxScore(Y, X, D, S, propensity1 = propensity1, p_s = propensity2,  scale_Y,
                                 discretization, cost_treatment,
                                 params, frontier_objective = frontier_objective, mu_hat11 = mu_hat11,
                                 mu_hat01 = mu_hat01, mu_hat00 = mu_hat00, mu_hat10 = mu_hat10,
                                 all_g_i, max_treated_units = max_treated_units, maxtime = maxtime1,
                                 warm_start = warm_start, alpha_seq = alpha_seq, noparity_constraint = no_parity_constraint,
                                 additional_fairness_constraint = additional_fairness_constraint,
                                 unique_values = unique_values,
                                 distance = distance, m0 = m0, m1 = m1, probabilistic = probabilistic,
                                 numcores = numcores, two_directions = two_directions, tolerance = tolerance_optimization)
  } else {
    result <- Est_fairnessMaxScore_threshold_probabilistic(Y, X, D, S, propensity1 = propensity1, p_s = propensity2,  scale_Y,
                                                 discretization, cost_treatment,
                                                 params, frontier_objective = frontier_objective, mu_hat11 = mu_hat11,
                                                 mu_hat01 = mu_hat01, mu_hat00 = mu_hat00, mu_hat10 = mu_hat10,
                                                 all_g_i, max_treated_units = max_treated_units, maxtime = maxtime1,
                                                 warm_start = warm_start, alpha_seq = alpha_seq, noparity_constraint = no_parity_constraint,
                                                 additional_fairness_constraint = additional_fairness_constraint,
                                                 unique_values = unique_values,
                                                 distance = distance, m0 = m0, m1 = m1, probabilistic = probabilistic,
                                                 numcores = numcores, two_directions = two_directions, tolerance = tolerance_optimization)
  }


  ## Note: equity is for the full welfare function (no relative improvement)
  return(list(result = result, frontier = frontier))

}

## Maximize empirical welfare with a probabilistic decision rule
## and with fairness constraints on the parameter space
## Unfairness_bound: maximal unfairness allowed
## Note: this computes the linear probabilistic allocaiton (linear program only)

welfare_fairness_constraints_probabilistic <- function(Y, X, D, S, propensity1, propensity2,
                                                       B=1, params = NA, model_only = F,
                                                       tolerance_constraint = 10**(-7),
                          scale_Y = F,
                          cost_treatment  = 0, g_i = NA,
                          max_treated_units, maxtime = 300, m1 = 0, m0 = 0,
                          additional_fairness_constraint = F,cores = 1, UnFairness_bound = Inf, distance,
                          mu_hat11 = NA, mu_hat01 = NA ,
                          mu_hat00 = NA, mu_hat10 = NA, two_directions = T, alpha = NA,
                          intercept = T) {

  if(is.na(alpha)[1]) alpha = propensity2
  if(distance == 'welfare') stop('distance= welfare not supported for probabilistic rule with constraints. Try relative_welfare')
  n <- dim(X)[1]
  p <- dim(X)[2]
  XX = X
  if(intercept) XX <- cbind(1, X)
  p = dim(XX)[2]
  g_i_S <- rep(0, n)
  g_i_S2 <- rep(0,n)
  if(is.na(g_i)[1]){
      ## Propensity for the sensitive attribute
      G_i1 <- (Y - m0) * (1-D) * S / ((1-propensity1)*propensity2) + m0 * S/propensity2
      G_i2 <- (Y - cost_treatment - m1) * D * S/ (propensity1*propensity2) + m1 * S/propensity2
      g_i_S = G_i2 - G_i1

      G_i12 <- (Y - m0) * (1-D) * (1 - S) / ((1-propensity1)*(1 - propensity2)) + m0 * (1 - S)/(1 - propensity2)
      G_i22 <- (Y - cost_treatment - m1) * D * (1 - S)/ (propensity1*(1 - propensity2)) + m1 * (1 - S)/(1 - propensity2)
      g_i_S2 = G_i22 - G_i12
      g_i = alpha * g_i_S + (1 - alpha) * g_i_S2 ## estimated welfare
    }


  constant1 <- constant2 <- 0
  if(distance == 'relative_welfare'){
    ## Unfairness is the sum of component1 and component2
    component1 <- g_i_S
    component2 <- -g_i_S2
  }

  if(distance == 'envy'){
    component1 <- mu_hat01*S/propensity2 -  mu_hat00*S/propensity2 - g_i_S ## policy 1
    component2 <- mu_hat11*(1 - S)/(1 - propensity2) -  mu_hat10*(1 - S)/(1 - propensity2) - g_i_S2 ## policy
    }

  if(distance == 'parity'){
    component1 <- S/propensity2
    component2 <- - (1 - S)/(1 - propensity2)
  }
  UnFairness_before <- component1 + component2
  unique_values <- unique(XX)
  index_unique <- apply(XX,1, function(y) which(apply(unique_values, 1, function(x) all(y == x))))
  g_i_before <- g_i
  g_i <- sapply(c(1:dim(unique_values)[1]), function(x) sum(g_i[which(x == index_unique)]))
  component1 <- sapply(c(1:dim(unique_values)[1]), function(x) sum(component1[which(x == index_unique)]))
  component2 <- sapply(c(1:dim(unique_values)[1]), function(x) sum(component2[which(x == index_unique)]))
  UnFairness <- component1 + component2

  XX <- unique_values
  n_original <- n
  n <- length(g_i)
  n_indexes <- sapply(c(1:n), function(x) sum(index_unique == x))
  model  <- list()
  AA <- rbind(cbind(diag(1, nrow = n), -XX), c(UnFairness, rep(0, p )),
              c(UnFairness, rep(0, p )), c(n_indexes, rep(0, p )))

  model$obj<- c(g_i, rep(0, p ))
  model$modelsense<-'max'
  if(two_directions){
    rhs <- c(rep(0, dim(AA)[1] - 3), UnFairness_bound, - UnFairness_bound, max_treated_units)
  } else {
    ## Only impose the constraint in one direction if two_directions = F
    rhs <- c(rep(0, dim(AA)[1] - 3), Inf, - UnFairness_bound, max_treated_units)
  }
  sense <- c(rep('=', dim(AA)[1] - 3), '<=', '>=', '<=')
  model$vtype<- c(rep('C', n), rep('C', p))
  # Put bounds on the parameter space (If commented, parameter space = real line)
  model$ub<- c(rep(1,n), rep(B,p))
  model$lb<- c(rep(0,n), rep(-B,p ))

  model$rhs <- rhs
  model$sense <- sense
  model$A <- as.matrix(AA)

  if(is.na(params)[1]){
    params<- list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9, TimeLimit = maxtime, BarConvTol = exp(-2),
                  Threads = cores, Disconnected=0, Heuristics = 0, NodefileStart = 0.5)}
  .require_gurobi()
  result<- .call_gurobi(model, params = params)
  beta_hat <- result$x[(n+1):(n+p)]
  if(intercept) {
    pi_est <- apply(cbind(1,X), 1, function(x) x%*%beta_hat)
  } else {
    pi_est <- apply(X, 1, function(x) x%*%beta_hat)
  }
  return(list(obj_est = sum(g_i_before*pi_est), g_i = g_i, result = result$x, pi = pi_est, beta = beta_hat,
              UnFairness = sum(UnFairness_before * pi_est), group1 = mean(g_i_before * pi_est * S/propensity2),
              group2 = mean(g_i_before * pi_est * (1 - S)/(1 - propensity2))))
}


## Maximize empirical welfare with a deterministic decision rule
## and with fairness constraints on the parameter space
## Unfairness_bound: maximal unfairness allowed
welfare_fairness_constraints_MS <- function(Y, X, D, S, propensity1, propensity2,
                                                       B=1, params = NA, model_only = F,
                                                       tolerance_constraint = 10**(-7),
                                                       scale_Y = F,
                                                       cost_treatment  = 0, alpha = 1/2, g_i = NA,
                                                       max_treated_units, maxtime = 300, m1 = 0, m0 = 0,
                                                       additional_fairness_constraint = F,
                                                       parity_constraint = '>=', cores = 1, UnFairness_bound = Inf,
                                            distance,
                                            mu_hat11 = NA, mu_hat01 = NA ,
                                            mu_hat00 = NA, mu_hat10 = NA, two_directions = T) {
  if(distance == 'welfare') stop('distance = welfare not supported for probabilistic rule with constraints. Try relative_welfare')
  n <- dim(X)[1]
  p <- dim(X)[2]
  XX <- cbind(1, X)
  g_i_S <- rep(0, n)
  g_i_S2 <- rep(0,n)
  if(is.na(g_i)[1]){
    ## Propensity for the sensitive attribute
    G_i1 <- (Y - m0) * (1-D) * S / ((1-propensity1)*propensity2) + m0 * S/propensity2
    G_i2 <- (Y - cost_treatment - m1) * D * S/ (propensity1*propensity2) + m1 * S/propensity2
    g_i_S = G_i2 - G_i1

    G_i1 <- (Y - m0) * (1-D) * (1 - S) / ((1-propensity1)*(1 - propensity2)) + m0 * (1 - S)/(1 - propensity2)
    G_i2 <- (Y - cost_treatment - m1) * D * (1 - S)/ (propensity1*(1 - propensity2)) + m1 * (1 - S)/(1 - propensity2)
    g_i_S2 = G_i2 - G_i1
    g_i = propensity2 * g_i_S + (1 - propensity2) * g_i_S2 ## estimated welfare
  }

  constant1 <- constant2 <- 0
  ## Different effects for distance being either welfare or envy
  if(distance == 'relative_welfare'){
    ## Unfairness is the sum of component1 and component2
    component1 <- g_i_S
    component2 <- -g_i_S2
  }

  if(distance == 'envy'){
    component1 <- mu_hat01*S/propensity2 -  mu_hat00*S/propensity2 - g_i_S ## policy 1
    component2 <- mu_hat11*(1 - S)/(1 - propensity2) -  mu_hat10*(1 - S)/(1 - propensity2) - g_i_S2 ## policy

    }

  if(distance == 'parity'){
    component1 <- S/propensity2
    component2 <- - (1 - S)/(1 - propensity2)
  }
  UnFairness_before <- component1 + component2
  unique_values <- unique(XX)
  index_unique <- apply(XX,1, function(y) which(apply(unique_values, 1, function(x) all(y == x))))
  g_i_before <- g_i
  g_i <- sapply(c(1:dim(unique_values)[1]), function(x) sum(g_i[which(x == index_unique)]))
  component1 <- sapply(c(1:dim(unique_values)[1]), function(x) sum(component1[which(x == index_unique)]))
  component2 <- sapply(c(1:dim(unique_values)[1]), function(x) sum(component2[which(x == index_unique)]))
  UnFairness <- component1 + component2

  maximum_X <- max(apply(XX, 1, function(h) max(abs(h))))
  XX <- XX/maximum_X

  XX <- unique_values
  n_original <- n
  n <- length(g_i)
  n_indexes <- sapply(c(1:n), function(x) sum(index_unique == x))
  C <- B*max(apply(XX, 1, function(x) sum(abs(x))))
  XX <- as.matrix(XX/C)


  model  <- list()
  AA <- rbind(cbind(diag(1, nrow = n), -XX), cbind(diag(1, nrow = n), -XX), c(n_indexes, rep(0, p +1)))
  sense <- c(rep('<=', dim(AA)[1]/2), rep('>', dim(AA)[1]/2), '<=')

  model$obj<- c(g_i, rep(0, p + 1))
  model$modelsense<-'max'
  rhs <- c(rep(1 - tolerance_constraint, dim(AA)[1]/2), rep(tolerance_constraint, dim(AA)[1]/2), max_treated_units)
  model$vtype<- c(rep('B', n), rep('C', p+1))
  # Put bounds on the parameter space (If commented, parameter space = real line)
  model$ub<- c(rep(1,n), rep(B,1+p))
  model$lb<- c(rep(0,n), rep(-B,p + 1))

  AA <-  rbind(AA, c(UnFairness, rep(0, p + 1)), c(UnFairness, rep(0, p + 1)))
  if(two_directions){
  rhs <- c(rhs, UnFairness_bound, - UnFairness_bound) } else{
    rhs <- c(rhs, Inf, - UnFairness_bound)
  }
  sense <- c(sense, '<=', '>=')
  model$rhs <- rhs
  model$sense <- sense
  model$A <- AA
  if(is.na(params)[1]){
    params<- list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9, TimeLimit = maxtime, BarConvTol = exp(-2),
                  Threads = cores, Disconnected=0, Heuristics = 0, NodefileStart = 0.5)}
  .require_gurobi()
  result<- .call_gurobi(model, params = params)
  beta_hat <- result$x[(n+1):(n+p+1)]
  pi_est <- apply(cbind(1,X), 1, function(x) ifelse(x%*%beta_hat > 0, 1, 0))
  return(list(obj_est = sum(g_i_before*pi_est), g_i = g_i, result = result$x, pi = pi_est, beta = beta_hat,
              UnFairness = sum(UnFairness_before * pi_est), group1 = mean(g_i_before * pi_est * S/propensity2),
              group2 = mean(g_i_before * pi_est * (1 - S)/(1 - propensity2))))
}


# ---- from tree_optimization.R ----

## Inputs: x: individual covariates
##         splits: points where to make the splits
##         var_vector: vector with corresponding covariates where to make the split
## return the label ('00', '01', etc) which indicates the position for a given split
## only support trees of length 1,2,3
helper_splits = function(x, splits, var_vector, k = 1){

  label = ifelse(x[unlist(var_vector[1])] < splits[1], 1, 0)

  if (length(var_vector) == 1) return(label)


  if(label == 1){
    label = paste0(label, ifelse(x[unlist(var_vector[2])] < splits[2], 1, 0))
  } else {
    label = paste0(label, ifelse(x[unlist(var_vector[3])] < splits[3], 1, 0))
  }

  if(length(var_vector) == 3) return(label)

  if(length(var_vector) == 7) {

    if(label == '11'){
      label = paste0(label, ifelse(x[unlist(var_vector[4])] < splits[4], 1, 0))
    }
    if(label == '10'){
      label = paste0(label, ifelse(x[unlist(var_vector[5])] < splits[5], 1, 0))
      }

    if(label == '01'){
      label = paste0(label, ifelse(x[unlist(var_vector[6])] < splits[6], 1, 0))
    }


    if(label == '00'){
      label = paste0(label, ifelse(x[unlist(var_vector[7])] < splits[7], 1, 0))
    }
    return(label)
  }

  if(length(var_vector) > 7){
    stop('Function only support trees of depth up to three.')
  }

}



estimate_simple_tree = function(Y, X, D, S, propensity1, propensity2, B=1, params = NA, model_only = F, tolerance_constraint = 10**(-7),
         scale_Y = F,
         cost_treatment  = 0, alpha = 1/2,
         max_treated_units, maxtime = 300, m1, m0,
         additional_fairness_constraint = F,
         parity_constraint = '>=', cores = 1) {


  my_results =  Est_max_score_probabilistic(Y, X, D, S, propensity1, propensity2, B=1, params, model_only = F, tolerance_constraint,
                                          scale_Y,
                                          cost_treatment, alpha, g_i = NA,
                                          max_treated_units, maxtime, m1, m0,
                                          additional_fairness_constraint,
                                          parity_constraint, cores, intercept = F)


  return(my_results)
}



## given a matrix of covariates, split points and variables where to split it construct
## a covariate matrix with dummies to assign observations into different classification buckets
## note: var_vector is formed such that (var1, var2, var3, ...) and splits contain the split across each variable
## var vector is first split is var1, second splits are var 2, var3 and so on

## the function return the matrix of covariates, the splits and var_vector passed to the function

generate_labels <- function(X, splits, var_vector){

  my_labels = apply(X, 1, function(x) helper_splits(x, splits, var_vector))
  if(length(var_vector) == 1) my_labels = factor(my_labels, levels = c('0', '1'))
  if(length(var_vector) == 3) my_labels = factor(my_labels, levels = c('00', '01', '10', '11'))
  if(length(var_vector) == 7) my_labels = factor(my_labels, levels = c('000', '001', '010', '011',
                                                                      '100', '111', '110', '101'))

  new_X = stats::model.matrix(~ my_labels - 1)
  return(list(new_X, splits, var_vector))

}

## Given a matrix of covariates X, number of splits considered and the vector of variables considered for the splits
## the function return a list of covariance matrix for classification (dummies for different classification buckets)
## return list of matrices

create_matrices_given_var <- function(X, num_splits, var_vector){



   n = dim(X)[1]
   seq_matrix = list()

   my_saved_list = list()

   acc = 1
   acc2 = 1
   for(j in 1:length(var_vector)){

     vv = unlist(var_vector[j])
     my_splits = unique(sort(X[, vv])[seq(from = n/num_splits, to = n - num_splits, length = num_splits  - 1)] + 0.000001)
     if(max(my_splits) > max(X[,vv])) my_splits = my_splits[-length(my_splits)]
     seq_matrix[[acc]] <- my_splits
     acc = acc + 1
   }



   my_grid = do.call(expand.grid, seq_matrix)
  # matrices = list()
  # for(j in 1:dim(my_grid)[1]){
  #   matrices[[j]] = generate_labels(X, my_grid[j,], var_vector)
  # }
   return(my_grid)


}

compute_predictions_tree <- function(my_tree_vector, X,
                                     length_tree = 2){
  num_variables = 2^length_tree - 1
  splits_variables = my_tree_vector[c(2:(num_variables + 1))]
  which_variables_to_split = my_tree_vector[c((num_variables + 2):(num_variables + 1 + num_variables))]
  coefficients = my_tree_vector[c((num_variables + 3 + num_variables):length(my_tree_vector))]

  matrix_covariates = generate_labels(X, splits_variables, which_variables_to_split)
  policy1 = apply(matrix_covariates[[1]], 1, function(x) x%*%unlist(coefficients))
  return(policy1)
}

## The function compute different measures of fairness for a classification tree
## my_tree_vector is a vector object returned by maximize_optimal_tree function
## no_parity_constraint indicates that there is no sensitive attribute S in the covariate matrix passed to the function
## two_directions indicates whether we take absolute values to compute fairness

compute_fairness_tree = function(my_tree_vector, Y, X, D, S, propensity1, propensity2, m1, m0, mu_hat11, mu_hat10, mu_hat01, mu_hat00,
                                 length_tree = 2, cost_treatment = 0, g_i = NA,
                                 no_parity_constraint = F, two_directions = T){

  num_variables = 2^length_tree - 1
  splits_variables = my_tree_vector[c(2:(num_variables + 1))]
  which_variables_to_split = my_tree_vector[c((num_variables + 2):(num_variables + 1 + num_variables))]
  coefficients = my_tree_vector[c((num_variables + 3 + num_variables):length(my_tree_vector))]
  if(no_parity_constraint){
    XX1 <- as.matrix(cbind(1, X[,-1]))
    XX0 <- as.matrix(cbind(0, X[,-1]))
  } else {
    XX1 <- X
    XX0 <- X
  }

  matrix_covariates1 = generate_labels(XX1, splits_variables, which_variables_to_split)
  matrix_covariates0 = generate_labels(XX0, splits_variables, which_variables_to_split)

  policy1 = apply(matrix_covariates1[[1]], 1, function(x) x%*%unlist(coefficients))
  policy0 = apply(matrix_covariates0[[1]], 1, function(x) x%*%unlist(coefficients))

  n <- dim(X)[1]
  p <- dim(X)[2]

  G_i1 <- (Y - m0) * (1-D) * S / ((1-propensity1)*propensity2) + m0 * S/propensity2
  G_i2 <- (Y - cost_treatment - m1) * D * S/ (propensity1*propensity2) + m1 * S/propensity2
  g_i_S = G_i2 - G_i1

  G_i12 <- (Y - m0) * (1-D) * (1 - S) / ((1-propensity1)*(1 - propensity2)) + m0 * (1 - S)/(1 - propensity2)
  G_i22 <- (Y - cost_treatment - m1) * D * (1 - S)/ (propensity1*(1 - propensity2)) + m1 * (1 - S)/(1 - propensity2)
  g_i_S2 = G_i22 - G_i12

  constant1 <- constant2 <- 0
  ## Different effects for distance being either welfare or envy

  welfare1 <- sum(g_i_S*policy1) + sum(G_i1)
  welfare0 <- sum(g_i_S2*policy0) + sum(G_i12)
  unfairness_welfare <-  welfare0 - welfare1
  if(two_directions) unfairness_welfare <- abs(welfare1 - welfare0)

  welfare1b <- sum(g_i_S*policy1)
  welfare0b <- sum(g_i_S2*policy0)
  unfairness_relative_welfare <-  welfare0b - welfare1b
  if(two_directions) unfairness_relative_welfare <- abs(welfare1b - welfare0b)

  w1 <- sum(S*policy1/propensity2)
  w0 <- sum((1 - S)*policy0/(1 - propensity2))
  unfairness_parity <- w0 - w1
  if(two_directions) unfairness_parity <- abs(w1 - w0)

  ## Compute the distance for the envy-based fairness

  welfare1 <- sum(g_i_S*policy1)
  welfare0 <- sum(g_i_S2*policy0)
  objective_warm_starts1 <- sum(mu_hat01*policy1*S/propensity2) +  sum(mu_hat00*(1 - policy1)*S/propensity2) -  welfare0
  objective_warm_starts2 <- sum(mu_hat11*policy0*(1 - S)/(1 - propensity2)) +  sum(mu_hat10*(1 - policy0)*(1 - S)/(1 - propensity2)) - welfare1
  unfairness_envy <- objective_warm_starts1 + objective_warm_starts2

  welfare1 <- sum(g_i_S*policy1) + sum(G_i1)
  welfare0 <- sum(g_i_S2*policy0)  + sum(G_i12)

  return(c(my_tree_vector, unfairness_parity = unfairness_parity, unfairness_envy = unfairness_envy,
              unfairness_relative_welfare = unfairness_relative_welfare, unfairness_welfare = unfairness_welfare,
           welfare1 = welfare1, welfare0 = welfare0))

}

## Find the welfare optimal tree on the frontier, for a given value of alpha
## ## additional fairness constraint, parity_constraint: add additional constraints
## cores: number of cores used for optimization
## num_splits: number of splits used
## length_tree: length of the tree
## epsilon_optimality: impose pareto optimality up to epsilon slackness

## return list with the best objective function and matrix with candidate trees in each row

maximize_welfare_tree <- function(Y, X, D, S, propensity1, propensity2, B=1, params = NA, model_only = F, tolerance_constraint = 10**(-7),
                                      scale_Y = F,
                                      cost_treatment  = 0, alpha = 1/2, g_i = NA,
                                      max_treated_units, maxtime = 300, m1 = 0, m0 = 0,
                                      additional_fairness_constraint = F,
                                      parity_constraint = '>=', cores = 1, num_splits = 4, length_tree = 2,
                                      epsilon_optimality = 0.001,
                                      constrained_policy_class = F,
                                      params_constrained_policy_class = list(UnFairness_bound = NA,
                                                                             distance = NA,
                                                                             two_directions = NA,
                                                                             mu_hat11 = NA,
                                                                             mu_hat10 = NA,
                                                                             mu_hat01 = NA,
                                                                             mu_hat00 = NA) ){
  if (!requireNamespace('foreach', quietly = TRUE)) .abort('The exhaustive tree optimizer requires the `foreach` package.')
  foreach <- foreach::foreach
  `%do%` <- foreach::`%do%`
  tmp = rep(list(seq_len(dim(X)[2])), 2^length_tree - 1)
  grid_1 = expand.grid(tmp)
  saved_res = c()
  best_object = -Inf
  list_of_trees = list()
  for(j in 1:dim(grid_1)[1]){


    var_vector = grid_1[j,]
    my_cov = create_matrices_given_var(X, num_splits, unlist(var_vector))
    acc = 1

    init_obj = -Inf


    my_results = foreach(k = 1:dim(my_cov)[1], .combine = rbind)%do%{

      aa = generate_labels(X, my_cov[k,], unlist(var_vector))
      if(constrained_policy_class == F){
      max_score_prob = tryCatch(Est_max_score_probabilistic(Y,
                                  as.matrix(aa[[1]]), D, S,
                                  propensity1, propensity2, B, params,
                                  model_only, tolerance_constraint,
                                  scale_Y,
                                  cost_treatment, alpha, g_i = NA,
                                  max_treated_units, maxtime, m1, m0,
                                  additional_fairness_constraint,
                                  parity_constraint, cores, intercept = F), error = function(e)
                                    list(obj_est = -Inf, beta = rep(NA, dim(as.matrix(aa[[1]]))[2] + 1 ))
                                )
      } else {

        UnFairness_bound = params_constrained_policy_class$UnFairness_bound
        distance = params_constrained_policy_class$distance
        two_directions = params_constrained_policy_class$two_directions
        mu_hat11 = params_constrained_policy_class$mu_hat11
        mu_hat01 = params_constrained_policy_class$mu_hat01
        mu_hat00 = params_constrained_policy_class$mu_hat00
        mu_hat10 = params_constrained_policy_class$mu_hat10

        max_score_prob = tryCatch(welfare_fairness_constraints_probabilistic(Y,
                                                                             X =  as.matrix(aa[[1]]), D, S,
                                                                             propensity1, propensity2, B, params,
                                                                             model_only, tolerance_constraint,
                                                                             scale_Y,
                                                                             cost_treatment, g_i,
                                                                             max_treated_units, maxtime, m1, m0,
                                                                             additional_fairness_constraint,
                                                                             mu_hat11 = mu_hat11, mu_hat01 = mu_hat01 ,
                                                                             mu_hat00 = mu_hat00, mu_hat10 = mu_hat10,
                                                                             alpha = alpha, UnFairness_bound = UnFairness_bound,
                                                                             distance = distance, two_directions = two_directions,
                                                                             intercept = F),
                                  error = function(e) list(obj_est = -Inf, beta = rep(NA, dim(as.matrix(aa[[1]]))[2] + 1 )))
      }



      # c(max_score_prob$obj_est, unlist(my_cov[[k]][2]), unlist(my_cov[[k]][3]), alpha, max_score_prob$beta)
      c(max_score_prob$obj_est,
        my_cov[k,], var_vector, alpha, max_score_prob$beta)
      }


    my_results = as.data.frame(my_results)
    best_results = my_results[which.max(my_results[,1]),]
    save_elements = which(as.numeric(my_results[,1]) > as.numeric(unlist(best_results[1])) - epsilon_optimality)
    if(j == 1){
      if(length(save_elements) == 1) {
        all_tree = matrix(my_results[save_elements,],nrow= 1)
        colnames(all_tree) = paste0('V', c(1:dim(all_tree)[2]))
      } else {
        all_tree = my_results[save_elements,]
      }
    } else {
      final_save = my_results[save_elements,]
      colnames(final_save) = colnames(all_tree)
      all_tree = rbind(all_tree, final_save)
    }
  }

  pruned_optimality = which(as.numeric(all_tree[,1]) > max(as.numeric(unlist(all_tree[,1]))) - epsilon_optimality)
  optimal_trees = all_tree[pruned_optimality,]
  best_object =  max(as.numeric(unlist(all_tree[,1])))

  return(list(best_object, optimal_trees))

}

### Estimate the pareto frontier for the tree
estimate_frontier_optimal_tree <- function(Y, X, D, S, propensity1, propensity2, B=1, params = NA, model_only = F,
                                           tolerance_constraint = 10**(-7),
                                           scale_Y = F,
                                           discretization = floor(sqrt(length(Y))),
                                           cost_treatment  = 0,
                                           alpha_seq = seq(from = 0, to = 1, length = discretization),
                                           g_i = NA,
                                           max_treated_units, maxtime = 300, m1 = 0, m0 = 0,
                                           additional_fairness_constraint = F,
                                           parity_constraint = '>=', cores = 12, num_splits = 4, length_tree = 2, parallel = FALSE,
                                           epsilon_optimality = 0.001) {
  if (!requireNamespace('foreach', quietly = TRUE)) .abort('The exhaustive tree optimizer requires the `foreach` package.')
  foreach <- foreach::foreach
  `%do%` <- foreach::`%do%`
  `%dopar%` <- foreach::`%dopar%`
  if(parallel){
  if (!requireNamespace('doParallel', quietly = TRUE)) .abort('Parallel tree optimization requires the `doParallel` package.')
  doParallel::registerDoParallel(cores)
  res <- foreach(i = alpha_seq, .combine = append, .export = c('Y', 'D', 'X', 'S', 'propensity1',
                                                              'propensity2', 'params', 'scale_Y',
                                                              'cost_treatment', 'max_treated_units',
                                                              'maxtime', 'm1', 'm0', 'additional_fairness_constraint',
                                                              'parity_constraint'))%dopar%{

                          aa =   maximize_welfare_tree(Y, X, D, S, propensity1, propensity2, B, params, model_only, tolerance_constraint,
                                                                                                      scale_Y,
                                                                                                      cost_treatment, alpha = i, g_i,
                                                                                                      max_treated_units, maxtime, m1, m0,
                                                                                                      additional_fairness_constraint,
                                                                                                      parity_constraint,
                                                       cores, num_splits, length_tree, epsilon_optimality)
                          list(aa)
                                                              }


  } else{
    res <- foreach(i = alpha_seq, .combine = append, .export = c('Y', 'D', 'X', 'S', 'propensity1',
                                                                'propensity2', 'params', 'scale_Y',
                                                                'cost_treatment', 'max_treated_units',
                                                                'maxtime', 'm1', 'm0', 'additional_fairness_constraint',
                                                                'parity_constraint'))%do%{

                                                        aa =          maximize_welfare_tree(Y, X, D, S, propensity1, propensity2, B, params, model_only, tolerance_constraint,
                                                                                       scale_Y,
                                                                                            cost_treatment, alpha = i, g_i,
                                                                                            max_treated_units, maxtime, m1, m0,
                                                                                            additional_fairness_constraint,
                                                                                            parity_constraint,
                                                                                            cores, num_splits, length_tree,
                                                                                       epsilon_optimality)
                                                    list(aa)

                                                          }
    }

  return(res)
}

## Construct matrix with all candidate optimal trees from the frontier

estimate_all_candidate_optimal_trees <- function(Y, X, D, S, propensity1, propensity2,  scale_Y = F,
                            discretization = floor(sqrt(length(Y))),
                            cost_treatment = 0, params=NA,
                            mu_hat11 = 0, mu_hat01 = 0, mu_hat00 = 0, mu_hat10 = 0,
                            max_treated_units, maxtime = 300, warm_start = NA, alpha_seq,
                            noparity_constraint = F,
                            additional_fairness_constraint = F,
                            m0 = 0, m1 = 0,
                            numcores = 12, two_directions = T, tolerance = 10**(-6),
                            num_splits = 4, length_tree = 2, epsilon_optimality = 0.001, parallel = FALSE, parity_constraint = '>='){



  frontier = estimate_frontier_optimal_tree(Y, X, D, S, propensity1, propensity2 =  propensity2, B=1, params = NA, model_only = F,
                                            tolerance_constraint = 10**(-7),
                                                        scale_Y,  discretization,
                                                        cost_treatment, alpha_seq, g_i = NA,
                                                        max_treated_units, maxtime, m1, m0,
                                                        additional_fairness_constraint,
                                                        parity_constraint, cores = numcores, num_splits, length_tree,
                                                        parallel = parallel, epsilon_optimality)


  ## Frontier contains all candidate trees


  all_res = frontier[[1]][[2]]
  ## remove duplicates
  if(is.null(dim(all_res))[1] == F) all_res = all_res[!duplicated(all_res[,1]),]


  if(is.null(dim(all_res))[1] == F){
    all_res = t(apply(all_res, 1, function(x) unlist(compute_fairness_tree(x, Y, X, D, S, propensity1, propensity2, m1, m0, mu_hat11, mu_hat10,
                                                                    mu_hat01, mu_hat00,
                                                        length_tree, cost_treatment, g_i,
                                                        noparity_constraint, two_directions))))
  } else{
    all_res = compute_fairness_tree(all_res, Y, X, D, S, propensity1, propensity2, m1, m0, mu_hat11, mu_hat10, mu_hat01, mu_hat00,
                                   length_tree, cost_treatment, g_i,
                                   noparity_constraint, two_directions)

  }

  for(k in frontier[-1]){
  ## remove duplicates
  without_duplicates = k[[2]]
  if(is.null(dim(k[[2]]))[1] == F) without_duplicates =  k[[2]][!duplicated(k[[2]][,1]),]

  if(is.null(dim(without_duplicates))[1]){
    without_duplicates = compute_fairness_tree(without_duplicates, Y, X, D, S, propensity1, propensity2, m1, m0, mu_hat11, mu_hat10,
                                               mu_hat01, mu_hat00,
                                    length_tree, cost_treatment, g_i,
                                    noparity_constraint, two_directions)
  } else{
    without_duplicates = t(apply(without_duplicates, 1, function(x) unlist(compute_fairness_tree(x, Y, X, D, S, propensity1, propensity2, m1, m0,
                                                                                          mu_hat11, mu_hat10, mu_hat01, mu_hat00,
                                               length_tree, cost_treatment, g_i,
                                               noparity_constraint, two_directions))))
  }
  all_res = rbind(all_res, without_duplicates)
  }

  ## return best policy and all the candidate policies from the frontier
  return(all_res)
}

## optimal the optimal tree in terms of fairness from the candidates on the estimated frontier

estimate_fairness_optimal_tree = function(Y, X, D, S, propensity1, propensity2,  scale_Y = F,
                                  discretization = floor(sqrt(length(Y))),
                                  cost_treatment = 0, params=NA,
                                  mu_hat11 = 0, mu_hat01 = 0, mu_hat00 = 0, mu_hat10 = 0,
                                  max_treated_units, maxtime = 300, alpha_seq,
                                  noparity_constraint = F,
                                  additional_fairness_constraint = F,
                                  m0 = 0, m1 = 0,
                                  numcores = 12, two_directions = T, tolerance = 10**(-6),
                                  num_splits = 6, length_tree = 2, epsilon_optimality = 0.001, parallel = FALSE, parity_constraint = '>=', warm_start = NA){


  results = estimate_all_candidate_optimal_trees(Y, X, D, S, propensity1, propensity2,  scale_Y,
                                 discretization,
                                 cost_treatment, params,
                                 mu_hat11, mu_hat01, mu_hat00, mu_hat10,
                                 max_treated_units, maxtime, warm_start, alpha_seq,
                                 noparity_constraint,
                                 additional_fairness_constraint,
                                 m0, m1,
                                 numcores, two_directions, tolerance,
                                 num_splits, length_tree, epsilon_optimality, parallel = parallel, parity_constraint = parity_constraint)
  results = as.data.frame(results)
  best_parity = results[which.min(results$unfairness_parity),]

  best_envy =  results[which.min(results$unfairness_envy),]

  best_relative_welfare = results[which.min(results$unfairness_relative_welfare),]

  best_welfare = results[which.min(results$unfairness_welfare),]

  rownames(results) = c(1:dim(results)[1])
  return(list(best_parity = best_parity, best_envy = best_envy, best_relative_welfare = best_relative_welfare,
             best_welfare = best_welfare, all_candidates = results))
}


