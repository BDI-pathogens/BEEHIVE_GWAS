##### DEFINE FUNCTIONS

solve.for.lik <- function(Sigma, VL, X){ # solve for likelihood given var-covar matrix Sigma, data VL, design matrix X (for fixed effects)
  # if n is the number of datapoints (length(VL)), Sigma must be of dimensions n*n and X have n rows
  n <- length(VL)
  stopifnot(nrow(Sigma) == n)
  stopifnot(nrow(X) == n)
  lik <- - 10^9
  beta <- NULL
  
  tryCatch( # catch error
    {
      tmp <- get.inv.det.sigma(Sigma = Sigma, n = n) # get inverse sigma and log determinant

      # if(all(Sigma[lower.tri(Sigma)] == 0, Sigma[upper.tri(Sigma)] == 0)){ # if matrix is diagonal:
      #   inverse.Sigma <- diag(n)
      #   diag(inverse.Sigma) <- 1/diag(Sigma)
      #   log.det.Sigma <- sum(log(diag(Sigma)))
      # } else { # non-diagonal matrix use slower functions
      #   inverse.Sigma <- chol2inv(chol(Sigma)) # inverse from cholesky decomposition is faster -> this is the slow calculation
      #   log.det.Sigma <- determinant(Sigma, logarithm=T)$modulus[1]
      # }
      tmp2 <- get.beta.lik.v2(inverse.Sigma = tmp$inverse.Sigma, log.det.Sigma = tmp$log.det.Sigma, X =  X, phenotype = VL, n = n)
      beta <- tmp2$beta
      least.square.mu <- tmp2$least.square.mu
      lik <- tmp2$lik
    }
    ,
    error = function(e) {
      print(e)
    }
  )
  return(list(beta = beta, loglik = -lik, log.det.Sigma = tmp$log.det.Sigma, inverse.Sigma = tmp$inverse.Sigma))
}
get.inv.det.sigma <- function(Sigma, n){ # this function gets the inverse matrix and log-determinant
  if(all(Sigma[lower.tri(Sigma)] == 0, Sigma[upper.tri(Sigma)] == 0)){ # if matrix is diagonal:
    inverse.Sigma <- diag(n)
    diag(inverse.Sigma) <- 1/diag(Sigma)
    log.det.Sigma <- sum(log(diag(Sigma)))
  } else { # non-diagonal matrix use slower functions
    inverse.Sigma <- fastMatrixInversion2(Sigma) # inverse from cholesky decomposition is faster -> this is the slow calculation
    log.det.Sigma <- determinant(Sigma, logarithm=T)$modulus[1]
  }
  return(list(inverse.Sigma = inverse.Sigma, log.det.Sigma = log.det.Sigma))
}

fastMatrixInversion1 <- function(mat){
  if(!is.matrix(mat)) stop('matrix to invert is not in the matrix format')
  return(chol2inv(chol(mat)))
}
fastMatrixInversion2 <- function(mat){
  if(!is.matrix(mat)) stop('matrix to invert is not in the matrix format')
  return(eigenMatInv(mat))
}
Rcpp::sourceCpp(code =
                  "
                // [[Rcpp::depends(RcppEigen)]]
                
                #include <RcppEigen.h>
                
                
                // [[Rcpp::export]]
                SEXP eigenMatInv(Eigen::MatrixXd A){
                Eigen::MatrixXd Ainv = A.inverse();
                
                return Rcpp::wrap(Ainv);
                }
                
                // [[Rcpp::export]]
                SEXP eigenSPDMatInv(Eigen::MatrixXd A, int N){
                Eigen::MatrixXd I = Eigen::MatrixXd::Identity(N, N); // define I as identity matrix of size N*N
                Eigen::MatrixXd Ainv = A.llt().solve(I);
                return Rcpp::wrap(Ainv);
                }
                
                "
)
get.beta.lik <- function(inverse.Sigma, log.det.Sigma, X, phenotype, n){ # given the var-covar matrix get best fixed effects
  beta <- solve(t(X) %*% inverse.Sigma %*% X)  %*% (t(X) %*% inverse.Sigma %*% phenotype) # gives you best fixed effects. e.g equation (2.26) in Christophe Lippert thesis
  least.square.mu <- X %*% beta # average at the tips according to least square regression
  lik <- -(1/2) * (n * log(2*pi) + log.det.Sigma + t(phenotype - X %*% beta) %*% inverse.Sigma %*% (phenotype - X %*% beta)) # see e.g. equation (2.16) in Christoph Lippert thesis
  return(list(beta = beta, least.square.mu = least.square.mu, lik = as.numeric(lik)))
}
get.beta.lik.v2 <- function(inverse.Sigma, log.det.Sigma, X, phenotype, n){ # same as get.beta.lik but faster
  
  prod1 <- crossprod(X, inverse.Sigma) # same as t(X) %*% inverse.Sigma but faster
  prod2 <- prod1 %*% X
  term1 <- solve(prod2)
  term2 <- prod1 %*% phenotype
  beta <- term1 %*% term2
  
  # and to compute the likelihood:
  prod3 <- X %*% beta
  term3 <- phenotype - prod3
  prod4 <- crossprod(term3, inverse.Sigma) # same as t(term3) %*% inverse.Sigma but faster
  term4 <- prod4 %*% term3
  #beta <- solve(t(X) %*% inverse.Sigma %*% X)  %*% (t(X) %*% inverse.Sigma %*% phenotype) # gives you best fixed effects. e.g equation (2.26) in Christophe Lippert thesis
  
  least.square.mu <- prod3 # average at the tips according to least square regression
  lik <- -(1/2) * (n * log(2*pi) + log.det.Sigma + term4) # see e.g. equation (2.16) in Christoph Lippert thesis
  return(list(beta = beta, least.square.mu = least.square.mu, lik = as.numeric(lik)))
}

optim.fun.repeated <- function(n.repeats, lik.fun, init.fun, verbose = F, cluster = NULL, lower = NULL, upper = NULL, ...){
  # a function to MINIMISE properly the negative log-likelihood function
  # using repeated use of function optim
  generate.init.fun <- match.fun(init.fun)
  lik.fun <- match.fun(lik.fun)
  all.opt <- list()
  
  opt_fun <- function(ii){
    init <- generate.init.fun()
    
    if(length(init)==1){ # one parameter function, use optimise
      if(is.null(upper) | is.null(lower)){ # if unspecified bounds for parameters, use (0,1)
        interval <- c(0, 1)
      } else {
        interval <- c(lower[1], upper[1])
      }
      tmp_opt <- optimise(lik.fun, interval = interval, ...)
      names(tmp_opt) <- c("par", "value")
      all.opt[[ii]]$message <- "Successful convergence"
    } else { # otherwise use nmk
      if(is.null(upper) | is.null(lower)){ # if unspecified bounds for parameters, use nmk
        tmp_opt <- nmk(init, lik.fun, ...)
      } else {
        tmp_opt <- nmkb(init, lik.fun, lower = lower, upper = upper, ...)
      }
    }  
    return(tmp_opt)
  }
  
  if(!is.null(cluster)) {
    all.opt <- parLapply(cluster, 1:n.repeats, opt_fun)
  } else {
    all.opt <- lapply(1:n.repeats, opt_fun)
  }
  
  all.opt <- matrix(unlist(lapply(all.opt, function(x) y <- c(x$par, x$value, x$feval, x$restarts, x$convergence))), nrow = n.repeats, byrow = T)
  colnames(all.opt) <-c(paste0("par", 1:length(generate.init.fun())), "value", "feval", "restarts", "convergence")
  return(data.frame(all.opt))
  # if(verbose) cat(sum(unlist(lapply(all.opt, function(x)x$message == "Successful convergence"))), "optimisations converge\n")
  # all.ll <- unlist(lapply(all.opt, function(x)x$value))
  # min(all.ll) -> minll.pow
  # if(verbose) cat(length(which(all.ll < minll.pow + 1 & all.ll > minll.pow)), " optim within 1 log lik of minimum likelihood\n")
  # output <- list(minll.pow, all.opt[[which.min(all.ll)]]$par)
  # names(output) <- c("lik", "pars")
  #return(output)
}

# two negative likelihood functions

# mynegativeloglik0 <- function(par, data, n){
# # just a function essentially constructing the design matrix X and variance-covariance matrix Sigma, from the vector of parameters par, then feeding that to the likelihood
#   sigma_e <- par[1]
#   sigma_w <- par[2]
#   Sigma <- sigma_e * diag(n) + sigma_w * diag(n) * 1 / data$N.SPVL
#   X <- as.matrix(rep(1, n), ncol = 1) # only an intercept for now
#   lik <- solve.for.lik(Sigma, data$VL, X)
#   return(- lik$loglik)
# }
#mynegativeloglik1 <- function(par, data = list(VL = suba$SPVL, N.SPVL = suba$N.SPVL, Kmat = K), n = nrow(suba)){
# mynegativeloglik1 <- function(par, data, X, n){
#   # just a function essentially constructing the variance-covariance matrix Sigma, from the vector of parameters par, then feeding that to the likelihood
#   sigma_e <- par[1]
#   sigma_w <- par[2]
#   sigma_g <- par[3]
#   stopifnot(nrow(X)==n)
#   Sigma <- sigma_e * diag(n) + sigma_w * diag(n) * 1 / data$N.SPVL + sigma_g * data$Kmat
#   lik <- solve.for.lik(Sigma, data$phenotype, X)
#   return(- lik$loglik)
# }
mynegativeloglik1 <- function(par, data, X, n, Xrand, npar, return_fixed = F){
  # a function essentially constructing the variance-covariance matrix Sigma, from the vector of parameters par, then feeding that to the likelihood
  par <- unlist(par)
  stopifnot(length(par) == npar)
  stopifnot(nrow(X) == n)
  stopifnot(length(Xrand) == npar) # one variance parameter for each element of the list
  Sigma <- generate_sigma1(par = par, n = n, Xrand = Xrand, npar = npar) # construct the random effect matrix as a sum of terms
  lik <- solve.for.lik(Sigma, data$phenotype, X)
  if(return_fixed) return(lik) else { return( lik$loglik) }
}
generate_sigma1 <- function(par, n, Xrand, npar, return_full = F){
  stopifnot(length(Xrand) == npar)
  stopifnot(length(par) == npar)
  # build SigmaG:
  SigmaG <- par[1] * Xrand[[1]]
  
  # build SigmaE:
  SigmaE <- matrix(0, ncol = n, nrow = n)
  
  for(i in 2:npar) SigmaE <- SigmaE + par[i] * Xrand[[i]]
  
  if(return_full) return(list(SigmaG = SigmaG, SigmaE = SigmaE)) else return(SigmaG + SigmaE)
}
mynegativeloglik2 <- function(par, data, X, n, Xrand, npar, return_fixed = F){
  # a function essentially constructing the variance-covariance matrix Sigma, from the vector of parameters par, then feeding that to the likelihood
  par <- unlist(par)
  stopifnot(nrow(X)==n)
  stopifnot(npar == 4 | npar == 5)
  stopifnot(length(Xrand) + 2 == npar)
  
  threshold <- par[npar] # last parameter is threshold for low/high error variance
  below_threshold_error <- data$phenotype < threshold
  stopifnot(length(below_threshold_error) == n)
  
  Sigma <- generate_sigma2(par, n, Xrand, below_threshold_error, npar) # construct the random effect matrix as a sum of terms
  lik <- solve.for.lik(Sigma, data$phenotype, X)
  if(return_fixed) return(lik) else { return( lik$loglik) }
}
generate_sigma2 <- function(par, n, Xrand, below_threshold_error, npar, return_full = F){
  
  # build SigmaG:
  SigmaG <- par[1] * Xrand[[1]]
  
  # build SigmaE:
  SigmaE <- matrix(0, ncol = n, nrow = n)
  
  if(npar == 5){
    SigmaE <- SigmaE + par[2] * Xrand[[2]] # add the identity matrix random effect which is in position 2, only when there are 5 parameters (i.e. we have values for the error variance)
  }
  vec_par <- rep(par[npar - 2], n) # the variance when phenotype is above threshold
  vec_par[below_threshold_error] <- par[npar - 1] # the (higher) variance when phenotype is below threshold
  SigmaE <- SigmaE + vec_par * Xrand[[npar - 2]] # add the error variance with that structure
  if(return_full) return(list(SigmaG = SigmaG, SigmaE = SigmaE)) else return(SigmaG + SigmaE)
}

clean_lik <- function(all_opt, npar){
  opt <- list(lik = min(all_opt$value), pars = all_opt[which.min(all_opt$value), paste0("par", 1:npar)])
  opt$pars <- unname(unlist(opt$pars))
  opt
}

get_h2 <- function(par, data, X, n, Xrand, npar, lik_fun_name, nsim = 500, return_full = F){ # computation of h2 from e and g
  
  lik_fun <- match.fun(lik_fun_name)
  
  # estimation of heritability based on simulations of ML model
  full_lik <- lik_fun(par = par, data = data, X = X, n = n, Xrand = Xrand, npar = npar, return_fixed = T)
  fixed_effects <- full_lik$beta
  loglik <- full_lik$loglik
  
  mu_e <- X %*% full_lik$beta # in that case X includes only environmental variables
  mu_g <- rep(0, n)           # genetic mean is 0
  
  if(lik_fun_name == "mynegativeloglik2"){ # when we have a phenotype threshold for the error
    threshold <- par[npar] # last parameter is threshold for low/high error variance
    below_threshold_error <- suba[, pheno_name] < threshold
    tmp <- generate_sigma2(par, n, Xrand, below_threshold_error = below_threshold_error, npar = npar, return_full = T)
  }
  if(lik_fun_name == "mynegativeloglik1"){ # when we don't
    tmp <- generate_sigma1(par, n, Xrand, npar = npar, return_full = T)
  }
  
  sigma_g <- tmp$SigmaG
  sigma_e <- tmp$SigmaE
  
  e <- MASS::mvrnorm(n = nsim, mu = as.matrix(mu_e), Sigma = as.matrix(sigma_e))
  g <- MASS::mvrnorm(n = nsim, mu = as.matrix(mu_g), Sigma = as.matrix(sigma_g))
  
  vare <- rowMeans(e^2) - rowMeans(e)^2
  varg <- rowMeans(g^2) - rowMeans(g)^2
  h2 <- mean(varg / (vare + varg)) # 16% for SPVL; 16% for normalised_adjusted_spvl; 14% for GSVL
  if(return_full) return(list(h2 = h2, e = e, g = g, loglik = loglik, fixed_effects = fixed_effects)) else return(h2)
}

get_h2_v2 <- function(par, data, X, n, Xrand, npar, lik_fun_name, nsim = 500, return_full = F){ # computation of h2 from e and z (useful when sigma_g is not positive definite)
  
  lik_fun <- match.fun(lik_fun_name)
  
  # estimation of heritability based on simulations of ML model
  full_lik <- lik_fun(par = par, data = data, X = X, n = n, Xrand = Xrand, npar = npar, return_fixed = T)
  fixed_effects <- full_lik$beta
  
  mu_e <- X %*% full_lik$beta # in that case X includes only environmental variables
  mu_g <- rep(0, n)           # genetic mean is 0
  
  if(lik_fun_name == "mynegativeloglik2"){
    threshold <- par[npar] # last parameter is threshold for low/high error variance
    below_threshold_error <- suba[, pheno_name] < threshold
    tmp <- generate_sigma2(par, n, Xrand, below_threshold_error = below_threshold_error, npar = npar, return_full = T)
  }
  if(lik_fun_name == "mynegativeloglik1"){
    tmp <- generate_sigma1(par, n, Xrand, npar = npar, return_full = T)
  }
  
  sigma_g <- tmp$SigmaG
  sigma_e <- tmp$SigmaE
  
  e <- MASS::mvrnorm(n = nsim, mu = mu_e, Sigma = sigma_e)
  z <- MASS::mvrnorm(n = nsim, mu = mu_g + mu_e, Sigma = sigma_g + sigma_e)
  
  vare <- rowMeans(e^2) - rowMeans(e)^2
  varz <- rowMeans(z^2) - rowMeans(z)^2
  h2 <- mean(1 - vare / varz) # 16% for SPVL; 16% for normalised_adjusted_spvl; 14% for GSVL
  if(return_full) return(list(h2 = h2, e = e, g = g)) else return(h2)
}

##### LOAD DATA AND OPTIMISE LIKELIHOOD (without genetic data here)
if(FALSE){
  a <- read.csv("~/DID/BEEHIVE_Hackathon/Data/SummaryData/BEEHIVE_summary05.06.2017.csv")
  to.keep <- !is.na(a$SPVL) & !is.na(a$N.SPVL)
  suba <- a[to.keep, ]
  
  # example function:
  mynegativeloglik0(c(0.5, 0.1))
  
  # optimising the function: (needs installing dfoptim package)
  require(dfoptim)
  opt <- optim.fun.repeated(n.repeats = 1, lik.fun = mynegativeloglik0, init.fun = function() runif(2, 0, 1), lower = c(0,0), upper = c(1,1))
  
  x <- seq(0.4, 0.6, 0.01)
  y <- sapply(x, function(ve) mynegativeloglik0(c(ve, 0.3568)))
  plot(x, y, type = "o", ylab = "log lik", xlab = "environmental variance")
  
  x <- seq(0.25, 0.45, 0.01)
  y <- sapply(x, function(vw) mynegativeloglik0(c(0.4728, vw)))
  plot(x, y, type = "o", ylab = "log lik", xlab = "within-host / error variance")
}

####                                    FUNCTION NMKB FROM DFOPTIM PACKAGE              ####

nmkb <- function (par, fn, lower = -Inf, upper = Inf, control = list(), ...) 
{
  ctrl <- list(tol = 1e-06, maxfeval = min(5000, max(1500, 
                                                     20 * length(par)^2)), regsimp = TRUE, maximize = FALSE, 
               restarts.max = 3, trace = FALSE)
  namc <- match.arg(names(control), choices = names(ctrl), 
                    several.ok = TRUE)
  if (!all(namc %in% names(ctrl))) 
    stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
  if (!is.null(names(control))) 
    ctrl[namc] <- control
  ftol <- ctrl$tol
  maxfeval <- ctrl$maxfeval
  regsimp <- ctrl$regsimp
  restarts.max <- ctrl$restarts.max
  maximize <- ctrl$maximize
  trace <- ctrl$trace
  n <- length(par)
  g <- function(x) {
    gx <- x
    gx[c1] <- atanh(2 * (x[c1] - lower[c1])/(upper[c1] - 
                                               lower[c1]) - 1)
    gx[c3] <- log(x[c3] - lower[c3])
    gx[c4] <- log(upper[c4] - x[c4])
    gx
  }
  ginv <- function(x) {
    gix <- x
    gix[c1] <- lower[c1] + (upper[c1] - lower[c1])/2 * (1 + 
                                                          tanh(x[c1]))
    gix[c3] <- lower[c3] + exp(x[c3])
    gix[c4] <- upper[c4] - exp(x[c4])
    gix
  }
  if (length(lower) == 1) 
    lower <- rep(lower, n)
  if (length(upper) == 1) 
    upper <- rep(upper, n)
  if (any(c(par < lower, upper < par))) 
    stop("Infeasible starting values!", call. = FALSE)
  low.finite <- is.finite(lower)
  upp.finite <- is.finite(upper)
  c1 <- low.finite & upp.finite
  c2 <- !(low.finite | upp.finite)
  c3 <- !(c1 | c2) & low.finite
  c4 <- !(c1 | c2) & upp.finite
  if (all(c2)) 
    stop("Use `nmk()' for unconstrained optimization!", call. = FALSE)
  if (maximize) 
    fnmb <- function(par) -fn(ginv(par), ...)
  else fnmb <- function(par) fn(ginv(par), ...)
  x0 <- g(par)
  if (n == 1) 
    stop(call. = FALSE, "Use `optimize' for univariate optimization")
  if (n > 30) 
    warning("Nelder-Mead should not be used for high-dimensional optimization")
  V <- cbind(rep(0, n), diag(n))
  f <- rep(0, n + 1)
  f[1] <- fnmb(x0)
  V[, 1] <- x0
  scale <- max(1, sqrt(sum(x0^2)))
  if (regsimp) {
    alpha <- scale/(n * sqrt(2)) * c(sqrt(n + 1) + n - 1, 
                                     sqrt(n + 1) - 1)
    V[, -1] <- (x0 + alpha[2])
    diag(V[, -1]) <- x0[1:n] + alpha[1]
    for (j in 2:ncol(V)) f[j] <- fnmb(V[, j])
  }
  else {
    V[, -1] <- x0 + scale * V[, -1]
    for (j in 2:ncol(V)) f[j] <- fnmb(V[, j])
  }
  f[is.nan(f)] <- Inf
  nf <- n + 1
  ord <- order(f)
  f <- f[ord]
  V <- V[, ord]
  rho <- 1
  gamma <- 0.5
  chi <- 2
  sigma <- 0.5
  conv <- 1
  oshrink <- 0
  restarts <- 0
  orth <- 0
  dist <- f[n + 1] - f[1]
  v <- V[, -1] - V[, 1]
  delf <- f[-1] - f[1]
  diam <- sqrt(colSums(v^2))
  sgrad <- c(crossprod(t(v), delf))
  alpha <- 1e-04 * max(diam)/sqrt(sum(sgrad^2))
  simplex.size <- sum(abs(V[, -1] - V[, 1]))/max(1, sum(abs(V[, 
                                                              1])))
  itc <- 0
  conv <- 0
  message <- "Succesful convergence"
  while (nf < maxfeval & restarts < restarts.max & dist > ftol & 
         simplex.size > 1e-06) {
    fbc <- mean(f)
    happy <- 0
    itc <- itc + 1
    xbar <- rowMeans(V[, 1:n])
    xr <- (1 + rho) * xbar - rho * V[, n + 1]
    fr <- fnmb(xr)
    nf <- nf + 1
    if (is.nan(fr)) 
      fr <- Inf
    if (fr >= f[1] & fr < f[n]) {
      happy <- 1
      xnew <- xr
      fnew <- fr
    }
    else if (fr < f[1]) {
      xe <- (1 + rho * chi) * xbar - rho * chi * V[, n + 
                                                     1]
      fe <- fnmb(xe)
      if (is.nan(fe)) 
        fe <- Inf
      nf <- nf + 1
      if (fe < fr) {
        xnew <- xe
        fnew <- fe
        happy <- 1
      }
      else {
        xnew <- xr
        fnew <- fr
        happy <- 1
      }
    }
    else if (fr >= f[n] & fr < f[n + 1]) {
      xc <- (1 + rho * gamma) * xbar - rho * gamma * V[, 
                                                       n + 1]
      fc <- fnmb(xc)
      if (is.nan(fc)) 
        fc <- Inf
      nf <- nf + 1
      if (fc <= fr) {
        xnew <- xc
        fnew <- fc
        happy <- 1
      }
    }
    else if (fr >= f[n + 1]) {
      xc <- (1 - gamma) * xbar + gamma * V[, n + 1]
      fc <- fnmb(xc)
      if (is.nan(fc)) 
        fc <- Inf
      nf <- nf + 1
      if (fc < f[n + 1]) {
        xnew <- xc
        fnew <- fc
        happy <- 1
      }
    }
    if (happy == 1 & oshrink == 1) {
      fbt <- mean(c(f[1:n], fnew))
      delfb <- fbt - fbc
      armtst <- alpha * sum(sgrad^2)
      if (delfb > -armtst/n) {
        if (trace) 
          cat("Trouble - restarting: \n")
        restarts <- restarts + 1
        orth <- 1
        diams <- min(diam)
        sx <- sign(0.5 * sign(sgrad))
        happy <- 0
        V[, -1] <- V[, 1]
        diag(V[, -1]) <- diag(V[, -1]) - diams * sx[1:n]
      }
    }
    if (happy == 1) {
      V[, n + 1] <- xnew
      f[n + 1] <- fnew
      ord <- order(f)
      V <- V[, ord]
      f <- f[ord]
    }
    else if (happy == 0 & restarts < restarts.max) {
      if (orth == 0) 
        orth <- 1
      V[, -1] <- V[, 1] - sigma * (V[, -1] - V[, 1])
      for (j in 2:ncol(V)) f[j] <- fnmb(V[, j])
      nf <- nf + n
      ord <- order(f)
      V <- V[, ord]
      f <- f[ord]
    }
    v <- V[, -1] - V[, 1]
    delf <- f[-1] - f[1]
    diam <- sqrt(colSums(v^2))
    simplex.size <- sum(abs(v))/max(1, sum(abs(V[, 1])))
    f[is.nan(f)] <- Inf
    dist <- f[n + 1] - f[1]
    sgrad <- c(crossprod(t(v), delf))
    if (trace & !(itc%%2)) 
      cat("iter: ", itc, "\n", "value: ", f[1], "\n")
  }
  if (dist <= ftol | simplex.size <= 1e-06) {
    conv <- 0
    message <- "Successful convergence"
  }
  else if (nf >= maxfeval) {
    conv <- 1
    message <- "Maximum number of fevals exceeded"
  }
  else if (restarts >= restarts.max) {
    conv <- 2
    message <- "Stagnation in Nelder-Mead"
  }
  return(list(par = ginv(V[, 1]), value = f[1] * (-1)^maximize, 
              feval = nf, restarts = restarts, convergence = conv, 
              message = message))
}

##### TEST OF THE MATRIX INVERSION ALGORITHMS #####
# t1 <- 0
# t2 <- 0
# for(i in 1:100){
#   dim <- 1000
#   Sigma <- matrix(runif(n = dim*dim, min = 0, max  = 1), ncol = dim); Sigma <- t(Sigma) %*% Sigma
#   stopifnot(isSymmetric.matrix(Sigma))
#   
#   t0 <- Sys.time()
#   inverse.Sigma1 <- fastMatrixInversion1(as.matrix(Sigma)) # inverse from cholesky decomposition is faster -> this is the slow calculation
#   t1 <- t1 + Sys.time() - t0
#   
#   t0 <- Sys.time()
#   inverse.Sigma2 <- fastMatrixInversion2(as.matrix(Sigma)) # TODO WRAPPER TO TEST SPD AND DIMENSION
#   t2 <- t2 + Sys.time() - t0
#   
#   
#   print(i)
#   print(summary(c(inverse.Sigma1) - c(inverse.Sigma2)))
#   
# }
# t1
# t2

