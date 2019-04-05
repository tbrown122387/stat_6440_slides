# p(y \mid mu, ss) = Normal(mu, ss)
# p(mu, ss) \propto 1/ss

# we know the posterior in this case, so we can compare
# p(\mu | ss, y) = Normal(mean(y), ss/n)
# p(ss | y) = Inv-Gamma((n-1)/2, (n-1)s^2/2)

# generate fake data
real_mu <- 5
real_ss <- 2
n <- 10
y <- rnorm(n, mean = real_mu, sd = sqrt(real_ss))
ssquared <- var(y)
ybar <- mean(y)

# find posterior mode with BFGS (a Quasi-Newton approach)
# Recall quasi-Newton only uses first derivatives of the log-likelihood
neg_log_unnorm_post <- function(params){
  mu <- params[1]
  theta <- params[2]
  chunk <- (n-1)*ssquared + n*(ybar - mu)^2
  n*theta + .5*exp(-2*theta + log(chunk)) 
}

neg_gradient <- function(params){
  mu <- params[1]
  theta <- params[2]
  c(
    -exp(-2*theta)*n*(ybar - mu),
    n - exp(-2*theta)*( (n-1)*ssquared + n*(ybar - mu)^2 )
    )
}

optim_results <- optim(par = c(5, 0), 
                       fn = neg_log_unnorm_post, 
                       gr = neg_gradient, # if left blank, finite difference approx. used
                       method = "BFGS", 
                       hessian = T)
if(optim_results$convergence==0){
  cat('successfully converged\n')
  cat('mode: mu = ', optim_results$par[1], ', and ss = ', exp(2*optim_results$par[2]), '\n')
}

# normal approximation to p(mu, theta | y
approx_mean <- optim_results$par
approx_cov_mat <- solve(optim_results$hessian)
approx_corr_mat <- cov2cor(approx_cov_mat)

unnorm_dens <- function(mu, log_sigma){
  params <- c(mu, log_sigma)
  lup <- exp(-neg_log_unnorm_post(params))
}
library(mvtnorm)
norm_approx_dens <- function(mu, log_sigma){
  dmvnorm(x = c(mu, log_sigma), mean = approx_mean, sigma = approx_cov_mat)
}
plotSurface <- function(lowerFirst, upperFirst, 
                        lowerSecond, upperSecond, 
                        numGridPointsOnEachAxis, f, contour = F, ...){
  A <- seq(lowerFirst, upperFirst, length.out = numGridPointsOnEachAxis)
  B <- seq(lowerSecond, upperSecond, length.out = numGridPointsOnEachAxis)
  args <- expand.grid(A,B)
  z <- mapply(f, args[,1], args[,2])
  dim(z) <- c(length(A), length(B))
  if(contour){
    contour(A, B, z)
  }else{
    persp(x=A, y=B, z=z, ...)
  }
}
# real but unnormalized
plotSurface(3,6,0,1.25,numGridPointsOnEachAxis = 50, f = unnorm_dens, contour = T)
title("Unnormalized true p(mu, theta | y)", font = 4)
plotSurface(3,6,0,1.25,numGridPointsOnEachAxis = 50, f = norm_approx_dens, contour = T)
title("Normal approx. p(mu, theta | y)", font = 4)


# compare posteriors
poss_mus <- seq(-2*abs(mu), 2*abs(mu), length.out = 100)
real_posterior <- dnorm(x = poss_mus, mean = mean(y), sd = 1/sqrt(length(y)))
bfgs_approx_mean <- optim_results$par
bfgs_approx_var <- 