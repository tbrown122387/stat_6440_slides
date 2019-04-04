
# data
y <- c(-0.59836613, -0.10194784, -0.06399318,  0.41957954, 
       -0.98565345, -0.09350700,  1.75391375,  0.79859181,  
       0.32095243, 0.72231544)
nu <- 5

#############
# Example 1 #
#############

# function that evaluates log posteriors
eval_log_unnormalized_posterior <- function(mu, ss){
  zsquareds <- (y - mu)*(y - mu)/ss
  lds <- -.5*(nu+1) * log( 1 + zsquareds/nu )
  sum(lds) - log(ss)
}

# graph parameters
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

plotSurface(-50, 50, 0.0001, 50, 20, eval_log_unnormalized_posterior, F,
            theta=-120, zlab = "log unnorm dens", xlab = "mu", ylab = "ss")

#############
# Example 2 #
#############

v <- seq(.001,.5,.01)
ss <- .01
yi_minus_mu <- .1
plot(v, dinvgamma(v, shape=(nu+1)/2, .5*(nu*ss + (yi_minus_mu)^2)), 
     type = "l", main = "p(v_1 | mu, ss=.01, y)", ylab = "density")

