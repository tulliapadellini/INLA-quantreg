
library(Qtools)
# vs
library(INLA)


# Assess crossing ---------------------------------------------------------

# This function checks whether or not there is at least one crossing in the domain x \in [0,1]

check.for.crossings <- function(beta1, beta2, intercept = T, tol= 1e-7, u.bound = 1, l.bound=0){
  
  A = rbind(beta1, beta2)
  
  if(intercept) A = A[,-1, drop = F]  
  
  
  ff = function(x) {
    if(intercept) x = c(1,x)
    abs(x %*% beta1 - x%*% beta2)
  }
  
  EE = optim(rep(u.bound/2, ncol(A)), fn = ff, lower = l.bound, upper =u.bound, method = "L-BFGS-B")
  return(EE$value < tol)
}


# take estimated curves and counts crossing among all quantile levels
count.crossings = function(betas, alpha.seq, intercept = TRUE, u.bound = 1, l.bound=0){
  
  pp = length(alpha.seq)
  error.mat = matrix(NA, pp,pp)
  for(x in 1:pp){
    for(y in x:pp){
      error.mat[x,y] = (check.for.crossings(beta1 = betas[x,], beta2 = betas[y,], intercept = intercept, u.bound = u.bound, l.bound = l.bound)  )
    }
  }  
  
  return(sum(error.mat[upper.tri(error.mat)]))
  
}


# Data generation + model estimation --------------------------------------

# each of these function returns the beta estimated for jittered and model based approach, using different ways to simulate x and y

gen.crossing = function(n, b=1, m = 10, intercept = TRUE, alpha.seq, covariates.gen, response.gen = "poisson"){
  
  p = length(b)
  
  x = switch(covariates.gen,
             "TN"         = matrix(truncnorm::rtruncnorm(n*p, a = 0, b = 1), nrow = n),
             "NT"         = matrix(truncnorm::rtruncnorm(n*p, a = 0, b = 1, mean = 1), nrow = n),
             "TN.outl"    = { 
               components <- sample(1:2,prob=c(0.95,0.05), size=n*p, replace=TRUE)
               means <- c(0,1)
               sds <- c(.5,.5)
               
               x <- truncnorm::rtruncnorm(n=n*p, a = 0, b = 1, mean = means[components], sd = sds[components])
               x <- matrix(x, nrow = n)
             },
             "Unif"      = matrix(runif(n*p, min = 0, max = 1), nrow = n),
             "Unif.outl" = {
               components <- sample(1:2,prob=c(0.95,0.05), size=n*p, replace=TRUE)
               l.bounds <- c(0,.85)
               u.bounds <- c(.85,1)
               
               x <- runif(n*p, min = l.bounds[components], max = u.bounds[components])
               x <- matrix(x, nrow = n)
             } 
  )
  
  
  colnames(x) <- paste0("x", 1:p)
  if(intercept) x[,1] <- 1
  eta = x%*%b
  lambda = exp(eta)
  
  
  y = numeric(n)
  
  y = switch(response.gen,
             "poisson" = rpois(n, lambda = lambda),
             "gaussian"= floor(truncnorm::rtruncnorm(n, a = 0, mean = lambda, sd = lambda))
  )
  
  
  formula = as.formula(paste("y~ -1 +", paste(colnames(x), collapse = "+")))
  
  df = data.frame(y = y, x)
  
  
  res.list = list()
  beta.q = matrix(NA, nrow = length(alpha.seq), ncol = p)
  
  for(alpha.idx in 1:length(alpha.seq)){
    
    
    temp = inla(formula, data = df, family = "poisson", 
                control.family = list(control.link = 
                                        list(model = "quantile", quantile = alpha.seq[alpha.idx])))
    res.list[[alpha.idx]] = temp
    beta.q[alpha.idx,] = temp$summary.fixed[,1]
    
    
  }
  
  
  machado.list = list()
  beta.j = matrix(NA, nrow = length(alpha.seq), ncol = p)
  for(alpha in 1:length(alpha.seq)){
    temp = Qtools::rq.counts(formula, data = df, tau = alpha.seq[alpha], M = m)
    machado.list[[alpha]] = temp
    beta.j[alpha,] = temp$coefficients
    
  }
  
  return(list(beta.mb = beta.q, beta.jitt = beta.j))
}






# Crossing count ----------------------------------------------------------

# add together -- generates the data + fit the model + counts the crossings for the 4 different models

cross.sim = function(n, beta.sim, alpha.seq, intercept=T, covariates.gen, response.gen="poisson"){
  
  x.sim = gen.crossing(n, b = beta.sim, intercept = intercept, alpha.seq = alpha.seq, covariates.gen = covariates.gen, response.gen = response.gen)
  x.mb  = count.crossings(x.sim$beta.mb, alpha.seq, intercept = intercept)
  x.jitt= count.crossings(x.sim$beta.jitt, alpha.seq, intercept = intercept)
  
  return(list(mb.cross = x.mb, jitt.cross = x.jitt))
  
}
