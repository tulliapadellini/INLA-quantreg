source("crossing-functions.R")

#install.packages("remotes")
remotes::install_github("grayclhn/dbframe-R-library")
library(dbframe)

scenario3 = c(1,1,1,0,0,0)

S3.pois.U.50 = RepParallel(1000, cross.sim(50, scenario3, alpha.seq = alpha.seq, intercept = T, covariates.gen = "Unif", response.gen = "poisson"))
apply(unlist(S3.pois.U.50 != 0), 1, mean)

S3.pois.TN.50 = RepParallel(1000, cross.sim(50, scenario3, alpha.seq = alpha.seq, intercept = T, covariates.gen = "TN", response.gen = "poisson"))
apply(unlist(S3.pois.TN.50 != 0), 1, mean)


S3.pois.U.100  = RepParallel(1000, cross.sim(100, scenario3, alpha.seq = alpha.seq, intercept = T, covariates.gen = "Unif", response.gen = "poisson"))
apply(unlist(S3.pois.U.100 != 0), 1, mean)

S3.pois.TN.100 = RepParallel(1000, cross.sim(100, scenario3, alpha.seq = alpha.seq, intercept = T, covariates.gen = "TN", response.gen = "poisson"))
apply(unlist(S3.pois.TN.100 != 0), 1, mean)


S3.pois.U.1000 = RepParallel(1000, cross.sim(1000, scenario3, alpha.seq = alpha.seq, intercept = T, covariates.gen = "Unif", response.gen = "poisson"))
apply(unlist(S3.pois.U.1000 != 0), 1, mean)

S3.pois.TN.1000 = RepParallel(1000, cross.sim(1000, scenario3, alpha.seq = alpha.seq, intercept = T, covariates.gen = "TN", response.gen = "poisson"))
apply(unlist(S3.pois.TN.1000 != 0), 1, mean)


S3.pois.U.1000 = RepParallel(1000, cross.sim(1000, scenario3, alpha.seq = alpha.seq, intercept = T, covariates.gen = "Unif", response.gen = "poisson"))
apply(unlist(S3.pois.U.1000 != 0), 1, mean)

S3.pois.TN.1000 = RepParallel(1000, cross.sim(1000, scenario3, alpha.seq = alpha.seq, intercept = T, covariates.gen = "TN", response.gen = "poisson"))
apply(unlist(S3.pois.TN.1000 != 0), 1, mean)





S3.gauss.U.50  = RepParallel(1000, cross.sim(50, scenario3, alpha.seq = alpha.seq, intercept = T, covariates.gen = "Unif", response.gen = "gaussian"))
apply(unlist(S3.gauss.U.50 != 0), 1, mean)

S3.gauss.TN.50 = RepParallel(1000, cross.sim(50, scenario3, alpha.seq = alpha.seq, intercept = T, covariates.gen = "TN", response.gen = "gaussian"))
apply(unlist(S3.gauss.TN.50 != 0), 1, mean)


S3.gauss.U.100  = RepParallel(1000, cross.sim(100, scenario3, alpha.seq = alpha.seq, intercept = T, covariates.gen = "Unif", response.gen = "gaussian"))
apply(unlist(S3.gauss.U.100 != 0), 1, mean)

S3.gauss.TN.100 = RepParallel(1000, cross.sim(100, scenario3, alpha.seq = alpha.seq, intercept = T, covariates.gen = "TN", response.gen = "gaussian"))
apply(unlist(S3.gauss.TN.100 != 0), 1, mean)


S3.gauss.U.1000  = RepParallel(1000, cross.sim(1000, scenario3, alpha.seq = alpha.seq, intercept = T, covariates.gen = "Unif", response.gen = "gaussian"))
apply(unlist(S3.gauss.U.1000 != 0), 1, mean)

S3.gauss.TN.1000 = RepParallel(1000, cross.sim(1000, scenario3, alpha.seq = alpha.seq, intercept = T, covariates.gen = "TN", response.gen = "gaussian"))
apply(unlist(S3.gauss.TN.1000 != 0), 1, mean)







