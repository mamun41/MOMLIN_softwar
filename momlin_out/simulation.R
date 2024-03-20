library(MASS) # mvrnorm
library(randcorr)
library(ggplot2)


generate_data <- function(n = 1000, n_cluster = 5, p = 100, 
                          trait = c('categ','conti')) {
  ## Generate clustered multi-view data
  set.seed(123)
  
  ## Correlation level
  pc = round(p * 0.1)
  r = runif(pc, 0.5, 0.8)
  cov = rbind(cbind(diag(pc), diag(r)), cbind(diag(r), diag(pc)))
  XY_c = mvrnorm(n = n, mu = c(rep(0,pc), rep(0,pc)), Sigma = cov, empirical = F)
  X_c = XY_c[,1:pc]
  Y_c = XY_c[,(pc+1):(2*pc)]
  
  ## Trait level (cluster)
  pt = round(p * 0.1)
  if(trait == 'conti'){
    rand = sample(1:n_cluster, size=n, replace=TRUE)
    mus <- seq(-1, 1, length.out = n_cluster)
    sds <- rep(1, n_cluster)
    trait = rnorm(n, mean=mus[rand], sd=sds[rand])
    r = runif(pt, 0.8, 1)
    cov = rbind(cbind(diag(pt), diag(r)), cbind(diag(r), diag(pt)))
    XY_t = mvrnorm(n = n, mu = rep(0,2*pt), Sigma = cov, empirical = F)
    X_t = rand * XY_t[,1:pt]
    Y_t = rand * XY_t[,pt+1:pt]
    clusters = rand
  }
  
  pct = pc + pt
  X_ct <- cbind(X_c, X_t)
  Y_ct <- cbind(Y_c, Y_t)
  
  # canonical weights
  A = matrix(runif(pct^2, min = -3, max = 3), ncol = pct)
  B = matrix(runif(pct^2, min = -3, max = 3), ncol = pct)
  
  # generate data features
  X_ct = X_ct %*% A
  Y_ct = Y_ct %*% B
  
  ## Individual level (noise)
  pxi = p - pct
  pyi = p - pct
  XY_i = matrix(rnorm(n*(pxi+pyi), mean = 0, sd = 1), nrow = n)
  X_i = XY_i[,1:pxi]
  Y_i = XY_i[,(pxi+1):(pxi+pyi)]
  
  X = cbind(X_ct, X_i)
  Y = cbind(Y_ct, Y_i)
  
  return(list(X = X, Y = Y, trait = trait, cluster = clusters))
}

#########################################
n_clust = 3
n_sample = 1000
n_feature = 2000 # sum of 2 views

data = generate_data(n=n_sample, n_cluster=n_clust, 
                     p=n_feature, trait='conti')
X = scale(data$X); Y = scale(data$Y); trait = data$trait; cluster = data$cluster
rownames(X) = rownames(Y) = names(trait) = paste0('s', 1:nrow(X))
colnames(X) = paste0('x', 1:ncol(X))
colnames(Y) = paste0('y', 1:ncol(Y))
# distribution of trait
hist(trait, breaks = 40)

