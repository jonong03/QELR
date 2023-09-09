#Table3: Quadratic Exponential Logistic Regression (Linear)
gc(); rm(list=ls())
library(dplyr)
library(data.table)
library(geepack)

beta <- c(-2.4, -2.0, -2.6)
gamma <- c(-1.2, -0.5)
n_seq <- c(100,300,500)
nmax <- max(n_seq)
p <- 15

rising <- function(n, V){
  a = c(0,1)
  m = dim(V)[1]
  config0 = expand.grid(rep(list(a),m)) %>% as.matrix
  
  pmf <- function(y) exp(sum((V%*%y)*y))
  pr <- sapply(1:nrow(config0), function(i) pmf(config0[i,]))
  cdf <- cumsum(pr)/sum(pr)
  Y <- sapply(1:n, function(x) config0[sum(cdf < runif(1)) + 1,])
  return(t(Y))
  
}

datagen <- function(beta, gamma, n, p){
  ngamma <- length(gamma)
  nbeta <- length(beta)
  X<- cbind(1, sapply(1:(nbeta-1), function(i) rnorm(n*p))) #
  Y<- matrix(0,n,p)
  
  w_f<- function(){
    seq <- sample(1:3, p, replace= T)
    V <- matrix(0,p,p)
    V[(seq %*% t(seq)) %in% seq^2] <- 1
    diag(V)<- 0
    return(V)
  }
  
  #U<- lapply(1:n, function(j) w_f())
  U<- lapply(1:n, function(i) lapply(1:ngamma, function(j) w_f()))
  
  for (i in 1:n){
    L <- lapply(1:ngamma, function(j) U[[i]][[j]] * gamma[j]/2) %>% Reduce(`+`,.)
    diag(L) <-  X[(i*p-p+1): (i*p),] %*% beta
    Y[i,]<- rising(1, L)
  }
  # Design matrix
  y <- c(t(Y))
  vec <- matrix(0, nrow=n*p, ncol= ngamma)
  for (k in 1:ngamma){
    vec[,k] <- sapply(1:n, function(j) U[[j]][[k]] %*% Y[j,]) ###
  }
  W <- cbind(X, vec)
  id <- rep(1:n, each=p)
  
  DAT <- data.frame(y=y, w=W, id=id) 
  
  return(list(Y=Y, X=X, DAT= DAT))  
  
}

QELR <- function(DATA){
  aa0 <- glm(y~.-1-id, family=binomial, data=DATA) %>% summary %>% coef
  col_name <- colnames(aa0)
  row_name <- rownames(aa0)
  rownames(aa0) <- paste0("glm_", row_name)
  aa1 <- geeglm(y~.-1-id, family=binomial, id=id, 
                data=DATA, corstr="independence") %>% summary %>% coef
  colnames(aa1) <- col_name
  rownames(aa1) <- paste0("geee_", row_name)
  
  aa2 <- geeglm(y~.-1-id, family=binomial, id=id, 
                data=DATA, corstr="exchangeable") %>% summary %>% coef
  colnames(aa2) <- col_name
  rownames(aa2) <- paste0("geecs_", row_name)
  
  out<- as.matrix(rbind(aa0, aa1, aa2))
  return(out)
  
}

#One Run
dat <- datagen(beta, gamma, n=max(n_seq), p)
QELR(dat$DAT)

#Simulation
simdriver<- function(iter, seedn, beta, gamma, p, n_seq){
  set.seed(seedn)
  true_par <- c(beta, gamma)
  np<- length(true_par)
  
  OUT <- array(0, dim=c(iter, length(n_seq), length(n_seq)*np, 4))
  for (i in 1:iter){
    dat <- datagen(beta, gamma, n=max(n_seq), p)
    for (n in 1:length(n_seq)){
      OUT[i,n,,] <-  QELR(dat$DAT[1:(n_seq[n]*p),])
    }
    cat("iter ",i, " is done at", date(),"\n")
    gc()
  }
  
  G <- OUT[1:iter,,,]
  se <- t(apply(G[,,1:np,1], c(2,3), sd))
  glm_b <- t(apply(G[,,1:np,1], c(2,3), mean)) - true_par
  glm_re <- t(apply(G[,,1:np,2], c(2,3), mean))/se
  geei_b <- t(apply(G[,,1:np+np,1], c(2,3), mean)) - true_par
  geei_re <- t(apply(G[,,1:np+np,2], c(2,3), mean))/se
  geee_b <- t(apply(G[,,1:np+(np*2),1], c(2,3), mean)) - true_par
  geee_re <- t(apply(G[,,1:np+(np*2),2], c(2,3), mean))/se
  K <- cbind(c(se), c(glm_b), c(glm_re), c(geei_b), c(geei_re), c(geee_b), c(geee_re))
  
  colnames(K)<- c("SE", "GLM_Bias", "GLM_RE", "GEEIND_Bias","GEEIND_RE", "GEECS_Bias", "GEECS_Re")
  rowname<- expand.grid(b= paste0("b",1:length(beta)-1),paste0("gamma",1:length(gamma)),n= paste0("n=",n_seq))
  K<- cbind(rowname, Truth= rep(true_par, by= length(n_seq)), K)

return(list(OUT= OUT, K= K))

}
sim2<- simdriver(iter=3, seedn=88, beta, gamma, p, n_seq)
sim2$K


