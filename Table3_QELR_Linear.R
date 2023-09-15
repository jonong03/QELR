#Table3: Quadratic Exponential Logistic Regression (Linear)

library(dplyr)
library(geepack)



# Create Functions --------------------------------------------------------

# Function: Random generation for QEBD
# n: Number of sample data to be generated
# V: A "mxm"square and diagonal matrix, with main effects placed at diagonal and interaction effects placed at off-diagonal position
RandomQEBD <- function(n, V){
  m = ncol(V)
  config0 = expand.grid(rep(list(c(0,1)),m)) %>% as.matrix
  
  pmf <- function(U, y) exp(sum((U%*%y)*y))
  pr <- sapply(1:nrow(config0), function(i) pmf(U= V, y= config0[i,]))
  cdf <- cumsum(pr)/sum(pr)
  Y <- sapply(1:n, function(x) config0[sum(cdf < runif(1)) + 1,])
  return(t(Y))
}

# Function: Generate simulation data for QELR (Linear Correlation Effect)
# beta: true parameter values for main effects. Length >= 1
# gamma: true parameter value for interaction effects. Length must be 1 since this is the simulation case for a common interaction effect. 
# n: sample size
# p: number of subjects 
GenData <- function(beta, gamma, n, p){
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
  
  U<- lapply(1:n, function(i) lapply(1:ngamma, function(j) w_f()))
  
  for (i in 1:n){
    L <- lapply(1:ngamma, function(j) U[[i]][[j]] * gamma[j]/2) %>% Reduce(`+`,.)
    diag(L) <-  X[(i*p-p+1): (i*p),] %*% beta
    Y[i,]<- RandomQEBD(1, L)
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

# Function: Estimate parameters for QELR Data
# DATA: `DAT` output from datagen().
EstQELR <- function(DATA){
  aa0 <- glm(y~.-1-id, family=binomial, data=DATA) %>% summary %>% coef
  col_name <- colnames(aa0)
  row_name <- rownames(aa0)
  rownames(aa0) <- paste0("glm_", row_name)
  
  aa1 <- geeglm(y~.-1-id, family=binomial, id=id, 
                data=DATA, corstr="independence") %>% summary %>% coef
  colnames(aa1) <- col_name
  rownames(aa1) <- paste0("geeind_", row_name)
  
  aa2 <- geeglm(y~.-1-id, family=binomial, id=id, 
                data=DATA, corstr="exchangeable") %>% summary %>% coef
  colnames(aa2) <- col_name
  rownames(aa2) <- paste0("geeexc_", row_name)
  
  out<- as.matrix(rbind(aa0, aa1, aa2))
  return(out)
  
}

# Run One Example ---------------------------------------------------------
# dat <- GenData(beta= c(-2.4, -2.0, -2.6), gamma= c(-1.2, -0.5), n=300, p= 15)
# EstQELR(dat$DAT) %>% round(.,3)

# Simulation --------------------------------------------------------------

SimDriver<- function(iter, seedn, beta, gamma, p, n_seq){
  set.seed(seedn)
  true_par <- c(beta, gamma)
  np<- length(true_par)
  
  OUT <- array(0, dim=c(iter, length(n_seq), length(n_seq)*np, 4))
  for (i in 1:iter){
    dat <- GenData(beta, gamma, n=max(n_seq), p)
    for (n in 1:length(n_seq)){
      OUT[i,n,,] <-  EstQELR(dat$DAT[1:(n_seq[n]*p),])
    }
    cat("iter ",i, " is done at", date(),"\n")
    gc()
  }
  
  G <- OUT[1:iter,,,]
  se <- t(apply(G[,,1:np,1], c(2,3), sd))
  glm <- t(apply(G[,,1:np,1], c(2,3), mean))
  glm_b <- glm - true_par
  glm_se <- t(apply(G[,,1:np,2], c(2,3), mean))
  glm_re <- glm_se/se
  geeind <- t(apply(G[,,1:np+np,1], c(2,3), mean))
  geeind_b <- geeind - true_par
  geeind_se <- t(apply(G[,,1:np+np,2], c(2,3), mean))
  geeind_re <- geeind_se/se
  geeexc <- t(apply(G[,,1:np+(np*2),1], c(2,3), mean))
  geeexc_b <- geeexc - true_par
  geeexc_se <- t(apply(G[,,1:np+(np*2),2], c(2,3), mean))
  geeexc_re <- geeexc_se/se
  K <- cbind(c(se), c(glm), c(glm_b), c(glm_se), c(glm_re), c(geeind), c(geeind_b), c(geeind_se), c(geeind_re), c(geeexc), c(geeexc_b), c(geeexc_se), c(geeexc_re))
  colnames(K)<- c("SE", "GLM","GLM Bias", "GLM SE","GLM R.E.", "GEEIND","GEEIND Bias","GEEIND SE","GEEIND R.E.", "GEEEXC","GEEEXC Bias","GEEEXC SE","GEEEXC R.E.")

  rowname<- expand.grid(b= c(paste0("beta",1:length(beta)-1),paste0("gamma",1:length(gamma))),n= paste0("n=",n_seq))
  K<- cbind(rowname, Truth= rep(true_par, by= length(n_seq)), K)

return(list(OUT= OUT, K= K))

}

# Parameter Setup
iter <- 500
seedn <- 88
beta <- c(-2.4, -2.0, -2.6)
gamma <- c(-1.2, -0.5)
n_seq <- c(100,300,500)
p <- 15

sim.table3<- SimDriver(iter, seedn, beta, gamma, p, n_seq)
sim.table3$K
