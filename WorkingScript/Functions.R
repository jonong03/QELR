# Functions only 

library(dplyr)
library(geepack)


# Create Functions --------------------------------------------------------

# Function 1: Density function
# v: a quantile vector with length m
# param: a set of parameters arranging in a `by row` order: theta11, theta12, ..., theta1k, theta22, theta23,..., theta2k, theta33, ...., thetakk.
# take_log: if TRUE, then returns a log-likelihood
DensityQEBD <- function(v, param, take_log=FALSE){
  m = dim(v)[2]
  m2 = m*(m+1)/2
  if(length(param) != m2) stop("wrong parameter dimension")
  
  config0 = expand.grid(rep(list(c(0,1)),m)) %>% as.matrix
  P <- matrix(0, m, m)
  P[lower.tri(P, diag=TRUE)] <- param
  
  pmf_inner <- function(y) exp(t(y)%*%P%*%y)
  den <- sapply(1:nrow(config0), function(i) pmf_inner(config0[i,])) %>% sum
  
  v <- as.matrix(v)
  if(dim(v)[2] == 1) v <- t(v)
  if(take_log){
    num <- (v%*%P%*%t(v)) %>% diag   #How to incorporate self-interaction term here?
    return(num - log(den))
  }else{
    num <- exp(v%*%P%*%t(v)) %>% diag
    return(num/den)
  }
}

# Function 2: Random generation for QEBD
# n: Number of sample data to be generated
# V: A "mxm"square and diagonal matrix, with main effects placed at diagonal and interaction effects placed at off-diagonal positions.
RandomQEBD <- function(n, V){
  m = ncol(V)
  config0 = expand.grid(rep(list(c(0,1)),m)) %>% as.matrix
  
  pmf <- function(U, y) exp(sum((U%*%y)*y))
  pr <- sapply(1:nrow(config0), function(i) pmf(U= V, y= config0[i,]))
  cdf <- cumsum(pr)/sum(pr)
  Y <- sapply(1:n, function(x) config0[sum(cdf < runif(1)) + 1,])
  return(t(Y))
}


# Function 3-1: Generate simulation data for QELR (Common Interaction Effect)
# beta: true parameter values for main effects. Length >= 1
# gamma: true parameter value for interaction effects. Length must be 1 since this is the simulation case for a common interaction effect. 
# n: sample size
# p: number of subjects 
GenData <- function(beta, gamma, n, p){
  nbeta <- length(beta)
  X<- cbind(1, sapply(1:(nbeta-1), function(i) rnorm(n*p))) #
  P<- lapply(1:n, function(i) matrix(gamma, p, p))
  Y<- matrix(0,n,p)
  
  for (i in 1:n){
    diag(P[[i]])<-  X[(i*p-p+1): (i*p),] %*% beta
    Y[i,]<- RandomQEBD(1, P[[i]])
  }
  
  # Design matrix
  y <- c(t(Y))
  U <- matrix(2 , p, p); diag(U) <- 0 ###
  vec <- sapply(1:n, function(j) U%*%Y[j,]) ###
  W <- cbind(X, c(vec))
  id <- rep(1:n, each=p)
  
  DAT <- data.frame(y=y, w=W, id=id) 
  
  return(list(Y=Y, X=X, DAT= DAT))  
}


# Function 3-2: Generate simulation data for QELR (Linear Correlation Effect)
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

# Function 3: Parameter estimation 
# Y: a "n x m" qebd data matrix.
# method: "mle" for maximum likelihood estimation, "logit" for logistic regression, "geeind" for GEE with independence working correlation.
EstQEBD <- function(Y, method){
  Y = as.matrix(Y)
  n = nrow(Y); m = ncol(Y)
  m2 = m*(m+1)/2
  
  Q = matrix(NA, m, m)
  Q[upper.tri(Q, diag=TRUE)] = 0
  loc<- which(Q==0, arr.ind=TRUE,useNames=TRUE) %>% as.data.frame %>% arrange(row) 
  var_name <- paste0(loc$row, ":",loc$col)
  
  # Design Matrix for estimating functions
  designX<- function(Y, varname){
    n = nrow(Y); m = ncol(Y)
    m2 = m*(m+1)/2
    
    X <- matrix(NA, nrow= n*m, ncol= m*m)
    I <- diag(m)
    for (i in 1:n){
      for (j in 1:m){
        y<- Y[i,]
        y[j]<- 1
        X[(j+m*(i-1)),] <- kronecker(y,I[,j])
      }
    }
    
    G<- matrix(0, nrow= m2, ncol= m*m)
    count<- 0
    for (i in 1:m){
      for (j in i:m){
        count <- count + 1
        if( i == j ) { 
          G[count,] <- kronecker(I[,i], I[,j])
        } else {
          G[count,] <- (kronecker(I[,i],I[,j]) + kronecker(I[,j], I[,i]))
        }
      }
    }
    
    X <- X %*% t(G)
    colnames(X)<- var_name
    
    return(X)
    
  }
  
  if(method=="mle"){
    
    loglik <- function(param) -sum(DensityQEBD(Y, param, take_log=TRUE))
    
    start= rep(0,m2) 
    est <- nlm(loglik, p= start, hessian=TRUE, gradtol = 1e-5)
    vcov_m <- est$hessian %>% solve
    se <- vcov_m %>% diag %>% sqrt
    
    output <- cbind(Estimate=est$estimate, se, zv= est$estimate/se, pv=pnorm(-abs(est$estimate/se)))
<<<<<<< HEAD
=======
    
>>>>>>> c4e24735a4bd4a98f8fbaa81bbe3fc59ceb63ab3
  }
  
  if(method=="geeind"){
    
    X<- designX(Y, var_name)
    ID= rep(1:n, each = m)
    
    gee <- geeglm(c(t(Y))~ 0+X, id=ID, family= binomial, corstr= "independence")
    est <- coef(gee)
    V <- vcov(gee)
    se <- diag(V) %>% sqrt
    
<<<<<<< HEAD
    o2<- QIC(gee)
    output<- cbind(est, se, est/se, pnorm(-abs(est/se)))

=======
    output<- cbind(est, se, est/se, pnorm(-abs(est/se)))
    
>>>>>>> c4e24735a4bd4a98f8fbaa81bbe3fc59ceb63ab3
  }
  
  if(method=="logit"){
    
    X<- designX(Y, var_name)
    
    fit <- glm(c(t(Y))~ 0+X, family= binomial)
    est <- coef(fit)
    V <- vcov(fit)
    se <- diag(V) %>% sqrt
    
    output<- cbind(est, se, est/se, pnorm(-abs(est/se)))
<<<<<<< HEAD

=======
    
>>>>>>> c4e24735a4bd4a98f8fbaa81bbe3fc59ceb63ab3
  }  
  
  # Reformat Output
  output<- as.data.frame(output)
  colnames(output)<- c("Estimate", "Std. Error", "z value", "Pr(|z|)")
  rownames(output)<- var_name
  
<<<<<<< HEAD
  if(method=="geeind") output <- list(output= output, qic = o2)
  
=======
>>>>>>> c4e24735a4bd4a98f8fbaa81bbe3fc59ceb63ab3
  return(output)
  
}


# Function 4: Estimate parameters for QELR Data
EstQELR <- function(DATA){
<<<<<<< HEAD
  #aa0 <- glm(y~.-1-id, family=binomial, data=DATA) %>% summary %>% coef
  #col_name <- colnames(aa0)
  #row_name <- rownames(aa0)
=======
  aa0 <- glm(y~.-1-id, family=binomial, data=DATA) %>% summary %>% coef
  col_name <- colnames(aa0)
  row_name <- rownames(aa0)
>>>>>>> c4e24735a4bd4a98f8fbaa81bbe3fc59ceb63ab3
  rownames(aa0) <- paste0("glm_", row_name)
  aa1 <- geeglm(y~.-1-id, family=binomial, id=id, 
                data=DATA, corstr="independence") %>% summary %>% coef
  colnames(aa1) <- col_name
  rownames(aa1) <- paste0("geeind_", row_name)
  
<<<<<<< HEAD
  #aa2 <- geeglm(y~.-1-id, family=binomial, id=id, 
  #              data=DATA, corstr="exchangeable") %>% summary %>% coef
  #colnames(aa2) <- col_name
  #rownames(aa2) <- paste0("geeexc_", row_name)
=======
  aa2 <- geeglm(y~.-1-id, family=binomial, id=id, 
                data=DATA, corstr="exchangeable") %>% summary %>% coef
  colnames(aa2) <- col_name
  rownames(aa2) <- paste0("geeexc_", row_name)
>>>>>>> c4e24735a4bd4a98f8fbaa81bbe3fc59ceb63ab3
  
  out<- as.matrix(rbind(aa0, aa1, aa2))
  return(out)
}
