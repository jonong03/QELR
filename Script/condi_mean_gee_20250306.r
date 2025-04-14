library(dplyr)
library(gee)
library(geepack)

ordi_GEE <- function(Y, X, W, a, maxite=20, tol=10^-6){
  n <- nrow(Y); m <- ncol(Y); pn <- dim(X)[2]; qn <- dim(W)[3]; d <- pn + qn
  
  ## prepare D matrix
  D <- array(0, dim=c(n, pn+qn, m))
  for(i in 1:n){
    for(j in 1:m){
      y_new <- Y[i,]; y_new[j] <- a
      D[i,,j] <- c(X[i,,j], W[i,j,,]%*%y_new)
    }
  }
  
  Z <- matrix(0, n*m, pn+qn)
  for(i in 1:n) Z[1:m+(i-1)*m, ] <- t(D[i,,])
  y <- c(t(Y))
  id <- rep(1:n, each=m)
  #out1 <- gee(y~-1+Z, id, family=binomial, corstr="independence", scale.fix=TRUE, scale.value=1)
  out2 <- geeglm(y~-1+Z, id = id, family = binomial, corstr = "independence")
  
  #return(list(gee=summary(out1)$coefficients,
  #            geepack=summary(out2)$coefficients))
  return(summary(out2)$coefficients[,1:2])
}

CondiMean_GEE <- function(Y, X, W, a, maxite=20, tol=10^-6){
  # Y: n X m
  # X: n X pn X m
  # W: n X m X qn X m
  n <- nrow(Y); m <- ncol(Y); pn <- dim(X)[2]; qn <- dim(W)[3]; d <- pn + qn
  
  ## prepare D matrix
  D <- array(0, dim=c(n, pn+qn, m))
  for(i in 1:n){
    for(j in 1:m){
      y_new <- Y[i,]; y_new[j] <- a
      D[i,,j] <- c(X[i,,j], W[i,j,,]%*%y_new)
    }
  }
  
  one_update <- function(theta){
    U <- rep(0, d); M <- matrix(0, d, d); B2 <- M; B3 <- M 
    for(i in 1:n){
      mu <- c(exp(t(D[i,,])%*%theta)); mu <- mu/(1+mu) 
      A <- diag(mu*(1-mu))
      Vinv <- solve(A)
    
      m_diff <- c(Y[i,] - mu)
      U <- U + D[i,,]%*%m_diff
      M <- M + D[i,,]%*%m_diff%*%t(m_diff)%*%t(D[i,,])
      B2 <- B2 + D[i,,]%*%A%*%t(D[i,,])

      Delta <- diag(m_diff*(1-2*mu))
      B3 <- B3 + D[i,,]%*%Delta%*%t(D[i,,])
    }
    return(list(U=U, M=M, B1=(B2-B3), B2=B2))
  }
  
  ## prepare D matrix
  Z <- matrix(0, n*m, pn+qn)
  for(i in 1:n) Z[1:m+(i-1)*m, ] <- t(D[i,,])
  y <- c(t(Y))
  theta <- glm(y~-1+Z, family=binomial) %>% coefficients
  for(i in 1:maxite){
    temp <- one_update(theta)
    theta_new <- theta + solve(temp$B1, temp$U)
    test <- max(abs(theta-theta_new))
    theta <- theta_new
    cat(i, test, theta, "\n")
    if(test < tol) break
  }
  
  temp <- one_update(theta)
  svar1 <- solve(temp$B1)%*%temp$M%*%solve(temp$B1) # w.o. taking expectation on the bum matrix
  svar2 <- solve(temp$B2)%*%temp$M%*%solve(temp$B2) # with taking expectation on the bum matrix
  
  return(list(theta=theta, V1=svar1, V2=svar2))
}



DGP <- function(n, m, pn, qn, beta, alpha){
  X <- array(rnorm(n*pn*m), dim=c(n, pn, m)); X[,1,] <- 1
  W <- array(rnorm(n*m*qn*m), dim=c(n, m, qn, m)) 
  Y <- matrix(0, n, m)
  Im <- diag(m)
  
  for(i in 1:n){
    mu <- t(X[i,,])%*%beta
    p <- mu[1]; p <- exp(p)/(1+exp(p))
    Y[i,1] <- rbinom(1, 1, p)
    W[i, 1, 1,] <- rep(0, m)
    for(j in 2:m){
      p <- mu[j] + alpha[1]*Y[i,j-1]
      p <- exp(p)/(1+exp(p))
      Y[i,j] <- rbinom(1, 1, p)
      W[i, j, ,] <- Im[j-1,]
    }
  }
  return(list(Y=Y, X=X, W=W))
}

n <- 200 # sample size
m <- 5 # repeated measures per subject
pn <- 3 # number of main effects
qn <- 1 # number of interaction effects

B <- 200
OUT <- array(0, dim=c(B, 4, 5))
for(i in 1:B){
  beta <- c(-1,0,2); alpha <- c(-1)
  data <- DGP(n, m, pn, qn, beta, alpha)
  Y <- data$Y; X <- data$X; W <- data$W

  out <- CondiMean_GEE(Y, X, W, 0, maxite=30)
  K <- cbind(out$theta, sqrt(diag(out$V1)), sqrt(diag(out$V2))) 

  outo <- ordi_GEE(Y, X, W, 0)
  OUT[i,,] <- as.matrix(cbind(K, outo))
  cat("Iteration ", i, "is done. \n")
}
OO <- cbind(c(beta,alpha), apply(OUT, c(2,3), mean), apply(OUT[,,1], 2, sd))
colnames(OO) <- c("true","gee_est", "gee_se_1", "gee_se_2", "pack_est", "pack_se", "emp_se")
OO

