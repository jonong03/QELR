# Table 1: Quadratic Exponential Distribution

rm(list=ls()); gc()
library(dplyr)
library(data.table)
library(geepack)

# Function 1: Density function
dising <- function(vv, m, theta, take_log=FALSE){
  if(length(theta) != m + (m*(m-1)/2)){
    cat("wrong dimension of parameter theta\n")
    return()
  }
  
  sign = 1 #Default as 1
  if(sign==1){a= c(0,1)}
  if(sign==2){a= c(-1,1)}
  
  config0 = expand.grid(rep(list(a),m)) %>% as.matrix
  P <- matrix(0, m, m)
  P[lower.tri(P, diag=TRUE)] <- theta
  #P[lower.tri(P)]<-  P[lower.tri(P)]/2
  #P[upper.tri(P)]<- P[lower.tri(P)]
  
  #pmf_inner <- function(y) exp(t(y)%*%P%*%y/2)
  pmf_inner <- function(y) exp(t(y)%*%P%*%y)
  den <- sapply(1:nrow(config0), function(i) pmf_inner(config0[i,])) %>% sum
  
  vv <- as.matrix(vv)
  if(dim(vv)[2] == 1) vv <- t(vv)
  if(take_log){
    num <- (vv%*%P%*%t(vv)) %>% diag   #How to incorporate self-interaction term here?
    return(num - log(den))
  }else{
    num <- exp(vv%*%P%*%t(vv)) %>% diag
    return(num/den)
  }
}

# Function 2: Random number generator
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

# Function 3: Max likelihood estimation for thetas. Methods include: newton, logit, gd (gradient descent)
MLE_ising <- function(Y, method, start= rep(0,ncol(Y)*(ncol(Y)+1)/2) ){
  options(dplyr.summarise.inform = FALSE)
  Y = as.matrix(Y)
  n = dim(Y)[1]
  m = dim(Y)[2]
  m2 = m*(m+1)/2
  
  colnames(Y)= paste0("x",1:m)
  vars = colnames(Y)
  
  Q = matrix(NA, m, m)
  Q[upper.tri(Q, diag=TRUE)] = 0
  loc<- which(Q==0, arr.ind=TRUE,useNames=TRUE) %>% as.data.frame %>% arrange(row) 
  t_list <- paste0(loc$row, ":",loc$col)
  
  if(method=="newton"){
    
    loglik <- function(param) -sum(dising(Y, m, param, take_log=TRUE))
    
    est <- nlm(loglik, p= start, hessian=TRUE, gradtol = 1e-5)
    vcov_m <- est$hessian %>% solve
    se <- vcov_m %>% diag %>% sqrt
    
    output = cbind(Estimate=est$estimate, se, zv= est$estimate/se, pv=pnorm(-abs(est$estimate/se)))

  }
  
  if(method=="QEBD"){
    
    # design matrix X ∈ Rnm×m2 where the j + m(k − 1)th row of X is y1k[j] ⊗ ej w
    # where k = 1, . . . , n and j = 1, . . . , m.
    # n: number of samples, m: number of people
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
    colnames(X)<- t_list
    
    ID= rep(1:n, each = m)
    
    gee <- geeglm(c(t(Y))~ 0+X, id=ID, family= binomial, corstr= "independence")
    est <- coef(gee)
    V <- vcov(gee)
    se <- diag(V) %>% sqrt
    
    output= cbind(est, se, est/se, pnorm(-abs(est/se)))

  }
  
  if(method=="logit"){
    
    # design matrix X ∈ Rnm×m2 where the j + m(k − 1)th row of X is y1k[j] ⊗ ej w
    # where k = 1, . . . , n and j = 1, . . . , m.
    # n: number of samples, m: number of people
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
    colnames(X)<- t_list
    
    fit <- glm(c(t(Y))~ 0+X, family= binomial)
    est <- coef(fit)
    V <- vcov(fit)
    se <- diag(V) %>% sqrt
    
    output= cbind(est, se, est/se, pnorm(-abs(est/se)))

    
  }  
  
  # Output reformat
  output= as.data.frame(output)
  output.names<- c("Estimate", "Std. Error", "z value", "Pr(|z|)")
  setnames(output, new= output.names)
  rownames(output)<- t_list
  
  output= output[match(t_list, row.names(output)),]
  
  return(output)

}

#Parameter setup
gamma<- 0.4
p<- 5
a<- rep(c(-1,1), length= p )
V<- (a %*% t(a)) * 0.4/2
diag(V) <- seq(-1.5, 1.5, length= p)

#One run
Y<- rising(n=300, V=V)
MLE_ising(Y, method="QEBD")
MLE_ising(Y, method="newton")
MLE_ising(Y, method="logit")

#Simulation
simdriver <- function(iter, seedn, n, V){
  set.seed(seedn)
  p<- dim(V)[1]
  m2 <- p*(p+1)/2
  
  Yn <- lapply(1:iter, function(i) rising(n, V))
  
  out1<- lapply(Yn, function(dat) MLE_ising(dat, method="QEBD")) %>% rbindlist(.) %>% as.matrix
  out2<- lapply(Yn, function(dat) MLE_ising(dat, method="logit")) %>% rbindlist(.) %>% as.matrix
  out3<- lapply(Yn, function(dat) MLE_ising(dat, method="newton")) %>% rbindlist(.) %>% as.matrix
  OUT<- list(qebd= out1, logit= out2, newton= out3)                            
  
  V[lower.tri(V)] <-   V[lower.tri(V)] * 2
  truth <-   V[lower.tri(V, diag=TRUE)]
  se<- apply(matrix(out3[,1], nrow = iter, ncol= m2, byrow = TRUE), 2, sd)
  mle_b <- apply(matrix(out3[,1], nrow = iter, ncol= m2, byrow = TRUE), 2, mean) - truth
  mle_re <- apply(matrix(out3[,2], nrow = iter, ncol= m2, byrow = TRUE), 2, mean)/ se
  logit_b <- apply(matrix(out2[,1], nrow = iter, ncol= m2, byrow = TRUE), 2, mean) - truth
  logit_re <- apply(matrix(out2[,2], nrow = iter, ncol= m2, byrow = TRUE), 2, mean)/ se
  qebd_b <- apply(matrix(out1[,1], nrow = iter, ncol= m2, byrow = TRUE), 2, mean) - truth
  qebd_re <- apply(matrix(out1[,2], nrow = iter, ncol= m2, byrow = TRUE), 2, mean)/ se
  
  K<- cbind(truth, se, mle_b, mle_re, logit_b, logit_re, qebd_b, qebd_re)

  return(list(OUT= OUT, K= K))
  
}
sim_table1<- simdriver(iter= 500, seedn= 88, n=300, V= V)
sim_table1$K %>% round(.,3)



