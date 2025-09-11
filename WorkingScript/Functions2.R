# Table 1: Quadratic Exponential Distribution
pacman::p_load(dplyr, geepack)

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


# Function 3: Parameter estimation 
# Y: a "n x m" qebd data matrix.
# method: "mle" for maximum likelihood estimation, "logit" for logistic regression, "geeind" for GEE with independence working correlation.
EstQEBD <- function(Y, method, ...){
  args <- list(...)
  
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
    fit <- NA  # this object is just to facilitate simdriver250413
  }
  if(method=="logit"){
    
    X<- designX(Y, var_name)
    
    fit <- glm(c(t(Y))~ 0+X, family= binomial)
    est <- coef(fit)
    V <- vcov(fit)
    se <- diag(V) %>% sqrt
    
    output<- cbind(est, se, est/se, pnorm(-abs(est/se)))
    
  }  
  
  if (method %in% c("geepack", "gee")) {
    X <- designX(Y, var_name)
    ID <- rep(1:n, each = m)
    GEEdata<- data.frame(cbind(y= c(t(Y)), X))
    if (!"var.type" %in% names(args)) {
      var.type <- "robust" 
    } else {
      var.type <- args$var.type
      args$var.type <- NULL
    }
    
    gee_args <- list( formula = y~ -1 + ., data= GEEdata, id = ID, family = binomial)
    gee_args <- append(gee_args, args)
    
    #if (method == "geepack") {
    #  geeobj <- do.call(geepack::geeglm, gee_args)
    #  V <- vcov(geeobj) } 
    #if (method == "gee") {
    #  geeobj <- do.call(gee::gee, gee_args)
    #  if (var.type == "naive") {
    #    V <- geeobj$naive.variance
    #  } else {
    #    V <- geeobj$robust.variance
    #  }
    #}
    
    #est <- coef(geeobj)
    #se <- sqrt(diag(V))
    if (method == "geepack") fit <- do.call(geepack::geeglm, gee_args)
    if (method == "gee") fit <- do.call(gee::gee, gee_args)
    
    output <- summary(fit)$coefficients
  }
  
  # Reformat Output
  output<- as.matrix(output)
  if(!(method %in% c("geepack", "gee"))) colnames(output)<- c("Estimate", "Std. Error", "z value", "Pr(|z|)")
  rownames(output)<- var_name
  
  if(exists("fit")) {
    return(list(output= output, fit= fit)) 
  } else {
    return(output)
  }
  
}


