#Table2: Quadratic Exponential Logistic Regression (Common)

rm(list=ls()); gc()
library(dplyr)
library(geepack)

beta <- c(-2.4, -2.0, -2.6)
gamma <- -1.4
n_seq <- c(100,300,500)
nmax <- max(n_seq)
p <- 15
iter <- 100
seedn <- 88

rqebd <- function(n, V){
  m = ncol(V)
  config0 = expand.grid(rep(list(c(0,1)),m)) %>% as.matrix
  
  pmf <- function(V, y) exp(sum((V%*%y)*y))
  pr <- sapply(1:nrow(config0), function(i) pmf(V, config0[i,]))
  cdf <- cumsum(pr)/sum(pr)
  Y <- sapply(1:n, function(x) config0[sum(cdf < runif(1)) + 1,])
  return(t(Y))
}

datagen <- function(beta, gamma, n, p){
  nbeta <- length(beta)
  X<- cbind(1, sapply(1:(nbeta-1), function(i) rnorm(n*p))) #
  P<- lapply(1:n, function(i) matrix(gamma, p, p))
  Y<- matrix(0,n,p)
  
  for (i in 1:n){
    diag(P[[i]])<-  X[(i*p-p+1): (i*p),] %*% beta
    Y[i,]<- rqebd(1, P[[i]])
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

simdriver<- function(iter, seedn, beta, gamma, p, n_seq){
  set.seed(seedn)
  
  OUT <- array(0, dim=c(iter, length(n_seq), 12, 4))
  for (i in 1:iter){
    gc()
    dat <- datagen(beta, gamma, n=max(n_seq), p)
    for (n in 1:length(n_seq)){
      OUT[i,n,,] <-  QELR(dat$DAT[1:(n_seq[n]*p),])
    }
    cat("iter ",i, " is done at", date(),"\n")
  }
  
  true_par <- c(beta, gamma)
  G <- OUT[1:iter,,,]
  se <- t(apply(G[,,1:4,1], c(2,3), sd))
  glm <- t(apply(G[,,1:4,1], c(2,3), mean))
  glm_b <- glm - true_par
  glm_se <- t(apply(G[,,1:4,2], c(2,3), mean))
  glm_re <- glm_se/se
  geeind <- t(apply(G[,,1:4+4,1], c(2,3), mean))
  geeind_b <- geeind - true_par
  geeind_se <- t(apply(G[,,1:4+4,2], c(2,3), mean))
  geeind_re <- geeind_se/se
  geeexc <- t(apply(G[,,1:4+8,1], c(2,3), mean))
  geeexc_b <- geeexc - true_par
  geeexc_se <- t(apply(G[,,1:4+8,2], c(2,3), mean))
  geeexc_re <- geeexc_se/se
  K <- cbind(c(se), c(glm), c(glm_b), c(glm_se), c(glm_re), c(geeind), c(geeind_b), c(geeind_se), c(geeind_re), c(geeexc), c(geeexc_b), c(geeexc_se), c(geeexc_re))
  colnames(K)<- c("SE", "GLM","GLM_Bias", "GLM_SE","GLM R.E.", "GEEIND","GEEIND_Bias","GEEIND_SE","GEEIND R.E.", "GEEEXC","GEEEXC_Bias","GEEEXC_SE","GEEEXC R.E.")
  
  rowname<- expand.grid(b= c(paste0("beta",1:length(beta)-1),"gamma" ),n= paste0("n=",n_seq))
  K<- cbind(rowname, Truth= rep(true_par, by= length(n_seq)), K)
  
  return(list(OUT= OUT, K= K))
  
}

sim1<- simdriver(iter, seedn, beta, gamma, p, n_seq)
sim1$K
table2<- sim1$K[,c(1:4,seq(6,16, by=2))]
table2



# Exclude Divergence Cases for GEE-EXC ------------------------------------
G<- sim1$OUT[1:iter,,,]

# Define and calculate divergence
threshold<- 1
max.se<- apply(G[,,1:4+8,2],c(1,2),max)  # Max SE in each iteration
divergence<- max.se > threshold   # If max SE of any sets of estimates > threshold, define as a divergence case.
divergence_rate <- apply(divergence,2, mean)
divergence_rate


# Recalculate statistics after excluding divergence cases

geeexc<- NULL; geeexc_se<- NULL
for (i in 1:length(n_seq)){
  new.est<- apply(G[!divergence[,i],i,1:4+8,1], 2, mean)
  new.se<- apply(G[!divergence[,i],i,1:4+8,2], 2, mean)
  geeexc<- cbind(geeexc, new.est)
  geeexc_se<- cbind(geeexc_se, new.se)
}

true_par <- c(beta, gamma)
geeexc_b <- geeexc - true_par
geeexc_b
se <- t(apply(G[,,1:4,1], c(2,3), sd))
geeexc_re <- geeexc_se/se
geeexc_re


# Update Table 2
table2$`GEEEXC_Bias`<- c(geeexc_b)
table2$`GEEEXC R.E.`<- c(geeexc_re)

table2





#One Run
#dat <- datagen(beta, gamma, n=max(n_seq), p)
#QELR(dat$DAT) %>% round(.,3)
