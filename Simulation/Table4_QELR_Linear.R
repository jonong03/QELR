#Table3: Quadratic Exponential Logistic Regression (Linear)
library(dplyr)
library(geepack)
source("Simulation/SimulationFunctions.R")

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
table3<- sim.table3$K[,c(1:4,seq(6,16, by=2))]
table3


# Table 3 Output ----------------------------------------------------------

xtable::xtable(table3, digits=3)

# Table 3 Output based on below setup:
# iter <- 500
# seedn <- 88
# beta <- c(-2.4, -2.0, -2.6)
# gamma <- c(-1.2, -0.5)
# n_seq <- c(100,300,500)
# p <- 15

#  \begin{table}[ht]
#  \centering
#  \begin{tabular}{rllrrrrrrrr}
#  \hline
#  & b & n & Truth & SE & GLM Bias & GLM R.E. & GEEIND Bias & GEEIND R.E. & GEEEXC Bias & GEEEXC R.E. \\ 
#  \hline
#  beta0 & n=100 & -2.400 & 0.242 & -0.008 & 0.882 & -0.008 & 1.000 & 0.127 & 0.997 \\ 
#  beta1 & n=100 & -2.000 & 0.145 & -0.035 & 1.083 & -0.035 & 1.084 & -0.018 & 1.076 \\ 
#  beta2 & n=100 & -2.600 & 0.178 & -0.031 & 1.032 & -0.031 & 1.037 & -0.009 & 1.032 \\ 
#  gamma1 & n=100 & -1.200 & 0.238 & -0.045 & 0.761 & -0.045 & 0.976 & -0.123 & 0.990 \\ 
#  gamma2 & n=100 & -0.500 & 0.193 & -0.023 & 0.807 & -0.023 & 1.011 & -0.095 & 1.028 \\ 
#  beta0 & n=300 & -2.400 & 0.142 & 0.002 & 0.860 & 0.002 & 0.981 & 0.135 & 0.976 \\ 
#  beta1 & n=300 & -2.000 & 0.087 & -0.012 & 1.030 & -0.012 & 1.038 & 0.004 & 1.031 \\ 
#  beta2 & n=300 & -2.600 & 0.104 & -0.011 & 1.009 & -0.011 & 1.018 & 0.009 & 1.012 \\ 
#  gamma1 & n=300 & -1.200 & 0.136 & -0.018 & 0.757 & -0.018 & 0.989 & -0.097 & 1.005 \\ 
#  gamma2 & n=300 & -0.500 & 0.111 & -0.008 & 0.802 & -0.008 & 1.030 & -0.082 & 1.049 \\ 
#  beta0 & n=500 & -2.400 & 0.113 & 0.003 & 0.833 & 0.003 & 0.953 & 0.135 & 0.948 \\ 
#  beta1 & n=500 & -2.000 & 0.067 & -0.008 & 1.026 & -0.008 & 1.033 & 0.008 & 1.027 \\ 
#  beta2 & n=500 & -2.600 & 0.082 & -0.010 & 0.992 & -0.010 & 1.003 & 0.010 & 0.997 \\ 
#  gamma1 & n=500 & -1.200 & 0.108 & -0.009 & 0.738 & -0.009 & 0.973 & -0.088 & 0.990 \\ 
#  gamma2 & n=500 & -0.500 & 0.089 & -0.011 & 0.771 & -0.011 & 0.988 & -0.084 & 1.007 \\ 
#  \hline
#  \end{tabular}
#  \end{table}
