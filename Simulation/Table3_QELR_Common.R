#Table2: Quadratic Exponential Logistic Regression (Common)

library(dplyr)
library(geepack)
source("Simulation/SimulationFunctions.R")


# Run One Example ---------------------------------------------------------
# dat <- GenData(beta= c(-2.4, -2.0, -2.6), gamma= -1.4, n=300, p= 15)
# EstQELR(dat$DAT) %>% round(.,3)


# Simulation --------------------------------------------------------------

SimDriver<- function(iter, seedn, beta, gamma, p, n_seq){
  set.seed(seedn)
  
  OUT <- array(0, dim=c(iter, length(n_seq), 12, 4))
  for (i in 1:iter){
    gc()
    dat <- GenData(beta, gamma, n=max(n_seq), p)
    for (n in 1:length(n_seq)){
      OUT[i,n,,] <-  EstQELR(dat$DAT[1:(n_seq[n]*p),])
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
  colnames(K)<- c("SE", "GLM","GLM Bias", "GLM SE","GLM R.E.", "GEEIND","GEEIND Bias","GEEIND SE","GEEIND R.E.", "GEEEXC","GEEEXC Bias","GEEEXC SE","GEEEXC R.E.")
  
  rowname<- expand.grid(b= c(paste0("beta",1:length(beta)-1),"gamma"),n= paste0("n=",n_seq))
  K<- cbind(rowname, Truth= rep(true_par, by= length(n_seq)), K)
  
  return(list(OUT= OUT, K= K))
  
}

# Parameter Setup
beta <- c(-2.4, -2.0, -2.6)
gamma <- -1.4
n_seq <- c(100,300,500)
p <- 15
iter <- 500
seedn <- 88

sim.table2<- SimDriver(iter, seedn, beta, gamma, p, n_seq)
sim.table2$K
table2<- sim.table2$K[,c(1:4,seq(6,16, by=2))]
table2


# Exclude Divergence Cases for GEE-EXC ------------------------------------
G<- sim.table2$OUT[1:iter,,,]

# Define and calculate divergence
threshold<- 2
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



# Table 2 Output ----------------------------------------------------------

xtable::xtable(table2, digits=3)

# Table 2 Output based on below setup:
# iter <- 500
# seedn <- 88
# beta <- c(-2.4, -2.0, -2.6)
# gamma <- -1.4
# n_seq <- c(100,300,500)
# p <- 15

# Output:
#   \begin{table}[ht]
#   \centering
#   \begin{tabular}{rrllrrrrrrrr}
#   \hline
#   Variable & Sample Size & Truth & Emp. S.D. & GLM Bias & GLM R.E. & GEE-IND Bias & GEE-IND R.E. & GEE-EXC Bias & GEE-EXC R.E.\\ 
#   \hline
#   beta0 & n=100 & -2.400 & 0.594 & 0.042 & 0.627 & 0.042 & 0.975 & 0.918 & 0.976 \\ 
#   beta1 & n=100 & -2.000 & 0.220 & -0.051 & 0.978 & -0.051 & 0.982 & 0.161 & 0.886 \\ 
#   beta2 & n=100 & -2.600 & 0.265 & -0.068 & 0.948 & -0.068 & 0.957 & 0.216 & 0.867 \\ 
#   gamma & n=100 & -1.400 & 0.298 & -0.082 & 0.555 & -0.082 & 0.937 & -0.438 & 0.999 \\ 
#   
#   beta0 & n=300 & -2.400 & 0.333 & 0.014 & 0.623 & 0.014 & 0.965 & -0.863 & 0.943 \\ 
#   beta1 & n=300 & -2.000 & 0.119 & -0.015 & 1.019 & -0.015 & 1.041 & 0.638 & 0.929 \\ 
#   beta2 & n=300 & -2.600 & 0.142 & -0.020 & 0.992 & -0.020 & 1.018 & -0.301 & 0.924 \\ 
#   gamma & n=300 & -1.400 & 0.167 & -0.023 & 0.546 & -0.023 & 0.938 & 0.969 & 1.026 \\ 
#   
#   beta0 & n=500 & -2.400 & 0.247 & 0.003 & 0.649 & 0.003 & 1.000 & -0.127 & 0.991 \\ 
#   beta1 & n=500 & -2.000 & 0.090 & -0.011 & 1.034 & -0.011 & 1.059 & -1.248 & 0.941 \\ 
#   beta2 & n=500 & -2.600 & 0.104 & -0.011 & 1.041 & -0.011 & 1.074 & 0.127 & 0.969 \\ 
#   gamma & n=500 & -1.400 & 0.125 & -0.013 & 0.560 & -0.013 & 0.963 & 0.371 & 1.061 \\ 
#   \hline
#   \end{tabular}
#   \end{table}

# divergence_rate <- c(0.530, 0.576, 0.598)

