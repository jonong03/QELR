# Table 1: Quadratic Exponential Distribution

library(dplyr)
library(geepack)
source("Simulation/SimulationFunctions.R")

# Run One Example ---------------------------------------------------------
# Y<- RandomQEBD(n=300, V=V)
# EstQEBD(Y, method="mle")
# EstQEBD(Y, method="logit")
# EstQEBD(Y, method="geeind")



# Simulation --------------------------------------------------------------

SimDriver <- function(iter, seedn, n, V){
  set.seed(seedn)
  p<- dim(V)[1]
  m2 <- p*(p+1)/2
  
  Yn <- lapply(1:iter, function(i) RandomQEBD(n, V))
  
  out1<- lapply(Yn, function(dat) EstQEBD(dat, method="geeind")) %>% do.call(rbind,.) %>% as.matrix
  out2<- lapply(Yn, function(dat) EstQEBD(dat, method="logit")) %>% do.call(rbind,.)  %>% as.matrix
  out3<- lapply(Yn, function(dat) EstQEBD(dat, method="mle")) %>% do.call(rbind,.) %>% as.matrix
  OUT<- list(geeind= out1, logit= out2, mle= out3)                            
  
  V[lower.tri(V)] <-   V[lower.tri(V)] * 2
  truth <-   V[lower.tri(V, diag=TRUE)]
  se<- apply(matrix(out3[,1], nrow = iter, ncol= m2, byrow = TRUE), 2, sd)
  mle <- apply(matrix(out3[,1], nrow = iter, ncol= m2, byrow = TRUE), 2, mean)
  mle_b <- mle - truth
  mle_se <- apply(matrix(out3[,2], nrow = iter, ncol= m2, byrow = TRUE), 2, mean)
  mle_re <- mle_se/ se
  logit <- apply(matrix(out2[,1], nrow = iter, ncol= m2, byrow = TRUE), 2, mean) 
  logit_b <- logit - truth
  logit_se <- apply(matrix(out2[,2], nrow = iter, ncol= m2, byrow = TRUE), 2, mean)
  logit_re <- logit_se/ se
  qebd <- apply(matrix(out1[,1], nrow = iter, ncol= m2, byrow = TRUE), 2, mean)
  qebd_b <- qebd - truth
  qebd_se <- apply(matrix(out1[,2], nrow = iter, ncol= m2, byrow = TRUE), 2, mean)
  qebd_re <- qebd_se/ se
  
  K<- cbind(truth, se, mle, mle_b, mle_se, mle_re, logit, logit_b, logit_se, logit_re, qebd, qebd_b, qebd_se, qebd_re)
  
  return(list(OUT= OUT, K= K))
  
}

# Parameter setup
iter<- 500
seedn<- 88
n<- 300
gamma<- 0.4
m<- 5
a<- rep(c(-1,1), length= m)
V<- (a %*% t(a)) * 0.4/2
diag(V) <- seq(-1.5, 1.5, length= m)

sim1<- SimDriver(iter, seedn, n, V)
table1<- sim1$K[,c(1,2,4,6,8,10,12,14)] %>% round(.,3)
table1

# Reformat table 1
Q<- matrix(NA, m, m)
Q[upper.tri(Q, diag=TRUE)] = 0
loc<- which(Q==0, arr.ind=TRUE,useNames=TRUE) %>% as.data.frame %>% arrange(row) 
var_name<- paste0(loc$row, ":",loc$col)
rownames(table1)<- var_name
colnames(table1)<- c("Truth", "Emp. S.D.", "MLE Bias", "MLE R.E.", "GLM Bias", "GLM R.E.", "GEE-IND Bias", "GEE-IND R.E.")

top<- sapply(1:length(var_name), function(x) unlist(strsplit(var_name[x],":"))[1] == unlist(strsplit(var_name[x],":"))[2] ) %>% which
bottom<- sapply(1:length(var_name), function(x) unlist(strsplit(var_name[x],":"))[1] != unlist(strsplit(var_name[x],":"))[2] ) %>% which
table1<- table1[c(top,bottom),]
table1

xtable::xtable(table1, digits=3)


# Table 1 Output ----------------------------------------------------------

# Table 1 Output based on below setup:
# iter<- 500
# seedn<- 88
# n<- 300
# gamma<- 0.4
# m<- 5
# a<- rep(c(-1,1), length= m)
# V<- (a %*% t(a)) * 0.4/2
# diag(V) <- seq(-1.5, 1.5, length= m)


# Output:
#   \begin{table}[ht]
#   \centering
#   \begin{tabular}{rrrrrrrrr}
#   \hline
#   & Truth & Emp. S.D. & MLE Bias & MLE R.E. & GLM Bias & GLM R.E. & GEE-IND Bias & GEE-IND R.E. \\ 
#   \hline
#   1:1 & -1.500 & 0.439 & -0.063 & 1.045 & -0.064 & 0.766 & -0.064 & 1.046 \\ 
#   2:2 & -0.750 & 0.381 & -0.015 & 0.989 & -0.015 & 0.737 & -0.015 & 0.988 \\ 
#   3:3 & 0.000 & 0.333 & 0.000 & 0.996 & 0.000 & 0.737 & 0.000 & 0.995 \\ 
#   4:4 & 0.750 & 0.328 & 0.018 & 0.987 & 0.018 & 0.746 & 0.018 & 0.985 \\ 
#   5:5 & 1.500 & 0.318 & 0.017 & 1.009 & 0.017 & 0.777 & 0.017 & 1.004 \\ 
#   1:2 & -0.400 & 0.385 & -0.025 & 0.962 & -0.024 & 0.676 & -0.024 & 0.961 \\ 
#   1:3 & 0.400 & 0.317 & -0.015 & 0.932 & -0.015 & 0.655 & -0.015 & 0.932 \\ 
#   1:4 & -0.400 & 0.309 & -0.012 & 0.940 & -0.012 & 0.661 & -0.012 & 0.941 \\ 
#   1:5 & 0.400 & 0.400 & 0.054 & 1.043 & 0.054 & 0.734 & 0.054 & 1.043 \\ 
#   2:3 & -0.400 & 0.280 & -0.014 & 1.004 & -0.014 & 0.705 & -0.014 & 1.004 \\ 
#   2:4 & 0.400 & 0.291 & 0.011 & 0.992 & 0.011 & 0.696 & 0.011 & 0.992 \\ 
#   2:5 & -0.400 & 0.343 & -0.009 & 0.971 & -0.009 & 0.682 & -0.009 & 0.970 \\ 
#   3:4 & -0.400 & 0.249 & -0.027 & 0.980 & -0.027 & 0.688 & -0.027 & 0.979 \\ 
#   3:5 & 0.400 & 0.316 & 0.027 & 0.976 & 0.027 & 0.686 & 0.027 & 0.974 \\ 
#   4:5 & -0.400 & 0.323 & 0.005 & 0.986 & 0.005 & 0.692 & 0.005 & 0.984 \\ 
#   \hline
#   \end{tabular}
#   \end{table}

