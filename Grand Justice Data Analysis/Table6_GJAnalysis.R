# Grand Justice Data Analysis

rm(list=ls()); gc()
pacman::p_load(geepack, dplyr, data.table, xtable)

# Define functions
# A Design Matrix Function, including 1st, 2nd and 3rd order effect for Time
# type=1: Common Interaction, type=2: Linear Interaction
designMatrix_GJ<- function(datalist, type=2, subject_effect=TRUE){
  Y <- data.frame(datalist$Y); S <- datalist$S; VE <- data.frame(datalist$VE); Time <- datalist$Time
  n <- nrow(Y); p <- ncol(Y); d <- ncol(VE)
  
  m2= length(S)
  EQ_subject <- function(X){
    p <- nrow(X); d <- ncol(X)
    E <- matrix(0, p, p)
    for(i in 1:p) E[i,] <- sapply(1:p, function(j) ifelse(sum(X[j,] == X[i,]) == d, 1, 0))
    return(E)
  }
  
  if (type==2) {W <- lapply(1:m2, function(i) EQ_subject(S[[i]])); names(W)<- names(S)}
  if (type==1) {W <- list(matrix(1,p,p)); names(W)<- "Common"}
  
  ID <- rep(1:n, each=p)
  y <- c(t(Y))
  time <- c(t(Time))
  time <- time/ max(time)
  
  ### predictors
  if(subject_effect== TRUE) {
    q = p
    VE2 <- matrix(0, nrow=n, ncol=q)
    colnames(VE2) <- colnames(Y)
    VE<- cbind(VE2, VE)
  } else {
    q = 0
  }
  
  X <- matrix(0, n*p, q+d+length(W) )
  
  for(i in 1:n){
    for(j in 1:p){
      y0 <- Y[i,]; y0[j] <- 0
      #X[(i-1)*p+j,] <- c(VE[i,], sum(y0*W[j,]))
      if(subject_effect==TRUE) {
        VE[i,(1:p)] <- rep(0, p)
        VE[i,j] <- 1 
      }
      X[(i-1)*p+j,] <- unlist(c(VE[i,], sapply(W, function(x) sum(x[j,]*y0))))
    }
  }
  
  colnames(X) <- c(colnames(VE), names(W) )
  X2<- cbind(y=y, X, time=time, time2= I(time^2), time3= I(time^3)) %>% as.data.frame
  
  return(list(X= X2, ID=ID))
}

getwd()
# Import Data
dt1<- read.csv("QELR/Grand Justice Data Analysis/GJVotes_P10.csv")
dt2<- read.csv("QELR/Grand Justice Data Analysis/GJBackground_P10.csv")

# Split Data
Y= dt1[,9:23]   # Responses (Votes) of each subject
VE= dt1[,5:8]   # Covariates describing each responses: Issue and Case Outcome
P1= dt2[,2:4]   # First correlated variable of each subject: Prior Occupation
P2= dt2[,5:7]   # Second correlated variable of each subject: 

p<- ncol(Y)   # number of subjects
Time = sapply(1:p, function(i) difftime(dt1$CaseDate, dt2$JoinDate[i], units="days")/3650)  # Compute 
colnames(Time) <- paste0("Time", colnames(Y))

datalist <- list(Y = Y, VE =VE, S = list(PriorOccupation=P1, Education= P2), Time= Time)

# Generate Design Matrix for model
dtX <- designMatrix_GJ(datalist)   

# Fit a full model
outgee_all<- geeglm(y~ -1+., id=dtX$ID, data= dtX$X, family=binomial, corstr ="independence")
summary(outgee_all)

# Variable Backward Selection ---------------------------------------------
source("/Grand Justice Data Analysis/Functions/stepGEE.R")
stepoo<- stepGEE(outgee_all, fixed_vec=1:15)

# xtable
outgee_all %>% summary %>% coef %>% xtable(digits=3)
stepoo %>% summary %>% coef %>% xtable(digits=3)
