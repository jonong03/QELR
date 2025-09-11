# StepGEE Function:
# A Backward Variable Selection Procedure for a GEE fitted model, based on QIC
# outgee: a gee fitted object
# fixed_vec: a numerical vector indicating which column should be retained

stepGEE<- function(outgee, fixed_vec){
  D<- outgee$model  # Design Matrix
  ID<- outgee$id
  m<- ncol(D)-length(fixed_vec)
  
  qic_vec<- rep(0, times= m)
  qic_vec[1]<- QIC(outgee)[1]  #step1: store starting model's QIC
  cat("Step 0 - Full Model QIC:", qic_vec[1], "\n")
  for(i in 1:m){
    cat("Step",i,"\n")
    drop_candidate<- (2+length(fixed_vec)):ncol(D)   # 1st column is 'y'
    # colnames(D)[drop_candidate]
    reduced_gee<- sapply(drop_candidate, function(j) QIC(geeglm(y~ -1+., id=ID, data= D[,-j], family=binomial, corstr ="independence")) [1] ) #step2
    cat(reduced_gee, "\n")
    drop<- which.min(reduced_gee) #step3: Best Model has the smallest QIC
    # drop_candidate[drop] #step3: Which variable to drop
    if(qic_vec[i]<= min(reduced_gee)) {  # If current QIC is better/ smaller, stop
      cat("Drop: NULL \n")
      break
    } else {    # Else: drop the variable that produced the smallest QIC
      cat("Drop: ", colnames(D)[drop_candidate[drop]],"; Updated QIC:",reduced_gee[drop] ,"\n")
      D<- D[,-drop_candidate[drop]]
    }
    qic_vec[i+1]<- reduced_gee[drop]
  }
  out<- geeglm(y~ -1+., id=ID, data= D, family=binomial, corstr ="independence")
  cat(colnames(D), "\n")
  cat("QIC: ", qic_vec[i], "\n")  
  return(out)
}

