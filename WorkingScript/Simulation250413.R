# Simulation 0413: Include QIC and redesign parameter values
getwd()
source("Script/Functions2.R")

# Simulation --------------------------------------------------------------
SimDriver <- function(iter, seedn, n, V){
  set.seed(seedn)
  p<- dim(V)[1]
  m2 <- p*(p+1)/2
  it<- 100
  Yn <- lapply(1:iter, function(i) RandomQEBD(n, V))  # Generate Data
  
  params <- list(
    logit = list(method= "logit"),
    gee = list(method= "geepack", corstr = "independence"),
    mle = list(method = "mle")
  )
  
  results <- lapply(params, function(param) {
    lapply(Yn, function(dat) {
      tryCatch({
        do.call(EstQEBD, c(list(Y = dat), param))
      }, error = function(e) NA)
    })
  })
  result_coef <- lapply(results, function(layer1) {
    lapply(layer1, function(layer2) layer2[[1]])
  })
  result_fitobj <- lapply(results, function(layer1) {
    lapply(layer1, function(layer2) layer2[[2]])
  })
  
  V[lower.tri(V)] <-   V[lower.tri(V)] * 2
  truth <-   V[lower.tri(V, diag=TRUE)]
  mleresults<- result_coef$mle %>% do.call(rbind,.)
  true_se<- apply(matrix(mleresults[,1], nrow = iter, ncol= m2, byrow = TRUE), 2, sd)
  
  #convergecase<- lapply(names(results), function(x) !sapply(results[[x]], anyNA))
  #names(convergecase)<- names(results)
  #useful_iter<- lapply(convergecase,sum)
  
  ## 1- Generate the regular summary table based on GEE full model
  ## 2- Apply step QIC to GEE full model, record: 1) reduced model, 2) reduced QIC value
  reorder_columns <- function(df) {
    # Extract column names
    y<- df$y
    df$y<- NULL
    cols <- colnames(df)
    
    # Create logical vector for identical-numbered columns
    identical_cols <- grepl("X(\\d+)\\.\\1$", cols)
    
    # Get column indices for identical and non-identical
    identical_indices <- which(identical_cols)
    non_identical_indices <- which(!identical_cols)
    
    # Reorder columns
    cbind(y, df[, c(identical_indices, non_identical_indices)])
  }
  
  one_step <- function(D, ID, fixed_vec){
    fixed <- length(fixed_vec)
    m <- ncol(D) - fixed - 1
    qic <- rep(0, m+1)
    #qic_names <- rep(NA, m+1)
    qic[1] <- QIC(out)[1]
    #qic_names[1]<- "full model"
    for(index in 1:m){
      drop_one <- index + fixed + 1
      qic[index+1] <- QIC(geeglm(y~-1+., data=D[,-drop_one], id=ID, family=binomial))[1]
      #qic_names[index+1] <- paste("Excl", colnames(D)[drop_one])
    }
    cat(qic,"\n")
    return(qic)
    #return(rbind(qic_names,qic))
  }
  
  step.GEE <- function(out, fixed_vec){
    D <- out$model %>% reorder_columns
    #m <- ncol(D) - length(fixed_vec) -1
    drop_list<- (2+length(fixed_vec)) : ncol(D)

    for (j in 1: length(drop_list)){
      temp <- one_step(D, out$id, fixed_vec)
      cat("Step", j, ":")
      drop <- which(temp == min(temp))[1] #drop indicates exclusion of which variable would reduce model QIC the most
      if(drop == 1){ #drop == 1 refers to the first model (full model)
        cat("choose drop=NULL\n")
        break
      }else{
        k<- (1+length(fixed_vec))+(drop-1)
        cat("choose drop =", colnames(D)[k], "\n")
        D <- D[,-k]  #exclude this variable that would reduce QIC the most after excluding it
      }
    }
    oo <- geeglm(y~-1+., data=D, id=out$id, family=binomial)
    qic_oo <- QIC(oo)
    return(list(fit=oo, qic= qic_oo))
  }
  
  out_reducedGEE<- lapply(result_fitobj$gee, function(q) step.GEE(q, fixed_vec=1:5))
  # prepare a table to record which variables are retained in reduced model
  # eg 1000 iterations x each variables matrix
  hname<- colnames(result_fitobj$gee[[1]]$model)
  qic_vote<- data.frame(matrix(0, nrow= iter, ncol= (1+m2)))
  colnames(qic_vote)<- hname
  qic_vote<- reorder_columns(qic_vote)
  qic_reduced<- rep(0, iter)
  qic_full<- rep(0, iter)
  
  for (row in 1:iter) {
    retainedX<- colnames(out_reducedGEE[[row]]$fit$model)
    qic_vote[row, retainedX] <- 1
    qic_reduced[row] <- out_reducedGEE[[row]]$qic[1]
    qic_full[row] <- QIC(result_fitobj$gee[[row]])[1]
  }
  qic_out<- cbind(qic_vote, qic_full, qic_reduced)
  
  summarytable<- function(dat){
    if(!is.null(unlist(dat))){
      q<- length(dat)
      dat <- dat %>% do.call(rbind,.)
      estmean <- apply(matrix(dat[,1], nrow = q, ncol= m2, byrow = TRUE), 2, mean)
      bias <- estmean - truth
      se <- apply(matrix(dat[,2], nrow = q, ncol= m2, byrow = TRUE), 2, mean)
      re <- se/ true_se
      cbind(estmean,bias,se,re)
    } else {
      NA
    }
  }
  
  K<- lapply(names(result_coef), function(x) summarytable(result_coef[[x]]))
  names(K)<- names(result_coef)
  for(i in 1: length(K)){
    if(!(anyNA(K[[i]]))){
      colnames(K[[i]])<- paste0(colnames(K[[i]]),"_", names(result_coef)[i])
      K[[i]]<- cbind(truth= truth, true_se= true_se,K[[i]])
    } else {
      K[[i]]<- K[[i]]
    }
  }
  
  #OUT: raw long matrix of estimates for all simulation
  #K: summary of estimated, a list object
  return(list(OUT= results, K= K, qic_out= qic_out))
  
}



# Parameter setup
iter<- 100
seedn<- 88
n<- 300
gamma<- 0.4
m<- 5
a<- rep(c(-1,1), length= m)
V<- (a %*% t(a)) * 0.4/2
diag(V) <- seq(-1.5, 1.5, length= m)
V[1,4]<- V[4,1]<- V[1,5] <- V[5,1] <- V[2,4] <- V[4,2] <- V[2,5] <- V[5,2] <- V[3,4] <- V[4,3] <- V[3,5] <- V[5,3]<- 0
V[3,1]<- V[1,3]<- 0.6

sim1<- SimDriver(iter, seedn, n, V)
