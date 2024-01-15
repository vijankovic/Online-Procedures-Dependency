##############################################################
### Adaptive procedures for potentially dependent p-values ###
##############################################################


# Continuous Adaptive-Spending #

continuous_spending <- function(z_scores, n_sample, gammai, lambda = 0.5, tau = 1, alpha = 0.05, closed = FALSE) {
  
  n_hyp <- length(z_scores)
  
  if (lambda < 0 || lambda > 1) {
    stop("lambda has to be between 0 and 1")
  }
  if (tau <= lambda || tau > 1) {
    stop("tau has to be between lambda and 1")
  }
  
  # Default value for gammai #
  if (missing(gammai)){
    gammai <- (6/pi^2)/(1:n_hyp)^2
  }
  # Linear interpolation of gammai #
  gam <- function(x,an) {
    k1 <- floor(x)
    k2 <- ceiling(x)
    y  <- x - k1
    return(an[k1] + y * (an[k2] - an[k1]))
  }
  
  p_values <- 1 - pnorm(z_scores)
  
  rejects <- alphai <- zetai <- vector("numeric",n_hyp)
  
  # Parametric bootstrap probabilities/ indicators #
  pseudo_resampling <- function(z_score,n_sample,lambda,tau,shrinkage = TRUE){
    if (shrinkage){
      shrinkage_estimate <- sqrt(sqrt(n_sample)/n_sample) * z_score 
    } else {
      shrinkage_estimate <- z_score 
    }
    prob_star <- pnorm(qnorm(1-lambda),shrinkage_estimate,1) - pnorm(qnorm(1-tau),shrinkage_estimate,1) # (parametric) bootstrap probability of p-value exceeding lambda
    return(prob_star)                                  
  }
  
  s <- 1 + gammai[1] * (tau - lambda - 0.5)  # scaling factor if tau - lambda != 0.5
  
  for (j in 1:n_hyp) {
    if (closed) {
      alphai[j]  <- alpha * (tau - lambda)/s * gam(1 + sum((1-rejects[seq_len(j-1)])*zetai[seq_len(j-1)]),gammai)
    } else {
      alphai[j]  <- alpha * (tau - lambda)/s * gam(1 + sum(zetai[seq_len(j-1)]),gammai)
    }
    rejects[j] <- p_values[j] <= alphai[j]
    zetai[j] <- pseudo_resampling(z_scores[j],n_sample,lambda,tau)
  } 
  
  results <- as.data.frame(cbind(p_values,alphai,zetai,rejects),row.names = 1:n_hyp)
  return(results)
}


# Continuous Adaptive-Graph #

continuous_adaptive_graph <- function(z_scores, n_sample, gammai, w, lambda = 0.5, tau = 1, alpha = 0.05, closed = FALSE){
  
  n_hyp <- length(z_scores)
  
  if (lambda < 0 || lambda > 1) {
    stop("lambda has to be between 0 and 1")
  }
  if (tau <= lambda || tau > 1) {
    stop("tau has to be between lambda and 1")
  }
  
  # Default value for gammai #
  if (missing(gammai)){
    gammai <- (6/pi^2)/(1:n_hyp)^2
  }
  
  # Default for w
  if (missing(w)){
    w=abs(matrix(1:n_hyp-1 , nrow = n_hyp, ncol = n_hyp, byrow = TRUE) - (1:n_hyp-1))
    w[w==0]=1
    w=matrix(gammai[w],n_hyp,n_hyp)
    w[upper.tri(w)==0]=0
  }
  
  p_values <- 1 - pnorm(z_scores)
  
  rejects <- alphai <- zetai <- vector("numeric",n_hyp)
  
  # Parametric bootstrap probabilities #
  pseudo_resampling <- function(z_score,n_sample,lambda,tau,shrinkage = TRUE){ 
    if (shrinkage){
      shrinkage_estimate <- sqrt(sqrt(n_sample)/n_sample)*z_score 
    } else {
      shrinkage_estimate <- z_score 
    }
    prob_star <- pnorm(qnorm(1-lambda),shrinkage_estimate,1) - pnorm(qnorm(1-tau),shrinkage_estimate,1) # (parametric) bootstrap probability of p-value exceeding lambda
    return(prob_star)                                  
  }
  
  #Adaptive-Graph with Bootstrap estimates:
  
  for (j in 1:n_hyp) {                     
    alphai[j]  <- (tau-lambda)*(alpha*gammai[j]+sum(alphai[seq_len(j-1)]*(1-zetai[seq_len(j-1)])*w[seq_len(j-1),j]/(tau-lambda))) # CHANGED
    rejects[j] <- p_values[j] <= alphai[j]
    prob_star <- pseudo_resampling(z_scores[j],n_sample,lambda,tau)
    if (closed){
      zetai[j]   <- (1 - rejects[j]) * prob_star
    } else{
      zetai[j]   <- prob_star
    }
  } 
  
  results <- as.data.frame(cbind(p_values,alphai,zetai,rejects),row.names = 1:n_hyp)
  return(results)
}
