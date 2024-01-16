#####################################################################################################
### Parallel Simulation framework for Online Testing FWER-Procedures in a platform trial scenario ###
#####################################################################################################

# Provides a parallelized simulation function for Online Testing Procedures estimating the FWER and Power
# in a setup of Gaussian Random Variables with different possible dependency conditions.

library(MASS)
library(onlineFDR)
library(foreach)
library(ggplot2)
library(ggpubr)

source("test_procedures.R")
source("plotting.R")

### Main Simulation parallelized ###

run_simulationPlatform <- function(n_hyp = 50,
                           n_runs = 20000,
                           n_sample = 100,
                           persons_per_unit = 10,
                           offset = 2,
                           procedures = c("Online-Fallback","Adaptive-Spending","Continuous-Spending","Continuous-Graph"),
                           pi1_vec = c(seq(0.05,0.15,by = 0.05),seq(0.2,0.9,by = 0.1)),
                           mu_0 = 0,
                           mu_1 = 5/sqrt(n_sample),    # shrinked effect of observation (so power does not increase with sample size for comparison)
                           mu_C = 0,    
                           sigma = 1,
                           alpha = 0.05,
                           lambda = 0.5,
                           tau = 1
){
  n_procs  <- length(procedures)
  fwer_mat <- pwr_mat <- matrix(0,nrow = n_procs,ncol = length(pi1_vec))
  
  # Default gamma sequence #
  gam   <- (6/pi^2)/seq_len(n_hyp)^2

  
  # Set up parallel backend #
  n_cores <- parallel::detectCores() - 2
  clus   <- parallel::makeCluster(n_cores,type = "PSOCK")
  doParallel::registerDoParallel(cl = clus)
  
  for (pi1 in pi1_vec) {
    time_start <- Sys.time()
    print(paste("Computing for pi1 =",pi1,"..."))
    cur_data <- foreach (j = 1:n_runs,.combine = 'rbind',.packages = c("MASS","onlineFDR"),.export = c("continuous_spending","continuous_adaptive_graph") ) %dopar% {
      set.seed(j) 
      fun <- function(vec){
        result <- vector(mode = "numeric",length = n_hyp)
        for (i in 1:n_hyp) {
          result[i] <- mean(vec[((i-1)*persons_per_unit*offset + 1):((i-1)*persons_per_unit*offset + n_sample)])
        }
        return(result)
      }
      # Data Generation #
      mixed_vec <- rbinom(n=n_hyp,size = 1,prob = pi1)     # Choosing Null vs Alternative
      mu_vec1   <- mu_1 * mixed_vec + mu_0 * !mixed_vec    # True parameters 
      x_data  <- mvrnorm(n=n_sample,mu=mu_vec1,Sigma = diag(n_hyp))  
      y_data <- rnorm(n_sample + persons_per_unit * offset * (n_hyp - 1),mean = mu_C,sd = 1)
      z_scores <- sqrt(n_sample/2) * (colMeans(x_data) - fun(y_data)) # Test statistics
      p_values  <- 1 - pnorm(z_scores)
      
      fwe_vec <- pwp_vec <- vector(mode = "numeric",length = n_procs)
      
      for (proc_index in 1:n_procs){
        proc <- procedures[proc_index]
        if (proc == "Alpha-spending"){
          rejects <- Alpha_spending(p_values, alpha = alpha, gammai = gam)$R
        } else if (proc == "Online-Fallback"){
          rejects <- online_fallback(p_values, alpha = alpha, gammai = gam)$R
        } else if (proc == "Adaptive-Spending"){
          rejects <- ADDIS_spending(p_values, alpha = alpha, lambda = lambda, tau = 1, gammai = gam)$R
        } else if (proc == "Continuous-Spending"){
          rejects <- continuous_spending(z_scores, n_sample = n_sample, alpha = alpha, lambda = lambda,tau = 1, gammai = gam, closed = TRUE)$rejects
        } else if (proc == "Continuous-Graph"){
          rejects <- continuous_adaptive_graph(z_scores, n_sample = n_sample, alpha = alpha, lambda = lambda, tau = 1, gammai = gam, closed = TRUE)$rejects
        } else {
          stop(paste("The specified procedure",proc, "doesn't exist!"),sep = " ")
        }
        false_rejects <- rejects * !mixed_vec
        true_rejects  <- rejects * mixed_vec
        fwe           <- sum(false_rejects) >= 1
        pwp           <- sum(true_rejects)/max(1,sum(mixed_vec))
        
        fwe_vec[proc_index] <- fwe
        pwp_vec[proc_index] <- pwp
      }
      c(fwe_vec,pwp_vec)
    }     
    
    fwe_mat <- cur_data[,1:n_procs]
    pwp_mat <- cur_data[,(n_procs + 1):(2 * n_procs)]
    
    if (n_procs == 1){
      cur_fwer <- mean(fwe_mat)
      cur_pwr  <- mean(pwp_mat)
    } else{
      cur_fwer <- colMeans(fwe_mat)
      cur_pwr  <- colMeans(pwp_mat)
    }
    fwer_mat[,match(pi1,pi1_vec)] <- cur_fwer
    pwr_mat[,match(pi1,pi1_vec)]  <- cur_pwr
    
    print(Sys.time()-time_start)
  }
  parallel::stopCluster(cl = clus)
  
  # save data #
  rownames(fwer_mat) <- rownames(pwr_mat) <- procedures
  colnames(fwer_mat) <- colnames(pwr_mat) <- pi1_vec
  
  parameters        <- c(n_hyp,n_runs,n_sample,alpha,lambda,tau,mu_0,mu_1,mu_C,sigma,persons_per_unit,offset)
  names(parameters) <- c("n_hyp","n_runs","n_sample","alpha","lambda","tau","mu0","mu1","muC","sigma","persons_per_unit","offset")
  
  results <- list(fwer_mat,pwr_mat,parameters)
  names(results) <- c("fwer_mat","pwr_mat","parameters")
  
  save(results,file = paste("data","n_hyp",n_hyp,"lambda",lambda,"mu0",mu_0,"mu1",mu_1,format(Sys.time(), "%Y-%m-%d_%H-%M"),sep = "_"))
  
  # generate plot #
  conv_res <- convert_results(results)
  plot <- plotting(conv_res)
  save_plot(name =  paste("plot","n_hyp",n_hyp,"lambda",lambda,"mu0",mu_0,"mu1",mu_1,format(Sys.time(), "%Y-%m-%d_%H-%M"),sep = "_"),plot = plot)
  
  
  return(results)
}
