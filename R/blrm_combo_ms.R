if(getRversion() >= "2.15.1") utils::globalVariables(c("seeds", "probability"))

# the core function: `blrm_combo_core'
blrm_combo_core <- function(prior, data){
  
  ################### the generic functions to be used ###################
  # probability of toxicity for the drug combo
  tox_prob <- function(d1, d2, param){  
    # d1: log(dose1/ref_dose1)
    # d2: log(dose2/ref_dose2)
    
    # drug 1
    logit_pi1 <- log(param[1]) + param[2] * d1
    pi1 <- exp(logit_pi1) / (1 + exp(logit_pi1))
    
    # drug 2
    logit_pi2 <- log(param[3]) + param[4] * d2
    pi2 <- exp(logit_pi2) / (1 + exp(logit_pi2))
    
    # odds of probability of toxicity of the drug combo
    odds_pi12 <- exp(param[5] * exp(d1) * exp(d2)) * (pi1  + pi2 - pi1 * pi2) / ((1 - pi1) * (1 - pi2))
    
    # probability of toxicity of the drug combo
    pi12 <- odds_pi12 / (1 + odds_pi12)
    
    return(list(pi1=pi1, pi2=pi2, pi12=pi12))
  }
  ##
  
  # prior for the parameters
  prior_summary <- function(prior){ 
    # mu: mean
    # se: standard deviation
    # corr: correlation
    # extract prior components
    prior.drug1 <- prior[[1]]  # drug 1
    prior.drug2 <- prior[[2]]  # drug 2
    prior.inter <- prior[[3]]  # their interaction
    
    mean1 <- prior.drug1[[1]]
    se1 <- prior.drug1[[2]]
    corr1 <- prior.drug1[[3]]
    
    mean2 <- prior.drug2[[1]]
    se2 <- prior.drug2[[2]]
    corr2 <- prior.drug2[[3]]
    
    mean3 <- prior.inter[[1]]
    se3 <- prior.inter[[2]]
    
    covariance_matrix1 <- matrix(c(se1[1]**2, se1[1]*se1[2]*corr1, se1[1]*se1[2]*corr1, se1[2]**2), ncol=2, byrow=T) # variance-covariance structure for drug 1
    covariance_matrix2 <- matrix(c(se2[1]**2, se2[1]*se2[2]*corr2, se2[1]*se2[2]*corr2, se2[2]**2), ncol=2, byrow=T) # variance-covaraicne structure for drug 2
    
    return(list(prior1=list(mean1, covariance_matrix1), 
                prior2=list(mean2, covariance_matrix2), 
                prior3=list(mean3, se3**2)))
  }
  ##
  
  # go through one observed cohort
  single_cohort <- function(){  
    
    # data for JAGS
    mydata <- list(nb_pat=sum(cohort_size), 
                   s=data_dlt, 
                   d1=data_sdose1, 
                   d2=data_sdose2, 
                   p1=prior_para$prior1[[1]], # prior mean for drug 1
                   p2=prior_para$prior1[[2]], # prior covariance matrix for drug 1
                   p3=prior_para$prior2[[1]], # prior mean for drug 2
                   p4=prior_para$prior2[[2]], # prior covariance matrix for drug 2
                   p5=prior_para$prior3[[1]], # prior mean for the interaction term, eta
                   p6=prior_para$prior3[[2]]) # prior variance for the interaction term, eta
    
    # number of iterations per chain, e.g., (1 + 0.2) * 10000/2, since the first 20% will be chopped off;
    # there are 2 parallel chains
    niters <- (1 + burn_in) * nsamples/2
    
    # justify the model for library(rjags)
    modelstring <- "model
    {
    Omega1[1:2,1:2] <- inverse(p2[,])
    log.alpha[1:2] ~ dmnorm(p1[], Omega1[,])
    
    Omega2[1:2,1:2] <- inverse(p4[,])
    log.alpha[3:4] ~ dmnorm(p3[], Omega2[,])
    
    tau <- 1/p6
    eta ~ dnorm(p5, tau)
    
    alpha[1] <- exp(log.alpha[1])
    alpha[2] <- exp(log.alpha[2])
    alpha[3] <- exp(log.alpha[3])
    alpha[4] <- exp(log.alpha[4])
    
    for(i in 1:nb_pat){
    logit(pi1[i]) <- log.alpha[1] + alpha[2] * d1[i]
    logit(pi2[i]) <- log.alpha[3] + alpha[4] * d2[i]
    odds_pi12[i] <- exp(eta * exp(d1[i]) * exp(d2[i])) * (pi1[i] + pi2[i] - pi1[i] * pi2[i])/((1 - pi1[i]) * (1 - pi2[i]))
    pi12[i] <- odds_pi12[i]/(1 + odds_pi12[i])
    s[i] ~ dbern(pi12[i])
    }
    }"
    
    # specify starting values for the parameters, random number generating (RNG) seeds, and RNG modules
    inits.list <- list(list(log.alpha=c(-3,0,-3,0), eta=0, .RNG.seed=seeds[1], .RNG.name="base::Wichmann-Hill"), 
                       list(log.alpha=c(-3,0,-3,0), eta=0, .RNG.seed=seeds[2], .RNG.name="base::Wichmann-Hill"))
    
    # rjags::jags.model: to create an object representing a Bayesian graphical model, specified 
    # with a BUGS-language description of the prior distribution, and a set of data; provide initial 
    # values for parameters with vague prior distributions
    jagsobj <- rjags::jags.model(textConnection(modelstring), 
                                 data=mydata, n.chains=2, quiet=TRUE, inits=inits.list) # n.chains: the number of parallel chains for the model
    update(jagsobj, n.iter=niters, progress.bar="none")
    res <- rjags::jags.samples(jagsobj, c("alpha", "eta"), n.iter=niters, progress.bar="none")
    
    # chop off the burn-in iterations, e.g., 20% of the number of targeted samples (e.g., 10000)
    # for each parameter; there are 2 chains generated in parallel
    alpha1 <- res$alpha[1,,][-c(1:(burn_in*nsamples/2)),]
    beta1 <- res$alpha[2,,][-c(1:(burn_in*nsamples/2)),]
    alpha2 <- res$alpha[3,,][-c(1:(burn_in*nsamples/2)),]
    beta2 <- res$alpha[4,,][-c(1:(burn_in*nsamples/2)),]
    eta <- res$eta[1,,][-c(1:(burn_in*nsamples/2)),]
    # the posterior samples of interest
    posterior_param <- cbind(c(alpha1), c(beta1), c(alpha2), c(beta2), c(eta)) # dimension of (nsamples * 5)
    
    # (posterior) parameter estimates: log(alpha1.hat), log(alpha2.hat), log(beta1.hat), log(beta2.hat) and eta.hat
    log_para <- apply(posterior_param[,1:4], 2, log)
    log_para <- cbind(log_para, posterior_param[,5])
    
    para_hat <- apply(log_para, 2, mean)  # MAP estimates
    para_sd <-  apply(log_para, 2, sd)    # standard deviation of the MAP estimates
    para_corr1 <- boot::corr(log_para[,1:2])    # correlation coefficient between log(alpha1.hat) and log(beta1.hat)
    para_corr2 <- boot::corr(log_para[,3:4])    # correlation coefficient between log(alpha2.hat) and log(beta2.hat)
    posterior_para_summary <- list(para_hat=para_hat, para_sd=para_sd, para_corr1=para_corr1, para_corr2=para_corr2)
    
    # collect posterior point estimates for pi12, i.e., DLT rate.hat | data: for each set 
    # of posterior parameter point estimates, we have a posterior point estimate of DLT rate
    samples_sdose <- matrix(0, nprov_dose, nsamples)
    for(i in 1:nprov_dose){
      for(j in 1:nsamples){
        samples_sdose[i,j] <- tox_prob(sprov_dose[i,1], sprov_dose[i,2], posterior_param[j,])$pi12
      }
    }
    
    # interval probabilities by dose: obtain the probability that DLT rate lies in each category interval
    posterior_prob_summary <- matrix(0, length(category_bound)+1, nprov_dose)
    for(i in 1:nprov_dose){
      posterior_prob_summary[,i] <- as.numeric(table(cut(samples_sdose[i,], breaks=c(0, category_bound[1], category_bound[2], 1), right=TRUE))/nsamples)
    }
    
    # P(pi12 | data): mean, standard deviation, 0.025, 0.5(median), and 0.975 quantile
    posterior_pi_summary <- matrix(0, 5, nprov_dose)
    for(i in 1:nprov_dose){
      posterior_pi_summary[,i] <- c(mean(samples_sdose[i,]), sd(samples_sdose[i,]), quantile(samples_sdose[i,], c(0.025, 0.5, 0.975)))
    }
    
    # the next dose: select the dose with the highest probability of targeted toxicity among all the safe doses
    ## justify safe doses
    safe_dose_range <- (posterior_prob_summary[3,] <= ewoc) # posterior_prob_summary[3,]: (0.33, 1]
    if(sum(safe_dose_range) != 0){
      
      ## convert matrix to data.frame
      posterior_prob_summary <- as.data.frame(posterior_prob_summary)
      posterior_pi_summary <- as.data.frame(posterior_pi_summary)
      names(posterior_prob_summary) <- paste(1:nrow(prov_dose))
      names(posterior_pi_summary) <- names(posterior_prob_summary)
      
      # justify the dose with the highest probability of targeted toxicity
      next_dose_index <- as.numeric(names(which.max(posterior_prob_summary[2,safe_dose_range]))) # justify the dose index
      next_dose_level <- prov_dose[next_dose_index,]                               # the recommended drug combo level
      next_dose_posterior_prob_summary <- posterior_prob_summary[,next_dose_index] # interval probabilities by dose
      next_dose_posterior_pi_summary <- posterior_pi_summary[,next_dose_index]     # posterior summary of DLT rate
      
      # combine them to a list
      next_dose <- list(index=next_dose_index,
                        dose=next_dose_level, 
                        posterior_prob_summary=list(next_dose_posterior_prob_summary), 
                        posterior_pi_summary=list(next_dose_posterior_pi_summary))     
    }else{
      next_dose <- NULL
    }
    
    # convert back to matrix
    posterior_prob_summary <- as.matrix(posterior_prob_summary)
    posterior_pi_summary <- as.matrix(posterior_pi_summary)
    
    # return what is needed
    return(list(posterior_prob_summary=posterior_prob_summary,
                posterior_para_summary=posterior_para_summary,
                posterior_pi_summary=posterior_pi_summary,
                next_dose=next_dose))
  }
  # end of a single cohort
  
  ########################## done with defining generic functions ##########################
  
  # extract prior
  prior_para <- prior_summary(prior)
  
  # extract data
  seeds <- data$seeds
  nsamples <- data$nsamples 
  burn_in <- data$burn_in
  n_pat <- data$n_pat
  dlt_test <- data$dlt
  ewoc <- data$ewoc
  category_bound <- data$category_bound
  category_name <- data$category_name
  
  # drug specific information
  ## drug 1
  drug1_name <- data$drug1_name
  prov_dose1 <- data$prov_dose1
  ref_dose1 <- data$ref_dose1
  dose1_unit <- data$dose1_unit
  dose1 <- data$dose1
  
  ## drug 2
  drug2_name <- data$drug2_name
  prov_dose2 <- data$prov_dose2
  ref_dose2 <- data$ref_dose2
  dose2_unit <- data$dose2_unit
  dose2 <- data$dose2
  
  # further process the data
  sdose1 <- log(dose1/ref_dose1) # log(standardized tested doses for drug 1)
  sdose2 <- log(dose2/ref_dose2) # log(standardized tested doses for drug 2)
  sdose <- cbind(sdose1, sdose2) # the tested dose combination set: we will go through each of them !!!
  ndose1 <- length(sdose1)       # the number of tested doses for drug 1
  ndose2 <- length(sdose2)       # the number of tested doses for drug 2
  ndose <- nrow(sdose)           # the total number of tested dose combinations
  
  prov_dose <- expand.grid(prov_dose1, prov_dose2)      # the provisional dose set
  sprov_dose <- expand.grid(log(prov_dose1/ref_dose1), log(prov_dose2/ref_dose2)) # the standardized provisional dose set      
  nprov_dose1 <- length(prov_dose1)                     # the number of provisional doses for drug 1          
  nprov_dose2 <- length(prov_dose2)                     # the number of provisinoal doses for drug 2
  nprov_dose <- nrow(prov_dose)                         # the number of provisional doses for the drug combination set
  
  # go through all the tested doses cumulatively
  current_dose_index <- 1
  dose_index <- cohort_size <- dlt <- NULL
  data_sdose1 <- data_sdose2 <- data_dlt <- NULL
  
  # save the results for each tested dose as a list
  cohort_all <- list()
  
  # go through each tested drug combo dose level
  # progress indicator
  t <- 1
  while(current_dose_index <= ndose){
    
    current_cohort_size <- n_pat[current_dose_index]
    current_dlt <- dlt_test[current_dose_index]
    
    # ... cumutatively
    dose_index <- c(dose_index, current_dose_index)
    cohort_size <- c(cohort_size, current_cohort_size)
    dlt <- c(dlt, current_dlt)
    
    # for each cohort, statistical analysis is based on all the accumulated observed cohorts information
    current_data_sdose1 <- rep(sdose[current_dose_index,1], current_cohort_size)
    current_data_sdose2 <- rep(sdose[current_dose_index,2], current_cohort_size)
    data_sdose1 <- c(data_sdose1, current_data_sdose1)
    data_sdose2 <- c(data_sdose2, current_data_sdose2)
    current_data_dlt <- ifelse(rep(current_dlt==0, current_cohort_size), rep(0, current_cohort_size), c(rep(1, current_dlt), rep(0, current_cohort_size-current_dlt)))
    data_dlt <- c(data_dlt, current_data_dlt)
    # start going through each observed cohort
    cohort <- single_cohort()
    cohort_all[[current_dose_index]] <- list(dose_index=dose_index,
                                             cohort_size=cohort_size,
                                             dlt=dlt,
                                             posterior_prob=cohort$posterior_prob_summary,
                                             posterior_para=cohort$posterior_para_summary,
                                             posterior_pi=cohort$posterior_pi_summary,
                                             next_dose=cohort$next_dose)
    if(is.null(cohort$next_dose)){
      # next dose can not be justified
      break
    }else{
      # continue to the next tested dose
      current_dose_index <- current_dose_index + 1
    }
    
    # progress meter
    cat("tested drug combination ", t, " is done.\n", sep="")
    t <- t + 1
  }
  # done with going through all the tested drug combinations
  
  # we only take interest to the final-round information!
  cohort_final <- cohort_all[[length(cohort_all)]]
  # extract posterior summary
  prob_posterior <- cohort_final$posterior_prob
  pi_posterior <- cohort_final$posterior_pi
  para_posterior <- cohort_final$posterior_para
  # extract the next dose level information
  next_dose <- cohort_final$next_dose
  
  # construct the framework for the simulation output
  current_summary <- matrix(0, 25+2*nprov_dose+3+5+3+1+1+5, max(nprov_dose1, nprov_dose2)+3)
  current_summary <- as.data.frame(current_summary)
  for(i in 1:nrow(current_summary)){
    for(j in 1:ncol(current_summary)){
      current_summary[i,j] <- ""
    }
  }
  
  # Header
  current_summary[1,1] <- "Header"
  current_summary[2,2:4] <- c("simulation running date", "drug name", "RNG seeds")
  current_summary[3,2:4] <- c(paste(Sys.Date()), paste(c(drug1_name, drug2_name), collapse=", "), paste(seeds, collapse=", "))
  
  # Input
  current_summary[4,1] <- "User-specified Input"
  ## prior
  current_summary[5,2] <- "prior"
  current_summary[5,4:8] <- c("log(alpha1)", "log(beta1)", "log(alpha2)", "log(beta2)", "eta")
  current_summary[6,3] <- "mean"
  current_summary[7,3] <- "standard deviation"
  current_summary[8,3] <- "correlation"
  current_summary[6,4:8] <- c(prior[[1]]$mean, prior[[2]]$mean, prior[[3]]$mean)
  current_summary[7,4:8] <- c(prior[[1]]$se, prior[[2]]$se, prior[[3]]$se)
  current_summary[8,c(4,6)] <- c(prior[[1]]$corr, prior[[2]]$corr)
  
  ## data
  current_summary[9,2] <- "data"
  current_summary[10:22,3] <- c("burn-in", "number of simulated samples", "dose unit", 
                                "provisional doses of drug 1", "provisional doses of drug 2", 
                                "reference dose", "tested doses of drug 1", "tested doses of drug 2", 
                                "number of patients", "number of DLTs", "category bounds", "category names", "EWOC")
  current_summary[10,4] <- burn_in
  current_summary[11,4] <- nsamples
  current_summary[12,4:5] <- c(paste0("drug1: ", dose1_unit), paste0("drug2: ", dose2_unit))
  current_summary[13,4:(4+nprov_dose1-1)] <- prov_dose1  # provisional doses of drug 1
  current_summary[14,4:(4+nprov_dose2-1)] <- prov_dose2  # provisional doses of drug 2
  current_summary[15,4:5] <- c(paste0("drug1: ", ref_dose1), (paste0("drug2: ", ref_dose2)))
  current_summary[16,4:(4+ndose1-1)] <- dose1
  current_summary[17,4:(4+ndose2-1)] <- dose2
  current_summary[18,4:(4+ndose-1)] <- n_pat
  current_summary[19,4:(4+ndose-1)] <- dlt
  current_summary[20,4:(4+length(category_bound)-1)] <- category_bound
  current_summary[21,4:(4+length(category_name)-1)] <- category_name
  current_summary[22,4] <- ewoc
  
  # output
  current_summary[23,1] <- "Simulation Output"
  ## posterior DLT rate estimates
  current_summary[24,2] <- "posterior DLT rate estimates"
  current_summary[25,5:9] <- c("mean", "standard deviation", "2.5% quantile", "median", "97.5% quantile")
  current_summary[25,3:4] <- paste("drug", 1:2)
  current_summary[26:(26+nprov_dose-1),3:4] <- prov_dose
  for(i in 1:nprov_dose){
    current_summary[26+i-1,5:9] <- round(pi_posterior[,i], 5)
  }
  
  ## interval probabilities by dose
  current_summary[25+nprov_dose+1,2] <- "interval probabilities by dose"
  current_summary[26+nprov_dose+1,5:7] <- c(paste0(category_name[1], ": (0, ", category_bound[1], "]"), 
                                            paste0(category_name[2], ": (", category_bound[1], ", ", category_bound[2], "]"), 
                                            paste0(category_name[3], ": (", category_bound[2], ", 1]"))
  
  current_summary[26+nprov_dose+1,3:4] <- paste("drug", 1:2)
  current_summary[(26+nprov_dose+1+1):(26+2*nprov_dose+2-1),3:4] <- prov_dose
  for(i in 1:nprov_dose){
    current_summary[26+nprov_dose+1+1+i-1,5:7] <- round(prob_posterior[,i], 5)
  }
  
  ## posterior parameter estimates
  current_summary[(25+2*nprov_dose+3),2] <- "posterior parameter estimates"
  current_summary[(26+2*nprov_dose+3+1):(26+2*nprov_dose+3+3),3] <- c("mean", "standard deviation", "correlation")
  current_summary[(26+2*nprov_dose+3),4:8] <- c("log(alpha1)", "log(beta1)", "log(alpha2)", "log(beta2)", "eta")
  current_summary[(26+2*nprov_dose+3+1),4:8] <- round(para_posterior$para_hat, 3)
  current_summary[(26+2*nprov_dose+3+2),4:8] <- round(para_posterior$para_sd, 3)
  current_summary[(26+2*nprov_dose+3+3),c(4,6)] <- round(c(para_posterior$para_corr1, para_posterior$para_corr2), 3)
  
  ## recommended next dose
  # dose combo
  current_summary[(25+2*nprov_dose+3+5),2] <- paste0("next dose: drug 1=", next_dose$dose[1], ", drug 2=", next_dose$dose[2])
  # interval probabilities by dose
  current_summary[(25+2*nprov_dose+3+5+1),3] <- "interval probabilities by dose"
  current_summary[(25+2*nprov_dose+3+5+1+1):(25+2*nprov_dose+3+5+3+1),4] <- c(paste0(category_name[1], ": (0, ", category_bound[1], "]"), 
                                                                              paste0(category_name[2], ": (", category_bound[1], ", ", category_bound[2], "]"), 
                                                                              paste0(category_name[3], ": (", category_bound[2], ", 1]"))
  current_summary[(25+2*nprov_dose+3+5+1+1):(25+2*nprov_dose+3+5+3+1),5] <- round(unlist(next_dose$posterior_prob_summary), 5) 
  # posterior DLT rate estimates
  current_summary[(25+2*nprov_dose+3+5+3+1+1),3] <- "posterior DLT rate estimates"
  current_summary[(25+2*nprov_dose+3+5+3+1+1+1):(25+2*nprov_dose+3+5+3+1+1+5),4] <- c("mean", "standard deviation", "2.5% quantile", "median", "97.5% quantile")
  current_summary[(25+2*nprov_dose+3+5+3+1+1+1):(25+2*nprov_dose+3+5+3+1+1+5),5] <- round(unlist(next_dose$posterior_pi_summary), 5)
  
  return(list(prob_posterior=prob_posterior,    # interval probabilities by dose
              para_posterior=para_posterior,    # posterior parameter estimates
              pi_posterior=pi_posterior,        # posterior DLT rate estimates
              next_dose=next_dose,              # recommended next dose
              current_summary=current_summary,  # summary of simulation outputs
              cohort_all=cohort_all))           # (accumulated) summary of each observed cohort
  }


# multi-scenarios
blrm_combo_ms <- function(prior, data, output_excel=FALSE, output_pdf=FALSE){
  
  # extract shared information among scenarios
  category_bound <- data[[1]]$category_bound
  category_name <- data[[1]]$category_name
  ewoc <- data[[1]]$ewoc
  burn_in <- data[[1]]$burn_in
  nsamples <- data[[1]]$nsamples
  
  drug1_name <- data[[1]]$drug1_name
  dose1_unit <- data[[1]]$dose1_unit
  ref_dose1 <- data[[1]]$ref_dose1
  
  drug2_name <- data[[1]]$drug2_name
  dose2_unit <- data[[1]]$dose2_unit
  ref_dose2 <- data[[1]]$ref_dose2
  
  # pre-allocate space
  all_prob_posterior <- all_para_posterior <- all_pi_posterior <- all_next_dose <- vector("list", length(data))
  simulation_summary <- NULL
  
  # save each scenario on a separate sheet of .xlsx file
  list_output <- vector("list", length(data))
  names(list_output) <- paste("Sheet", 1:length(data))
  
  # justify the maximum number of provisional doses: for purpose of simulation output
  nprov_max <- 12
  
  # Let us go through each scenario!
  for(k in 1:length(data)){
    
    # provisional doses
    prov_dose1 <- data[[k]]$prov_dose1
    prov_dose2 <- data[[k]]$prov_dose2
    nprov_dose1 <- length(prov_dose1)
    nprov_dose2 <- length(prov_dose2)
    prov_dose <- expand.grid(prov_dose1, prov_dose2)
    nprov_dose <- nrow(prov_dose)
    
    # tested doses
    dose1 <- data[[k]]$dose1
    dose2 <- data[[k]]$dose2
    dose <- cbind(dose1, dose2)
    ndose1 <- ndose2 <- ndose <- length(dose1) # = length(dose2) = nrow(dose)
    
    # observed cohorts: number of patients and number of DLTs
    n_pat <- data[[k]]$n_pat
    dlt <- data[[k]]$dlt
    
    # go for the specified scenario!!!
    scenario <- blrm_combo_ss(prior=prior, data=data[[k]])
    all_prob_posterior[[k]] <- scenario$prob_posterior
    all_para_posterior[[k]] <- scenario$para_posterior
    all_pi_posterior[[k]] <- scenario$pi_posterior
    all_next_dose[[k]] <- scenario$next_dose
    
    # construct the framework
    current_summary <- matrix(0, 25+2*nprov_dose+3+5+3+1+1+5, nprov_max+4)
    current_summary <- as.data.frame(current_summary)
    for(i in 1:nrow(current_summary)){
      for(j in 1:ncol(current_summary)){
        current_summary[i,j] <- ""
      }
    }
    
    # the kth scenario
    current_summary[1,1] <- paste0("Scenario", k)
    
    # Header
    current_summary[1,2] <- "Header"
    current_summary[2,3:5] <- c("simulation running date", "drug name", "RNG seeds")
    current_summary[3,3:5] <- c(paste(Sys.Date()), paste(c(drug1_name, drug2_name), collapse=", "), paste(seeds, collapse=", "))
    
    # Input
    current_summary[4,2] <- "User-specified Input"
    ## prior
    current_summary[5,2] <- "prior"
    current_summary[5,4:8] <- c("log(alpha1)", "log(beta1)", "log(alpha2)", "log(beta2)", "eta")
    current_summary[6,3] <- "mean"
    current_summary[7,3] <- "standard deviation"
    current_summary[8,3] <- "correlation"
    current_summary[6,4:8] <- c(prior[[1]]$mean, prior[[2]]$mean, prior[[3]]$mean)
    current_summary[7,4:8] <- c(prior[[1]]$se, prior[[2]]$se, prior[[3]]$se)
    current_summary[8,c(4,6)] <- c(prior[[1]]$corr, prior[[2]]$corr)
    
    ## data
    current_summary[9,2] <- "data"
    current_summary[10:22,3] <- c("burn-in", "number of simulation samples", "dose unit", 
                                  "provisional doses of drug 1", "provisional doses of drug 2", 
                                  "reference dose", "tested doses of drug 1", "tested doses of drug 2", 
                                  "number of patients", "number of DLTs", "category bounds", "category names", "EWOC")
    current_summary[10,4] <- burn_in
    current_summary[11,4] <- nsamples
    current_summary[12,4:5] <- c(paste0("drug1: ", dose1_unit), paste0("drug2: ", dose2_unit))
    current_summary[13,4:(4+nprov_dose1-1)] <- prov_dose1  # provisional doses of drug 1
    current_summary[14,4:(4+nprov_dose2-1)] <- prov_dose2  # provisional doses of drug 2
    current_summary[15,4:5] <- c(paste0("drug1: ", ref_dose1), (paste0("drug2: ", ref_dose2)))
    current_summary[16,4:(4+ndose1-1)] <- dose1
    current_summary[17,4:(4+ndose2-1)] <- dose2
    current_summary[18,4:(4+ndose-1)] <- n_pat
    current_summary[19,4:(4+ndose-1)] <- dlt
    current_summary[20,4:(4+length(category_bound)-1)] <- category_bound
    current_summary[21,4:(4+length(category_name)-1)] <- category_name
    current_summary[22,4] <- ewoc
    
    # output
    current_summary[23,2] <- "Simulation Output"
    
    ## posterior DLT rate estimates
    current_summary[24,3] <- "posterior DLT rate estimates"
    current_summary[25,6:10] <- c("mean", "standard deviation", "2.5% quantile", "median", "97.5% quantile")
    current_summary[25,4:5] <- paste("drug", 1:2)
    current_summary[26:(26+nprov_dose-1), 4:5] <- prov_dose
    for(i in 1:nprov_dose){
      current_summary[26+i-1,6:10] <- round(all_pi_posterior[[k]][,i], 5)
    }
    
    ## interval probabilities by dose
    current_summary[25+nprov_dose+1,3] <- "interval probabilities by dose"
    current_summary[26+nprov_dose+1,6:8] <- c(paste0(category_name[1], ": (0, ", category_bound[1], "]"), 
                                              paste0(category_name[2], ": (", category_bound[1], ", ", category_bound[2], "]"), 
                                              paste0(category_name[3], ": (", category_bound[2], ", 1]"))
    
    current_summary[26+nprov_dose+1,4:5] <- paste("drug", 1:2)
    current_summary[(26+nprov_dose+1+1):(26+2*nprov_dose+2-1),4:5] <- prov_dose
    for(i in 1:nprov_dose){
      current_summary[26+nprov_dose+1+1+i-1,6:8] <- round(all_prob_posterior[[k]][,i], 5)
    }
    
    ## posterior parameter estimates
    current_summary[(25+2*nprov_dose+3),3] <- "posterior parameter estimates"
    current_summary[(26+2*nprov_dose+3+1):(26+2*nprov_dose+3+3),4] <- c("mean", "standard deviation", "correlation")
    current_summary[(26+2*nprov_dose+3),5:9] <- c("log(alpha1)", "log(beta1)", "log(alpha2)", "log(beta2)", "eta")
    current_summary[(26+2*nprov_dose+3+1),5:9] <- round(all_para_posterior[[k]]$para_hat, 3)
    current_summary[(26+2*nprov_dose+3+2),5:9] <- round(all_para_posterior[[k]]$para_sd, 3)
    current_summary[(26+2*nprov_dose+3+3),c(5,7)] <- round(c(all_para_posterior[[k]]$para_corr1, all_para_posterior[[k]]$para_corr2), 3)
    
    ## recommended next dose
    # dose level
    current_summary[(25+2*nprov_dose+3+5),3] <- paste0("next dose: c(drug 1=", all_next_dose[[k]]$dose[1], ", drug 2=", all_next_dose[[k]]$dose[2], ")")
    # interval probabilities by dose
    current_summary[(25+2*nprov_dose+3+5+1),4] <- "interval probabilities by dose"
    current_summary[(25+2*nprov_dose+3+5+1+1):(25+2*nprov_dose+3+5+3+1),5] <- c(paste0(category_name[1], ": (0, ", category_bound[1], "]"), 
                                                                                paste0(category_name[2], ": (", category_bound[1], ", ", category_bound[2], "]"), 
                                                                                paste0(category_name[3], ": (", category_bound[2], ", 1]"))
    current_summary[(25+2*nprov_dose+3+5+1+1):(25+2*nprov_dose+3+5+3+1),6] <- round(unlist(all_next_dose[[k]]$posterior_prob_summary), 5)
    # posterior DLT rate estimates
    current_summary[(25+2*nprov_dose+3+5+3+1+1),4] <- "posterior DLT rate estimates"
    current_summary[(25+2*nprov_dose+3+5+3+1+1+1):(25+2*nprov_dose+3+5+3+1+1+5),5] <- c("mean", "standard deviation", "2.5% quantile", "median", "97.5% quantile")
    current_summary[(25+2*nprov_dose+3+5+3+1+1+1):(25+2*nprov_dose+3+5+3+1+1+5),6] <- round(unlist(all_next_dose[[k]]$posterior_pi_summary), 5)
    
    list_output[[paste("Sheet", k)]] <- current_summary # save each scenario on a separate sheet of the .xlsx file
    
    simulation_summary <- rbind(simulation_summary, current_summary) # combine each scenario iteratively since we need it as a returned argument
    
    cat("Scenario ", k, " is done.\n", sep="")
  }
  
  if(output_excel==TRUE){
     suppressMessages(openxlsx::write.xlsx(list_output, file=paste0(drug1_name, "-", drug2_name, "_BLRM_ms_", Sys.Date(), ".xlsx"), colNames=FALSE, rowNames=FALSE))
  }
  
  if(output_pdf==TRUE){
    
    # all the scenarios will be outputted in ONE file
    pdf(file=paste0(drug1_name, "-", drug2_name, "_BLRM_ms_", Sys.Date(), ".pdf"), width=10, height=10)
    
    for(k in 1:length(data)){
      
      # provisional doses
      prov_dose1 <- data[[k]]$prov_dose1
      prov_dose2 <- data[[k]]$prov_dose2
      prov_dose <- expand.grid(prov_dose1, prov_dose2)
      nprov_dose <- nrow(prov_dose)
      
      ###### Interval Probabilities by Dose: `(0.33, 1]' is the interval we take interest to ######
      prob_posterior_3 <- cbind(prov_dose, all_prob_posterior[[k]][3,]) # dim(prob_posterior_3): 3 * nrow(prov_dose)
      names(prob_posterior_3) <- c("drug1", "drug2", "probability")
      prob_posterior_3 <- transform(prob_posterior_3, level=ifelse(probability > ewoc, 1, 0))
      for(i in 1:nrow(prob_posterior_3)){
        if(as.numeric(rownames(prob_posterior_3[i,])) == as.numeric(all_next_dose[[k]]$index)){
          prob_posterior_3[i,4] <- 2
        }
      }
      prob_posterior_3_wide <- reshape2::dcast(prob_posterior_3[,-3], drug2 ~ drug1, value.var="level")
      prob_posterior_3_wide <- prob_posterior_3_wide[,-1]
      rownames(prob_posterior_3_wide) <- paste(prov_dose2)
      colnames(prob_posterior_3_wide) <- paste(prov_dose1)
      cols <- matrix(NA, nrow=length(prov_dose2), ncol=length(prov_dose1))
      for(i in 1:nrow(prob_posterior_3_wide)){
        for(j in 1:ncol(prob_posterior_3_wide)){
          if(prob_posterior_3_wide[i,j]==0){
            cols[i,j] <- "green"
          }else if(prob_posterior_3_wide[i,j]==1){
            cols[i,j] <- "red"
          }else{
            cols[i,j] <- "green" # "blue"
          }
        }
      }
      ## generate the plot
      plot(NA, NA, type='n', xaxt='n', yaxt='n', cex.lab=1.5, cex.main=2,
           xlim=range(1:length(prov_dose1)), 
           ylim=range(1:length(prov_dose2)),
           xlab=paste0(drug1_name, "(", dose1_unit, ")"), 
           ylab=paste0(drug2_name, "(", dose2_unit, ")"),
           main=paste0("Scenario ", k, ": ", "Dose combo Categorization"))
      abline(h=1:length(prov_dose2), v=1:length(prov_dose1), lty=2, lwd=1, col="gray")
      axis(1, at=1:length(prov_dose1), labels=paste(prov_dose1), las=0) # 1: bottom
      axis(2, at=1:length(prov_dose2), labels=paste(prov_dose2), las=2) # 2: left
      # add dose combos falling within different categories
      for(i in 1:length(prov_dose2)){
        for(j in 1:length(prov_dose1)){
          points(j, i, pch=19, col=cols[i,j], cex=4)
        }
      }
      legend("topright", c(" <= EWOC", " > EWOC"),# "Recommended Next Dose"), 
             cex=1.1, pch=rep(19, 2), col=c("green", "red"),# "blue"), 
             pt.cex=2, xpd=TRUE, horiz=TRUE, inset=c(0, -0.045), bty='n')
      
      ###### Posterior Distribution of DLT Rate ######
      par(mfrow=c(1,1))
      labels <- apply(prov_dose, 1, paste, collapse=",") # labels of x axis
      plot(1:nrow(prov_dose), all_pi_posterior[[k]][4,], type='p', 
           pch=20, xlab="drug combo", xaxt='n', cex.lab=1.5,
           ylab="DLT rate", ylim=c(0, max(all_pi_posterior[[k]])), 
           main=paste0("Scenario ", k, ": ", "Posterior Distribution of DLT Rate"), 
           cex.main=2.0, bty='n')
      axis(1, at=1:nrow(prov_dose), labels=FALSE)
      text(x=1:nrow(prov_dose), par("usr")[3]-0.03, labels=labels, srt=90, pos=1, xpd=TRUE, cex=0.5)
      # 95% credible interval
      arrows(1:nrow(prov_dose), all_pi_posterior[[k]][3,],
             1:nrow(prov_dose), all_pi_posterior[[k]][5,],
             code=3, angle=90, length=0.1, lwd=1.5, col=1) 
      if(max(all_pi_posterior[[k]][5,]) >= category_bound[2]){
        abline(h=category_bound, lty=2, col=c(rgb(0,1,0,alpha=0.8), rgb(1,0,0,alpha=0.8))) 
        legend("top", c(paste(category_bound), "median", "95% credible interval"), 
               lty=c(2,2,NA,1), lwd=c(1,1,NA,1.5), pch=c(NA,NA,20,NA), 
               col=c(rgb(0,1,0,alpha=0.8), rgb(1,0,0,alpha=0.8), 1, 1), 
               xpd=TRUE, horiz=TRUE, inset=c(0, -0.035), bty='n')
      }else if((max(all_pi_posterior[[k]][5,]) >= category_bound[1]) && (max(all_pi_posterior[[k]][5,]) < category_bound[2])){
        abline(h=category_bound[1], lty=2, col=rgb(0,1,0,alpha=0.8)) 
        legend("top", c(paste(category_bound[1]), "median", "95% credible interval"), 
               lty=c(2,NA,1), lwd=c(1,NA,1.5), pch=c(NA,20,NA), 
               col=c(rgb(0,1,0,alpha=0.8), 1, 1), bty='n',
               xpd=TRUE, horiz=TRUE, inset=c(0, -0.035))
      }else{
        legend("top", c("median", "95% credible interval"), 
               lty=c(NA,1), lwd=c(NA,1.5), pch=c(20,NA), col=c(1,1),
               xpd=TRUE, horiz=TRUE, inset=c(0, -0.035), bty='n')
      }
      #
    }
    dev.off()
  }
  
  return(list(prob_posterior=all_prob_posterior,
              para_posterior=all_para_posterior,
              pi_posterior=all_pi_posterior,
              next_dose=all_next_dose, 
              simulation_summary=simulation_summary))
}
