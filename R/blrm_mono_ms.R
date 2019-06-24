if(getRversion() >= "2.15.1") utils::globalVariables(c("Omega1", "inverse", "p2", "log.alpha", "dmnorm", "p1", "alpha", "d", "logit", "dbern", "s", "nb_pat", "logit<-"))

## the core function: ``blrm_mono_core"
blrm_mono_core <- function(prior, data){
  
  ################### generic functions to be used
  # BLRM for a single drug
  blrm_mono <- function(){  # p1: mean structure of the prior
    # p2: variance-covariance structure of the prior
    # log.alpha: log(alpha) & log(beta)
    # d: log(dose/ref_dose); should be a numeric vector of size nb_pat
    
    Omega1[1:2,1:2] <- inverse(p2[,])
    log.alpha[1:2] ~ dmnorm(p1[], Omega1[,])
    
    # alpha & beta: the original form of parameters
    alpha[1] <- exp(log.alpha[1])
    alpha[2] <- exp(log.alpha[2])
    
    for(i in 1:nb_pat){     # nb_pat: number of patients so far
        logit(pi1[i]) <- log.alpha[1] + alpha[2] * d[i]
        s[i] ~ dbern(pi1[i])
    }
  }
  
  # probability of toxicity
  tox_prob <- function(d, para){  # d: a numeric scalar
    # para: the original form of alpha & beta
    # Once we have (posterior) estimates of parameters, we are able to calculate probability of toxicity for each provisional dose
    
    logit_pi1 <- log(para[1]) + para[2] * d
    pi1 <- exp(logit_pi1)/(1 + exp(logit_pi1))
    return(list(pi1=pi1))
  }
  
  # prior for log(alpha) and log(beta)
  prior_summary <- function(prior){ # `prior` is a list composed of three elements, and 
    # each element is a numeric vector containing information 
    # for log(alpha) and log(beta), respectively, and their correlation
    mu <- prior[[1]]
    se <- prior[[2]]
    corr <- prior[[3]]
    covariance_matrix <- matrix(c(se[1]**2, se[1]*se[2]*corr, se[1]*se[2]*corr, se[2]**2), ncol=2, byrow=T)
    return(list(mean=mu, var=covariance_matrix))
  }
  
  # go through cumulated observed cohorts
  single_cohort <- function(prior_para){
    
    # data for JAGS
    mydata <- list(nb_pat=sum(cohort_size), 
                   s=data_dlt,    
                   d=data_sdose,
                   p1=prior_para$mean,
                   p2=prior_para$var)
    
    niters <- (1+burn_in)*nsamples/2 # number of iterations, e.g., (1+0.2)*10000, since 20% of 10000 will be chopped off
    
    # path.model <- file.path(tempdir(), "model.file.txt") # define model path
    # R2WinBUGS::write.model(blrm_mono, path.model)
    
    modelstring <- "model
    {
      Omega1[1:2, 1:2] <- inverse(p2[, ])
      log.alpha[1:2] ~ dmnorm(p1[], Omega1[, ])
      alpha[1] <- exp(log.alpha[1])
      alpha[2] <- exp(log.alpha[2])
      for (i in 1:nb_pat) {
         logit(pi1[i]) <- log.alpha[1] + alpha[2] * d[i]
         s[i] ~ dbern(pi1[i])
      }
    }"
    
    inits.list <- list(list(log.alpha=c(-3,0), .RNG.seed=seeds[1], .RNG.name="base::Wichmann-Hill"), # initialize the parameters in the log form
                       list(log.alpha=c(-3,0), .RNG.seed=seeds[2], .RNG.name="base::Wichmann-Hill"))
    # jagsobj <- rjags::jags.model(path.model, data=mydata, n.chains=2, quiet=TRUE, inits=inits.list)
    jagsobj <- rjags::jags.model(textConnection(modelstring), data=mydata, n.chains=2, quiet=TRUE, inits=inits.list)
    update(jagsobj, n.iter=niters, progress.bar="none")
    res <- rjags::jags.samples(jagsobj, "alpha", n.iter=niters, progress.bar="none") # `alpha`: alpha & beta
    # chop off the burn-in iterations, e.g., 20% of the targeted number of samples (e.g., 10000)
    # for each parameter; there are 2 chains generated in parallel
    alpha1 <- res$alpha[1,,][-c(1:(burn_in*nsamples/2)),]
    beta1 <- res$alpha[2,,][-c(1:(burn_in*nsamples/2)),]
    # the posterior samples of interest
    posterior_param <- cbind(c(alpha1), c(beta1)) 
    
    # (posterior) parameter estimates: log(alpha.hat) & log(beta.hat)
    log_para <- log(posterior_param)
    para_hat <- apply(log_para, 2, mean) # point estimates: alpha.hat and beta.hat
    para_sd <- apply(log_para, 2, sd)    # standard deviation of the point estimates
    para_corr <- corr(log_para)          # correlation coefficient between alpha.hat and beta.hat
    posterior_para_summary <- list(para_hat=para_hat, para_sd=para_sd, para_corr=para_corr)
    
    # collect posterior estimates for pi1, i.e., Pr(DLT | data): for each pair of posterior parameter estimates, 
    # we have a posterior estimate of DLT rate (toxic probability)
    samples_sdose <- matrix(0, nprov_dose, nsamples) # go through each standardized provisional dose
    for(i in 1:nprov_dose){
      for(j in 1:nsamples){
        samples_sdose[i,j] <- tox_prob(sprov_dose[i], posterior_param[j,])$pi1  
      }
    }
    
    # interval probabilities by dose: obtain the probability that Pr(DLT | data) lies in each category interval
    posterior_prob_summary <- matrix(0, length(category_bound)+1, nprov_dose)
    for(i in 1:nprov_dose){
      posterior_prob_summary[,i] <- as.numeric(table(cut(samples_sdose[i,], breaks=c(0, category_bound[1], category_bound[2], 1), right=TRUE))/nsamples)
    }
    
    # (posterior) pi1 estimates: mean & sd and 0.025, 0.5 & 0.975 quantile
    posterior_pi_summary <- matrix(0, 5, nprov_dose)
    for(i in 1:nprov_dose){
      posterior_pi_summary[,i] <- c(mean(samples_sdose[i,]), sd(samples_sdose[i,]), quantile(samples_sdose[i,], c(0.025, 0.5, 0.975)))
    }
    
    # next dose
    safe_dose_range <- (posterior_prob_summary[3,] <= ewoc)
    
    if(sum(safe_dose_range) != 0){
      # select the dose with the highest probability of targeted toxicity among safe doses
      next_dose_index <- which.max(posterior_prob_summary[2, safe_dose_range])
      next_dose_level <- prov_dose[next_dose_index] # dose level
      next_dose_posterior_prob_summary <- posterior_prob_summary[,next_dose_index] # interval probabilities by dose
      next_dose_posterior_pi_summary <- posterior_pi_summary[,next_dose_index]     # posterior summary of DLT rate
      next_dose <- list(dose=next_dose_level, 
                        posterior_prob_summary=list(next_dose_posterior_prob_summary), 
                        posterior_pi_summary=list(next_dose_posterior_pi_summary))
    }else next_dose <- NULL
    
    # return the posterior summary & the recommended next dose
    return(list(posterior_pi_summary=posterior_pi_summary,
                posterior_prob_summary=posterior_prob_summary,
                posterior_para_summary=posterior_para_summary, 
                next_dose=next_dose))
  }
  ################################################################
  
  # extract prior
  prior_para <- prior_summary(prior)    # prior for the parameters
  
  # extract data
  seeds <- data$seeds
  nsamples <- data$nsamples 
  burn_in <- data$burn_in
  drug_name <- data$drug_name
  prov_dose <- data$prov_dose
  ref_dose <- data$ref_dose
  dose_unit <- data$dose_unit
  #
  dose <- data$dose
  n_pat <- data$n_pat
  dlt_test <- data$dlt
  #
  category_bound <- data$category_bound
  category_name <- data$category_name
  ewoc <- data$ewoc
  
  # further process data
  sdose <- log(dose/ref_dose)           # log(standardized tested doses)
  ndose <- length(sdose)                # number of tested doses
  
  sprov_dose <- log(prov_dose/ref_dose) # log(standardized provisional doses)
  nprov_dose <- length(sprov_dose)      # number of provisional doses
  
  # go through all the tested doses cumulatively
  current_dose_index <- 1
  dose_index <- cohort_size <- dlt <- NULL
  data_sdose <- data_dlt <- NULL
  
  # save the results for each tested dose as a list
  cohort_all <- list()
  
  # go through each tested dose
  while(current_dose_index <= ndose){
    
    current_cohort_size <- n_pat[current_dose_index]
    current_dlt <- dlt_test[current_dose_index]
    
    # ... cumutatively
    dose_index <- c(dose_index, current_dose_index)
    cohort_size <- c(cohort_size, current_cohort_size)
    dlt <- c(dlt, current_dlt)
    
    # for each cohort, statistical analysis is based on all the cumulated observed cohorts information
    current_data_sdose <- rep(sdose[current_dose_index], current_cohort_size)
    data_sdose <- c(data_sdose, current_data_sdose)
    current_data_dlt <- ifelse(rep(current_dlt==0, current_cohort_size), rep(0, current_cohort_size), c(rep(1, current_dlt), rep(0, current_cohort_size-current_dlt)))
    data_dlt <- c(data_dlt, current_data_dlt)
    
    # go through the observed cohorts cumutatively
    cohort <- single_cohort(prior_para)
    cohort_all[[current_dose_index]] <- list(posterior_prob=cohort$posterior_prob_summary,
                                             posterior_para=cohort$posterior_para_summary,
                                             posterior_pi=cohort$posterior_pi_summary,
                                             next_dose=cohort$next_dose)
    
    # next dose cannot be determined
    if(is.null(cohort$next_dose)){
      break
    }else{
      # continue to the next tested dose
      current_dose_index <- current_dose_index + 1
    }
  }
  # done with going through all the tested doses
  
  # we only take interest to the final-round information
  cohort_final <- cohort_all[[length(cohort_all)]]
  
  # extract posterior summary
  prob_posterior <- cohort_final$posterior_prob
  pi_posterior <- cohort_final$posterior_pi
  para_posterior <- cohort_final$posterior_para
  next_dose <- cohort_final$next_dose
  
  # construct the framework for the simulation output
  current_summary <- matrix(0, 23+2*nprov_dose+3+5+3+1+1+5, nprov_dose+3)
  current_summary <- as.data.frame(current_summary)
  for(i in 1:nrow(current_summary)){
    for(j in 1:ncol(current_summary)){
        current_summary[i,j] <- ""
    }
  }
  
  # Header
  current_summary[1,1] <- "Header"
  current_summary[2,2:4] <- c("algorithm running date", "drug name", "RNG seeds")
  current_summary[3,2:4] <- c(paste(Sys.Date()), drug_name, paste(seeds, collapse=", "))
  
  # Input
  current_summary[4,1] <- "User-specified Input"
  ## prior
  current_summary[5,2] <- "prior"
  current_summary[5,4:5] <- c("log(alpha)", "log(beta)")
  current_summary[6,3] <- "mean"
  current_summary[7,3] <- "standard deviation"
  current_summary[8,3] <- "correlation"
  current_summary[6,4:5] <- prior$mean
  current_summary[7,4:5] <- prior$se
  current_summary[8,4] <- prior$corr
  
  ## data
  current_summary[9,2] <- "data"
  current_summary[10:20,3] <- c("number of simulated samples", "burn-in", "dose unit", "provisional doses", "reference dose", 
                                "tested doses", "number of patients", "number of DLTs", "category bounds", "category names", "EWOC")
  current_summary[10,4] <- nsamples
  current_summary[11,4] <- burn_in
  current_summary[12,4] <- dose_unit
  current_summary[13,4:(4+nprov_dose-1)] <- prov_dose
  current_summary[14,4] <- ref_dose
  current_summary[15,4:(4+ndose-1)] <- dose
  current_summary[16,4:(4+ndose-1)] <- n_pat
  current_summary[17,4:(4+ndose-1)] <- dlt
  current_summary[18,4:(4+length(category_bound)-1)] <- category_bound
  current_summary[19,4:(4+length(category_name)-1)] <- category_name
  current_summary[20,4] <- ewoc
  
  # output
  current_summary[21,1] <- "Simulation Output"
  
  ## posterior DLT rate estimates
  current_summary[22,2] <- "posterior DLT rate estimates"
  current_summary[23,4:8] <- c("mean", "standard deviation", "2.5% quantile", "median", "97.5% quantile")
  current_summary[23,3] <- "dose"
  current_summary[24:(24+nprov_dose-1), 3] <- prov_dose
  for(i in 1:nprov_dose){
    current_summary[24+i-1,4:8] <- round(pi_posterior[,i], 5)
  }
  
  ## interval probabilities by dose
  current_summary[23+nprov_dose+1,2] <- "interval probabilities by dose"
  current_summary[23+nprov_dose+2,4:6] <- c(paste0(category_name[1], ": (0, ", category_bound[1], "]"), 
                                            paste0(category_name[2], ": (", category_bound[1], ", ", category_bound[2], "]"), 
                                            paste0(category_name[3], ": (", category_bound[2], ", 1]"))
  current_summary[23+nprov_dose+2,3] <- "dose"
  current_summary[(23+nprov_dose+2+1):(23+2*nprov_dose+2+1-1),3] <- prov_dose
  for(i in 1:nprov_dose){
    current_summary[23+nprov_dose+2+i,4:6] <- round(prob_posterior[,i], 5)
  }
  
  ## posterior parameter estimates
  current_summary[(23+2*nprov_dose+3),2] <- "posterior parameter estimates"
  current_summary[(23+2*nprov_dose+3+1+1):(23+2*nprov_dose+3+1+2+1),3] <- c("mean", "standard deviation", "correlation")
  current_summary[(23+2*nprov_dose+3+1),4:5] <- c("log(alpha)", "log(beta)")
  current_summary[(23+2*nprov_dose+3+1+1),4:5] <- round(para_posterior$para_hat, 3)
  current_summary[(23+2*nprov_dose+3+2+1),4:5] <- round(para_posterior$para_sd, 3)
  current_summary[(23+2*nprov_dose+3+3+1),4] <- round(para_posterior$para_corr, 3)
  
  ## recommended next dose
  # dose level
  current_summary[(23+2*nprov_dose+3+5),2] <- paste0("next dose: ", next_dose$dose)
  
  # interval probabilities by dose
  current_summary[(23+2*nprov_dose+3+5+1),3] <- "interval probabilities by dose"
  current_summary[(23+2*nprov_dose+3+5+1+1):(23+2*nprov_dose+3+5+3+1),4] <- c(paste0(category_name[1], ": (0, ", category_bound[1], "]"), 
                                                                              paste0(category_name[2], ": (", category_bound[1], ", ", category_bound[2], "]"), 
                                                                              paste0(category_name[3], ": (", category_bound[2], ", 1]"))
  current_summary[(23+2*nprov_dose+3+5+1+1):(23+2*nprov_dose+3+5+3+1),5] <- round(unlist(next_dose$posterior_prob_summary), 5)
  
  # posterior DLT rate estimates
  current_summary[(23+2*nprov_dose+3+5+3+1+1),3] <- "posterior DLT rate estimates"
  current_summary[(23+2*nprov_dose+3+5+3+1+1+1):(23+2*nprov_dose+3+5+3+1+1+5),4] <- c("mean", "standard deviation", "2.5% quantile", "median", "97.5% quantile")
  current_summary[(23+2*nprov_dose+3+5+3+1+1+1):(23+2*nprov_dose+3+5+3+1+1+5),5] <- round(unlist(next_dose$posterior_pi_summary), 5)
  
  return(list(prob_posterior=prob_posterior,    # interval probabilities by dose
              para_posterior=para_posterior,    # posterior parameter estimates
              pi_posterior=pi_posterior,        # posterior DLT rate estimates
              next_dose=next_dose,              # recommended next dose
              current_summary=current_summary,  # summary of simulation outputs
              cohort_all=cohort_all))           # (accumulated) summary of each observed cohort
}

## multiple scenarios
blrm_mono_ms <- function(prior, data, output_excel=FALSE, output_pdf=FALSE){
  
     # extract common scenario information
     seeds <- data[[1]]$seeds
     drug_name <- data[[1]]$drug_name
     dose_unit <- data[[1]]$dose_unit
     prov_dose <- data[[1]]$prov_dose
     ref_dose <- data[[1]]$ref_dose
     category_bound <- data[[1]]$category_bound
     category_name <- data[[1]]$category_name
     ewoc <- data[[1]]$ewoc
     burn_in <- data[[1]]$burn_in
     nsamples <- data[[1]]$nsamples
     
     # obtain the number of provisional dose levels
     nprov_dose <- length(prov_dose)
     
     # # extract the version number of the algorithm: e.g., "_v1.R"
     # files <- list.files(path=".", pattern="\\v[0-9]+.R$")
     # vn <- regmatches(files, regexpr("\\v[0-9]+", files))[1]
     # 
     # # algorithm release date (year-month)
     # rd <- format(Sys.Date(), "%Y-%m")
     
     # pre-allocate space
     all_prob_posterior <- all_para_posterior <- all_pi_posterior <- all_next_dose <- all_cohort_all <- vector("list", length(data))
     simulation_summary <- NULL
     
     # save each scenario on a separate sheet of .xlsx file
     list_output <- vector("list", length(data))
     names(list_output) <- paste("Sheet", 1:length(data))
     
     for(k in 1:length(data)){
         
         # extract scenario specific information
         dose <- data[[k]]$dose
         n_pat <- data[[k]]$n_pat
         dlt <- data[[k]]$dlt
         
         # calculate the number of dose levels for each scenario
         ndose <- length(dose)
         
         # go for the specified scenario!!!
         scenario <- blrm_mono_core(prior=prior, data=data[[k]])  # for each scenario, data is different, while prior is the same
         all_prob_posterior[[k]] <- scenario$prob_posterior
         all_para_posterior[[k]] <- scenario$para_posterior
         all_pi_posterior[[k]] <- scenario$pi_posterior
         all_next_dose[[k]] <- scenario$next_dose
         all_cohort_all[[k]] <- scenario$cohort_all
         
         # construct the framework of simulation output
         current_summary <- matrix(0, 24+2*nprov_dose+3+5+1+3+1+5, nprov_dose+4)
         current_summary <- as.data.frame(current_summary)
         for(i in 1:nrow(current_summary)){
             for(j in 1:ncol(current_summary)){
                 current_summary[i,j] <- ""
             }
         }
         
         # specified scenario
         current_summary[1,1] <- paste0("Scenario ", k)
         # header
         current_summary[2,2] <- "Header"
         current_summary[3,3:5] <- c("algorithm running date", "drug name", "RNG seeds")
         current_summary[4,3:5] <- c(paste(Sys.Date()), drug_name, paste(seeds, collapse=", "))
         
         # input
         current_summary[5,2] <- "User-specified Input"
         ## prior
         current_summary[6,3] <- "prior"
         current_summary[6,5:6] <- c("log(alpha)", "log(beta)")
         current_summary[7:9,4] <- c("mean", "standard deviation", "correlation")
         current_summary[7,5:6] <- prior$mean
         current_summary[8,5:6] <- prior$se
         current_summary[9,5] <- prior$corr
         
         ## data
         current_summary[10,3] <- "data"
         current_summary[11:21,4] <- c("number of simulated samples", "burn-in", "dose unit", "provisional doses", "reference dose", 
                                       "tested doses", "number of patients", "number of DLTs", "category bounds", "category names", "EWOC")
         current_summary[11,5] <- nsamples
         current_summary[12,5] <- burn_in
         current_summary[13,5] <- dose_unit
         current_summary[14,5:(5+nprov_dose-1)] <- prov_dose
         current_summary[15,5] <- ref_dose
         current_summary[16,5:(5+ndose-1)] <- dose
         current_summary[17,5:(5+ndose-1)] <- n_pat
         current_summary[18,5:(5+ndose-1)] <- dlt
         current_summary[19,5:(5+length(category_bound)-1)] <- category_bound
         current_summary[20,5:(5+length(category_name)-1)] <- category_name
         current_summary[21,5] <- ewoc
         
         # output
         current_summary[22,2] <- "Simulation Output"
         
         ## posterior DLT rate estimates
         current_summary[23,3] <- "posterior DLT rate estimates"
         current_summary[24,5:9] <- c("mean", "sd", "2.5% quantile", "median", "97.5% quantile")
         current_summary[24,4] <- "dose"
         current_summary[25:(25+nprov_dose-1),4] <- prov_dose
         for(i in 1:nprov_dose){
             current_summary[25+i-1,5:9] <- round(all_pi_posterior[[k]][,i], 5)
         }
         
         ## interval probabilities by dose
         current_summary[24+nprov_dose+1,3] <- "Interval probabilities by dose"
         current_summary[24+nprov_dose+1+1,5:7] <- c(paste0(category_name[1], ": (0, ", category_bound[1], "]"), 
                                                   paste0(category_name[2], ": (", category_bound[1], ", ", category_bound[2], "]"), 
                                                   paste0(category_name[3], ": (", category_bound[2], ", 1]"))
         current_summary[24+nprov_dose+1+1,4] <- "dose"
         current_summary[(24+nprov_dose+1+1+1):(24+2*nprov_dose+1+1+1-1),4] <- prov_dose
         for(i in 1:nprov_dose){
             current_summary[24+nprov_dose+3+i-1,5:7] <- round(all_prob_posterior[[k]][,i], 5)
         }
         
         ## posterior parameter estimates
         current_summary[(24+2*nprov_dose+3+1-1),3] <- "posterior parameter estimates"
         current_summary[(24+2*nprov_dose+3+1+1):(24+2*nprov_dose+3+4),4] <- c("mean", "standard deviation", "correlation")
         current_summary[(24+2*nprov_dose+3+1),5:6] <- c("log(alpha)", "log(beta)")
         current_summary[(24+2*nprov_dose+3+1+1),5:6] <- round(all_para_posterior[[k]]$para_hat, 3)
         current_summary[(24+2*nprov_dose+3+2+1),5:6] <- round(all_para_posterior[[k]]$para_sd, 3)
         current_summary[(24+2*nprov_dose+3+3+1),5] <- round(all_para_posterior[[k]]$para_corr,3)
         
         ## recommended next dose
         # dose level
         current_summary[(24+2*nprov_dose+3+5),3] <- paste0("Next dose: ", all_next_dose[[k]]$dose)
         # interval probabilities by dose
         current_summary[(24+2*nprov_dose+3+5+1),4] <- "interval probabilities by dose"
         current_summary[(24+2*nprov_dose+3+5+1+1):(24+2*nprov_dose+3+5+1+3),5] <- c(paste0(category_name[1], ": (0, ", category_bound[1], "]"), 
                                                                                     paste0(category_name[2], ": (", category_bound[1], ", ", category_bound[2], "]"), 
                                                                                     paste0(category_name[3], ": (", category_bound[2], ", 1]"))
         current_summary[(24+2*nprov_dose+3+5+1+1):(24+2*nprov_dose+3+5+1+3),6] <- round(as.numeric(unlist(all_next_dose[[k]]$posterior_prob_summary)), 5)
         # posterior DLT rate estimates
         current_summary[(24+2*nprov_dose+3+5+1+3+1),4] <- "posterior DLT rate estimates"
         current_summary[(24+2*nprov_dose+3+5+1+3+1+1):(24+2*nprov_dose+3+5+1+3+1+5),5] <- c("mean", "sd", "2.5% quantile", "median", "97.5% quantile")
         current_summary[(24+2*nprov_dose+3+5+1+3+1+1):(24+2*nprov_dose+3+5+1+3+1+5),6] <- round(as.numeric(unlist(all_next_dose[[k]]$posterior_pi_summary)), 5)
         
         # save each scenario on a separate sheet of the .xlsx file
         list_output[[paste("Sheet", k)]] <- current_summary
         
         # combine each scenario iteratively since we need it as a returned argument
         simulation_summary <- rbind(simulation_summary, current_summary)
         
         # progress meter
         cat("Scenario ", k, "is done.\n", sep="")
     }
     
     if(output_excel==TRUE){
        suppressMessages(openxlsx::write.xlsx(list_output, file=paste0(drug_name, "_BLRM_ms_", Sys.Date(), ".xlsx"), 
                                              colNames=FALSE, rowNames=FALSE))
     }
     
     if(output_pdf==TRUE){
        
        # all the scenarios will be printed out in ONE file
        pdf(file=paste0(drug_name, "_BLRM_ms_", Sys.Date(), ".pdf"), width=8, height=8)
        for(k in 1:length(data)){
            
            # extract scenario specific information
            dose <- data[[k]]$dose
            n_pat <- data[[k]]$n_pat
            dlt <- data[[k]]$dlt
            ndose <- length(dose)
            
            # interval probabilities by dose
            cols <- c("green", "red")[(all_prob_posterior[[k]][3,] > ewoc)+1] # red: stop; green: pass
            yrange <- function(max.y) yrange.final <- ifelse(rep((1.1*max.y <= 1), 2), c(0, 1.1*max.y), c(0, max.y)) # max.y: the maximum of the response
            layout(matrix(c(1,2,3), 3, 1, byrow=TRUE))
            ## e.g., (0.33, 1]
            barplot(all_prob_posterior[[k]][3,], xlab=paste0("dose", "(", dose_unit, ")"), ylab="probability",
                    ylim=yrange(max(all_prob_posterior[[k]][3,])), names.arg=paste(prov_dose), col=cols, 
                    main=paste0(category_name[3], ": (", category_bound[2], ", 1]"), cex.main=1.3, font.main=4)
            if(max(all_prob_posterior[[k]][3,]) >= ewoc){
               abline(h=ewoc, lty=2, col="red")
               text(x=0.45, y=1.15*ewoc, labels=paste0("EWOC=", ewoc), col="red", cex=1.2, font.main=4)
            }else{
               text(x=0.9, y=max(all_prob_posterior[[k]][3,]), labels=paste0("EWOC=", ewoc), col="red", cex=1.2, font.main=4)
            }
            ## e.g., (0.16, 0.33]
            barplot(all_prob_posterior[[k]][2,], xlab=paste0("dose", "(", dose_unit, ")"), ylab="probability", 
                    ylim=yrange(max(all_prob_posterior[[k]][2,])), names.arg=paste(prov_dose), col="green", 
                    main=paste0(category_name[2], ": (", category_bound[1], ", ", category_bound[2], "]"), cex.main=1.3, font.main=4)
            ## e.g., (0, 0.16]
            barplot(all_prob_posterior[[k]][1,], xlab=paste0("dose", "(", dose_unit, ")"), ylab="probability",  
                    ylim=yrange(max(all_prob_posterior[[k]][1,])), names.arg=paste(prov_dose), col="green", 
                    main=paste0(category_name[1], ": (0, ", category_bound[1], "]"), cex.main=1.3, font.main=4)
            ## add a main title to the three barplots together
            mtext(paste0("Scenario", k, ": Interval Probabilities by Dose"), side=3, outer=TRUE, 
                  line=-1.5, at=par("usr")[1]+0.039*diff(par("usr")[1:2]), cex=1, font=2)
            
            # posterior distribution of DLT rate
            par(mfrow=c(1,1))
            plot(prov_dose, all_pi_posterior[[k]][4,], type="p", pch=20, xlab=paste0("dose", "(", dose_unit, ")"), ylab="DLT rate", 
                 xlim=range(prov_dose), ylim=c(0, max(all_pi_posterior[[k]])), main=paste0("Scenario", k, ": Posterior Distribution of DLT Rate"), bty="n")
            arrows(prov_dose, all_pi_posterior[[k]][3,], prov_dose, all_pi_posterior[[k]][5,], code=3, angle=90, length=0.1, lwd=1.5, col=1) # 95% credible interval
            if(max(all_pi_posterior[[k]][5,]) >= category_bound[2]){
               abline(h=category_bound, lty=2, col=c(rgb(0,1,0,alpha=0.8), rgb(1,0,0,alpha=0.8))) 
               legend("topleft", c(paste(category_bound), "median", "95% credible interval"), 
                      lty=c(2,2,NA,1), lwd=c(1,1,NA,1.5), pch=c(NA,NA,20,NA), col=c(rgb(0,1,0,alpha=0.8), rgb(1,0,0,alpha=0.8), 1, 1), bty="n")
            }else if((max(all_pi_posterior[[k]][5,]) >= category_bound[1]) && (max(all_pi_posterior[[k]][5,]) < category_bound[2])){
              abline(h=category_bound[1], lty=2, col=rgb(0,1,0,alpha=0.8)) 
              legend("topleft", c(paste(category_bound[1]), "median", "95% credible interval"), 
                     lty=c(2,NA,1), lwd=c(1,NA,1.5), pch=c(NA,20,NA), col=c(rgb(0,1,0,alpha=0.8), 1, 1), bty="n")
            }else{
              legend("topleft", c("median", "95% credible interval"), 
                     lty=c(NA,1), lwd=c(NA,1.5), pch=c(20,NA), col=c(1,1), bty="n")
            }
        }
        dev.off()
     }
     
     return(list(prob_posterior=all_prob_posterior,
                 para_posterior=all_para_posterior,
                 pi_posterior=all_pi_posterior,
                 next_dose=all_next_dose,
                 cohort_all=all_cohort_all,
                 simulation_summary=simulation_summary))
}
