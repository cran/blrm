blrm_mono_ss <- function(prior, data, output_excel=FALSE, output_pdf=FALSE){
  
  ################### generic functions to be used ###################
  # probability of toxicity
  tox_prob <- function(d, para){  
    # d: a numeric scalar
    # para: the original form of alpha & beta
    # Once we have (posterior) estimates of parameters, we are able to calculate probability of toxicity for each provisional dose
    
    logit_pi1 <- log(para[1]) + para[2] * d
    pi1 <- exp(logit_pi1)/(1 + exp(logit_pi1))
    return(list(pi1=pi1))
  }
  
  # prior for the log(alpha) and log(beta)
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
    
    # justify the model for library(rjags)
    modelstring <- "model
    {
    Omega1[1:2, 1:2] <- inverse(p2[, ])
    log.alpha[1:2] ~ dmnorm(p1[], Omega1[, ])
    alpha[1] <- exp(log.alpha[1])
    alpha[2] <- exp(log.alpha[2])
    for(i in 1:nb_pat){
    logit(pi1[i]) <- log.alpha[1] + alpha[2] * d[i]
    s[i] ~ dbern(pi1[i])
    }
    }"
    # initialize the parameters in the log form, and set up the RNG seeds
    inits.list <- list(list(log.alpha=c(-3,0), .RNG.seed=seeds[1], .RNG.name="base::Wichmann-Hill"),
                       list(log.alpha=c(-3,0), .RNG.seed=seeds[2], .RNG.name="base::Wichmann-Hill"))
    jagsobj <- rjags::jags.model(textConnection(modelstring), data=mydata, n.chains=2, quiet=TRUE, inits=inits.list)
    update(jagsobj, n.iter=niters, progress.bar="none")
    res <- rjags::jags.samples(jagsobj, "alpha", n.iter=niters, progress.bar="none") # `alpha`: alpha & beta
    # chop off the burn-in iterations, e.g., 20% of the target number of samples (e.g., 10000)
    # for each parameter; there are 2 chains generated in parallel, and thus, 10% of each will be chopped off
    alpha1 <- res$alpha[1,,][-c(1:(burn_in*nsamples/2)),]
    beta1 <- res$alpha[2,,][-c(1:(burn_in*nsamples/2)),]
    # the posterior samples of interest
    posterior_param <- cbind(c(alpha1), c(beta1)) 
    
    # (posterior) parameter estimates: log(alpha.hat) & log(beta.hat)
    log_para <- log(posterior_param)
    para_hat <- apply(log_para, 2, mean) # point estimates: alpha.hat and beta.hat
    para_sd <- apply(log_para, 2, sd)    # standard deviation of the point estimates
    para_corr <- boot::corr(log_para)          # correlation coefficient between alpha.hat and beta.hat
    posterior_para_summary <- list(para_hat=para_hat, para_sd=para_sd, para_corr=para_corr)
    
    # collect posterior estimates for pi1, i.e., Pr(DLT | data): for each pair of 
    # posterior parameter estimates, we have a posterior estimate of DLT rate (probability of toxicity)
    samples_sdose <- matrix(NA, nrow=nprov_dose, ncol=nsamples) # go through each standardized provisional dose
    for(i in 1:nprov_dose){
      for(j in 1:nsamples){
        samples_sdose[i,j] <- tox_prob(sprov_dose[i], posterior_param[j,])$pi1  
      }
    }
    
    # interval probabilities by dose: obtain the probability that Pr(DLT | data) lies in each category interval
    posterior_prob_summary <- matrix(NA, nrow=length(category_bound)+1, ncol=nprov_dose)
    for(i in 1:nprov_dose){
      posterior_prob_summary[,i] <- as.numeric(table(cut(samples_sdose[i,], breaks=c(0, category_bound[1], category_bound[2], 1), right=TRUE))/nsamples)
    }
    
    # (posterior) pi1 estimates: mean & sd and 0.025, 0.5 & 0.975 quantile
    posterior_pi_summary <- matrix(0, 5, nprov_dose)
    for(i in 1:nprov_dose){
      posterior_pi_summary[,i] <- c(mean(samples_sdose[i,]), sd(samples_sdose[i,]), quantile(samples_sdose[i,], c(0.025, 0.5, 0.975)))
    }
    
    # the next dose: select the dose with the highest probability of targeted-toxicity among safe doses
    ## justify safe doses first
    safe_dose_range <- (posterior_prob_summary[3,] <= ewoc) # posterior_prob_summary[3,]: (0.33, 1]
    if(sum(safe_dose_range) != 0){
      
      ## subset safe doses
      safe_dose <- prov_dose[safe_dose_range]
      safe_dose_prob <- posterior_prob_summary[,safe_dose_range]
      safe_dose_pi <- posterior_pi_summary[,safe_dose_range]
      
      # justify the dose with the highest probability of targeted toxicity
      next_dose_index <- which.max(safe_dose_prob[2,])  
      next_dose_level <- safe_dose[next_dose_index]                        # the recommended dose level
      next_dose_posterior_prob_summary <- safe_dose_prob[,next_dose_index] # interval probabilities by dose
      next_dose_posterior_pi_summary <- safe_dose_pi[,next_dose_index]     # posterior summary of DLT rate
      
      # combine them to a list
      next_dose <- list(dose=next_dose_level, 
                        posterior_prob_summary=list(next_dose_posterior_prob_summary), 
                        posterior_pi_summary=list(next_dose_posterior_pi_summary))     
    }else{
      next_dose <- NULL
    }
    
    # return the posterior summary & the recommended next dose
    return(list(posterior_pi_summary=posterior_pi_summary,
                posterior_prob_summary=posterior_prob_summary,
                posterior_para_summary=posterior_para_summary, 
                next_dose=next_dose))
  }
  ## end of single cohort
  
  # extract prior
  prior_para <- prior_summary(prior)
  
  # extract parameters
  seeds <- data$seeds
  nsamples <- data$nsamples 
  burn_in <- data$burn_in
  category_bound <- data$category_bound
  category_name <- data$category_name
  ewoc <- data$ewoc
  
  # extract drug information
  drug_name <- data$drug_name
  prov_dose <- data$prov_dose
  # order provisional doses from low to high (for purpose of index tracking)
  prov_dose <- sort(prov_dose, decreasing=FALSE)
  ref_dose <- data$ref_dose
  dose_unit <- data$dose_unit
  
  # extract observed cohorts
  dose <- data$dose
  n_pat <- data$n_pat
  dlt_test <- data$dlt
  
  # further process drug-related information
  sdose <- log(dose/ref_dose)           # log(scaled tested doses)
  ndose <- length(sdose)                # number of tested doses
  
  sprov_dose <- log(prov_dose/ref_dose) # log(scaled provisional doses)
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
    
    # cumulatively
    dose_index <- c(dose_index, current_dose_index)
    cohort_size <- c(cohort_size, current_cohort_size)
    dlt <- c(dlt, current_dlt)
    
    # for each cohort, statistical analysis is based on all the cumulated observed cohorts information
    current_data_sdose <- rep(sdose[current_dose_index], current_cohort_size)
    data_sdose <- c(data_sdose, current_data_sdose)
    current_data_dlt <- ifelse(rep(current_dlt==0, current_cohort_size), rep(0, current_cohort_size), c(rep(1, current_dlt), rep(0, current_cohort_size-current_dlt)))
    data_dlt <- c(data_dlt, current_data_dlt)
    
    # go through the observed cohorts cumulatively
    cohort <- single_cohort(prior_para)
    cohort_all[[current_dose_index]] <- list(posterior_prob=cohort$posterior_prob_summary,
                                             posterior_para=cohort$posterior_para_summary,
                                             posterior_pi=cohort$posterior_pi_summary,
                                             next_dose=cohort$next_dose)
    
    if(is.null(cohort$next_dose)){ # no recommended next dose
      if(current_dose_index < ndose){
        stop("There is no RND even not all the tested doses are gone through!")
      }else{  # current_dose_index == ndose
        break
      }
    }else{ # continue to the next tested dose
      current_dose_index <- current_dose_index + 1
    }
  }
  # done with going through all the tested doses
  
  # we only take interest to the final-round information
  cohort_final <- cohort_all[[length(cohort_all)]]
  
  # extract posterior summary and next dose
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
  current_summary[2,2:4] <- c("simulation running date", "drug name", "RNG seeds")
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
  
  if(output_excel==TRUE){  
    
    # write the simulation output to .xlsx
    output_1 <- current_summary[1:20,]
    output_2 <- current_summary[21:(23+2*nprov_dose+3+5+3+1+1+5),]
    list_output <- list("Sheet 1"=output_1, "Sheet 2"=output_2)
    suppressMessages(openxlsx::write.xlsx(list_output, file=paste0(drug_name, "_BLRM_ss_", Sys.Date(), ".xlsx"), 
                                          colNames=FALSE, rowNames=FALSE))
  }
  
  if(output_pdf==TRUE){
    
    pdf(file=paste0(drug_name, "_BLRM_ss_", Sys.Date(), ".pdf"), width=8, height=8)
    
    # interval probabilities by dose
    cols <- c("green", "red")[((prob_posterior[3,]) > ewoc)+1]
    yrange <- function(max.y){
      yrange.final <- ifelse(rep((1.1*max.y <= 1), 2), c(0, 1.1*max.y), c(0, max.y)) # max.y: the maximum of the responses
    }
    layout(matrix(c(1,2,3), 3, 1, byrow=TRUE))
    ## e.g., (0.33, 1]
    barplot(prob_posterior[3,], xlab=paste0("dose", "(", dose_unit, ")"), ylab="probability",
            ylim=yrange(max(prob_posterior[3,])), names.arg=paste(prov_dose), col=cols,
            main=paste0(category_name[3], ": (", category_bound[2], ", 1]"), 
            cex.main=1.5, font.main=4, cex.axis=1.5, cex.names=1.5, cex.lab=1.5)
    if(max(prob_posterior[3,]) >= ewoc){
      abline(h=ewoc, lty=3, col="black")
      text(x=0.4, y=1.15*ewoc, labels=paste0("EWOC=", ewoc), col="red", cex=1.5, font.main=4)
    }else{
      text(x=0.8, y=max(prob_posterior[3,]), labels=paste0("EWOC=", ewoc), col="red", cex=1.5, font.main=4)
    }
    ## e.g., (0.16, 0.33]
    barplot(prob_posterior[2,], xlab=paste0("dose", "(", dose_unit, ")"), ylab="probability",
            ylim=yrange(max(prob_posterior[2,])), names.arg=paste(prov_dose), col="green",
            main=paste0(category_name[2], ": (", category_bound[1], ", ", category_bound[2], "]"), 
            cex.main=1.5, font.main=4, cex.axis=1.5, cex.names=1.5, cex.lab=1.5)
    ## e.g., (0, 0.16]
    barplot(prob_posterior[1,], xlab=paste0("dose", "(", dose_unit, ")"), ylab="probability",
            ylim=yrange(max(prob_posterior[1,])), names.arg=paste(prov_dose), col="green",
            main=paste0(category_name[1], ": (0, ", category_bound[1], "]"), 
            cex.main=1.5, font.main=4, cex.axis=1.5, cex.names=1.5, cex.lab=1.5)
    ## add a main title to the three barplots together
    mtext("Interval Probabilities by Dose", side=3, outer=TRUE, line=-2,
          at=par("usr")[1]+0.035*diff(par("usr")[1:2]), cex=1.2, font=2)
    
    # posterior distribution of DLT rate
    par(mfrow=c(1,1))
    plot(prov_dose, pi_posterior[4,], type='p', pch=20, xlab=paste0("dose", "(", dose_unit, ")"), ylab="DLT rate", cex.lab=1.5,
         xlim=range(prov_dose), ylim=c(0, max(pi_posterior)), main="Posterior Distribution of DLT Rate", bty='n') 
    arrows(prov_dose, pi_posterior[3,], prov_dose, pi_posterior[5,], code=3, angle=90, length=0.1, lwd=1.5, col=1) # 95% credible interval
    if(max(pi_posterior[5,]) >= category_bound[2]){
      abline(h=category_bound, lty=2, col=c(rgb(0,1,0,alpha=0.8), rgb(1,0,0,alpha=0.8)))
      legend("topleft", c(paste(category_bound), "median", "95% credible interval"),
             lty=c(2,2,NA,1), lwd=c(1,1,NA,1.5), pch=c(NA,NA,20,NA),
             col=c(rgb(0,1,0,alpha=0.8), rgb(1,0,0,alpha=0.8), 1, 1), bty='n')
    }else if((max(pi_posterior[5,]) >= category_bound[1]) && (max(pi_posterior[5,]) < category_bound[2])){
      abline(h=category_bound[1], lty=2, col=rgb(0,1,0,alpha=0.8))
      legend("topleft", c(paste(category_bound[1]), "median", "95% credible interval"),
             lty=c(2,NA,1), lwd=c(1,NA,1.5), pch=c(NA,20,NA), col=c(rgb(0,1,0,alpha=0.8), 1, 1), bty='n')
    }else{
      legend("topleft", c("median", "95% credible interval"), lty=c(NA,1), lwd=c(NA,1.5), pch=c(20,NA), col=c(1,1), bty='n')
    }
    
    dev.off()
  }
  
  return(list(prob_posterior=prob_posterior,    # interval probabilities by dose
              para_posterior=para_posterior,    # posterior parameter estimates
              pi_posterior=pi_posterior,        # posterior DLT rate estimates
              next_dose=next_dose,              # recommended next dose
              current_summary=current_summary,  # summary of simulation outputs
              cohort_all=cohort_all))           # (accumulated) summary of each observed cohort
}
