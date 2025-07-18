## Edit DMC 20240215: get_sims has new parameter 'tmax' for maximum time used to compute AUC.
## Defaults to max of tseq2 (i.e. final point on plot), however in general should be set so that S(tmax) ~=0
## Tried to edit get_auc to calculate the upper limit of the integral as the 0.9999 quantile of the relevant distribution
## This does not work since some distributions plateau!


library(flexsurv)
library(mvtnorm)
library(Hmisc)
library(GGally)
library(tidyverse)

## Need to update to take multiple timepoints
## get_sims and is_surv_viz_gg need to be updated at a minimum
## Not sure about is_surv and IS 

## Might also be useful to have a function that calculates ex_info automatically given tstar, mu, sd etc.


##---
## Overall structure ----
##---

## main functions

## -- ex_info --
## defines loss function for external
## output includes location of timepoints, loss function and associated parameters 
## (e.g., dnorm, dmvnorm; mu.t, sigma.t, etc, ...)
## for now we will just do this manually

## -- surv_fit -- 
## likelihood term for loss function 
## e.g., flexsurvreg object, or output from this -- just going to use this for now
##  needs dlist$inv.transforms, surv_fit$dfns$p, surv_fit$coefficients, surv_fit$cov as inputs
## N.B., dlist$inv.transforms and dfns$p can be defined in advance 

##-- is_surv --
## this takes external information and fitted model and returns outputs from modified survival function
## internally calls IS and get_alpha 


## Edit DMC: this is no longer working; replaced with is_surv_viz_gg
##-- is_surv_viz --
## visualises output from is_surv
## diagnostic plots to assess alpha; plot parameters old vs new; plot survival curves

## internal functions -- i.e., called internally by main functions 

## transforms_fun --
## s_fun --
## IS --
## get_alpha --
## N.B., get_alpha calls IS

## New functions DMC:
## get_sims: extract parameter simulations and corresponding survival curves and AUC estimates
## get_auc: calculates AUC
## inv_trans: implements inverse transformations


##---
## Define functions ----
##---





# The code below creates named lists of cumulative density functions, parameter transforms and inverse transforms used for each distribution

# List of standard distributions
sdists<-c("exp","weibull","gompertz","lnorm",
                   "llogis","gengamma")

inv_transforms<-list()
transforms<-list()
d_funs<-list()

for(d in sdists)
{
  inv_transforms[[d]]<-flexsurv::flexsurv.dists[[d]][["inv.transforms"]]
  transforms[[d]]<-flexsurv::flexsurv.dists[[d]][["transforms"]]
  d_funs[[d]]<- get(paste0("p",d))
}

# Longer distribution names for compatibility
dists<-c("exponential","weibull","gompertz","lognormal",
         "llogis","gengamma")

# Full names of distributions (named vector) for printing
distributions<-c("exponential"="Exponential",
                 "weibull"="Weibull",
                 "lognormal"="Log-normal",
                 "llogis"="Log-logistic",
                 "gompertz"="Gompertz",
                 "gengamma"="Gen. Gamma")

names(inv_transforms)<-dists
names(transforms)<-dists
names(d_funs)<-dists


## s_fun
## Returns survival at time t for distribution dist given parameters trans_pars

s_fun <- function(t, dist, trans_pars){
  
  if(NCOL(trans_pars) == 1){
    ret <- 1 - d_funs[[dist]](t, trans_pars)  
  } 
  if(NCOL(trans_pars) == 2){
    ret <- 1 - d_funs[[dist]](t, trans_pars[,1], trans_pars[,2])  
  } 
  if(NCOL(trans_pars) == 3){
    ret <- 1 - d_funs[[dist]](t, trans_pars[,1], trans_pars[,2], trans_pars[,3])  
  }
  return(ret)
} 


## Function for importance sampling 
## Need to set alpha and run this repeatedly
## Inputs: MLE parameters, distribution, and prior
## Outputs: importance sampling diagnostics (used to choose alpha), posterior parameter estimates, covariance and samples




IS <- function(alpha = 10, 
               #parameter estimates from flexsurv
               coeff,
               # covariance matrix
               cov,
               # transformed to MVN scale already?
               transformed=TRUE,
               # Distribution: one of ("exponential","weibull","gompertz","lognormal","llogis","gengamma")
               dist,
               # Prior/loss function incorporating external information 
               ex_info,
               iter = 5000){
  
  if (!transformed)
  {
    # Need some code here that will do the transformations
    stop("Sorry, this functionality is not available yet. Please transform the parameters to the MVN scale.")
  }
  
 
  tstar. <- ex_info$tstar
  
  ## Generate samples  - re-weight cov matrix using alpha. 
  q_samples <- rmvnorm(iter, mean = coeff, sigma = alpha * cov) 
  
  ## For each sample compute survival at time(s) tstar
  
  if(length(tstar.) == 1){
    if(ncol(q_samples) == 1){
      
      q_samples_trans <- inv_transforms[[1]][[1]](q_samples)
      
      q_S_tstar <- s_fun(q_samples_trans, dist=dist, t = tstar.)
    } else{
      ## transform parameter samples to natural scales to compute survival
      q_samples_trans <- q_samples 
      for(nc in 1:ncol(q_samples)){
        
        q_samples_trans[, nc] <- inv_transforms[[dist]][[nc]](q_samples[, nc])
      }
      ## evaluate survival 
      
      q_S_tstar <- s_fun(q_samples_trans, dist=dist, t = tstar.)
    }
  } else{
    len_t <- length(tstar.)
    q_S_tstar <- matrix(0, nrow = nrow(q_samples), ncol = len_t)
    if(ncol(q_samples) == 1){
   
      
      q_samples_trans <- inv_transforms[[1]][[1]](q_samples)
      for(t in 1:len_t){
        
   
        q_S_tstar[, t] <- s_fun(q_samples_trans, dist=dist, t = tstar.[t])  
      }
      
    } else{
      ## transform parameter samples to natural scales to compute survival
      q_samples_trans <- q_samples 
      for(nc in 1:ncol(q_samples)){
        
   
        q_samples_trans[, nc] <- inv_transforms[[dist]][[nc]](q_samples[, nc])
      }
      ## evaluate survival 
      for(t in 1:len_t){
        
        q_S_tstar[, t] <- s_fun(q_samples_trans, dist=dist, t = tstar.[t])  
      } 
    }
  }
  
  
  ## need to compute weights on log scale for stability 

  # Edit DMC: need to handle multiple timepoints with 'apply'
  if(length(tstar.) == 1){
     log_p <- ex_info$loss(q_S_tstar) + 
     dmvnorm(q_samples, coeff, cov, log = TRUE)
  } else {
    log_p <- apply(q_S_tstar,1,ex_info$loss) + 
      dmvnorm(q_samples, coeff, cov, log = TRUE)
    }
  
  log_q <- dmvnorm(q_samples, mean = coeff, 
                   sigma = alpha * cov, log = TRUE)
  
  log_w <- log_p - log_q
  log_w <- log_w - max(log_w)   
  
  w0 <- exp(log_w)  
  wbar <- mean(w0)
  cov_calc <- cov.wt(q_samples, wt = w0)
  
  post_mean <- cov_calc$center
  post_cov <- cov_calc$cov
  cv_star <- var(w0)/mean(w0)^2
  ESS1 = iter/(1 + cv_star)
  ESS2 = sqrt(mean((w0/wbar - 1)^2))
  
  return(list(diag = c(cv_star, ESS1, ESS2), weights = w0, pars = list(post_mean, post_cov), samples=q_samples))
} # IS


## Function to identify good choice of alpha 

get_alpha <- function(len = 2*1e2, alpha_min = .5, alpha.max = 10, ...){
  
  alpha_vec <- seq(alpha_min, alpha.max, length.out=len)
  store_mat <- matrix(0, nrow = len, ncol = 3) 
  for(i in 1:len) store_mat[i, ] <- IS(alpha = alpha_vec[i], ...)$diag
  
  ## choose best alpha and compute required number of iterations based on effective sample size
  alpha_star <- alpha_vec[which.max(store_mat[, 2])]
  iter_star <- round(5000 * (5000/max(store_mat[, 2])))
  
  return(list(alpha_star = alpha_star, iter_star = iter_star, ESS_mat = store_mat, alpha_vec = alpha_vec))
}


## main function 
## Inputs: MLE parameter estimates, covariance, distribution, prior
## Outputs: posterior mean and covariance, samples, MLE parameter estimates, prior, importance sampling diagnostics
##

is_surv <- function(
    #parameter estimates from flexsurv
    coeff,
    # covariance matrix
    cov,
    # transformed to MVN scale already?
    transformed=TRUE,
    # Distribution: one of ("exponential","weibull","gompertz","lognormal","llogis","gengamma")
    dist,
    # Loss function/prior
    ex_info){
  
    ## find best alpha 
    cat("Choosing alpha...\n")
  
    # Choose alpha for importance sampling
    alpha_list <- get_alpha(coeff=coeff,cov=cov,dist=dist, ex_info = ex_info)
  
  ## estimate pars of restricted model 
  cat("fitting restricted model ...\n")
  
  # Run importance sampling for optimal alpha and get posterior parameter estimates
  new_coeffs <- IS(alpha = alpha_list$alpha_star, iter = alpha_list$iter_star, 
                   #surv_fit = surv_fit, 
                   coeff=coeff,
                   cov=cov,
                   dist=dist,
                   ex_info = ex_info)$pars
  
  post_mean <-  new_coeffs[[1]]
  post_cov <- new_coeffs[[2]]
  par_new <- rmvnorm(5000, mean = post_mean, sigma = post_cov)
  
  cat("Done.\n")
  
  output <- list()
  output$post_mean <- post_mean
  output$post_cov <- post_cov
  output$par_new <- par_new
  
  output$orig <- list("coefficients"=coeff,
                      "cov"=cov,
                      "dist"=dist)
  output$ex_info <- ex_info
  output$alpha_star <- alpha_list$alpha_star
  output$iter_star <- alpha_list$iter_star
  output$ESS_mat <- alpha_list$ESS_mat
  output$alpha_vec <- alpha_list$alpha_vec
  output
}

## compute updated coefficients and survival terms etc.  
## Updated version with ggplot


is_surv_viz_gg <- function(is_surv, times=tseq2, 
                           #tstar,
                           dist,what=1:3,...){
  ind <- rep(FALSE, 3)
  ind[what] <- TRUE
  ## assess alpha -- diagnostic plots
  
  ## list to save plots
  outplots<-list()
  
  if(ind[1]) {
    
    df<-data.frame("alpha"=is_surv$alpha_vec,
                   #"COV"=is_surv$ESS_mat[,1],
                   "ESS"=is_surv$ESS_mat[,2]#,
                   #"ESS2"=is_surv$ESS_mat[,3]
                   ) #%>%
     # pivot_longer(cols=2:4)
    
    g<-ggplot(data=df,aes(x=alpha,y=ESS)) + 
     geom_line() + 
     geom_vline(xintercept=is_surv$alpha_star) + 
    labs(title=paste0("Effective sample size - ",distributions[dist]))
    # + facet_wrap(name~.,ncol=1,scales="free_y")
      
    
    print(g)
    
    cat('\n\n') 
  outplots[["diagnostics"]]<-g
  
  }
  
  ## plot parameters old vs new
  
  par_old <- rmvnorm(5000, mean = is_surv$orig$coefficients, sigma = is_surv$orig$cov)
  # par_new <- rmvnorm(5000, mean = is_surv$post_mean, sigma = is_surv$post_cov)
  if(ind[2]) {
    
    if(dist == "exponential"){
      par_old <- data.frame(log_rate = par_old)
      par_old["Method"]<-"Trial data only"
      
      par_new<-data.frame(log_rate = is.models$exponential$par_new)
      par_new["Method"]<-"Trial data with external information"
      
      df<-bind_rows(par_old,par_new)
      
      g2<- ggplot(df, aes(x = log_rate, colour = Method, fill = Method)) + 
        geom_density(alpha=0.4 ) +
        labs(title=paste0("Parameters - ",distributions[dist]))
    } else{
      par_old<-data.frame(par_old)
      par_old["Method"]<-"Trial data only"
      
      par_new<-data.frame(is_surv$par_new)
      par_new["Method"]<-"Trial data with external information"
      
      df<-bind_rows(par_old,par_new)
      
      g2<-ggpairs(df, aes(colour = Method, alpha = 0.4),
                  columns = 1:(ncol(df)-1),
                  title=paste0("Parameters - ",distributions[dist]))
    }
    
        
    print(g2)

    cat('\n\n') 
    
    outplots[["parameters"]]<-g2
  }
  

 
  ## Plot survival curves
  
  if(ind[3]){
    
    # Extract survival quantiles
    oldcurve<-get_sims(dist=dist,coeff = is_surv$orig$coefficients, cov = is_surv$orig$cov, tst=is_surv$ex_info$tstar, times=times, tmax=max(times))$survsummary
    newcurve<-get_sims(dist=dist,coeff = is_surv$post_mean, cov = is_surv$post_cov, tst=is_surv$ex_info$tstar, times=times, tmax=max(times))$survsummary
    
    oldcurve["Method"]="Trial data only"
    newcurve["Method"]="Trial data with external information"
    plotcurves<-bind_rows(oldcurve,newcurve)
    
    # Lower and upper limits to plot prior
    
    #lower_prob<-is_surv$ex_info_loss
    
    ## Data frame of confidence intervals for priors
    prior.df<-data.frame("tst"=is_surv[["ex_info"]]$tstar,
                         "lower"=is_surv[["ex_info"]]$lower.probs,
                         "upper"=is_surv[["ex_info"]]$upper.probs
                         )
    
    g3<-ggplot(data=plotcurves,aes(x=time))+
      geom_line(aes(y=S_median,colour=Method),lwd=1)+
      geom_ribbon(aes(ymin=S_lower,ymax=S_upper,colour=Method,fill=Method),alpha=0.1,linetype="dashed")+
      #geom_segment(x=is_surv$ex_info$tstar, xend=is_surv$ex_info$tstar,y=is_surv[["ex_info"]]$lower.probs,yend=is_surv[["ex_info"]]$upper.probs,colour="black",lwd=1)+
      geom_segment(data=prior.df,
                   aes(x=tst,xend=tst,y=lower,yend=upper),
                   colour="black",
                   lwd=1)+
      labs(y="S(t)",x="time (t)",title=paste0("Survival - ",distributions[dist]))
    
    print(g3)
    
    cat('\n\n') 
    
    outplots[["survival"]]<-g3
    }
  
  return(outplots)
}

## Additional functions added by DMC

## get_auc

# This function takes as input a data frame of parameter simulations, one per row
# and on the multivariate normal scale, and integrates over the range of
# interest to calculate AUC, returning a vector of length equal to the
# number of simulations

get_auc<-function(params,
                  dist, # Distribution
                  lwr=0,
                  upr=max(tseq2))
{
  
  integrate_auc<-function(x,dist)
  {
    if(dist==1|dist=="exponential")
    {
      #upr<-qexp(p=1-1e-2,rate = exp(x[1]))
      return(integrate(pexp,lower=lwr,upper=upr,lower.tail=F,
                      rate = exp(x[1]))$value)}
    if(dist==2|dist=="weibull")
    {
      #upr<-qweibull(p=1-1e-2,shape = exp(x[1]), scale = exp(x[2]))
      return(integrate(pweibull,lower=lwr,upper=upr,lower.tail=F,
                      shape = exp(x[1]), scale = exp(x[2]))$value)}
    if(dist==3|dist=="gompertz")
    {
      #upr<-qgompertz(p=1-1e-2,shape = x[1], rate = exp(x[2]))
      return(integrate(pgompertz,lower=lwr,upper=upr,lower.tail=F,
                      shape = x[1], rate = exp(x[2]))$value)}
    if(dist==4|dist=="lognormal")
    {
     # upr<-qlnorm(p=1-1e-2, meanlog = x[1], sdlog = exp(x[2]))
      return(integrate(plnorm,lower=lwr,upper=upr,lower.tail=F,
                      meanlog = x[1], sdlog = exp(x[2]))$value)}
    if(dist==5|dist=="llogis")
    {
    #  upr<-qllogis(p=1-1e-2,shape = exp(x[1]), scale = exp(x[2]))
      return(integrate(pllogis,lower=lwr,upper=upr,lower.tail=F,
                      shape = exp(x[1]), scale = exp(x[2]))$value)}
    if(dist==6|dist=="gengamma")
    {
     # upr<-qgengamma(p=1-1e-2, mu = x[1],sigma = exp(x[2]),Q = x[3])
      return(integrate(pgengamma,lower=lwr,upper=upr,lower.tail=F,
                      mu = x[1],sigma = exp(x[2]),Q = x[3])$value)}
    
  }
  
  num_cols <- unlist(lapply(data.frame(params), is.numeric))  
  
  return(apply(data.frame(params[,num_cols]),1,integrate_auc,dist=dist))
  
}

## inv_trans
## Transforms parameters back from MVN scale to natural scale
## Takes as input a matrix or DF of parameter simulations, 1 per row

inv_trans<-function(dist,sims)
{
  if(ncol(sims)==1)
  {
    return(inv_transforms[["exponential"]][[1]](sims))
  } else
  {
    out<-matrix(NA,nrow=nrow(sims),ncol=ncol(sims))
    for(nc in 1:ncol(sims)){
      out[, nc] <- inv_transforms[[dist]][[nc]](sims[, nc])
    }
    return(out)
  }
  
}

## trans()
## Transforms parameters to MVN scale
## Takes as input a matrix or DF of parameter simulations, 1 per row

trans<-function(dist,sims)
{
  if(ncol(sims)==1)
  {
    return(transforms[["exponential"]][[1]](sims))
  } else
  {
    out<-matrix(NA,nrow=nrow(sims),ncol=ncol(sims))
    for(nc in 1:ncol(sims)){
      out[, nc] <- transforms[[dist]][[nc]](sims[, nc])
    }
    return(out)
  }
  
}


## get_sims

## This function takes survival curve parameter mean and covariance estimates on the MVN scales
## Returns: 
## 1. simulations on both the MVN and natural scales
## 2. simulated survival predictions at time tstar
## 3. simulated AUC 

## Edit DMC 20240215: new parameter 'tmax' for maximum time used to compute AUC.
## Default to max of tseq2

## Edit DMC 20240605: get rid of 'survtimes' output to improve speed

get_sims<-function(dist,coeff,cov,nsim=5000,tst=tstar,times=tseq2,tmax=max(tseq2))
{
  out<-list()
  
  out[["coeff"]]<-coeff
  out[["cov"]]<-cov
  # Simulated parameter draws on the MVN scale
  out[["sims.mvn"]]<-as.matrix(rmvnorm(n=nsim,mean=coeff,sigma=cov))
  # flexsurv doesn't give the exponential parameter a 'name' so add one manually
  if(dist=="exponential")
  {
    colnames(out[["sims.mvn"]])<-"rate"
  }  
  
  # Simulated parameter draws on the natural scale
  out[["sims.nat"]]<-inv_trans(dist=dist,sims=out[["sims.mvn"]])
  
  # Simulated survival estimates at time tstar
  # This doesn't need cases but I've arranged it this way to avoid breaking anything for the length 1 case
  if (length(tst==1))
  {out[["s.tstar"]]<-data.frame("S.tstar"=as.numeric(s_fun(t=tst,dist=dist,trans_pars=out[["sims.nat"]])))}
  if (length(tst>1)){
    out[["s.tstar"]]<-data.frame(sapply(tst,function(x){s_fun(t=x,dist=dist,trans_pars=out[["sims.nat"]])}))
    names(out[["s.tstar"]])<-paste0("S(",tst,")")
     }
  
  
  # Simulated survival estimates across all timepoints
  
  # Compute matrix of survival time estimates
  S_mat <- matrix(0, nrow = length(times), ncol = nsim)
  S_mat<-t(sapply(times,
                  function(x) {s_fun(t=x,dist=dist,trans_pars=out[["sims.nat"]])}
                  ))
  #for(t1 in 1:length(times))
  #{S_mat[t1, ] <- s_fun(trans_pars=out[["sims.nat"]], 
  #                      dist=dist, 
  #                      t = times[t1])}
  
  #out[["survtimes"]]<-bind_cols(data.frame("time"=times),data.frame(S_mat))
  
  # Quantiles of survival esimtates over time
  
  S_upper <- apply(S_mat, 1, quantile, probs = 0.975)
  S_lower <- apply(S_mat, 1, quantile, probs = 0.025)
  S_median <- apply(S_mat, 1, median)
  
  S_mat_summary <- cbind(S_lower,S_median,S_upper)
  
  out[["survsummary"]]<-data.frame("time"=times) %>%
    bind_cols(S_mat_summary)
  
  
  # AUC simulations
  out[["AUC"]]<-get_auc(dist=dist,params = out[["sims.mvn"]],upr = tmax)
  
  return(out)
}

