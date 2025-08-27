## idea of doing this separately is that individual models can tuned as needed

ES_exp  <- fit.models.expert(
  formula=surv.form,
  data=surv.IPD,
  distr=c("exp"),
  method="bayes",
  iter = 5000,
  pool_type = "linear pool", 
  opinion_type = "survival",
  times_expert = tstar, 
  param_expert = expert_op)

## N.B., use exArgs to pass commands to rstan and JAGS (for c("gomp", "gga"))
## work horse function is expertsurv:::runBAYES

ES_wei  <- fit.models.expert(
  formula=surv.form,
  data=surv.IPD,
  distr=c("weibull"),
  method="bayes",
  iter = 5000,
  pool_type = "linear pool", 
  opinion_type = "survival",
  times_expert = tstar, 
  param_expert = expert_op)

ES_gomp  <- fit.models.expert(
  formula=surv.form,
  data=surv.IPD,
  distr=c("gompertz"),
  method="bayes",
  iter = 5000,
  pool_type = "linear pool", 
  opinion_type = "survival",
  times_expert = tstar, 
  param_expert = expert_op)

ES_ln  <- fit.models.expert(
  formula=surv.form,
  data=surv.IPD,
  distr=c("lognormal"),
  method="bayes",
  iter = 5000,
  pool_type = "linear pool", 
  opinion_type = "survival",
  times_expert = tstar, 
  param_expert = expert_op)

ES_logl  <- fit.models.expert(
  formula=surv.form,
  data=surv.IPD,
  distr=c("loglogistic"),
  method="bayes",
  iter = 5000,
  pool_type = "linear pool", 
  opinion_type = "survival",
  times_expert = tstar, 
  param_expert = expert_op)

ES_gengamma  <- fit.models.expert(
  formula=surv.form,
  data=surv.IPD,
  distr=c("gengamma"),
  method="bayes",
  iter = 5000,
  pool_type = "linear pool", 
  opinion_type = "survival",
  times_expert = tstar, 
  param_expert = expert_op)

## discrepancies mainly with gen Gamma and log normal so double check the performance here 
## Discrepancies should not be down to differences in Monte Carlo behaviour (i.e., poor MCMC performance) if at all possible

## Check MCMC performance ----

## R hat etc.
print(ES_exp$models[[1]])
print(ES_wei$models[[1]])
print(ES_gomp$models[[1]])

print(ES_ln$models[[1]])
print(ES_logl$models[[1]])
print(ES_gengamma$models[[1]])


##trace plots 
plot(as.mcmc(as.matrix(ES_exp$models[[1]])))
plot(as.mcmc(as.matrix(ES_wei$models[[1]])))
plot(as.mcmc(ES_gomp$models[[1]])) 

plot(as.mcmc(as.matrix(ES_ln$models[[1]])))
plot(as.mcmc(as.matrix(ES_logl$models[[1]])))
plot(as.mcmc(ES_gengamma$models[[1]])) 

