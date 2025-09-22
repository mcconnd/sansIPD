source("../functions.R")
library(tidyverse)
library(flexsurv)
library(expertsurv)
set.seed(12345)




### Basic input data (required for all methods)

# Load digitised IPD
ipd.list <- readRDS("ipd.list.RDS")
# Interim analysis 1, OS, pembro+axi arm only
surv.IPD <- ipd.list$ia1.OS.PA$IPD

#Formula for survival data
surv.form <- Surv(time, status)~1

# Timepoints for plotting (note: time in months)

tseq2 <- seq(0, 12*20, len = 1e4)

# Distribution names
dists<-c("exponential","weibull","gompertz","lognormal",
         "llogis","gengamma")

# Long names of distributions (named vector)
distributions<-c("exponential"="Exponential",
                 "weibull"="Weibull",
                 "gompertz"="Gompertz",
                 "lognormal"="Log-normal",
                 "llogis"="Log-logistic",
                 "gengamma"="Gen. Gamma")

# Maximum possible OS for computing AUC
# Set this to 100 years - actually some distributions never converge so we need to make a choice here
t.max<-12*100


### Specification of prior information

# This is the time point at which the prior for survival is to be specified.
tstar<-60 # months in this example

## Expert survival based on Clinical Opinion doc Expert 2 Table 1
## Also accounting for NICE slides https://www.nice.org.uk/guidance/ta650/documents/1 
lower_prob <- 0.20
upper_prob <- 0.55
mu_t <- (lower_prob+upper_prob)/2


## Further we assume estimate is normally distributed with 95% confidence between limits

sigma_t <- (upper_prob - mu_t) / (1.96)


# Set up prior for IS method
ex.info<-list("tstar"=tstar,
                  "loss"=function(x){dnorm(x,mean=mu_t,sd=sigma_t,log=TRUE)},
                  "lower.probs"=lower_prob,
                  "upper.probs"=upper_prob
)

# Set up prior for expertsurv method
expert_op <- list()
expert_op[[1]] <- data.frame(dist = c("norm"),
                             wi = c(1),
                             param1 = c(mu_t),
                             param2 = c(sigma_t),
                             param3 = c(NA))

# Colourblind palette for plots
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

