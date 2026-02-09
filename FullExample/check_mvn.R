## Assuming that run_all.R has been run already 
## Assuming that you are in correct directory (".../SurvISArticle/FullExample/")
# source("qqmvn.R")

lapply(is.models, check_mvn_wrapper)
## returns a weird NULL output for some reason?? 
## But gives us the plots we want 

## Looking at this on log-log scale (for MVN) seems to be more helpful??
is.models2 <- NULL
is.models2[[1]] <- is.models[[2]]
is.models2[[2]] <- is.models[[3]]
is.models2[[3]] <- is.models[[4]]
is.models2[[4]] <- is.models[[5]]
is.models2[[5]] <- is.models[[6]]
lapply(is.models2, check_mvn_wrapper , log="xy"      )  



## Other diagnostic plots (for two param models only)

is.models3<-is.models

for(i in 1:6)
{
  sample.w<-sample(1:length(is.models3[[i]]$sample_weights),
                   replace=TRUE,
                   size=5000,
                   prob = is.models3[[i]]$sample_weights)
  
  is.models3[[i]]$mvn_resamples<-is.models3[[i]][["raw_samples_mvn"]][sample.w,]
  
}

# Exponential

plot(density(is.models3[[1]]$mvn_resamples))

# 3d density plots for two parameter models

library(MVN)

res2<-mvn(is.models3[[2]]$mvn_resamples)
res3<-mvn(is.models3[[3]]$mvn_resamples)
res4<-mvn(is.models3[[4]]$mvn_resamples)
res5<-mvn(is.models3[[5]]$mvn_resamples)

plot(res2, diagnostic = "multivariate", type = "persp")
plot(res2, diagnostic = "multivariate", type = "contour")
plot(res2, diagnostic = "multivariate", type = "qq")

plot(res3, diagnostic = "multivariate", type = "persp")
plot(res3, diagnostic = "multivariate", type = "contour")
plot(res3, diagnostic = "multivariate", type = "qq")


plot(res4, diagnostic = "multivariate", type = "persp")
plot(res4, diagnostic = "multivariate", type = "contour")
plot(res4, diagnostic = "multivariate", type = "qq")


plot(res5, diagnostic = "multivariate", type = "persp")
plot(res5, diagnostic = "multivariate", type = "contour")
plot(res5, diagnostic = "multivariate", type = "qq")

# Gen gamma - I think only QQ plots are sensible here

res6<-mvn(is.models3[[6]]$mvn_resamples)

plot(res6, diagnostic = "multivariate", type = "qq")



