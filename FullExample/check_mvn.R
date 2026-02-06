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
lapply(is.models2, check_mvn_wrapper, log = "xy")  
