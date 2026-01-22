# This code summarises AUC simulations from the model and plots them

#AUC for the models


# Extract AUC simulations and combine into single data frame
auc.is<-data.frame(matrix(NA,nrow=5000,ncol=6))
names(auc.is)<-dists
auc.mle<-auc.is

for( dist in dists)
{
  auc.is[dist]<-is.sims[[dist]][["AUC"]]
  auc.mle[dist]<-mle.sims[[dist]][["AUC"]]
}


auc.is["Method"]="IS"
auc.mle["Method"]="MLE"

auc.df<-bind_rows(auc.is,auc.mle) %>%
  pivot_longer(1:6,
               names_to="Distribution",
               values_to="AUC") %>%
  mutate(Distribution=str_replace_all(Distribution,distributions))


### AUC (actually 100 year RMST) summary

auc.summary<-auc.df %>%
  group_by(Distribution,Method) %>%
  summarise(mean=mean(AUC),
            sd=sd(AUC),
            median=median(AUC),
            lwr.95=quantile(AUC,0.025),
            upr.95=quantile(AUC,0.975)) 

#print(auc.summary)



# AUC density plot
# No truncation of AUC values but truncation of graph
library(ggridges)

g.auc<-ggplot(data=auc.df,aes(x=AUC,y=Distribution,colour=Method,fill=Method))+
  geom_density_ridges(alpha=0.5,scale=0.95)+
  scale_fill_manual(values=cbPalette[2:3])+
  scale_colour_manual(values=cbPalette[2:3])+
  labs(x="Mean OS, Months",y="Density")+
  xlim(0,360)

##g.auc