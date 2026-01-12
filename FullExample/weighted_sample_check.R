### Compute 5-year OS and mean OS using weighted samples from IS method

wt.summary<-data.frame()

for(dist in dists)
{
  # Extract parameter sims and weights
  raw.samples<-is.models[[dist]]$raw_samples_nat
  raw.samples.mvn<-is.models[[dist]]$raw_samples_mvn
  sample.weights<-is.models[[dist]]$sample_weights
  
  # Predict survival
  # From raw samples
  st.raw<- data.frame("S.tstar"=s_fun(t=tstar,dist=dist,trans_pars = raw.samples),
                      "S.mean"=get_auc(params=raw.samples.mvn,dist=dist,lwr=0,upr = t.max),
                      "Method"="Samples",
                      "Weight"=sample.weights,
                      "Distribution"=dist)
  
  # From summary MVN distribution
  st.mvn<-data.frame("S.tstar"=unname(is.sims[[dist]][["s.tstar"]]),
                     "S.mean"=is.sims[[dist]][["AUC"]],
                     "Method"="MVN",
                     "Weight"=1,
                     "Distribution"=dist)
 
 
 
  # Combine into existing datafram
  
  wt.summary<-rbind.data.frame(wt.summary,st.raw,st.mvn)
  
}


##Summarise
wt.summary.long<-wt.summary %>% pivot_longer(1:2) 

wt.summary.long%>%
  group_by(Distribution,name,Method) %>%
  summarise(mean=weighted.mean(value,w=Weight),
            median=wtd.quantile(value,weights=Weight,probs=0.5),
            lwr.95=wtd.quantile(value,weights=Weight,probs=0.025),
            upr.95=wtd.quantile(value,weights=Weight,probs=0.975)
            ) %>%
  pivot_wider(values_from=mean:upr.95,names_from = Method)
  
## Boxplot

ggplot(data=wt.summary.long,aes(y=value,x=Distribution,fill = Method,weight=Weight))+
  geom_boxplot(outliers=FALSE)+
  facet_wrap(.~name,ncol=2,scales="free_y")+
    theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank())

## 