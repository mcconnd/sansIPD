## This code summarises and plots probability density of survival estimates at time 10 years and 20 years from the prior, likelihood, and 'posterior'

surv.tstar<-data.frame()

for(dist in dists)
{
  # Extract parameter sims
  param.sims.mle<-mle.sims[[dist]][["sims.nat"]]
  param.sims.is<-is.sims[[dist]][["sims.nat"]]
  
  # Predict survival
  
  df.mle<-data.frame(sapply(c(120,240),function(x){s_fun(t=x,dist=dist,trans_pars=param.sims.mle)}))
  names(df.mle)<-c("S.10","S.20")
  df.mle[,"method"]<-"MLE"
  
  df.is<-data.frame(sapply(c(120,240),function(x){s_fun(t=x,dist=dist,trans_pars=param.sims.is)}))
  names(df.is)<-c("S.10","S.20")
  df.is[,"method"]<-"IS"
  
  df.st<-rbind.data.frame(df.mle,df.is)
  df.st[,"distribution"]<-dist
  
  # Combine into existing datafram
  
  surv.tstar<-rbind.data.frame(surv.tstar,df.st)
  
}

surv.tstar.long<-pivot_longer(surv.tstar,cols=1:2,values_to = "S.tstar") %>%
  mutate(Distribution=str_replace_all(distribution,distributions),
         Timepoint=str_replace_all(name,c("S.10"="10 years","S.20"="20 years")),
         Output=str_replace_all(method,c("MLE"="Trial data only","IS"="Trial data with external information")))


###


surv.s10s20.summary<-surv.tstar.long %>% 
  group_by(Distribution,Timepoint,method) %>%
  summarise(mean=mean(S.tstar),
            median=median(S.tstar),
            sd=sd(S.tstar),
            lwr.95=quantile(S.tstar,0.025),
            upr.95=quantile(S.tstar,0.975),
            ) %>%
  ungroup()%>%
  pivot_wider(names_from=method,
              values_from=4:8) %>%
  mutate(mean_diff=mean_IS-mean_MLE,
         var_ratio=sd_IS^2/sd_MLE^2,
         across(where(is.numeric),~format(round(.x,2),nsmall=2)),
         CrI_IS=paste0("(",lwr.95_IS," to ",upr.95_IS,")"),
         CrI_MLE=paste0("(",lwr.95_MLE," to ",upr.95_MLE,")"),
  ) %>%
  select(Distribution,Timepoint,mean_MLE,median_MLE,CrI_MLE,mean_IS,median_IS,CrI_IS,mean_diff,var_ratio)
  
#write.csv(surv.s10s20.summary,"./output/summary.s10s20.csv")


print(surv.s10s20.summary)


# Version in use
g.s10s20<-ggplot(data=surv.tstar.long,aes(x=S.tstar,fill=Output,height=after_stat(scaled)))+
  geom_density(alpha=0.5,#scale=1,
                      #rel_min_height=0.01,
                       trim=TRUE,
                      #stat="scaled",
                      bounds=c(0,1)
  )+
  #geom_histogram()+
  scale_fill_manual(values=cbPalette[2:3])+
  scale_colour_manual(values=cbPalette[2:3])+
  labs(x="OS probability",y="Density")+
  facet_wrap(Distribution~Timepoint,ncol=2,scales="free_y",
             labeller = 
               label_wrap_gen(multi_line=FALSE),
             strip.position="top"
               #labeller()
             )+
 # facet_grid(Distribution~Timepoint,scales="free",shrink=TRUE)+
  theme(legend.position = "bottom",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())





#ggsave("./output/plot_s10s20.png",plot=g.s10s20,width=7,height=7,units="in")