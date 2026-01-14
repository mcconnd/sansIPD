## Data frame of fitted survival curves (mean values only)
mle.points<-data.frame("time"=tseq2)

# Plot deterministic estimates of survival
for(dist in dists)
{
  # Transform parameters to natural scale
  
  pars.mvn<-is.models[[dist]][["orig"]][["coefficients"]]
  pars.nat<-inv_trans(dist=dist,
                      sims=as.matrix(t(pars.mvn),nrow=1,byrow=T))
  
  # Calculate survival times
  
  mle.points[,dist]=s_fun(tseq2,dist=dist,pars.nat)
}

# Extract KM estimates for plotting

KM.df<-data.frame("time"=summary(survfit(formula = surv.form,data=surv.IPD))[["time"]],
                  "OS"=summary(survfit(formula = surv.form,data=surv.IPD))[["surv"]],
                  "OS.lwr"=summary(survfit(formula = surv.form,data=surv.IPD))[["lower"]],
                  "OS.upr"=summary(survfit(formula = surv.form,data=surv.IPD))[["upper"]]
)


# Reshape for plotting

mle.points.df<-pivot_longer(mle.points,cols=2:7,names_to="dist",values_to = "OS")

## Extract deterministic IS method curve estimates

is.mean.df<-data.frame("time"=tseq2)

for (dist in dists)
{
  # Transform parameters to natural scale
  
  
  pars.nat<-inv_trans(dist=dist,
                      sims=as.matrix(t(is.models[[dist]][["post_mean"]]),nrow=1,byrow=T))
  
  
  # Calculate survival times
  
  is.mean.df[,dist]=s_fun(tseq2,dist=dist,pars.nat)
}


is.mean.df<-pivot_longer(is.mean.df,cols=2:7,names_to="dist",values_to = "OS") %>%
  mutate(Method="IS")


mle.points.df<-mutate(mle.points.df,Method="MLE")

is.compare.df<-bind_rows(is.mean.df,mle.points.df) %>%
  mutate(Distribution=dist)

# Rename distributions

is.compare.df$Distribution<-str_replace_all(is.compare.df$Distribution,
                                            distributions)


outplot<-ggplot(data=is.compare.df,aes(x=time,y=OS,colour=Distribution))+
  geom_line(lwd=1,alpha=0.8)+
  facet_wrap(.~factor(Method,
                      levels=c("MLE","IS"),
                      labels=c("Trial data only","Trial data with external information")),
             nrow=2)+
  geom_segment(x=60,xend=60,y=mu_t-1.96*sigma_t,yend=mu_t+1.96*sigma_t,lwd=1,colour="grey40")+
  geom_segment(x=59,xend=61,y=mu_t-1.96*sigma_t,yend=mu_t-1.96*sigma_t,lwd=1,colour="grey40")+
  geom_segment(x=59,xend=61,y=mu_t+1.96*sigma_t,yend=mu_t+1.96*sigma_t,lwd=1,colour="grey40")+
  theme_classic()+
  theme(legend.position = "bottom",
        panel.grid.major = element_line(colour="gray90"),
        axis.line.x.bottom = element_line(),
        plot.caption = element_text(hjust = 0))+
  labs(x="Months",y="Overall Survival")+
  scale_x_continuous(limits=c(0,120),breaks = seq(0,120,by=12),expand = expansion(mult=c(0.00,0.02)))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),expand = expansion(mult=c(0.001,0.05)))+
  geom_step(data=KM.df,aes(x=time,y=OS),colour="black")+  
  geom_step(data=KM.df,aes(x=time,y=OS.lwr),colour="black",linetype="dashed",lwd=0.5)+
  geom_step(data=KM.df,aes(x=time,y=OS.upr),colour="black",linetype="dashed",lwd=0.5)

plot(outplot)

ggsave("./output/plot_mle_is.png",plot=outplot,width=7, height=7,units="in",dpi=300)

outplot_wide<-ggplot(data=is.compare.df,aes(x=time,y=OS,colour=Distribution))+
  geom_line(lwd=1.2,alpha=0.8)+
  facet_wrap(factor(Method,
                    levels=c("MLE","IS"),
                    labels=c("Trial data only","Trial data with external information"))~.,
             ncol=2)+
  geom_segment(x=60,xend=60,y=mu_t-1.96*sigma_t,yend=mu_t+1.96*sigma_t,lwd=1,colour="grey40")+
  geom_segment(x=59,xend=61,y=mu_t-1.96*sigma_t,yend=mu_t-1.96*sigma_t,lwd=1,colour="grey40")+
  geom_segment(x=59,xend=61,y=mu_t+1.96*sigma_t,yend=mu_t+1.96*sigma_t,lwd=1,colour="grey40")+
  theme_classic()+
  theme(legend.position = "bottom",
        panel.grid.major = element_line(colour="gray90"),
        axis.line.x.bottom = element_line(),
        plot.caption = element_text(hjust = 0),
        text=element_text(size=15))+
  labs(x="Months",y="Overall Survival")+
  scale_x_continuous(limits=c(0,120),breaks = seq(0,120,by=12),expand = expansion(mult=c(0.00,0.02)))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),expand = expansion(mult=c(0.001,0.05)))+
  geom_step(data=KM.df,aes(x=time,y=OS),colour="black")+  
  geom_step(data=KM.df,aes(x=time,y=OS.lwr),colour="black",linetype="dashed",lwd=0.5)+
  geom_step(data=KM.df,aes(x=time,y=OS.upr),colour="black",linetype="dashed",lwd=0.5)

plot(outplot_wide)
ggsave("./output/plot_mle_is_wide.png",plot=outplot_wide,width=5, height=6,units="in")