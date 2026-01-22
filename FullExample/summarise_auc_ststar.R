# Combined AUC and S(t*) density plots for article
# Need to run auc_summary.R and s_tstar_summary.R first

# Revised to include % agreement with prior

# New code (grobs)
g.auc2<-g.auc+theme(legend.position = "none")
g.ststar2<-g.ststar+
  labs(y=element_blank())+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  #geom_vline(xintercept=s5.obs,linetype="dashed",colour="gray50")+
  geom_segment(aes(y=Distribution,yend=as.numeric(factor(Distribution))+0.9
                   ),
               x=s5.obs,xend=s5.obs,linetype="dashed",colour="gray50")

# Extract legend

library(ggpubr)
plot.auc.ststar<- ggarrange(g.auc,g.ststar2,
                            nrow=1,common.legend = TRUE,
                            legend = "bottom",
                            widths=c(1.3,1),
                            legend.grob = get_legend(g.ststar2,position="bottom"))

print(plot.auc.ststar)

ggsave("./output/plot_auc_s5.png",plot.auc.ststar,width=7,height=7,units="in",dpi=300)

### Tables for article

library(scales)

# Summary table of survival at time t_star
# Calculate mean difference and variance ratio 
# All other code is just formatting

s.tstar.print<-surv.tstar.summary %>%
  pivot_wider(names_from = name,values_from = 3:10) %>% 
  mutate(mean_diff=mean_IS-mean_MLE,
         var.ratio=sd_IS^2/sd_MLE^2) %>%
  select(Distribution,
         mean_MLE,lwr.95_MLE,upr.95_MLE,
         mean_IS,lwr.95_IS,upr.95_IS,
         mean_diff, 
         var.ratio,
         MSE_MLE,MSE_IS,
         MAE_MLE,MAE_IS,
         Agreement_MLE,Agreement_IS
  )%>%
  ungroup()%>%
  inner_join(surv.tstar.det,by=join_by(Distribution)) %>%
  mutate(
    Det_diff=Det_IS-Det_MLE,
    (across(c(mean_MLE:mean_diff,Det_MLE:Det_diff),~percent(.x,accuracy=0.1))),
         var.ratio=format(round(var.ratio,2),nsmall=2),
         across(MSE_MLE:MAE_IS,~format(round(.x,3),nsmall=3)),
         across(Agreement_MLE:Agreement_IS,~percent(.x,accuracy=0.1))
         ) %>%
  mutate(MLE=paste0(mean_MLE," (",lwr.95_MLE,", ",upr.95_MLE,")"),
         IS=paste0(mean_IS," (",lwr.95_IS,", ",upr.95_IS,")")
  ) %>%
  ungroup() %>%
  select(Distribution,
         Det_MLE,Det_IS,Det_diff,
         Prob_MLE=MLE,Agreement_MLE,Prob_IS=IS,Agreement_IS,
         mean_diff,var.ratio,
         MSE_MLE,MSE_IS,MAE_MLE,MAE_IS) 


print(s.tstar.print)

write.csv(s.tstar.print,"./output/summary.ststar.csv")

# Same for AUC

auc.print <- auc.summary %>%
  pivot_wider(names_from = Method,values_from = 3:7) %>% 
  mutate(mean_diff=mean_IS-mean_MLE,
         var.ratio=sd_IS^2/sd_MLE^2) %>% 
  mutate(across(where(is.numeric),~format(round(.x,2),nsmall=2)))  %>%
  mutate(MLE=paste0(mean_MLE," (",lwr.95_MLE,", ",upr.95_MLE,")"),
         IS=paste0(mean_IS," (",lwr.95_IS,", ",upr.95_IS,")"),
         Var.Ratio=var.ratio
  ) %>%
  ungroup() %>%
  select(Distribution,MLE,IS,Mean.Diff=mean_diff,Var.Ratio) 


print(auc.print)
write.csv(auc.print,"./output/summary.auc.csv")