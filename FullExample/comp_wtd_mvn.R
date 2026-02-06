# Generate summary survival curves (posterior median and 95% CrI) from the weighted samples and corresponding multivariate normal approximation
# Note: code is slow to run, currently sampling a reduced grid of 1000 timepoints. Best run outside of Rmarkdown.

wt.mvn.curves<-list()

for (dist in dists)
{
  
  wt.mvn.curves[[dist]]<-compare_weighted_sims(sims_raw=is.models[[dist]]$raw_samples_nat,
                                           sims_wts = is.models[[dist]]$sample_weights,
                                           sims_mvn = is.models[[dist]]$par_new,
                                           times=sample(tseq2,1000),
                                           dist=dist
                                           )
  
}


saveRDS(wt.mvn.curves,"wt.mvn.curves.RDS")