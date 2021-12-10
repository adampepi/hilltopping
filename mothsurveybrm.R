library(ggplot2)
library(AICcmodavg)
library(brms)
library(loo)
library(tidybayes)
library(tidyverse)
library(ggplot2)
library(modelr)
library(cowplot)
library(performance)
setwd("~/Documents/Research/dissertation/12 sites analyses/12sitesanalyses")
moths<-read.csv("mothsurvey2.csv")
str(moths)

head(moths)
mothsCI<-read.csv("mothsCI.csv")


#####
#Moths Connectivity
str(mothsCI)
mothsCI$CIS<-scale(log(mothsCI$connectivityalpha1k))

prior1<-c(set_prior(class="b",prior="normal(0,10)"),set_prior(class="Intercept",prior="normal(0,10)"))

prior2<-c(set_prior(class="Intercept",prior="normal(0,10)"))



mci1<-brm(count~CIS+(1|site),data=mothsCI,family='poisson',prior=prior1)

summary(mci1)


mci2<-brm(count~CIS+elevation+(1|site),data=mothsCI,family='poisson',prior=prior1)
summary(mci2)
mci3<-brm(count~CIS*elevation+(1|site),data=mothsCI,family='poisson',prior=prior1)
summary(mci3)
mci4<-brm(count~elevation+(1|site),data=mothsCI,family='poisson',prior=prior1)
summary(mci4)
mci5<-brm(count~(1|site),data=mothsCI,family='poisson',prior=prior2)
summary(mci4)

brms::pp_check(mci2,nsamples=50)
brms::pp_check(mci1,nsamples=50)
brms::pp_check(mci3,nsamples=50)

mci1<-add_criterion(mci1,"waic")
mci2<-add_criterion(mci2,"waic")
mci3<-add_criterion(mci3,"waic")
mci4<-add_criterion(mci4,"waic")
mci5<-add_criterion(mci5,"waic")

mcicomp<-loo_compare(mci1,mci2,mci3,mci4,mci5,criterion = "waic")
print(mcicomp,simplify = F)

write.table(mcicomp,'mcicomp.csv')

str(mothsCI)
connplot<-mothsCI %>%
  data_grid(CIS = seq_range(CIS, n = 101),elevation=mean(elevation)) %>%
  add_fitted_draws(mci2, n = 100,re_formula = NA) %>%
  ggplot(aes(x = CIS, y =count )) +
  geom_line(aes(y = .value,group = paste(.draw)), alpha = .1) +
  stat_lineribbon(aes(y = .value),.width = c(0, 0, 0),show.legend = F)+
  geom_point(data = mothsCI) +theme_classic()+scale_y_log10()+xlab('Connectivity')+ylab("Moth count")

elevationplot<-mothsCI %>%
  data_grid(elevation= seq_range(elevation, n = 101),CIS=mean(CIS)) %>%
  add_fitted_draws(mci2, n = 100,re_formula = NA) %>%
  ggplot(aes(x = elevation, y =count )) +
  geom_line(aes(y = .value,group = paste(.draw)), alpha = .1)  +
  stat_lineribbon(aes(y = .value),.width = c(0, 0, 0),show.legend = F)+
  geom_point(data = mothsCI) +
  scale_color_brewer(palette = "Dark2")+theme_classic()+scale_y_log10()+xlab('Elevation')+ylab("Moth count")


plot_grid(connplot,elevationplot,nrow=2,labels=c("A","B"))

###Cate4pillars
###
catCI<-read.csv("catmothCI.csv")
catCI
str(catCI)
catCI$logcount<-log(catCI$catcount*15+1)
catCI$count<-as.integer(catCI$catcount*15)
catCI$CI[catCI$CI==0]<-4.941077e-313
catCI$CIS<-scale(log(catCI$CI))

mcci1<-brm(logcount~CIS+(1|site),data=catCI,iter=10000, prior=prior1)
summary(mcci1)
posterior_interval(mcci1,prob=0.95)
posterior_interval(mcci1,prob=0.90)
bayes_R2(mcci1)

mcci1.1<-brm(logcount~(1|site),data=catCI,iter=10000, prior=prior2)
summary(mcci1.1)


mcci1<-add_criterion(mcci1,"waic")
mcci1.1<-add_criterion(mcci1.1,"waic")
mccicomp<-loo_compare(mcci1,mcci1.1,criterion = "waic")
str(mcicomp)
mccitable<-as.data.frame(mccicomp)
mccitable$deltawaic<-mccitable$waic-min(mccitable$waic)
mccitable

connplot2<-catCI %>%
  data_grid(CIS = seq_range(CIS, n = 101)) %>%
  add_fitted_draws(mcci1, n = 100,re_formula = NA) %>%
  ggplot(aes(x = CIS, y =logcount )) +
  geom_line(aes(y = .value,group = paste(.draw)), alpha = .1) +
  stat_lineribbon(aes(y = .value),.width = c(0, 0, 0),show.legend = F)+
  geom_point(data = catCI) +theme_classic()+xlab('Connectivity')+ylab("Log caterpillar density")
connplot2
