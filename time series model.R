rm(list=ls())
library (rstan)
library(loo)
library(tidyverse)
library(tidybayes)
library(bayesplot)
library(cowplot)
library(loo)


setwd("~/Documents/Research/dissertation/12 sites analyses/hilltopping")

##Data processing
connectivity<-read.csv('connectivity2.csv')
CImat1<-connectivity
CImat<-CImat1[,2:16]
rownames(CImat)<-CImat1[,1]
connectivity2<-read.csv('hanskiconnectivity.csv')
CImat2<-connectivity2
CImat2
HCImat<-CImat2[,2:16]
rownames(HCImat)<-CImat2[,1]
HCImat
str(HCImat)
HCImat<-log(as.matrix(HCImat)+7.015732e-321)##add small constant and log
str(HCImat)
HCImat
connectivity3<-read.csv('connectivityelevation.csv')
CImat3<-connectivity3
ECImat<-CImat3[,2:16]
rownames(ECImat)<-CImat3[,1]
HCImat<-as.matrix(HCImat)
CImat<-as.matrix(CImat)
ECImat<-as.matrix(ECImat)
HCImatscaled<-(HCImat-mean(HCImat))/sd(HCImat)
ECImatscaled<-(ECImat-mean(ECImat))/sd(ECImat)
CImatscaled<-CImat
CImat<-as.matrix(CImatscaled)
CImat<-t(CImat)
CImat<-t(CImat)
HCImat<-as.matrix(HCImatscaled)
ECImat<-as.matrix(ECImatscaled)
CImat
wetdry2<-read.csv('wetdry.csv')
wetdry1<-wetdry2[,2]
weather<-read.csv('weather.csv')
str(weather)
tempscaled<-as.vector(scale(weather$avgtm1))
precipscaled<-as.vector(scale(weather$preciptm1))
datamat<-read.csv("sitesdatamat.csv")
y<-as.matrix(datamat[,2:16])
y
Z.model4<-seq(1,12,by=1)
states<-Z.model4
obsVariances = rep(1, nrow(y))
N = ncol(y)
M = nrow(y)
y<-t(y)
y<-t(y)
N
M
row_indx_pos = matrix((rep(1:M, N)), M, N)[which(!is.na(y))]
col_indx_pos = matrix(sort(rep(1:N, M)), M, N)[which(!is.na(y))]
n_pos = length(row_indx_pos)
y = y[which(!is.na(y))]
y
W=2
wetdry<-wetdry1[1:12]
length(wetdry)

##Load text models
y
fit1<-rstan::stan(file='models/all same.stan',
                   data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=CImat,
                               "states"=states, "S" = max(states), "W"=W,
                               "wetdry"=wetdry,
                               "obsVariances"=obsVariances,
                               "n_obsvar" = max(obsVariances),  
                               "n_pos" = n_pos,
                               "col_indx_pos" = col_indx_pos,
                               "row_indx_pos" = row_indx_pos,
                               "y_int"=round(y)),
                   pars = c("pred","log_lik", "alpha_0","alpha_1","alpha_2","beta_precip","beta_CI","beta_temp","beta_pt","sigma_obs"), chains = 3,
                   iter = 20000, thin = 3, cores=3)
print(fit1)

fit2<-rstan::stan(file='models/full dd temp precip.stan',
                   data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=CImat,
                               "states"=states, "S" = max(states), "W"=W,
                               "wetdry"=wetdry,
                               "obsVariances"=obsVariances,
                               "n_obsvar" = max(obsVariances),  
                               "n_pos" = n_pos,
                               "col_indx_pos" = col_indx_pos,
                               "row_indx_pos" = row_indx_pos,
                               "y_int"=round(y)),
                   pars = c("pred","log_lik", "alpha_0","alpha_1","alpha_2","beta_precip","beta_temp"), chains = 3,
                   iter = 20000, thin = 3, cores=3)
fit3<-rstan::stan(file='models/full model no pt interaction.stan',
                   data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=CImat,
                               "states"=states, "S" = max(states), "W"=W,
                               "wetdry"=wetdry,
                               "obsVariances"=obsVariances,
                               "n_obsvar" = max(obsVariances),  
                               "n_pos" = n_pos,
                               "col_indx_pos" = col_indx_pos,
                               "row_indx_pos" = row_indx_pos,
                               "y_int"=round(y)),
                   pars = c("pred","log_lik", "alpha_0","alpha_1","alpha_2","beta_precip","beta_CI","beta_temp"), chains = 3,
                   iter = 20000, thin = 3, cores=3)

fit4<-rstan::stan(file='models/full dd temp.stan',
                   data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=CImat,
                               "states"=states, "S" = max(states), "W"=W,
                               "wetdry"=wetdry,
                               "obsVariances"=obsVariances,
                               "n_obsvar" = max(obsVariances),  
                               "n_pos" = n_pos,
                               "col_indx_pos" = col_indx_pos,
                               "row_indx_pos" = row_indx_pos,
                               "y_int"=round(y)),
                   pars = c("pred","log_lik", "alpha_0","alpha_1","alpha_2","beta_temp"), chains = 3,
                   iter = 20000, thin = 3, cores=3)

fit5<-rstan::stan(file='models/full dd precip.stan',
                   data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=CImat,
                               "states"=states, "S" = max(states), "W"=W,
                               "wetdry"=wetdry,
                               "obsVariances"=obsVariances,
                               "n_obsvar" = max(obsVariances),  
                               "n_pos" = n_pos,
                               "col_indx_pos" = col_indx_pos,
                               "row_indx_pos" = row_indx_pos,
                               "y_int"=round(y)),
                   pars = c("pred","log_lik", "alpha_0","alpha_1","alpha_2","beta_precip"), chains = 3,
                   iter = 20000, thin =5, cores=3)
fit6<-rstan::stan(file='models/full model no dd2.stan',
                   data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=CImat,
                               "states"=states, "S" = max(states), "W"=W,
                               "wetdry"=wetdry,
                               "obsVariances"=obsVariances,
                               "n_obsvar" = max(obsVariances),  
                               "n_pos" = n_pos,
                               "col_indx_pos" = col_indx_pos,
                               "row_indx_pos" = row_indx_pos,
                               "y_int"=round(y)),
                   pars = c("pred","log_lik", "alpha_0","alpha_1","beta_precip","beta_CI","beta_temp","beta_pt"), chains = 3,
                   iter = 20000, thin = 3, cores=3)

fit7<-rstan::stan(file='models/full dd CI.stan',
                   data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=CImat,
                               "states"=states, "S" = max(states), "W"=W,
                               "wetdry"=wetdry,
                               "obsVariances"=obsVariances,
                               "n_obsvar" = max(obsVariances),  
                               "n_pos" = n_pos,
                               "col_indx_pos" = col_indx_pos,
                               "row_indx_pos" = row_indx_pos,
                               "y_int"=round(y)),
                   pars = c("pred","log_lik", "alpha_0","alpha_1","beta_CI"), chains = 3,
                   iter = 20000, thin = 3, cores=3)
fit8<-rstan::stan(file='models/full dd precip CI.stan',
                   data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=CImat,
                               "states"=states, "S" = max(states), "W"=W,
                               "wetdry"=wetdry,
                               "obsVariances"=obsVariances,
                               "n_obsvar" = max(obsVariances),  
                               "n_pos" = n_pos,
                               "col_indx_pos" = col_indx_pos,
                               "row_indx_pos" = row_indx_pos,
                               "y_int"=round(y)),
                   pars = c("pred","log_lik", "alpha_0","alpha_1","beta_precip","beta_CI"), chains = 3,
                   iter = 20000, thin = 3, cores=3)
fit9<-rstan::stan(file='models/full dd temp CI.stan',
                   data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=CImat,
                               "states"=states, "S" = max(states), "W"=W,
                               "wetdry"=wetdry,
                               "obsVariances"=obsVariances,
                               "n_obsvar" = max(obsVariances),  
                               "n_pos" = n_pos,
                               "col_indx_pos" = col_indx_pos,
                               "row_indx_pos" = row_indx_pos,
                               "y_int"=round(y)),
                   pars = c("pred","log_lik", "alpha_0","alpha_1","beta_CI","beta_temp"), chains = 3,
                   iter = 20000, thin = 3, cores=3)


fit1H<-rstan::stan(file='models/all same.stan',
                  data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=HCImat,
                              "states"=states, "S" = max(states), "W"=W,
                              "wetdry"=wetdry,
                              "obsVariances"=obsVariances,
                              "n_obsvar" = max(obsVariances),  
                              "n_pos" = n_pos,
                              "col_indx_pos" = col_indx_pos,
                              "row_indx_pos" = row_indx_pos,
                              "y_int"=round(y)),
                  pars = c("pred","log_lik", "alpha_0","alpha_1","alpha_2","beta_precip","beta_CI","beta_temp","beta_pt","sigma_obs"), chains = 3,
                  iter = 20000, thin = 3, cores=3)

fit2H<-rstan::stan(file='models/full dd temp precip.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=HCImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","alpha_2","beta_precip","beta_temp"), chains = 3,
                    iter = 20000, thin = 3, cores=3)
fit3H<-rstan::stan(file='models/full model no pt interaction.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=HCImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","alpha_2","beta_precip","beta_CI","beta_temp"), chains = 3,
                    iter = 20000, thin = 3, cores=3)

fit4H<-rstan::stan(file='models/full dd temp.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=HCImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","alpha_2","beta_temp"), chains = 3,
                    iter = 20000, thin = 3, cores=3)

fit5H<-rstan::stan(file='models/full dd precip.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=HCImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","alpha_2","beta_precip"), chains = 3,
                    iter = 20000, thin =5, cores=3)
fit6H<-rstan::stan(file='models/full model no dd2.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=HCImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","beta_precip","beta_CI","beta_temp","beta_pt"), chains = 3,
                    iter = 20000, thin = 3, cores=3)

fit7H<-rstan::stan(file='models/full dd CI.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=HCImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","beta_CI"), chains = 3,
                    iter = 20000, thin = 3, cores=3)
fit8H<-rstan::stan(file='models/full dd precip CI.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=HCImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","beta_precip","beta_CI"), chains = 3,
                    iter = 20000, thin = 3, cores=3)
fit9H<-rstan::stan(file='models/full dd temp CI.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=HCImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","beta_CI","beta_temp"), chains = 3,
                    iter = 20000, thin = 3, cores=3)
fit1E<-rstan::stan(file='models/all same.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=ECImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","alpha_2","beta_precip","beta_CI","beta_temp","beta_pt"), chains = 3,
                    iter = 20000, thin = 3, cores=3)

fit2E<-rstan::stan(file='models/full dd temp precip.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=ECImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","alpha_2","beta_precip","beta_temp"), chains = 3,
                    iter = 20000, thin = 3, cores=3)
fit3E<-rstan::stan(file='models/full model no pt interaction.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=ECImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","alpha_2","beta_precip","beta_CI","beta_temp"), chains = 3,
                    iter = 20000, thin = 3, cores=3)

fit4E<-rstan::stan(file='models/full dd temp.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=ECImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","alpha_2","beta_temp"), chains = 3,
                    iter = 20000, thin = 3, cores=3)

fit5E<-rstan::stan(file='models/full dd precip.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=ECImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","alpha_2","beta_precip"), chains = 3,
                    iter = 20000, thin =5, cores=3)
fit6E<-rstan::stan(file='models/full model no dd2.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=ECImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","beta_precip","beta_CI","beta_temp","beta_pt"), chains = 3,
                    iter = 20000, thin = 3, cores=3)

fit7E<-rstan::stan(file='models/full dd CI.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=ECImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","beta_CI"), chains = 3,
                    iter = 20000, thin = 3, cores=3)
fit8E<-rstan::stan(file='models/full dd precip CI.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=ECImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","beta_precip","beta_CI"), chains = 3,
                    iter = 20000, thin = 3, cores=3)
fit9E<-rstan::stan(file='models/full dd temp CI.stan',
                    data = list("N"=N,"M"=M, "y"=y, 'precip'=precipscaled,'temp'=tempscaled, 'CI'=ECImat,
                                "states"=states, "S" = max(states), "W"=W,
                                "wetdry"=wetdry,
                                "obsVariances"=obsVariances,
                                "n_obsvar" = max(obsVariances),  
                                "n_pos" = n_pos,
                                "col_indx_pos" = col_indx_pos,
                                "row_indx_pos" = row_indx_pos,
                                "y_int"=round(y)),
                    pars = c("pred","log_lik", "alpha_0","alpha_1","beta_CI","beta_temp"), chains = 3,
                    iter = 20000, thin = 3, cores=3)





posterior1 <- as.matrix(fit1)
mcmc_areas(posterior1,
           pars = c("alpha_0","alpha_1",'alpha_2',"beta_precip",'beta_temp',"beta_pt",'beta_CI'),
           prob = 0.90, prob_outer = 1,) +scale_y_discrete(labels=c(expression(alpha[0]),expression(alpha[1]),expression(alpha[2]),expression(beta[precip]),expression(beta[temp]),expression(beta[temp%*%precip]),expression(beta[connectivity])))+theme(plot.title = element_text(hjust = 0.5,),text=element_text(size=20))





log_lik1 <-extract_log_lik(fit1, merge_chains = FALSE)
nrow(log_lik1)
waic1<-loo::waic(log_lik1)
rel_n_eff1 <- relative_eff(exp(log_lik1), chain_id = 1:3)
loo1<-loo(log_lik1, r_eff = rel_n_eff1, cores = 2)

log_lik2 <-extract_log_lik(fit2, merge_chains = FALSE)
waic2<-loo::waic(log_lik2)
rel_n_eff2 <- relative_eff(exp(log_lik2))
loo2<-loo(log_lik2, r_eff = rel_n_eff2, cores = 2)

log_lik3 <-extract_log_lik(fit3, merge_chains = FALSE)
waic3<-loo::waic(log_lik3)
rel_n_eff3 <- relative_eff(exp(log_lik3))
loo3<-loo(log_lik3, r_eff = rel_n_eff3, cores = 2)

log_lik4 <-extract_log_lik(fit4, merge_chains = FALSE)
waic4<-loo::waic(log_lik4)
rel_n_eff4 <- relative_eff(exp(log_lik4))
loo4<-loo(log_lik4, r_eff = rel_n_eff4, cores = 2)

log_lik5 <-extract_log_lik(fit5, merge_chains = FALSE)
waic5<-loo::waic(log_lik5)
rel_n_eff5 <- relative_eff(exp(log_lik5))
loo5<-loo(log_lik5, r_eff = rel_n_eff5, cores = 2)

log_lik6 <-extract_log_lik(fit6, merge_chains = FALSE)
waic6<-loo::waic(log_lik6)
rel_n_eff6 <- relative_eff(exp(log_lik6))
loo6<-loo(log_lik6, r_eff = rel_n_eff6, cores = 2)

log_lik7 <-extract_log_lik(fit7, merge_chains = FALSE)
waic7<-loo::waic(log_lik7)
rel_n_eff7 <- relative_eff(exp(log_lik7))
loo7<-loo(log_lik7, r_eff = rel_n_eff7, cores = 2)

log_lik8 <-extract_log_lik(fit8, merge_chains = FALSE)
waic8<-loo::waic(log_lik8)
rel_n_eff8 <- relative_eff(exp(log_lik8))
loo8<-loo(log_lik8, r_eff = rel_n_eff8, cores = 2,save_psis = T)


log_lik9 <-extract_log_lik(fit9, merge_chains = FALSE)
waic9<-loo::waic(log_lik9)
rel_n_eff9 <- relative_eff(exp(log_lik9))
loo9<-loo(log_lik9, r_eff = rel_n_eff9, cores = 2,save_psis = T)



log_lik1H <-extract_log_lik(fit1H, merge_chains = FALSE)
waic1H<-loo::waic(log_lik1H)
rel_n_eff1H <- relative_eff(exp(log_lik1H), chain_id = 1:3)
loo1H<-loo(log_lik1H, r_eff = rel_n_eff1H, cores = 2)

log_lik2H <-extract_log_lik(fit2H, merge_chains = FALSE)
waic2H<-loo::waic(log_lik2H)
rel_n_eff2H <- relative_eff(exp(log_lik2H))
loo2H<-loo(log_lik2H, r_eff = rel_n_eff2H, cores = 2)

log_lik3H <-extract_log_lik(fit3H, merge_chains = FALSE)
waic3H<-loo::waic(log_lik3H)
rel_n_eff3H <- relative_eff(exp(log_lik3H))
loo3H<-loo(log_lik3H, r_eff = rel_n_eff3H, cores = 2)

log_lik4H <-extract_log_lik(fit4H, merge_chains = FALSE)
waic4H<-loo::waic(log_lik4H)
rel_n_eff4H <- relative_eff(exp(log_lik4H))
loo4H<-loo(log_lik4H, r_eff = rel_n_eff4H, cores = 2)

log_lik5H <-extract_log_lik(fit5H, merge_chains = FALSE)
waic5H<-loo::waic(log_lik5H)
rel_n_eff5H <- relative_eff(exp(log_lik5H))
loo5H<-loo(log_lik5H, r_eff = rel_n_eff5H, cores = 2)

log_lik6H <-extract_log_lik(fit6H, merge_chains = FALSE)
waic6H<-loo::waic(log_lik6H)
rel_n_eff6H <- relative_eff(exp(log_lik6H))
loo6H<-loo(log_lik6H, r_eff = rel_n_eff6H, cores = 2)

log_lik7H <-extract_log_lik(fit7H, merge_chains = FALSE)
waic7H<-loo::waic(log_lik7H)
rel_n_eff7H <- relative_eff(exp(log_lik7H))
loo7H<-loo(log_lik7H, r_eff = rel_n_eff7H, cores = 2)

log_lik8H <-extract_log_lik(fit8H, merge_chains = FALSE)
waic8H<-loo::waic(log_lik8H)
rel_n_eff8H <- relative_eff(exp(log_lik8H))
loo8H<-loo(log_lik8H, r_eff = rel_n_eff8H, cores = 2)

log_lik9H <-extract_log_lik(fit9H, merge_chains = FALSE)
waic9H<-loo::waic(log_lik9H)
rel_n_eff9H <- relative_eff(exp(log_lik9H))
loo9H<-loo(log_lik9H, r_eff = rel_n_eff9H, cores = 2)

log_lik1E <-extract_log_lik(fit1E, merge_chains = FALSE)
waic1E<-loo::waic(log_lik1E)
rel_n_eff1E <- relative_eff(exp(log_lik1E), chain_id = 1:3)
loo1E<-loo(log_lik1E, r_eff = rel_n_eff1E, cores = 2)

log_lik2E <-extract_log_lik(fit2E, merge_chains = FALSE)
waic2E<-loo::waic(log_lik2E)
rel_n_eff2E <- relative_eff(exp(log_lik2E))
loo2E<-loo(log_lik2E, r_eff = rel_n_eff2E, cores = 2)

log_lik3E <-extract_log_lik(fit3E, merge_chains = FALSE)
waic3E<-loo::waic(log_lik3E)
rel_n_eff3E <- relative_eff(exp(log_lik3E))
loo3E<-loo(log_lik3E, r_eff = rel_n_eff3E, cores = 2)

log_lik4E <-extract_log_lik(fit4E, merge_chains = FALSE)
waic4E<-loo::waic(log_lik4E)
rel_n_eff4E <- relative_eff(exp(log_lik4E))
loo4E<-loo(log_lik4E, r_eff = rel_n_eff4E, cores = 2)

log_lik5E <-extract_log_lik(fit5E, merge_chains = FALSE)
waic5E<-loo::waic(log_lik5E)
rel_n_eff5E <- relative_eff(exp(log_lik5E))
loo5E<-loo(log_lik5E, r_eff = rel_n_eff5E, cores = 2)

log_lik6E <-extract_log_lik(fit6E, merge_chains = FALSE)
waic6E<-loo::waic(log_lik6E)
rel_n_eff6E <- relative_eff(exp(log_lik6E))
loo6E<-loo(log_lik6E, r_eff = rel_n_eff6E, cores = 2)

log_lik7E <-extract_log_lik(fit7E, merge_chains = FALSE)
waic7E<-loo::waic(log_lik7E)
rel_n_eff7E <- relative_eff(exp(log_lik7E))
loo7E<-loo(log_lik7E, r_eff = rel_n_eff7E, cores = 2)

log_lik8E <-extract_log_lik(fit8E, merge_chains = FALSE)
waic8E<-loo::waic(log_lik8E)
rel_n_eff8E <- relative_eff(exp(log_lik8E))
loo8E<-loo(log_lik8E, r_eff = rel_n_eff8E, cores = 2)

log_lik9E <-extract_log_lik(fit9E, merge_chains = FALSE)
waic9E<-loo::waic(log_lik9E)
rel_n_eff9E <- relative_eff(exp(log_lik9E))
loo9E<-loo(log_lik9E, r_eff = rel_n_eff9E, cores = 2)



###Selection of base model


comp1.1 <- loo_compare(waic1,waic2,waic3,waic4,waic5,waic6,waic7,waic8,waic9)

comp2.1 <- loo_compare(loo1,loo2,loo3,loo4,loo5,loo6,loo6,loo8,loo9)

print(comp1.1, digits = 2)

print(comp2.1, digits = 2)

comptable1<-data.frame(waic=comp1.1[,7])
comptable1$delta_waic<-comptable1$waic-comptable1$waic[1]



comp1.1H <- loo_compare(waic1H,waic2H,waic3H,waic4H,waic5H,waic6H,waic7H,waic8H,waic9H)
comp1.1H.1<-as.data.frame(comp1.1H)
str(comp1.1H.1)
comp1.1H.2<-comp1.1H.1[match(rownames(comptable1),rownames(comp1.1H.1)),]
comp1.1H.2

comp1.1E <- loo_compare(waic1E,waic2E,waic3E,waic4E,waic5E,waic6E,waic7E,waic8E,waic9E)


comp1.1E.1<-as.data.frame(comp1.1E)
str(comp1.1E.1)
comp1.1E.2<-comp1.1E.1[match(rownames(comptable1),rownames(comp1.1E.1)),]
comp1.1E.2


comptable1$waicH<-comp1.1H.2[,7]
comptable1$dH<-comptable1$waicH-comptable1$waic
comptable1$waicE<-comp1.1E.2[,7]
comptable1$dE<-comptable1$waicE-comptable1$waic
comptable1$model<-rownames(comptable1)
print(comptable1)

write.table(comptable1,"modelcompstruct3.csv")


#### Posterior predicitive simulation of interaction between temp and precip

##Fit 13
pfit13<-rstan::extract(fit1)
precips<-seq(from=min(precipscaled),to=max(precipscaled),by=0.01)
ppred1<-data.frame(Precipitation=precips,temp='Mean',Predicted=0,Draw=0)
ppred2<-data.frame(Precipitation=precips,temp='Min',Predicted=0,Draw=0)
ppred3<-data.frame(Precipitation=precips,temp='Max',Predicted=0,Draw=0)
temps<-seq(from=min(tempscaled),to=max(precipscaled),by=0.01)
tpred1<-data.frame(Temperature=temps,precip='Mean',Predicted=0)
tpred2<-data.frame(Temperature=temps,precip='Min',Predicted=0)
tpred3<-data.frame(Temperature=temps,precip='Max',Predicted=0)



preciplines<-function(
  alpha_0=sample(pfit13$alpha_0,size=1),
  alpha_1=sample(pfit13$alpha_1,size=1),
  alpha_2=sample(pfit13$alpha_2,size=1),
  beta_precip=sample(pfit13$beta_precip,size=1),
  beta_temp=sample(pfit13$beta_temp,size=1),
  beta_pt=sample(pfit13$beta_pt,size=1),
  beta_CI=sample(pfit13$beta_CI,size=1),
  tbar=mean(tempscaled),
  tmin=min(tempscaled),
  tmax=max(tempscaled),
  CIbar=mean(CImat),
  xbar=mean(y),
  precips=seq(from=min(precipscaled),to=max(precipscaled),by=0.01),
  draw=1
)
{
  for(i in 1:length(precips)){
    ppred1[i,3] = alpha_0 + alpha_1*xbar + alpha_2*xbar + beta_precip*precips[i] + beta_temp*tbar + beta_CI*CIbar + beta_pt*precips[i]*tbar
    ppred2[i,3] = alpha_0 + alpha_1*xbar + alpha_2*xbar + beta_precip*precips[i] + beta_temp*tmin + beta_CI*CIbar + beta_pt*precips[i]*tmin
    ppred3[i,3] = alpha_0 + alpha_1*xbar + alpha_2*xbar + beta_precip*precips[i] + beta_temp*tmax + beta_CI*CIbar + beta_pt*precips[i]*tmax
  }
  ppred<-rbind(ppred1,ppred2,ppred3)
  ppred<-as.data.frame(ppred)
  ppred$Draw<-draw
  return(ppred)
}
l1<-preciplines()
str(l1)
reps<-50
lps<-length(precips)*3 
lps
length(ppred1)
str(ppred1)


output<-data.frame(Precipitation=rep(precips,reps*3),temp='Max',Predicted=0,Draw=0)
str(output)
for(i in 1:reps){
  start<-(i-1)*lps+1
  end<-i*lps
  output[start:end,]<-preciplines(draw=i)
  
}

precipintplot<-ggplot(output,aes(x=Precipitation,y=Predicted))+geom_line(aes(group=Draw),alpha=0.2)+facet_grid(.~temp)+theme_classic()+ylab("Predicted Density")+stat_smooth(method='lm',se=F,color="black",size=2)

##with temp
templines<-function(
  alpha_0=sample(pfit13$alpha_0,size=1),
  alpha_1=sample(pfit13$alpha_1,size=1),
  alpha_2=sample(pfit13$alpha_2,size=1),
  beta_precip=sample(pfit13$beta_precip,size=1),
  beta_temp=sample(pfit13$beta_temp,size=1),
  beta_pt=sample(pfit13$beta_pt,size=1),
  beta_CI=sample(pfit13$beta_CI,size=1),
  pbar=mean(precipscaled),
  pmin=min(precipscaled),
  pmax=max(precipscaled),
  CIbar=mean(CImat),
  xbar=mean(y),
  temps=seq(from=min(tempscaled),to=max(tempscaled),by=0.01),
  draw=1
)
{
  for(i in 1:length(temps)){
    tpred1[i,3] = alpha_0 + alpha_1*xbar + alpha_2*xbar + beta_precip*pbar + beta_temp*temps[i] + beta_CI*CIbar + beta_pt*pbar*temps[i]
    tpred2[i,3] = alpha_0 + alpha_1*xbar + alpha_2*xbar + beta_precip*pmin + beta_temp*temps[i] + beta_CI*CIbar + beta_pt*pmin*temps[i]
    tpred3[i,3] = alpha_0 + alpha_1*xbar + alpha_2*xbar + beta_precip*pmin + beta_temp*temps[i] + beta_CI*CIbar + beta_pt*pmax*temps[i]
  }
  
  tpred<-rbind(tpred1,tpred2,tpred3)
  tpred<-as.data.frame(tpred)
  tpred$Draw<-draw
  return(tpred)
}

temps=seq(from=min(tempscaled),to=max(tempscaled),by=0.01)
lps<-length(temps)*3
output2<-data.frame(Temperature=rep(temps,reps*3),precip='Max',Predicted=0,Draw=0)
str(output2)
for(i in 1:reps){
  start<-(i-1)*lps+1
  end<-i*lps
  output2[start:end,]<-templines(draw=i)
  
}


tempintplot<-ggplot(na.omit(output2),aes(x=Temperature,y=Predicted))+geom_line(aes(group=Draw),alpha=0.2)+facet_grid(.~precip)+theme_classic()+ylab("Predicted Density")+stat_smooth(method='lm',se=F,color="black",size=2)

plot_grid(precipintplot,tempintplot,nrow=2,labels=c("A","B"))

