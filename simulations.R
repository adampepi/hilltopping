rm(list=ls())
###BMR metapopulation dynamics simulation
library(ncf)
library(ggplot2)
library(tidyverse)
library(cowplot)
setwd("~/Documents/Research/dissertation/12 sites analyses/hilltopping")


sites<-read.csv("12sitesseriesdensity.csv")
str(sites)
sites
sitesdata<-log(sites[,5:18]*15+1)##Convert to  count scale
sitesdata
sites2<-read.csv("12sitesserieslong4.csv")##Sites are organized different in this file than
str(sites2)
sites2

sitesdatamat<-read.csv('sitesdatamat.csv')

weather<-read.csv('weather.csv')

wetdry<-read.csv('wetdry.csv')

distmat<-read.csv('distmat.csv')
distmat
ncol(distmat)
distmat<-distmat[,2:21]
distmat2<-distmat

####Version with hilltop connectivity  


hillsim<-function(Hsd=183.9366,##Scaling parameters for hilltop connectivity
                  Hmean=-191.8851,##Scaling parameters for hilltop connectivity
                  CImean=-107.9676,##Scaling parameters for connectivity
                  CIsd=98.34425,##Scaling parameters for connectivity
                  period=14,
                  temp=scale(weather$avgtm1),
                  precip=scale(weather$preciptm1),
                  alpha_0=2.52,
                  alpha_1=0.17,
                  alpha_2=c(-0.63,-0.63),
                  beta_precip=0.7,
                  beta_temp=0.08,
                  beta_pt=0.68,
                  beta_CI=0.44,
                  nsites=12,
                  wetdry=c(2,2,2,2,2,2,2,2,1,1,1,1,1),
                  init1=sitesdata[1:12,1],
                  init2=sitesdata[1:12,2],
                  alpha=1,
                  distmat=distmat2
                  
){
  pops<-matrix(0,ncol=period,nrow=nsites)
  CIpatch<-matrix(nrow=nsites,ncol=period)
  CIhills<-matrix(nrow=8,ncol=period)
  pops[,1]<-init1
  pops[,2]<-init2
  for(t in 3:period){
    
    for(i in 1:8){
      ##Loop to calculate connectivity for each hilltop
      exps<-rep(0,12)
      for(q in 1:nsites){
        ##Each patch
        dist<-distmat[12+i,q]
        area<-pops[q,t-1]
        exps[q]<-exp(-alpha*dist)*area
        
      }
      CIhills[i,t-1]<-sum(exps)
    }
    CIhills[,t-1]<-log(CIhills[,t-1])
    CIhills[is.nan(CIhills[,t-1]),t-1]<-min(CIhills)
    CIhills[,t-1]<- (CIhills[,t-1]-Hmean)/Hsd+2
    for(r in 1:nsites){
      exps<-rep(0,8)
      for(j in 1:8){
        ##Using each hilltop
        dist<-distmat[r,12+j]
        area<-CIhills[j,t-1]
        exps[j]<-exp(-alpha*dist)*area
      }
      CIpatch[r,t-1]<-sum(exps)
    }
    CIpatch[,t-1]<-log(CIpatch[,t-1]+7.020758e-229)
    CIpatch[is.nan(CIpatch[,t-1]),t-1]<--2.5
    CIpatch[,t-1]<-(CIpatch[,t-1]-CImean)/CIsd
    
    
    for(s in 1:nsites){
      
      r.mean<- alpha_0 + alpha_1*pops[s,t-1] + alpha_2[wetdry[s]]*pops[s,t-2] + beta_precip*precip[t-1] + 
        beta_temp*temp[t-1] + beta_CI*CIpatch[s,t-1] + beta_pt*precip[t-1]*temp[t-1] ##Simulation from time series model
      
      pops[s,t]<- log(rpois(1,exp(r.mean))+1)  #Data generation, using Poisson process but logged
    } 
  }
  return(pops)
}
  CIpatch
 
pop1<-hillsim()
matplot(exp(t(pop1)),type="l",log='y')
pop1


####Version with hilltop & elevation connectivity  


elevationsim<-function(Hsd=183.9366,##Scaling paratemeters for hilltop connectivity
                  Hmean=-191.8851,##Scaling paratemeters for hilltop connectivity
                  CImean=-107.9676,##Scaling paratemeters for connectivity
                  CIsd=98.34425,##Scaling paratemeters for connectivity
                  period=15,
                  temp=scale(weather$avgtm1),
                  precip=scale(weather$preciptm1),
                  alpha_0=2.52,
                  alpha_1=0.17,
                  alpha_2=c(-0.63,-0.63),
                  beta_precip=0.7,
                  beta_temp=0.08,
                  beta_pt=0.68,
                  beta_CI=0.44,
                  nsites=12,
                  wetdry=c(2,2,2,2,2,2,2,2,1,1,1,1,1),
                  init1=sitesdata[1:12,1],
                  init2=sitesdata[1:12,2],
                  alpha=1,
                  distmat=distmat2
                  
){
  pops<-matrix(0,ncol=period,nrow=nsites)
  CIpatch<-matrix(nrow=nsites,ncol=period)
  CIhills<-matrix(nrow=8,ncol=period)
  elevations<-c(40,40,90,80,40,140,190,40)
  pops[,1]<-init1
  pops[,2]<-init2
  for(t in 3:period){
    
    for(i in 1:8){
      ##Loop to calculate connectivity for each hilltop
      exps<-rep(0,12)
      for(q in 1:nsites){
        ##Each patch
        dist<-distmat[12+i,q]
        area<-pops[q,t-1]
        exps[q]<-exp(-alpha*dist)*area
        
      }
      CIhills[i,t-1]<-sum(exps)
    }
    CIhills[,t-1]<-log(CIhills[,t-1])
    CIhills[is.nan(CIhills[,t-1]),t-1]<-min(CIhills)
    CIhills[,t-1]<- (CIhills[,t-1]-Hmean)/Hsd+2
    for(r in 1:nsites){
      exps<-rep(0,8)
      for(j in 1:8){
        ##Using each hilltop
        dist<-distmat[r,12+j]
        area<--0.89+1.27*CIhills[j,t-1]+0.02*elevations[j]
        exps[j]<-exp(-alpha*dist)*area
      }
      CIpatch[r,t-1]<-sum(exps)
    }
    CIpatch[,t-1]<-log(CIpatch[,t-1]+7.020758e-229)
    CIpatch[is.nan(CIpatch[,t-1]),t-1]<--2.5
    CIpatch[,t-1]<-(CIpatch[,t-1]-CImean)/CIsd
    
    
    for(s in 1:nsites){
      
      r.mean<- alpha_0 + alpha_1*pops[s,t-1] + alpha_2[wetdry[s]]*pops[s,t-2] + beta_precip*precip[t-1] + 
        beta_temp*temp[t-1] + beta_CI*CIpatch[s,t-1] + beta_pt*precip[t-1]*temp[t-1]
      
      pops[s,t]<- log(rpois(1,exp(r.mean))+1)
    } 
  }
  return(pops)
}
CIpatch

pop1<-elevationsim()
matplot(exp(t(pop1)),type="l",log='y')
pop1
####Version with Hanski connectivity  
##Use the same values or not?


hanskisim<-function(
  HCIsd=180.345, ##Scaling paratemeters for hilltop connectivity
  HCImean=-83.99448, ##Scaling paratemeters for hilltop connectivity
  period=14,
  temp=scale(weather$avgtm1),
  precip=scale(weather$preciptm1),
  alpha_0=2.52,
  alpha_1=0.17,
  alpha_2=c(-0.63,-0.63),
  beta_precip=0.7,
  beta_temp=0.08,
  beta_pt=0.68,
  beta_CI=0.44,
  nsites=12,
  wetdry=c(2,2,2,2,2,2,2,2,1,1,1,1,1),
  alpha=1,
  init1=sitesdata[1:12,1],
  init2=sitesdata[1:12,2],
  distmat=distmat2
){
  pops2<-matrix(0,ncol=period,nrow=nsites)
  
  CIpatch<-matrix(nrow=nsites,ncol=period)
  pops2[,1]<-init1
  pops2[,2]<-init2
 
  for(t in 3:period){
    
    for(r in 1:nsites){
      ##Loop to calculate connectivity for each patch
      
      exps<-rep(0,8)
      for(j in 1:nsites){
        ##Each patch
        dist<-distmat[j,r]
        area<-pops2[j,t-1]
        exps[j]<-exp(-alpha*dist)*area*20
        
      }
      CIpatch[r,t-1]<-sum(exps)
    }
    CIpatch[,t-1]<-log(CIpatch[,t-1]+7.015732e-321)
    CIpatch[,t-1]<-(CIpatch[,t-1]-HCImean)/HCIsd
    
    for(s in 1:nsites){
      
      r.mean<- alpha_0 + alpha_1*pops2[s,t-1] + alpha_2[wetdry[s]]*pops2[s,t-2] + beta_precip*precip[t-1] + beta_temp*temp[t-1] + beta_CI*CIpatch[s,t-1] + beta_pt*precip[t-1]*temp[t-1]
      
      pops2[s,t]<- log(rpois(1,exp(r.mean))+1)
    }
  }
  return(pops2)
}
pop1<-hanskisim()
matplot(exp(t(pop1)),type="l",log='y')
pop1

##Running simulations
reps<-10000
V1<-vector(length=reps)
V2<-vector(length=reps)
V3<-vector(length=reps)
for(i in 1:reps){
  pop2<-hillsim()
  V1[i]<-sum(pop2==0)
  
  
  empty<-which(pop2[,1:ncol(pop2)-1]==0,arr.ind = T)
  empty2<-empty
  empty3<-empty
  empty2[,2]<-empty[,2]-1
  empty3[,2]<-empty[,2]+1
  
  ntm1exct<-pop2[empty2]
  ntp1col<-pop2[empty3]
  extinct<-sum(pop2==0)-sum(ntm1exct==0)
  col<-sum(ntp1col>0)
  V2[i]<- extinct
  V3[i]<- col
}
mean(V1)
mean(V2)
mean(V3)

sum(pop2==0)

V1.1<-vector(length=reps)

V2.1<-vector(length=reps)
V3.1<-vector(length=reps)
for(i in 1:reps){
  pop2<-elevationsim()
  V1.1[i]<-sum(pop2==0)
  
  empty<-which(pop2[,1:ncol(pop2)-1]==0,arr.ind = T)
  empty2<-empty
  empty3<-empty
  empty2[,2]<-empty[,2]-1
  empty3[,2]<-empty[,2]+1
  
  ntm1exct<-pop2[empty2]
  ntp1col<-pop2[empty3]
  extinct<-sum(pop2==0)-sum(ntm1exct==0)
  col<-sum(ntp1col>0)
  V2.1[i]<- extinct
  V3.1[i]<- col
}
mean(V1.1)
mean(V2.1)
mean(V3.1)



reps<-10000
V4<-vector(length=reps)
V5<-vector(length=reps)
V6<-vector(length=reps)
for(i in 1:reps){
  pop1<-hanskisim()
  
  empty<-which(pop1[,1:ncol(pop1)-1]==0,arr.ind = T)
  empty2<-empty
  empty3<-empty
  empty2[,2]<-empty[,2]-1
  empty3[,2]<-empty[,2]+1
  
  ntm1exct<-pop1[empty2]
  ntp1col<-pop1[empty3]
  extinct<-sum(pop1==0)-sum(ntm1exct==0)
  col<-sum(ntp1col>0)
  V5[i]<- extinct
  V6[i]<- col
  V4[i]<-sum(pop1==0)
}
mean(V4)
mean(V5)
mean(V6)
reps<-10000
V7<-vector(length=reps)
for(i in 1:reps){
  V7[i]<- mSynch(x=hanskisim())$real[1]
}
V7
mean(V7)

reps<-10000
V8<-vector(length=reps)
for(i in 1:reps){
  V8[i]<-mSynch(x=hillsim())$real[1]
}
mean(V8)
V8 
reps<-10000
V9<-vector(length=reps)
for(i in 1:reps){
  V9[i]<-mSynch(x=elevationsim())$real[1]
}
mean(V9)
V9

##Collect results into data frame

hanskidat<-data.frame(occupancy=1-V4/(12*15),extinction=V5,colonization=V6,synchrony=V7,model='Hanski')#Occupancy calculated as 1-unoccupied/total sites*years
hilldat<-data.frame(occupancy=1-V1/(12*15),extinction=V2,colonization=V3,synchrony=V8,model='hill')
elevdat<-data.frame(occupancy=1-V1.1/(12*15),extinction=V2.1,colonization=V3.1,synchrony=V9,model='elevation')

str(hanskidat)
str(hilldat)
simdat<-rbind(hanskidat,hilldat,elevdat)
simdat<-as.data.frame(simdat)
str(simdat)
write.csv(simdat,"simdat.csv")
simdat<-read.csv('simdat.csv')
str(simdat)
levels(simdat$model)
colMeans(simdat[simdat$model=="hill",2:5])
quantile(simdat$occupancy[simdat$model=="hill"],c(0.25,0.5,0.75))
quantile(simdat$extinction[simdat$model=="hill"],c(0.25,0.5,0.75))
quantile(simdat$colonization[simdat$model=="hill"],c(0.25,0.5,0.75))
quantile(simdat$synchrony[simdat$model=="hill"],c(0.25,0.5,0.75))

quantile(simdat$occupancy[simdat$model=="Hanski"],c(0.25,0.5,0.75))
quantile(simdat$extinction[simdat$model=="Hanski"],c(0.25,0.5,0.75))
quantile(simdat$colonization[simdat$model=="Hanski"],c(0.25,0.5,0.75))
quantile(simdat$synchrony[simdat$model=="Hanski"],c(0.25,0.5,0.75))

quantile(simdat$occupancy[simdat$model=="elevation"],c(0.25,0.5,0.75))
quantile(simdat$extinction[simdat$model=="elevation"],c(0.25,0.5,0.75))
quantile(simdat$colonization[simdat$model=="elevation"],c(0.25,0.5,0.75))
quantile(simdat$synchrony[simdat$model=="elevation"],c(0.25,0.5,0.75))

medians<-simdat%>%group_by(model)%>%summarise(occupancy=median(occupancy), colonization=median(colonization),extinction=median(extinction),synchrony=median(synchrony))
medians

###Get true values
sitesdata
popreal<-sitesdata[1:12,]
occupancyr<-sum(popreal==0)
occupancyrr<-1-occupancyr/(12*15)


empty<-which(popreal[,1:ncol(popreal)-1]==0,arr.ind = T)
empty2<-empty
empty3<-empty
empty2[,2]<-empty[,2]-1
empty3[,2]<-empty[,2]+1

ntm1exct<-popreal[empty2]
ntp1col<-popreal[empty3]
extinctr<-sum(popreal==0)-sum(ntm1exct==0)
colr<-sum(ntp1col>0)
synchr<- mSynch(x=as.matrix(popreal))$real[1]
real.values<-data.frame(occupancy=occupancyrr,colonization=colr,extinction=extinctr,synchrony=synchr)
real.values
medians[4,1]<-"actual"
medians[4,2:5]<-real.values


##Make plots

occupancyplot<-ggplot(data=simdat,aes(x=occupancy,color=model))+geom_density(adjust=1.5,lwd=1)+geom_vline(data=medians,aes(xintercept=occupancy,color=model),lty=2,lwd=1)+scale_color_viridis_d()+theme_classic()+xlab("Mean occupancy (proportion)")+ylab('Kernel density')+theme(legend.position = 'top' )+labs(color="Model:")+xlim(0.7,0.99)
occupancyplot

medians
colonizationplot<-ggplot(data=simdat,aes(x=colonization,color=model))+geom_density(adjust=1.5,lwd=1)+geom_vline(data=medians,aes(xintercept=colonization,color=model),lty=2,lwd=1)+scale_color_viridis_d()+theme_classic()+xlab("Mean colonization events")+ylab('Kernel density')+theme(legend.position = 'none')+xlim(2,23)

extinctionplot<-ggplot(data=simdat,aes(x=extinction,color=model))+geom_density(adjust=1.5,lwd=1)+geom_vline(data=medians,aes(xintercept=extinction,color=model),lty=2,lwd=1)+scale_color_viridis_d()+theme_classic()+xlab("Mean extinction events")+ylab('Kernel density')+theme(legend.position = 'none')+xlim(2,23)

extinctionplot


synchronyplot<-ggplot(data=simdat,aes(x=synchrony,color=model))+geom_density(adjust=1.5,lwd=1)+geom_vline(data=medians,aes(xintercept=synchrony,color=model),lty=2,lwd=1)+scale_color_viridis_d()+theme_classic()+xlab("Mean synchrony")+ylab('Kernel density')+theme(legend.position = 'none')+xlim(0.45,.7)

saveRDS(simdat, file = "simdat.rds")


plot_grid(occupancyplot,colonizationplot,extinctionplot,synchronyplot,labels=c("A","B",'C',"D"),nrow=4,rel_heights = c(1.25,1,1,1))



