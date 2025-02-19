# Set path
Root_Path <- 'C:/zhao088/Desk20210304/Manu.cusal abun.trait/3.Trout Lake84-/CCM_AL fish'
OPath <- paste(Root_Path,sep='')
setwd(OPath)



#####################################################
library('rEDM')
packageVersion('rEDM') #[1] '1.2.3'
library('dplyr')
library('tidyr')
library('igraph')
library('MASS')
library('NetIndices')
library('Cairo')
library('vegan')
library('Kendall') # Kendall's tau test for the convergence of CCM

library('doParallel')
library('parallel')
library('foreach')
library('glmnet')
library('psych') #trace
library('pracma')
library('matlib') #maritx

####Function for time series standardization (normalization, detrend and deseason)
# x: time series data (vector form)
# normalization: normalization to zero mean & unit variance 
# dseason: logic argument for Deseasonalization by monthly mean
# season_sd: Deseasonalization by both monthly mean and monthly standard deviation
# sea: The period of seasonality (e.g. 12 months for monthly data)
# dtrend: logic argument for detrend
# dTtype: Type of detrend by first difference (first) or linear regression (linear)
nomz=function(x, normalization=T){
  x=as.numeric(x)
  xt=x
  # Normalization (zero mean & unity variance)
  if(normalization==T){xt=(xt-mean(xt,na.rm=T))/sd(xt,na.rm=T)}
  return(xt)
  
}

####### Function for generating lag time series 
laf=function(x,y,lagf){
  n <- NROW(x)
  x.t=x;y.t=y
  if(lagf<=0){x.t=x.t[(1-lagf):n];y.t=y.t[1:(n+lagf)]} # if lagf<0, y is leading
  if(lagf>0){x.t=x.t[1:(n-lagf)];y.t=y.t[(1+lagf):n]}  # if lagf>0, x is leading           
  return(cbind(x.t,y.t))
}

###Part 1 <80% zoo
new_abund<-read.csv("b.AL.fish.caught_per_effort.csv")
head(new_abund)
str(new_abund)
new_abund<-new_abund[,3:ncol(new_abund)]#has date
target_row=NULL
ii=2
for (ii in 2:ncol(new_abund)) { # no date
  if(sum(new_abund[,ii]==0)/nrow(new_abund)<0.6) {target_row=c(target_row,ii)}
}
print(target_row)

# then use these colnames to select biomass
new_biomass<-read.csv("b.AL.fish.length.csv")
new_biomass_matched=new_biomass[,c( colnames(new_abund[,target_row]) )]# find matched column names
# remove the rows that has NA values 
row_no_NA=NULL
for (i2 in 1:ncol(new_biomass_matched)){ if(!anyNA(new_biomass_matched[,i2])){row_no_NA=c(row_no_NA,i2)} }
row_no_NA

# bind
finals=as.data.frame( cbind(new_abund$year,new_abund[,target_row[row_no_NA]],new_biomass_matched[,row_no_NA]) )
#rename
colnames(finals)<-c("year",colnames(finals)[2:ncol(finals)])


# bind all
alga_zo=as.data.frame( cbind(new_abund$year,
                            new_abund[,target_row[row_no_NA]], # zoo abund
                            new_biomass_matched[,row_no_NA] # zoo bioma
                            )) # zoo biomass

#rename
colnames(alga_zo)<-c("sampledate",colnames(alga_zo)[2:ncol(alga_zo)])



# part 3: combine
alga_zo<-data.frame(alga_zo[,2:ncol(alga_zo)]) # no year
dim(alga_zo) 
colnames(alga_zo)
head(alga_zo)
#
n=nrow(alga_zo)
seed=25647

##################
# Part 1: trait affects abundance
##################

# Detrend + deseason time series
sdat=data.frame(apply(alga_zo,2,nomz))
##The index for testing causal links (a total of 12 links)
indmat=matrix(0,ncol(sdat)/2,2);colnames(indmat)=c('Effect','Cause')
indmat[,1]=colnames(sdat)[1:(ncol(sdat)/2)]
indmat[,2]=colnames(sdat)[(ncol(sdat)/2+1):ncol(sdat)]
indmat



# Determine the embedding dimensions (En) in CCM by try-and-error with best hindcast (tp=-1) skill
Emax=6
En=NULL
for(i in 1:nrow(indmat)){
  E.test=NULL
  for(E.t in 2:Emax){
    cmxy.t <- ccm(sdat, E = E.t,
                  lib_column = indmat[i,1], target_column = indmat[i,2],
                  lib_sizes = n, tp=-1,random_libs = F) 
    E.test=c(E.test,mean(cmxy.t$rho))
  }
  # Select the embedding dimension that makes the model with the highest hindcast (tp=-1) predictive skill
  En=c(En,which.max(E.test)+1) 
}

i=1
j=0
################################################################################
### CCM analysis for all causality testlinks
lib_siz=sort(c(seq(11,n,5),n)) # a sequence of library size
ccmda=NULL
for(i in 1:nrow(indmat)){
  ccmda.t=NULL
  for(j in 0:-2){
    da.t=laf(sdat[,indmat[i,1]],sdat[,indmat[i,2]],lagf=j) # Varying time lags
    colnames(da.t)=indmat[i,]
    # CCM analysis cross-mapping from an effect variable to its cause
    x_xmap_y <- ccm(da.t, E = En[i], # The embedding dimension E for each link were determined in previous step
                    lib_column = indmat[i,'Effect'], target_column = indmat[i,'Cause'],
                    lib_sizes = lib_siz, tp=0,RNGseed = seed,
                    num_samples = 100,replace=F)
    
    # Take average for the predictive skill under each library size
    aveg=cbind(unique(x_xmap_y$lib_size),aggregate(x_xmap_y[,c('rho')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'],
               aggregate(x_xmap_y[,c('mae')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'],
               aggregate(x_xmap_y[,c('rmse')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'])
    ccm_mean=data.frame(lag=rep(j,nrow(aveg)),x_xmap_y[1:nrow(aveg),]);
    ccm_mean[,c('lib_size','rho','mae','rmse')]=aveg
    ccm_mean[ccm_mean[,'rho']<0,'rho']=0
    
    ###########################
    # Convergence test in CCM
    ###########################
    # Fisher's delta rho Z test
    rho.Lmax=ccm_mean$rho[which.max(ccm_mean$lib_size)]
    rho.Lmin=ccm_mean$rho[1]
    ns=min(sum(!is.na(sdat[,indmat[i,1]])),sum(!is.na(sdat[,indmat[i,2]])))
    delta_rho=rho.Lmax-rho.Lmin
    z=abs(0.5*(log((1+rho.Lmax)/(1-rho.Lmax))-log((1+rho.Lmin)/(1-rho.Lmin)))*(2/(ns-3))^-0.5)
    z.p=(1-pnorm(z))
    # Kendall's tau test
    if(length(ccm_mean$rho)>3){
      kend=MannKendall(ccm_mean$rho)
      kend.tau=kend$tau[1]
      kend.p=kend$sl[[1]]
    }else{
      kend.tau=NA
      kend.p=NA
    }
    # Compile all the testing results
    ccmda.t=rbind(ccmda.t,
                  unlist(c(ccm_mean[1,1:5],rho_Lmax=rho.Lmax,rho_Lmin=rho.Lmin,
                           Z=z,p_Z=z.p,Kendall_tau=kend.tau,Kendall_p=kend.p)))
  }
  # Select the CCM results based on predictive skills 
  ccmda=rbind(ccmda,ccmda.t[which.max(ccmda.t[,'rho_Lmax']),])
}

ccmda=data.frame(indmat,ccmda)
Convergence=ccmda$p_Z<0.05 & ccmda$Kendall_p<=0.05 & ccmda$Kendall_tau>0
(ccmda=data.frame(ccmda,Convergence))

# selecting convergence.by me
#ccmda<-ccmda[ccmda$Convergence==1,]

# Standardized linkage strength by dividing the maximal linkage strength within the system
istd.1=data.frame(system=rep('Lake Allequash'),ccmda[,c('Cause','Effect','rho_Lmax','p_Z','Kendall_tau','Kendall_p','Convergence')],
                Standardized_linkage_strength=ccmda$rho_Lmax/max(ccmda$rho_Lmax[1:length(ccmda$rho_Lmax)]))

#
istd.1$direction<-c("trait affects abundance")


#### add mean abundance and mean biomass
j=3
istd.1$abundances=istd.1$biomasss=NA
for (j in 1:nrow(istd.1)) {
  istd.1[j, "abundances"]<- mean(alga_zo[,istd.1[j,"Effect"]], na.rm = TRUE) 
  istd.1[j, "biomasss"]<- mean(0.02*alga_zo[,paste(istd.1[j,"Effect"],".1",sep="") ]^3/10^2, na.rm = TRUE) 
}
#### End to add mean abundance and mean biomass



##################
# Part 2:  abundance affects trait
##################

# Detrend + deseason time series
sdat2=data.frame(apply(alga_zo,2,nomz))
##The index for testing causal links (a total of 12 links)
indmat2=matrix(0,ncol(sdat2)/2,2);colnames(indmat2)=c('Effect','Cause')
indmat2[,1]=colnames(sdat2)[(ncol(sdat2)/2+1):ncol(sdat2)] 
indmat2[,2]=colnames(sdat2)[1:(ncol(sdat2)/2)]
indmat2



# Determine the embedding dimensions (En) in CCM by try-and-error with best hindcast (tp=-1) skill
Emax=6
En=NULL
for(i in 1:nrow(indmat2)){
  E.test=NULL
  for(E.t in 2:Emax){
    cmxy.t <- ccm(sdat2, E = E.t,
                  lib_column = indmat2[i,1], target_column = indmat2[i,2],
                  lib_sizes = n, tp=-1,random_libs = F) 
    E.test=c(E.test,mean(cmxy.t$rho))
  }
  # Select the embedding dimension that makes the model with the highest hindcast (tp=-1) predictive skill
  En=c(En,which.max(E.test)+1) 
}

i=1

################################################################################
### CCM analysis for all causality testlinks
lib_siz=sort(c(seq(11,n,5),n)) # a sequence of library size
ccmda=NULL
for(i in 1:nrow(indmat2)){
  ccmda.t2=NULL
  for(j in 0:-2){
    da.t=laf(sdat2[,indmat2[i,1]],sdat2[,indmat2[i,2]],lagf=j) # Varying time lags
    colnames(da.t)=indmat2[i,]
    # CCM analysis cross-mapping from an effect variable to its cause
    x_xmap_y <- ccm(da.t, E = En[i], # The embedding dimension E for each link were determined in previous step
                    lib_column = indmat2[i,'Effect'], target_column = indmat2[i,'Cause'],
                    lib_sizes = lib_siz, tp=0,RNGseed = seed,
                    num_samples = 100,replace=F)
    
    # Take average for the predictive skill under each library size
    aveg=cbind(unique(x_xmap_y$lib_size),aggregate(x_xmap_y[,c('rho')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'],
               aggregate(x_xmap_y[,c('mae')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'],
               aggregate(x_xmap_y[,c('rmse')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'])
    ccm_mean=data.frame(lag=rep(j,nrow(aveg)),x_xmap_y[1:nrow(aveg),]);
    ccm_mean[,c('lib_size','rho','mae','rmse')]=aveg
    ccm_mean[ccm_mean[,'rho']<0,'rho']=0
    
    ###########################
    # Convergence test in CCM
    ###########################
    # Fisher's delta rho Z test
    rho.Lmax=ccm_mean$rho[which.max(ccm_mean$lib_size)]
    rho.Lmin=ccm_mean$rho[1]
    ns=min(sum(!is.na(sdat2[,indmat2[i,1]])),sum(!is.na(sdat2[,indmat2[i,2]])))
    delta_rho=rho.Lmax-rho.Lmin
    z=abs(0.5*(log((1+rho.Lmax)/(1-rho.Lmax))-log((1+rho.Lmin)/(1-rho.Lmin)))*(2/(ns-3))^-0.5)
    z.p=(1-pnorm(z))
    # Kendall's tau test
    if(length(ccm_mean$rho)>3){
      kend=MannKendall(ccm_mean$rho)
      kend.tau=kend$tau[1]
      kend.p=kend$sl[[1]]
    }else{
      kend.tau=NA
      kend.p=NA
    }
    # Compile all the testing results
    ccmda.t2=rbind(ccmda.t2,
                   unlist(c(ccm_mean[1,1:5],rho_Lmax=rho.Lmax,rho_Lmin=rho.Lmin,
                            Z=z,p_Z=z.p,Kendall_tau=kend.tau,Kendall_p=kend.p)))
  }
  # Select the CCM results based on predictive skills 
  ccmda=rbind(ccmda,ccmda.t2[which.max(ccmda.t2[,'rho_Lmax']),])
}

ccmda=data.frame(indmat2,ccmda)
Convergence=ccmda$p_Z<0.05 & ccmda$Kendall_p<=0.05 & ccmda$Kendall_tau>0
(ccmda=data.frame(ccmda,Convergence))

print(which(Convergence=="FALSE"))

# selecting convergence.by me
#ccmda<-ccmda[ccmda$Convergence==1,] # only using convergence,not using non-convergence

# Standardized linkage strength by dividing the maximal linkage strength within the system
istd.2=data.frame(system=rep('Lake Allequash'),ccmda[,c('Cause','Effect','rho_Lmax','p_Z','Kendall_tau','Kendall_p','Convergence')],
                  Standardized_linkage_strength=ccmda$rho_Lmax/max(ccmda$rho_Lmax[1:length(ccmda$rho_Lmax)]))


#
istd.2$direction<-c("abundance affects trait")



#### add mean abundance and mean biomass
j=3
istd.2$abundances=istd.2$biomasss=NA
for (j in 1:nrow(istd.2)) {
  istd.2[j, "abundances"]<- mean(alga_zo[,istd.2[j,"Cause"]], na.rm = TRUE) 
  istd.2[j, "biomasss"]<- mean(0.02*alga_zo[,paste(istd.2[j,"Cause"],".1",sep="") ]^3/10^2, na.rm = TRUE) 
}
#### End to add mean abundance and mean biomass



# merge istd.1 istd.2
finals<-data.frame(rbind(istd.1, istd.2))
nrow(finals)
# re-Standardized linkage strength by dividing the maximal linkage strength within the system
finals$Standardized_linkage_strength=finals$rho_Lmax/max(finals$rho_Lmax[1:length(finals$rho_Lmax)])

#finals=finals[finals$Standardized_linkage_strength>=0.3,]
finals

# total
finals$trait.effect.abundance.total= sum(finals$Convergence[finals$direction=="trait affects abundance"]=="FALSE")/(nrow(finals)/2)
finals$abundance.effect.trait.total= sum(finals$Convergence[finals$direction=="abundance affects trait"]=="FALSE")/(nrow(finals)/2)


##################
# Part 3: The computation of warming rate and Oxygenation, PH, mean TDP and TDN
##################
temp.oxge<-read.csv("4.physical.seasonal.AL-temp-O2yearly.csv")
#temp.oxge=temp.oxge[112:151,]#
#
head(temp.oxge)

#wtemp,o2,O2..µmol.l.,X.PO4.3...µmol.l.,total_Dissolved_N

# temperature
library(zyp)# Theil-Sen median based estimator
temp.mon=data.frame(temp=temp.oxge$wtemp[!is.na(temp.oxge$wtemp)],
                    month=c(1:length(temp.oxge$wtemp[!is.na(temp.oxge$wtemp)])))
fit<- zyp.sen(temp~month,temp.mon)
(Warming_rate.t=fit$coefficients[2]*1)
finals$Warming_rate<-Warming_rate.t
finals$temp_mean<-mean(temp.oxge$wtemp,na.rm=TRUE)
#skewness
library(moments)
# Calculate skewness
finals$skewness_temp <- as.numeric(skewness(temp.oxge[,"wtemp"], na.rm = TRUE))
#IQR(): Computes the interquartile range by subtracting the 25th percentile (Q1) from the 75th percentile (Q3).
finals$IQR_temp <- as.numeric(IQR(temp.oxge[,"wtemp"], na.rm = TRUE))
##median value
finals$median_temp <- median(temp.oxge[,"wtemp"], na.rm = TRUE)


# oxgen
library(zyp)# Theil-Sen median based estimator
temp.mon=data.frame(temp=temp.oxge$o2[!is.na(temp.oxge$o2)],
                    month=c(1:length(temp.oxge$o2[!is.na(temp.oxge$o2)])))
fit<- zyp.sen(temp~month,temp.mon)
(Oxgen_rate.t=fit$coefficients[2]*1)
finals$Oxgen_rate<-Oxgen_rate.t
finals$Oxgen_mean<-mean(temp.oxge$o2,na.rm=TRUE)
#skewness
library(moments)
# Calculate skewness
finals$skewness_Oxgen <- as.numeric(skewness(temp.oxge[,"o2"], na.rm = TRUE))
#IQR(): Computes the interquartile range by subtracting the 25th percentile (Q1) from the 75th percentile (Q3).
finals$IQR_Oxgen <- as.numeric(IQR(temp.oxge[,"o2"], na.rm = TRUE))
##median value
finals$median_Oxgen <- median(temp.oxge[,"o2"], na.rm = TRUE)



##################
temp.oxge<-read.csv("5.physical.seasonal.AL-PH nutreintyearly.csv")
#temp.oxge=temp.oxge[111:150,]#
# pH
library(zyp)# Theil-Sen median based estimator
temp.mon=data.frame(temp=temp.oxge$ph[!is.na(temp.oxge$ph)],
                    month=c(1:length(temp.oxge$ph[!is.na(temp.oxge$ph)])))
fit<- zyp.sen(temp~month,temp.mon)
(pH_rate.t=fit$coefficients[2]*1)
finals$pH_rate.t<-pH_rate.t
finals$pH_mean<-mean(temp.oxge$ph,na.rm=TRUE)
#skewness
library(moments)
# Calculate skewness
finals$skewness_pH <- as.numeric(skewness(temp.oxge[,"ph"], na.rm = TRUE))
#IQR(): Computes the interquartile range by subtracting the 25th percentile (Q1) from the 75th percentile (Q3).
finals$IQR_pH <- as.numeric(IQR(temp.oxge[,"ph"], na.rm = TRUE))
##median value
finals$median_pH <- median(temp.oxge[,"ph"], na.rm = TRUE)



# X.PO4.3...µmol.l.
library(zyp)# Theil-Sen median based estimator
temp.mon=data.frame(temp=temp.oxge$TDP[!is.na(temp.oxge$TDP)],
                    month=c(1:length(temp.oxge$TDP[!is.na(temp.oxge$TDP)])))
fit<- zyp.sen(temp~month,temp.mon)
(PO4.3_rate.t=fit$coefficients[2]*1)
finals$PO4.3_rate.t<-PO4.3_rate.t
finals$PO4.3_mean<-mean(temp.oxge$TDP,na.rm=TRUE)
#skewness
library(moments)
# Calculate skewness
finals$skewness_PO4.3 <- as.numeric(skewness(temp.oxge[,"TDP"], na.rm = TRUE))
#IQR(): Computes the interquartile range by subtracting the 25th percentile (Q1) from the 75th percentile (Q3).
finals$IQR_PO4.3 <- as.numeric(IQR(temp.oxge[,"TDP"], na.rm = TRUE))
##median value
finals$median_PO4.3 <- median(temp.oxge[,"TDP"], na.rm = TRUE)



# total_Dissolved_N
library(zyp)# Theil-Sen median based estimator
temp.mon=data.frame(temp=temp.oxge$NO3N[!is.na(temp.oxge$NO3N)],
                    month=c(1:length(temp.oxge$NO3N[!is.na(temp.oxge$NO3N)])))
#fit<- zyp.sen(temp~month,temp.mon)
#(total_Dissolved_N_rate.t=fit$coefficients[2]*4)
#finals$total_Dissolved_N_rate.t<-total_Dissolved_N_rate.t
#finals$otal_Dissolved_N_mean<-mean(temp.oxge[,"total_Dissolved_N"],na.rm=TRUE)
finals$NO3_mean<-mean(temp.oxge$NO3N,na.rm=TRUE)
#skewness
library(moments)
# Calculate skewness
finals$skewness_NO3 <- as.numeric(skewness(temp.oxge[,"NO3N"], na.rm = TRUE))
#IQR(): Computes the interquartile range by subtracting the 25th percentile (Q1) from the 75th percentile (Q3).
finals$IQR_NO3 <- as.numeric(IQR(temp.oxge[,"NO3N"], na.rm = TRUE))
##median value
finals$median_NO3 <- median(temp.oxge[,"NO3N"], na.rm = TRUE)



# Chl.a
temp.oxge<-read.csv("2.Chlorolla.seasonal.ALyearly.csv")
#temp.oxge=temp.oxge[111:150,]#
library(zyp)# Theil-Sen median based estimator
temp.mon=data.frame(temp=temp.oxge$chlor[!is.na(temp.oxge$chlor)],
                    month=c(1:length(temp.oxge$chlor[!is.na(temp.oxge$chlor)])))
finals$Chlorophyll.a<-mean(temp.oxge$chlor,na.rm=TRUE)
#skewness
library(moments)
# Calculate skewness
finals$skewness_Chlorophyll.a <- as.numeric(skewness(temp.oxge[,"chlor"], na.rm = TRUE))
#IQR(): Computes the interquartile range by subtracting the 25th percentile (Q1) from the 75th percentile (Q3).
finals$IQR_Chlorophyll.a <- as.numeric(IQR(temp.oxge[,"chlor"], na.rm = TRUE))
##median value
finals$median_Chlorophyll.a <- median(temp.oxge[,"chlor"], na.rm = TRUE)


write.csv(finals,"finals.csv")
