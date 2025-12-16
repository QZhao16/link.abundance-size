# Set path

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
#STL (Seasonal-Trend decomposition using Loess)
stl_detrend <- function(df, value_col, start_year, start_period,
                        frequency = 4, robust = TRUE) {
  # Extract the numeric vector
  y <- df[[value_col]]
  # Convert to time-series object
  ts_y <- ts( y,
              start = c(start_year, start_period),
              frequency = frequency # 4 is for mothly data
  )
  # STL decomposition
  stl_fit <- stl(
    ts_y,
    s.window = "periodic",
    robust = robust
  )
  # extract remainder to dataframe
  df$remainder <- as.numeric(stl_fit$time.series[, "remainder"])
  return(df)
}

# linear detrending
nomz=function(x, normalization=T){
  xt=data.frame(value=as.numeric(x),time=1:length(x))
  # Normalization (zero mean & unity variance)
  if(normalization==T){ lm_fit <- lm(value ~ time,data=xt); resid.lm_fit<-resid(lm_fit) }
  return(resid.lm_fit)
  
}

####### Function for generating lag time series 
laf=function(x,y,lagf){
  n <- NROW(x)
  x.t=x;y.t=y
  if(lagf<=0){x.t=x.t[(1-lagf):n];y.t=y.t[1:(n+lagf)]} # if lagf<0, y is leading
  if(lagf>0){x.t=x.t[1:(n-lagf)];y.t=y.t[(1+lagf):n]}  # if lagf>0, x is leading           
  return(cbind(x.t,y.t))
}

###Part 1 
new_abund<-read.csv("b.AL.fish.caught_per_effort.csv")
head(new_abund)
new_abund<-new_abund[,3:ncol(new_abund)]#has date
dim(new_abund)
target_row=2:ncol(new_abund)#only species

# then use these colnames to select biomass
new_biomass<-read.csv("b.AL.fish.length.csv") 
new_biomass<-new_biomass[,3:ncol(new_biomass)]#has date
dim(new_biomass)
new_biomass_matched=new_biomass[,c(target_row)]# find matched column names

# bind
finals=as.data.frame( cbind(new_abund$year,new_abund[,target_row],new_biomass_matched) )
#rename
colnames(finals)<-c("year",colnames(finals)[2:ncol(finals)])


# bind all
alga_zo=as.data.frame( cbind(new_abund$year,
                            new_abund[,target_row], # zoo abund
                            new_biomass_matched # zoo bioma
                            )) # zoo biomass

#rename
colnames(alga_zo)<-c("sampledate",colnames(alga_zo)[2:ncol(alga_zo)])



# part 3: combine
alga_zo<-data.frame(alga_zo[,2:ncol(alga_zo)]) # no year
dim(alga_zo) 
colnames(alga_zo)[12]
head(alga_zo)
#
n=nrow(alga_zo)
seed=25647

##################
# Part 1: trait affects abundance
##################

# Detrend + deseason time series
sdat=as.data.frame(apply(alga_zo,2,nomz))
##The index for testing causal links (a total of 12 links)
indmat=matrix(0,ncol(sdat)/2,2);colnames(indmat)=c('Effect','Cause')
indmat[,1]=colnames(sdat)[1:(ncol(sdat)/2)]
indmat[,2]=colnames(sdat)[(ncol(sdat)/2+1):ncol(sdat)]
indmat



# Determine the embedding dimensions (En) in CCM by try-and-error with best hindcast (tp=-1) skill
Emax=7
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
lib_siz=sort(c(seq(7,n,5),n)) # a sequence of library size
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
Emax=7
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
lib_siz=sort(c(seq(7,n,5),n)) # a sequence of library size
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
# Part 3: The computation of warming rate 
##################
temp.oxge<-read.csv("4.physical.seasonal.AL-temp-O2yearly.csv")
# temperature
library(trend)# Theil-Sen median based estimator
(Warming_rate.t=sens.slope(temp.oxge$wtemp)$estimates*1)#warming/year
finals$Warming_rate<-Warming_rate.t
finals$temp_mean<-mean(temp.oxge$wtemp,na.rm=TRUE)
