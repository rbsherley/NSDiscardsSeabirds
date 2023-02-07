require(jagsUI)
require(plotrix)
require(denstrip)
require(xts)
require(raster)
options(scipen = 999)
set.seed(42) # to reproduce the same results on randon numbers generated in R each time -- remove for truly random numer generation
MainDir = "~/Google Drive/Main_File_Store/Papers/My Papers/In review/Discards North Sea/Model 1990_offal_sensitivity"
input.dir = paste0(MainDir,"/Input data 1990")
output.dir = paste0(MainDir,"/Output")
#dir.create(input.dir,showWarnings = FALSE)
dir.create(output.dir,showWarnings = FALSE)
setwd(input.dir)
################################################
#            1. DATA IMPORTATION               #
################################################

## 1. seabird body data (data from the literature),
## 2. seabird proportions of the year spent breeding (data from the literature),
## 3. how much of each species' diet is made up of discards (data from the literature), 
## 4. seabird relative abundance proportions for 2010 and 2011 (from ESAS database), 
## 5. observed proportions of each discards type that gets consumed (based on Garthe et al. 1996), 
## 6. the observed energetic content for each discard group (not inverts): 2.5, 50, 97.5 percentile and mean (from Heath and Cook 2015)
## 7. the observed assimilation efficiency on fish-based diet from 21 seabird species in 14 studies

Species <- read.table('seabird_weights.txt', header=T, sep='\t') ## 1. seabird body mass: min, max, mean and SD (See supplementary material for sources),
Sp.names=names(Species)
nspecies <- ncol(Species)
breedtimes <- read.table('seabird_breed_times.txt', header=T, sep='\t') ## 2. seabird days and proportions of the year spent breeding (from Schreiber and Burger 2002 or Cramp et al. 1983),
discard_diet<-read.table('species_discard.txt', header=T, sep='\t') ## 3. how much of each species' diet is made up of discards (data from the literature)
abundance_proportion<-read.table("abundance_1990.txt", header=T, sep='\t') ## 4. seabird relative abundance proportions for 2010 and 2011 (from ESAS database), 
discard_percent<-read.table("discard_percent.txt", header=T, sep='\t') ## 5. observed proportions of each discards type (other than invertebrates) that gets consumed (based on Garthe et al. 1996), 
energy_content<-read.table("energy_content_1990_down10.txt", header=T, sep='\t') ## 6. the observed energetic content for each discard group (not inverts): 2.5, 50, 97.5 percentile and mean (from Heath and Cook 2015)
energy_content <- as.matrix(energy_content[1:5])
assim_eff <- read.table("assimilation_efficiency.txt", header=T, sep='\t') ## 7. the observed assimilation efficiency on fish-based diet from 21 seabird species in 14 studies
assim_eff <- as.matrix(assim_eff[2])
ndiscards<-ncol(energy_content)
discard.names=c("roundfish", "elasmobranch", "flatfish", "offal","other")
FMR_dat <- read.table('FMR_Dunn.txt', header=T, sep='\t')
################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################


########################################################################
#                  2. SEABIRD FMR DATA PREPARATION                     #
########################################################################

# Convert mean, upper and lower confidence intervals of species specific FMR from Dunn et al. 2018 into shape and rate parameters for gamma priors and combine the three distributions into one using simulations of a gamma distribution
# nsims = 333333
# Isd=Bsd=Csd=FMR_shape_N=FMR_rate_N=0
# for(i in 1:dim(FMR_dat)[1]){
#   Isd[i] = (((FMR_dat[i,3]-FMR_dat[i,4])/1.96)+((FMR_dat[i,5]-FMR_dat[i,3])/1.96))/2 # standard deviation based on the quantiles of a normal distribution, i.e principal that 95% falls within	1.96 sds of the mean
#   Bsd[i] = (((FMR_dat[i,6]-FMR_dat[i,7])/1.96)+((FMR_dat[i,8]-FMR_dat[i,6])/1.96))/2
#   Csd[i] = (((FMR_dat[i,9]-FMR_dat[i,10])/1.96)+((FMR_dat[i,11]-FMR_dat[i,9])/1.96))/2
#   FMR_shape_N[i] <- (FMR_dat[i,3]^2/Isd[i]^2)
#   FMR_rate_N[i] <- (FMR_dat[i,3]/Isd[i]^2)
# }
# 
# sampI=sampB=sampC=matrix(data=0,ncol=length(Isd),nrow=nsims)
# for(i in 1:length(Isd)){
#   for(j in 1:nsims){
#     sampI[j,i]=rgamma(n=1,shape=(FMR_dat[i,3]^2/Isd[i]^2),rate=(FMR_dat[i,3]/Isd[i]^2))
#     sampB[j,i]=rgamma(n=1,shape=(FMR_dat[i,6]^2/Bsd[i]^2),rate=(FMR_dat[i,6]/Bsd[i]^2))
#     sampC[j,i]=rgamma(n=1,shape=(FMR_dat[i,9]^2/Csd[i]^2),rate=(FMR_dat[i,9]/Csd[i]^2))
#   }}
# 
# FMR.sims = rbind(sampI,sampB,sampC)
# 
# FMR_shape=FMR_rate=0
# for(i in 1:dim(FMR_dat)[1]){
#   FMR_shape[i]=(mean(FMR.sims[,i])/sd(FMR.sims[,i]))^2
#   FMR_rate[i]=mean(FMR.sims[,i])/sd(FMR.sims[,i])^2
# }
# save(FMR_shape,file=paste0(output.dir,"/FMR_shape.rdata"))
# save(FMR_shape_N,file=paste0(output.dir,"/FMR_shape_N.rdata"))
# save(FMR_rate,file=paste0(output.dir,"/FMR_rate.rdata"))
# save(FMR_rate_N,file=paste0(output.dir,"/FMR_rate_N.rdata"))

load("~/Google Drive/Main_File_Store/Papers/My Papers/In review/Discards North Sea/Model 1990/Output/FMR_shape.rdata")
load("~/Google Drive/Main_File_Store/Papers/My Papers/In review/Discards North Sea/Model 1990/Output/FMR_rate.rdata")
load("~/Google Drive/Main_File_Store/Papers/My Papers/In review/Discards North Sea/Model 1990/Output/FMR_shape_N.rdata")
load("~/Google Drive/Main_File_Store/Papers/My Papers/In review/Discards North Sea/Model 1990/Output/FMR_rate_N.rdata")

########################################################################
#     3. SEABIRD ABUNDANCE, BODY MASS AND DIETS: DATA PREPARATION      #
########################################################################

wm<-wsd<-wmax<-wmin<-bd<-nd<-pb<-pn<-diet.mn<-diet.max<-diet.min<-diet.n<-abundance.NB<-abundance.B<-0
for(i in 1:nspecies){
  ## Relative abundances of seabirds in the North Sea:
  abundance.B[i]<-(abundance_proportion[i,1])   #the abundance average for the breeding season (may to august)
  abundance.NB[i]<-(abundance_proportion[i,2])   #the abundance average for the non-breeding season (september to april)
  
  ## Species-specific times and proportion of the year seabirds in the North Sea spend breeding:
  bd[i] <- breedtimes[1,i]
  nd[i] <- (365-breedtimes[1,i])
  pb[i] <- breedtimes[2,i]
  pn[i] <- (1-breedtimes[2,i])
  
  ## Body mass values - to use to produce a informative normal prior:
  wm[i] <- Species[3,i]  # mean weight
  wmax[i] <- Species[2,i]   # max weight
  wmin[i] <- Species[1,i]   # max weight
  wsd[i] <- Species[4,i]   # standard deviation of weight
  
  ## Observed proportions of discards in the diet of each seabird species:
  diet.mn[i]<-mean(discard_diet[discard_diet[,3]==i,1])   # the mean proportion of consumed discards for each species - used for diagnostic plots at the end
  diet.n[i]<-length(discard_diet[discard_diet[,3]==i,1])  # the number of estimates for each species - used for diagnostic plots at the end
  diet.max[i]<-max(discard_diet[discard_diet[,3]==i,1])   # the max proportion of consumed discards for each species - used to bound a uniform prior
  diet.min[i]<-min(discard_diet[discard_diet[,3]==i,1])   # the max proportion of consumed discards for each species - used to bound a uniform prior
} # for

################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################

########################################################################
#               4. INVERTEBRATE ENERGY: DATA PREPARATION               #
########################################################################

# The energy available from discarded invertebrates is based on an energetic equivalent factor of 2.0107 and between 4 and 11 kg of inverts for every kg of sole landed.
# Here the sole landings (tonnes) are read in, converted to kJ and multiplied by the energetic equivalent factor of 2.0107. The multiplication for every kg of sole caught happens in the JAGS
sole_landings<-matrix(data=c(21850,26110,30560,26694.23347), nrow=1, ncol=4, byrow=T) #Sole landings from Heath and Cook 2010.
invert_energy<-matrix(data = c(sole_landings*(2.0107*1000000)), nrow=4, ncol=1, byrow=FALSE, dimnames = list(c("2.5", "50", "97.5", "mean"), c("Sole energy")))
# Convert the mean, low and upper credible intervals into a mean and sd to use as an informative prior.
invert.energy.mn <- invert_energy[4]
invert.energy.sd <- (((invert_energy[3]-invert_energy[4])/1.96)+((invert_energy[4]-invert_energy[1])/1.96))/2 # # derive the variance based on the central limit theorm
invert.energy.shape <- ((invert.energy.mn)^2/(invert.energy.sd)^2)
invert.energy.rate <- invert.energy.mn/(invert.energy.sd)^2
### parameters to specify a gamma prior assuming the parameterization gamma(shape,rate): shape = mean^2/sd^2, rate = mean/sd^2 (Kruschke 2015)

################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################

########################################################################
#              5. OTHER DISCARDS ENERGY: DATA PREPARATION              #
########################################################################

# Convert the mean, low and upper credible intervals into parameters for a gamma distribution to use as an informative prior.
non.invert.energy.mn <- energy_content[4,]
non.invert.energy.sd <- (((energy_content[3,]-energy_content[4,])/1.96)+((energy_content[4,]-energy_content[1,])/1.96))/2
non.invert.energy.rate <- non.invert.energy.mn/(non.invert.energy.sd)^2
non.invert.energy.shape <- ((non.invert.energy.mn)^2/(non.invert.energy.sd)^2)
### parameters to specify a gamma prior assuming the parameterization gamma(shape,rate): shape = mean^2/sd^2, rate = mean/sd^2 (Kruschke 2015)

################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################

########################################################################
#                         6. JAGS MODEL CODE                           #
########################################################################

sink("supported90.jags")
cat("
    model { 

#################################################################################  
#         SEABIRD BODY MASS, BMR, FMR AND ANNUAL ENERGY REQUIREMENTS            #
################################################################################# 

    for (i in 1:S){
      FMR.N[i] ~ dgamma(FMR.shape.N[i],FMR.rate.N[i]) # Gamma priors for FMR during the breeding season (FMR.B) and non-breeding season (FMR.N)
      FMR.B[i] ~ dgamma(FMR.shape[i],FMR.rate[i])
      B.ER[i] <- (FMR.B[i]*bd[i])       # Energy requirements during the breeding season (kJ)
      NB.ER[i] <- (FMR.N[i]*nd[i])      # Energy requirements during the non-breeding season (kJ)

#################################################################################  
#                       SEABIRD DISCARDS IN THE DIET                            #
################################################################################# 
  ### uninformative priors for the parameters to map observed proportions of discards in the diet of the 8 seabird species to distrubtions ### 
       Ps[i] ~ dunif(min(diet.min),max(diet.max))
       t0[i] ~ dnorm(0,0.001) 
       tau.s[i] <- exp(t0[i])
       u.s[i] ~ dnorm(0,tau.u.s)      ### hyperprior for random effect for discard type
  #################################################################################
    } # close i loop

  ###  Uninformative prior for the variance for the random effect (of seabird species)
       tau.u.s ~ dgamma(0.001,0.001)    
  ################################################################################# 

  ### likelihood to map observed discard proportions in the seabird diet to distrubtions (Ps) ###
    for (j in 1:length(d.p.s)){
      d.p.s[j] ~ dbeta(p.s[j], q.s[j])
      p.s[j] <- mu2.s[j] * tau.s[species[j]]
      q.s[j] <- (1 - mu2.s[j]) * tau.s[species[j]]
      mu2.s[j] <- Ps[species[j]] + u.s[species[j]]
    } # close j loop
  ################################################################################# 


#################################################################################  
#                               FISH DISCARD GROUPS                             #
################################################################################# 

  ### Uninformative priors for the parameters to map observed discard proportions of fish, elasmobranches and other discards to distrubtions ### 
    for(k in 1:max(discard.type)){     ### (one level per discard type)
       Ct[k] ~ dunif(0,1)
       t2[k] ~ dnorm(0,0.001) 
       tau.n[k] <- exp(t2[k])
       u.d[k] ~ dnorm(0,tau.u.d)      ### hyperprior for random effect for discard type
  #################################################################################
    }
  ###  Uninformative prior for the variance for the random effect (of discard type)
       tau.u.d ~ dgamma(0.001,0.001)
    #################################################################################

  ### Uninformative priors for the parameters to map observed assimilation efficiency to a distrubtion ###
       Amod ~ dunif(0,1)
       t3 ~ dnorm(0,0.001) 
       tau.ad <- exp(t3)
  ################################################################################# 
    
  ### Likelihood to map observed assimilation efficiency to a distrubtion (Amod) ###
    for(m in 1:length(ae)){
      ae[m] ~ dbeta(p.ae[m], q.ae[m])
      p.ae[m] <- mu2.ae[m] * tau.ad
      q.ae[m] <- (1 - mu2.ae[m]) * tau.ad
      mu2.ae[m] <- Amod
    } # close m loop
  #################################################################################
  
  ### likelihood to map observed discard proportions of fish, elasmobranches and other discards to distrubtions (Ct) ###
    for (j in 1:N){
      d.p.n[j] ~ dbeta(p.n[j], q.n[j])
      p.n[j] <- mu2.n[j] * tau.n[discard.type[j]]
      q.n[j] <- (1 - mu2.n[j]) * tau.n[discard.type[j]]
      mu2.n[j] <- Ct[discard.type[j]] + u.d[discard.type[j]]
    } # close j loop

   ###  Gamma prior for the energy content of fish and offal discarded in the North Sea (to avoid negative values) ###
    for(l in 1:ndiscards){     ### (one level per fish discard type)
      Eft[l] ~ dgamma(non.invert.energy.shape[l],non.invert.energy.rate[l]) 
   #################################################################################
    
   ### Energy that can be assimilated from fish, elasmobranches and other discards (Aft) ####
      Aft[l] <- ((Ct[l]*Eft[l])*Amod)
    } # close l loop
   #################################################################################
   
#################################################################################  
#                             INVERTEBRATE DISCARDS                             #
#################################################################################     
    
  ### Energy content of invertebrates discarded in the North Sea (read in as kJ from Garthe et al. 1996) ###
    Ei ~ dgamma(invert.energy.shape,invert.energy.rate)
  #################################################################################
    
  ### Uniform prior for the weight of invertebrates discarded in the North Sea (based on a range of 4 to 11 kg of invertebrates for every kg of sole caught) ###
    Wi ~ dunif(4, 11)
  #################################################################################
    
  ### Estimated distribution for the energy content of invertebrates discarded in the North Sea
    T.Ei <- Ei*Wi
  #################################################################################
    
  ### Energy that can be assimilated from invertebrate discards ####
    Ai <- ((Wi*Ei)*Amod)*Ct[6]
  ### energy assimilated from benthic inverts = distribution for the prop. of invert. discards consumed (Ci) 
  ### x distribution for the energy content of invertebrates discarded in the North Sea (Ei) 
  ### x assimilation efficiency (assumed to be 0.65 for inverts.)
  ### x distribution for conversion factor from weight of sole caught to weight of inverts. caught (Wi)
  #################################################################################

################################################################################# 
#                   DERIVED PARAMETERS AND MODEL CALCULATIONS                   #
#################################################################################

   ### Total energy available to seabirds as discards in the North Sea in 2010 ###
    Discarded.Biomass <- sum(Eft)+T.Ei
    Discard.Energy <- sum(Aft)+Ai  
   ################################################################################# 
    
   ### Total energy available to each seabird species in discards:
     for (i in 1:S){
       B.EA[i]<-((Discard.Energy*pb[i])*abundance.B[i])     # In the breeding season
       NB.EA[i]<-((Discard.Energy*pn[i])*abundance.NB[i])   # In the non-breeding season
     } # close i loop  
       Esa <- B.EA+NB.EA                         # Across the year

   ###  Annual energy requirements of each seabird species across the year (kJ)   
       Rsa <- B.ER+NB.ER         

   ### The amount of energy needed in the form of discards by each species across the year (Energy requirements * discard proportions in the seabird diet (Ps))
       Dsa <- Rsa*Ps
       
   ### The number of seabirds in each species that can be supported by discards.
       Ts<-(Esa/Dsa)               # Total energy available to each seabird species in discards/The amount of energy needed in the form of discards by each species
      
   ### The total number of seabirds that can be supported by discards.
       T <- sum(Ts)
   #####################################################
    
   ##### Split estimates by breeding and non-breeding periods ############
   ### Total energy available to each seabird species in discards:
    for (i in 1:S){
      sB.EA[i]<-((Discard.Energy*pb[i])*(abundance.B[i]*pb[i]))     # In the breeding season
      sNB.EA[i]<-((Discard.Energy*pn[i])*(abundance.NB[i]*pn[i]))   # In the non-breeding season
      } # close i loop  

      BRE <- sB.EA/(B.ER*Ps)
      NON <- sNB.EA/(NB.ER*Ps)
      T.BRE <- sum(BRE)
      T.NON <- sum(NON)

   # End

########################################################################################################################
########################################################################################################################
    }
    ",fill = TRUE)
sink()

########################################################################
#                 6. DATA AND PARAMETERS TO PASS TO JAGS               #
########################################################################
  
### bundle the data to pass to JAGS
jags.data <- list(
  # Seabirds
  S=nspecies, bd=bd, nd=nd, pb=pb, pn=pn,
  abundance.NB=abundance.NB, abundance.B=abundance.B,  
  species=discard_diet$species.n, d.p.s=discard_diet$dtp,
  diet.min=diet.min,diet.max=diet.max,
  FMR.shape=FMR_shape,FMR.rate=FMR_rate,FMR.shape.N=FMR_shape_N,FMR.rate.N=FMR_rate_N,
  # Discards
  d.p.n=discard_percent$dpt, N=nrow(discard_percent),
  discard.type=discard_percent$type_n,
  non.invert.energy.shape=non.invert.energy.shape,non.invert.energy.rate=non.invert.energy.rate,
  ndiscards=ndiscards, ae=assim_eff[,1],
  # Inverts
  invert.energy.shape=invert.energy.shape,invert.energy.rate=invert.energy.rate
) # close list

inits <- function() list(Ms = rnorm(nspecies,wm,wsd))

### parameters to trace  
parameters <- c("Ts","T","Discarded.Biomass","Discard.Energy","Rsa","Esa","Eft","T.Ei","Ei","Ps","Ct","Ai","Ai1","Amod","Ms","BRE","NON","T.BRE","T.NON","BMR","FMR.B","FMR.N")

### MCMC settings
ni <- 1100000
nt <- 10
nb <- 100000
nc <- 3

### Call JAGS from R
support90down<- jags(jags.data,parallel = T, inits=inits, parameters, "supported90.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC = F)
   
### Print model output
print(support90down, 3)

### Save model object
save(support90down,file=paste0(output.dir,"/1990_model_down10.rdata"))

### Write model output to file
write.csv(support90down$summary, file=paste0(output.dir,"/1990_Output_down10.csv"))


##################################################################################################
#                                    7. Diagnostics Plots                                        #
##################################################################################################

######## PLOT TO CHECK CONVERGENCE OF FINAL ESTIMATED NUMBER OF SEABIRDS SUPPORTED ########
pdf(paste0(output.dir,"/TracePlot_TOT_90.pdf"),width=6,height=3) 
par(mfrow=c(1,1),mar=c(4,4,1,1))
traceplot(support90, parameters = "T",main="")
dev.off()
################################################################################################
################################################################################################
load("~/Google Drive/Main_File_Store/Papers/My Papers/In review/Discards North Sea/Model 1990/Output/1990_model.RData")
density(support$sims.list$T)

#Image <- function(){
quartz(width=3.15,height=3.15)   
par(mfrow=c(1,1),mar=c(3.5, 3.5, 0.1, 0.1),mgp = c(2, 0.6, 0))
plotCI(x=c(1:3),y=c(0,0,0),ui=c(0,0,0),li=c(0,0,0),xaxt='n',yaxt='n',xlim=c(0,1),
       ylim=c(0,10*10^6),
       xlab='',ylab=' ',main='',pch=16, col='white',lty=2)
abline(h=support90$mean$T, lty=2, col="grey70")

denstrip(density(support90down$sims.list$T)$x,density(support90down$sims.list$T)$y, at=0.2, horiz=FALSE,colmax='grey20',
         mticks=c(support90down$mean$T),mlen=1.5,mwd=1.5,mcol='black',
         ticks=c(support90down$q2.5$T,support90down$q97.5$T),
         tcol=c('grey70'),tlen=1, twd=1.5,scale=0.8, colmin = 'grey95')

denstrip(density(support90$sims.list$T)$x,density(support90$sims.list$T)$y, at=0.5, horiz=FALSE,colmax='grey20',
         mticks=c(support90$mean$T),mlen=1.5,mwd=1.5,mcol='black',
         ticks=c(support90$q2.5$T,support90$q97.5$T),
         tcol=c('grey70'),tlen=1, twd=1.5,scale=0.8, colmin = 'grey95')

denstrip(density(support90up$sims.list$T)$x,density(support90up$sims.list$T)$y, at=0.8, horiz=FALSE,colmax='grey20',
         mticks=c(support90up$mean$T),mlen=1.5,mwd=1.5,mcol='black',
         ticks=c(support90up$q2.5$T,support90up$q97.5$T),
         tcol=c('grey70'),tlen=1, twd=1.5,scale=0.8, colmin = 'grey95')
mtext(~"Individuals x "*10^6*"", side=2, line=1.9,cex=1.2)
mtext("Offal", side=1, line=2.5,cex=1.2)
axis(2,at=seq(0,25*10^6,5*10^6),labels=seq(0,25,5))
axis(2,at=seq(0,25*10^6,2.5*10^6),labels=NA)
axis(1,at=c(0.2,0.5,0.8),labels=c("\u219310%","Base","\u219110%"))
box()
quartz.save(file=paste0(output.dir,"/Figure_S9.jpg"),type='jpeg',dpi=800) ## replace with savePlot('Figure_2.pdf',type='pdf') if using a computer running windows



}
pdf(file=paste0(output.dir,"/Figure_S9.pdf"), width = 3.15, height = 3.15, onefile = TRUE,pagecentre=TRUE,colormodel="rgb")
Image()
dev.off()
jpeg(file=paste0(output.dir,"/Figure_S9.jpg"), width = 3.15, height = 3.15, units='in',res=350)
Image()
dev.off()
