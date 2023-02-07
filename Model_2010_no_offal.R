require(jagsUI)
require(plotrix)
require(denstrip)
require(xts)
require(raster)
options(scipen = 999)
set.seed(42) # to reproduce the same results on randon numbers generated in R each time -- remove for truly random numer generation
MainDir = "~/Google Drive/Main_File_Store/Papers/My Papers/In review/Discards North Sea/Model 2010_no_offal"
input.dir = paste0(MainDir,"/Input data_no_offal")
output.dir = paste0(MainDir,"/Output")
dir.create(input.dir,showWarnings = FALSE)
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
abundance_proportion<-read.table("abundance_2010.txt", header=T, sep='\t') ## 4. seabird relative abundance proportions for 2010 and 2011 (from ESAS database), 
discard_percent<-read.table("discard_percent.txt", header=T, sep='\t') ## 5. observed proportions of each discards type (other than invertebrates) that gets consumed (based on Garthe et al. 1996), 
energy_content<-read.table("energy_content_2010.txt", header=T, sep='\t') ## 6. the observed energetic content for each discard group (not inverts): 2.5, 50, 97.5 percentile and mean (from Heath and Cook 2015)
energy_content <- as.matrix(energy_content[,1:4])
assim_eff <- read.table("assimilation_efficiency.txt", header=T, sep='\t') ## 7. the observed assimilation efficiency on fish-based diet from 21 seabird species in 14 studies
assim_eff <- as.matrix(assim_eff[2])
ndiscards<-ncol(energy_content)
discard.names=c("roundfish", "elasmobranch", "flatfish","other")
FMR_dat <- read.table('FMR_Dunn.txt', header=T, sep='\t')
################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################

########################################################################
#                  2. SEABIRD FMR DATA PREPARATION                     #
########################################################################
Dir1990s = "~/Google Drive/Main_File_Store/Papers/My Papers/In review/Discards North Sea/Model 1990/Output/"
load(paste0(Dir1990s,"FMR_shape.rdata"))
load(paste0(Dir1990s,"FMR_shape_N.rdata"))
load(paste0(Dir1990s,"FMR_rate.rdata"))
load(paste0(Dir1990s,"FMR_rate_N.rdata"))

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
sole_landings<-matrix(data=c(10360,12420,14740,12561.2), nrow=1, ncol=4, byrow=T) #Sole landings from Heath and Cook 2010.
invert_energy<-matrix(data = c(sole_landings*(2.0107*1000000)), nrow=4, ncol=1, byrow=FALSE, dimnames = list(c("2.5", "50", "97.5", "mean"), c("Sole energy")))
# Convert the mean, low and upper credible intervals into a mean and sd to use as an informative prior.
invert.energy.mn <- invert_energy[4]
invert.energy.sd <- (((invert_energy[3]-invert_energy[4])/1.96)+((invert_energy[4]-invert_energy[1])/1.96))/2 # derive the variance based on the central limit theorm
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
non.invert.energy.shape <- ((non.invert.energy.mn)^2/(non.invert.energy.sd)^2)
non.invert.energy.rate <- non.invert.energy.mn/(non.invert.energy.sd)^2
### parameters to specify a gamma prior assuming the parameterization gamma(shape,rate): shape = mean^2/sd^2, rate = mean/sd^2 (Kruschke 2015)

################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################

########################################################################
#                         6. JAGS MODEL CODE                           #
########################################################################

sink("supported.jags")
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
       u.s[i] ~ dnorm(0,tau.u.s)      ### hyperprior for random effect seabird species
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
    
  ### Gamma prior for the energy content of invertebrates discarded in the North Sea (read in as kJ based on landings of sole * energetic equivalent of 2.0107) ###
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
  ### x assimilation efficiency (Amod)
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

### initials (not strictly needed)
inits <- function() list(Ms = rnorm(nspecies,wm,wsd))

### parameters to trace  
parameters <- c("Ts","T","Discarded.Biomass","Discard.Energy","Rsa","Esa","Eft","T.Ei","Ei","Ps","Ct","Ai","Ai1","Amod","Ms","BRE","NON","T.BRE","T.NON","BMR","FMR.B","FMR.N")

### MCMC settings
ni <- 1100000
nt <- 10
nb <- 100000
nc <- 3

### Call JAGS from R
support10no<- jags(jags.data,parallel = T, inits=inits, parameters, "supported.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC = F)
   
### Print model output
print(support10no, 3)

### Save model object
save(support10no,file=paste0(output.dir,"/2010_model_no_offal.rdata"))

### Write model output to file
write.csv(support10no$summary, file=paste0(output.dir,"/2010_Output_no_offal.csv"))


################################################################################################
################################################################################################

load("~/Google Drive/Main_File_Store/Papers/My Papers/In review/Discards North Sea/Model 2010/Output/2010_model.RData")
load("~/Google Drive/Main_File_Store/Papers/My Papers/In review/Discards North Sea/Model 1990/Output/1990_model.RData")
d10 <- density(support$sims.list$T)
b10 <- density(support$sims.list$T.BRE)
n10 <- density(support$sims.list$T.NON)

d <- density(support10no$sims.list$T)
b <- density(support10no$sims.list$T.BRE)
n <- density(support10no$sims.list$T.NON)

d90 <- density(support90$sims.list$T)
b90 <- density(support90$sims.list$T.BRE)
n90 <- density(support90$sims.list$T.NON)

quartz(width=3.15,height=9) ## replace with windows(width=6,height=6) if using a computer running windows
par(mfrow=c(3,1),mar=c(4, 4.5, 0.1, 0.1),mgp = c(2, 0.7, 0), family='sans')
# Panel 1
plot(d90$x,(d90$y/max(d90$y)),main='',xlab='',xaxt='n',yaxt='n',ylab='',ylim=c(0,1),xlim=c(0,12*10^6),type='l',col='white')
lines(b90$x,(b90$y/max(b90$y)))
polygon(c(b90$x, rev(b90$x)), c(b90$y/max(b90$y), seq(0,0,length.out = length(b90$y))),col=adjustcolor("green",alpha.f=0.4), border = 'NA')
lines(n90$x,(n90$y/max(n90$y)))
polygon(c(n90$x, rev(n90$x)), c(n90$y/max(n90$y), seq(0,0,length.out = length(n90$y))),col=adjustcolor("purple",alpha.f=0.4), border = 'NA')
abline(v=support90$mean$T.BRE,lty=1,col='green')
abline(v=support90$q2.5$T.BRE,lty=2,col='green')
abline(v=support90$q97.5$T.BRE,lty=2,col='green')
abline(v=support90$mean$T.NON,lty=1,col="purple")
abline(v=support90$q2.5$T.NON,lty=2,col="purple")
abline(v=support90$q97.5$T.NON,lty=2,col="purple")
axis(1,at=seq(0,12*10^6,2*10^6),lab=seq(0,12,2),las=1,cex.axis=1.5)
axis(2,at=c(0,0.25,0.5,0.75,1),lab=c("0","0.25","0.5","0.75","1"),las=2,cex.axis=1.5)
mtext(~"Individuals x "*10^6*"", side=1, line=2.4,cex=1.2)
mtext("Probability density", side=2, line=3.1,cex=1.2)
mtext("(A)", at = 1,side=2, line=2.5,cex=1.1,las=1)
text(11.1*10^6,0.99,"1990", cex=2)
box()

# Panel 2
plot(d10$x,(d10$y/max(d10$y)),main='',xlab='',xaxt='n',yaxt='n',ylab='',ylim=c(0,1),xlim=c(0,12*10^6),type='l',col='white')
lines(b10$x,(b10$y/max(b10$y)))
polygon(c(b10$x, rev(b10$x)), c(b10$y/max(b10$y), seq(0,0,length.out = length(b10$y))),col=adjustcolor("green",alpha.f=0.4), border = 'NA')
lines(n10$x,(n10$y/max(n10$y)))
polygon(c(n10$x, rev(n10$x)), c(n10$y/max(n10$y), seq(0,0,length.out = length(n10$y))),col=adjustcolor("purple",alpha.f=0.4), border = 'NA')
abline(v=support$mean$T.BRE,lty=1,col='green')
abline(v=support$q2.5$T.BRE,lty=2,col='green')
abline(v=support$q97.5$T.BRE,lty=2,col='green')
abline(v=support$mean$T.NON,lty=1,col="purple")
abline(v=support$q2.5$T.NON,lty=2,col="purple")
abline(v=support$q97.5$T.NON,lty=2,col="purple")
axis(1,at=seq(0,12*10^6,2*10^6),lab=seq(0,12,2),las=1,cex.axis=1.5)
axis(2,at=c(0,0.25,0.5,0.75,1),lab=c("0","0.25","0.5","0.75","1"),las=2,cex.axis=1.5)
mtext(~"Individuals x "*10^6*"", side=1, line=2.4,cex=1.2)
mtext("Probability density", side=2, line=3.1,cex=1.2)
mtext("(B)", at = 1,side=2, line=2.5,cex=1.1,las=1)
text(11.1*10^6,0.99,"2010", cex=2)
box()

# Panel 3
plot(d$x,(d$y/max(d$y)),main='',xlab='',xaxt='n',yaxt='n',ylab='',ylim=c(0,1),xlim=c(0,12*10^6),type='l',col='white')
lines(b$x,(b$y/max(b$y)))
polygon(c(b$x, rev(b$x)), c(b$y/max(b$y), seq(0,0,length.out = length(b$y))),col=adjustcolor("green",alpha.f=0.4), border = 'NA')
lines(n$x,(n$y/max(n$y)))
polygon(c(n$x, rev(n$x)), c(n$y/max(n$y), seq(0,0,length.out = length(n$y))),col=adjustcolor("purple",alpha.f=0.4), border = 'NA')
abline(v=support10no$mean$T.BRE,lty=1,col='green')
abline(v=support10no$q2.5$T.BRE,lty=2,col='green')
abline(v=support10no$q97.5$T.BRE,lty=2,col='green')
abline(v=support10no$mean$T.NON,lty=1,col="purple")
abline(v=support10no$q2.5$T.NON,lty=2,col="purple")
abline(v=support10no$q97.5$T.NON,lty=2,col="purple")
axis(1,at=seq(0,12*10^6,2*10^6),lab=seq(0,12,2),las=1,cex.axis=1.5)
axis(2,at=c(0,0.25,0.5,0.75,1),lab=c("0","0.25","0.5","0.75","1"),las=2,cex.axis=1.5)
mtext(~"Individuals x "*10^6*"", side=1, line=2.4,cex=1.2)
mtext("Probability density", side=2, line=3.1,cex=1.2)
mtext("(C)", at = 1,side=2, line=2.5,cex=1.1,las=1)
text(9.1*10^6,0.99,"2010 no offal", cex=2)
box()
quartz.save(file=paste0(output.dir,"/Figure_2_3panel.pdf"),type='pdf') ## replace with savePlot('Figure_2.pdf',type='pdf') if using a computer running windows
dev.off()

####### OR Figure S9 2 panel.

quartz(width=3.15,height=6) ## replace with windows(width=6,height=6) if using a computer running windows
par(mfrow=c(2,1),mar=c(4, 4, 0.1, 0.1),mgp = c(2, 0.6, 0), family='sans')
# Panel 1
plot(d10$x,(d10$y/max(d10$y)),main='',xlab='',xaxt='n',yaxt='n',ylab='',ylim=c(0,1),xlim=c(0,12*10^6),type='l',col='white')
lines(b10$x,(b10$y/max(b10$y)))
polygon(c(b10$x, rev(b10$x)), c(b10$y/max(b10$y), seq(0,0,length.out = length(b10$y))),col=adjustcolor("green",alpha.f=0.4), border = 'NA')
lines(n10$x,(n10$y/max(n10$y)))
polygon(c(n10$x, rev(n10$x)), c(n10$y/max(n10$y), seq(0,0,length.out = length(n10$y))),col=adjustcolor("purple",alpha.f=0.4), border = 'NA')
abline(v=support$mean$T.BRE,lty=1,col='green')
abline(v=support$q2.5$T.BRE,lty=2,col='green')
abline(v=support$q97.5$T.BRE,lty=2,col='green')
abline(v=support$mean$T.NON,lty=1,col="purple")
abline(v=support$q2.5$T.NON,lty=2,col="purple")
abline(v=support$q97.5$T.NON,lty=2,col="purple")
axis(1,at=seq(0,12*10^6,2*10^6),lab=seq(0,12,2),las=1,cex.axis=1.1)
axis(2,at=c(0,0.25,0.5,0.75,1),lab=c("0","0.25","0.5","0.75","1"),las=2,cex.axis=1.1)
mtext(~"Individuals x "*10^6*"", side=1, line=2.4,cex=1.2)
mtext("Probability density", side=2, line=3.1,cex=1.2)
mtext("(A)", at = 1,side=2, line=2.5,cex=1.1,las=1)
text(11.1*10^6,0.99,"2010", cex=1.3)
box()

# Panel 2
plot(d$x,(d$y/max(d$y)),main='',xlab='',xaxt='n',yaxt='n',ylab='',ylim=c(0,1),xlim=c(0,12*10^6),type='l',col='white')
lines(b$x,(b$y/max(b$y)))
polygon(c(b$x, rev(b$x)), c(b$y/max(b$y), seq(0,0,length.out = length(b$y))),col=adjustcolor("green",alpha.f=0.4), border = 'NA')
lines(n$x,(n$y/max(n$y)))
polygon(c(n$x, rev(n$x)), c(n$y/max(n$y), seq(0,0,length.out = length(n$y))),col=adjustcolor("purple",alpha.f=0.4), border = 'NA')
abline(v=support10no$mean$T.BRE,lty=1,col='green')
abline(v=support10no$q2.5$T.BRE,lty=2,col='green')
abline(v=support10no$q97.5$T.BRE,lty=2,col='green')
abline(v=support10no$mean$T.NON,lty=1,col="purple")
abline(v=support10no$q2.5$T.NON,lty=2,col="purple")
abline(v=support10no$q97.5$T.NON,lty=2,col="purple")
axis(1,at=seq(0,12*10^6,2*10^6),lab=seq(0,12,2),las=1,cex.axis=1.1)
axis(2,at=c(0,0.25,0.5,0.75,1),lab=c("0","0.25","0.5","0.75","1"),las=2,cex.axis=1.1)
mtext(~"Individuals x "*10^6*"", side=1, line=2.4,cex=1.2)
mtext("Probability density", side=2, line=3.1,cex=1.2)
mtext("(B)", at = 1,side=2, line=2.5,cex=1.1,las=1)
text(9*10^6,0.99,"2010 no offal", cex=1.3)
box()
quartz.save(file=paste0(output.dir,"/Figure_2_3panel.pdf"),type='pdf') ## replace with savePlot('Figure_2.pdf',type='pdf') if using a computer running windows
dev.off()



