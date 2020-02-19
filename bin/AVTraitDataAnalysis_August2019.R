#August 2, 2019
#Renee H. Petipas
#AVTraitDataAnalysis_Revisions1
#Re-run trait analysis using all the samples including liquid inoculum that was 
#not included in the first submission. Made a new clean folder in my Dropbox so I can
# get the code and everything ready for submission
#I used the same trait data file that we used to create the phyloseq object


rm(list=ls(all=TRUE)) # clear R environment

#set working directory
setwd("~/Dropbox/PanicumMicrobiome_CleanFolder/data/Traits")
#load .csv file
metadata.E <-read.csv("AV_TraitDataE_Aug2019.csv", row.names = 1)
head(metadata.E)
metadata.E[metadata.E == "na"] <- NA
write.csv(metadata.E, "metadata.E2.csv")

metadata.E <-read.csv("metadata.E2.csv", row.names = 1)
str(metadata.E)




#change so R recognizes factors as factors
metadata.E$G5Block<-as.factor(metadata.E$G5Block)
metadata.E$InoculumVolume<-as.factor(metadata.E$InoculumVolume)
metadata.E$PresenceofMoldAug31<-as.factor(metadata.E$PresenceofMoldAug31) 
str(metadata.E)


#include inoculum type in the models below


#load packages
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)

library(emmeans)
library(ggsignif)


#Calculate R.S and RTD
metadata.E$Root_Tissue_Density <- metadata.E$CoreRootMass /metadata.E$CoreRootVolume
metadata.E$R.S <- metadata.E$CoreRootMass /metadata.E$TotalShootMass
str(metadata.E)

#remove samples that have NA for root traits

metadata.E <- dplyr:: filter(metadata.E, CoreRootLength != "NA")
str(metadata.E)
#Type I sum of squares are “sequential.”  In essence the factors are tested in the order they are listed in the model.  
#Type I is not appropriate for unbalanced design
#Type III are “partial.”In essence, every term in the model is tested in light of every other term in the model.  


#do I need to set contrasts because it is a type III anova?
#Use sum contrasts to compare each group against grand mean.
options(contrasts=c('contr.sum','contr.poly'))

##########Plot residuals with a lmer model##########

M<-lmer (TotalShootMass ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)



##plot 
plot(M, add.smooth = FALSE)

##################Residuals vs. fitted understand normality distribution##################### 
#Step 1 Pull out residuals
E<-resid(M)

## Step 2 plot a histogram of residual frequency-
hist(E, xlab = "Residuals", main = "")

##################QQPlots to understand linearity of data############
##plotting quantiles (i.e. percentiles) from our distribution against a theoretical distribution
## QQplots of sample vs. theoretical quantiles also indicate non-normality
qqnorm(metadata.E$TotalShootMass)
qqline(metadata.E$TotalShootMass)

##########Plot residuals with a lmer model now log transformed##########


M<-lmer (log(TotalShootMass) ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
##plot 
plot(M, add.smooth = FALSE)
#better with a log transformation

#Step 1 Pull out residuals
E<-resid(M)

## Step 2 plot a histogram of residual frequency-
hist(E, xlab = "Residuals", main = "")



################################################################################################
###Shoot mass-log transformed data looks good (sort of...)
####################################################################################################
#the samples sizes are slightly uneven-type 3 appropriate?
#Autoclaved, LowN, Liquid (9) and solid (10)=19
#Autoclaved, HighN, Liquid (11) and solid (11)=22
#Not autoclaved, LowN, Liquid (10) and solid (12)=22
#Not autoclaved, HighN, Liquid (11) and solid (11)=22

#Full Model shoot mass

#log transformation improves model fit
Mod.shootmass<-lmer (log( TotalShootMass) ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType + (1|G5Block), data=metadata.E)
summary(Mod.shootmass)
anova(Mod.shootmass)
#Analysis of Variance Table of type III  with  Satterthwaite 
#approximation for degrees of freedom

#                      Sum Sq Mean Sq NumDF  DenDF F.value  Pr(>F)    
#Nlevel              114.297 114.297     1 75.140 165.462 < 2e-16 ***
#Sterility             0.310   0.310     1 75.199   0.449 0.50466    
#PresenceofMoldAug31   0.149   0.149     1 76.385   0.216 0.64333    
#InoculumType          2.012   2.012     1 75.625   2.913 0.09196 .  
#Nlevel:Sterility      0.235   0.235     1 75.130   0.341 0.56109    
#remove inoculum type
  
#re-run emmeans or change to lsmeans (use other computer...)
emmeans(Mod.shootmass, list(pairwise~Nlevel), adjust = "tukey")
#High and low N are sig dif
#$`pairwise differences of Nlevel`
#contrast   estimate    SE   df t.ratio p.value
#High - Low     2.41 0.191 35.1 12.636  <.0001 


############################################################################################
##Assumptions of linear regression R.S 
############################################################################################


str(data)

M<-lmer (R.S  ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
plot(M, add.smooth = FALSE)
E<-resid(M)
E <- resid(M)
## Histogram of residual frequency-
hist(E, xlab = "Residuals", main = "")

#### log might fit a little better-tradeoff skewness vs. homoegeneity
Mlog<-lmer (log(R.S) ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
plot(Mlog, add.smooth = FALSE)
E<-resid(Mlog)
E <- resid(Mlog)
## Histogram of residual frequency
hist(E, xlab = "Residuals", main = "")

## QQplots of sample vs. theoretical quantiles 
qqnorm(metadata.E$R.S )
qqline(metadata.E$R.S )

#Ratios never behave nicely


################################################################################################
###R.S 
####################################################################################################
#Full Model R.S 

Mod.R.S <-lmer (log(R.S ) ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType + (1|G5Block) , data=metadata.E)
summary(Mod.R.S )
anova(Mod.R.S )


#Analysis of Variance Table of type III  with  Satterthwaite 
#approximation for degrees of freedom
#                      Sum Sq Mean Sq NumDF DenDF F.value    Pr(>F)    
#Nlevel              6.4863  6.4863     1    78 13.6059 0.0004152 ***
#Sterility           0.0371  0.0371     1    78  0.0778 0.7810581    
#PresenceofMoldAug31 0.0650  0.0650     1    78  0.1363 0.7130290    
#InoculumType        0.0740  0.0740     1    78  0.1553 0.6946046    
#Nlevel:Sterility    0.0388  0.0388     1    78  0.0813 0.7763161 



emmeans(Mod.R.S, list(pairwise~Nlevel*Sterility), adjust = "tukey")

#sig difs are between low and high N level


################################################################################################
###Belowground Traits
####################################################################################################






################################################################################################
###Core root length
####################################################################################################
M<-lmer (CoreRootLength  ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType + (1|G5Block), data=metadata.E)

plot(M, add.smooth = FALSE)
E<-resid(M)
E <- resid(M)
## Histogram of residual frequency-
hist(E, xlab = "Residuals", main = "")

#### log might fit a little better-tradeoff skewness vs. homoegeneity
Mlog<-lmer (log(CoreRootLength) ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
plot(Mlog, add.smooth = FALSE)
E<-resid(Mlog)
E <- resid(Mlog)
## Histogram of residual frequency
hist(E, xlab = "Residuals", main = "")

## QQplots of sample vs. theoretical quantiles 
qqnorm(metadata.E$CoreRootLength )
qqline(metadata.E$CoreRootLength )

#log transformation is an improvement
#####Full Model
Mod.CRL<-lmer (log(CoreRootLength)~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
summary(Mod.CRL)
anova(Mod.CRL)

########## N level is the important determinant of CRL 
#Analysis of Variance Table of type III  with  Satterthwaite 
#approximation for degrees of freedom
#Sum Sq Mean Sq NumDF  DenDF F.value    Pr(>F)    
#Nlevel              42.970  42.970     1 75.244  56.199 1.059e-10 ***
#  Sterility            0.000   0.000     1 75.363   0.000    0.9998    
#PresenceofMoldAug31  0.171   0.171     1 77.352   0.224    0.6376    
#InoculumType         0.616   0.616     1 76.154   0.806    0.3722    
#Nlevel:Sterility     0.451   0.451     1 75.225   0.590    0.4449    




emmeans(Mod.CRL, list(pairwise~Nlevel*Sterility), adjust = "tukey")

IC.emm <- emmeans(Mod.CRL, ~ Nlevel*Sterility)

contrast(IC.emm, "consec", simple = "each", combine = TRUE, adjust = "mvt")

pairs(IC.emm, simple = "Sterility", adjust = "mvt")

#Tukey's test=high nitrogen treatments between perturbed (sterile) and unperturbed (not sterile)
# are not significantly different
################################################################################################
###Core_Mean_Root_Diameter
####################################################################################################
M<-lmer(CoreMeanRootDiameter  ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
plot(M, add.smooth = FALSE)
E<-resid(M)
E <- resid(M)
## Histogram of residual frequency-
hist(E, xlab = "Residuals", main = "")

#### log might fit a little better-tradeoff skewness vs. homoegeneity
Mlog<-lmer (log(CoreMeanRootDiameter) ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
plot(Mlog, add.smooth = FALSE)
E<-resid(Mlog)
E <- resid(Mlog)
## Histogram of residual frequency
hist(E, xlab = "Residuals", main = "")

## QQplots of sample vs. theoretical quantiles 
qqnorm(metadata.E$CoreMeanRootDiameter)
qqline(metadata.E$CoreMeanRootDiameter)
#so nicely behaved :)

#####Full Model
Mod.CRD<-lmer (CoreMeanRootDiameter ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
summary(Mod.CRD)
anova(Mod.CRD)


##############################################

#Analysis of Variance Table of type III  with  Satterthwaite 
#approximation for degrees of freedom
#                      Sum Sq    Mean Sq NumDF  DenDF F.value   Pr(>F)   
#Nlevel              9.8397e-05 9.8397e-05     1 75.388  9.3104 0.003145 **
#Sterility           1.2377e-05 1.2377e-05     1 75.631  1.1711 0.282615   
#PresenceofMoldAug31 8.9100e-07 8.9100e-07     1 78.164  0.0843 0.772324   
#InoculumType        2.4100e-06 2.4100e-06     1 77.014  0.2281 0.634315   
#Nlevel:Sterility    1.4045e-05 1.4045e-05     1 75.349  1.3289 0.252641   


################################################################################################
###Core_Specific_Root_Length
####################################################################################################
M<-lmer(CoreSpecificRootLength  ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
plot(M, add.smooth = FALSE)
E<-resid(M)
E <- resid(M)
## Histogram of residual frequency-
hist(E, xlab = "Residuals", main = "")

#### log might fit a little better-tradeoff skewness vs. homoegeneity
Mlog<-lmer (log(CoreSpecificRootLength) ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
plot(Mlog, add.smooth = FALSE)
E<-resid(Mlog)
E <- resid(Mlog)
## Histogram of residual frequency
hist(E, xlab = "Residuals", main = "")

## QQplots of sample vs. theoretical quantiles 
qqnorm(metadata.E$CoreSpecificRootLength)
qqline(metadata.E$CoreSpecificRootLength)

#also very well behaved
#####Full Model
Mod.SRL<-lmer (CoreSpecificRootLength ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType + (1|G5Block), data=metadata.E)
summary(Mod.SRL)
anova(Mod.SRL)



#                       Sum Sq Mean Sq NumDF  DenDF F.value   Pr(>F)   
#Nlevel              1649.34 1649.34     1 74.873  7.6891 0.007008 **
#Sterility            300.12  300.12     1 75.175  1.3991 0.240595   
#PresenceofMoldAug31    8.53    8.53     1 77.944  0.0397 0.842491   
#InoculumType         329.31  329.31     1 76.834  1.5352 0.219107   
#Nlevel:Sterility      21.84   21.84     1 74.826  0.1018 0.750534 

################################################################################################
##Root Tissue Density
####################################################################################################
M<-lmer(Root_Tissue_Density   ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
plot(M, add.smooth = FALSE)
E<-resid(M)
E <- resid(M)
## Histogram of residual frequency-
hist(E, xlab = "Residuals", main = "")

#### log might fit a little better-tradeoff skewness vs. homoegeneity
Mlog<-lmer (log(Root_Tissue_Density )  ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
plot(Mlog, add.smooth = FALSE)
E<-resid(Mlog)
E <- resid(Mlog)
## Histogram of residual frequency
hist(E, xlab = "Residuals", main = "")

## QQplots of sample vs. theoretical quantiles 
qqnorm(metadata.E$Root_Tissue_Density)
qqline(metadata.E$Root_Tissue_Density)


#####Full Model
Mod.RTD<-lmer (Root_Tissue_Density  ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
summary(Mod.RTD)
anova(Mod.RTD)
########## Nlevel
#                       Sum Sq    Mean Sq NumDF  DenDF F.value Pr(>F)
#Nlevel              0.00007870 0.00007870     1 75.025 0.22256 0.6385
#Sterility           0.00066817 0.00066817     1 75.175 1.88953 0.1733
#PresenceofMoldAug31 0.00005224 0.00005224     1 77.529 0.14772 0.7018
# InoculumType        0.00064481 0.00064481     1 76.148 1.82346 0.1809
# Nlevel:Sterility    0.00000684 0.00000684     1 75.001 0.01936 0.8897



################################################################################################
###CoreRootMass 
####################################################################################################
M<-lmer(CoreRootMass  ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
plot(M, add.smooth = FALSE)
E<-resid(M)
E <- resid(M)
## Histogram of residual frequency-
hist(E, xlab = "Residuals", main = "")

#### log might fit a little better-tradeoff skewness vs. homoegeneity
Mlog<-lmer (log(CoreRootMass ) ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
plot(Mlog, add.smooth = FALSE)
E<-resid(Mlog)
E <- resid(Mlog)
## Histogram of residual frequency
hist(E, xlab = "Residuals", main = "")

## QQplots of sample vs. theoretical quantiles 
qqnorm(metadata.E$CoreRootMass)
qqline(metadata.E$CoreRootMass)

#####Full Model
Mod.CRM<-lmer (log(CoreRootMass) ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType +(1|G5Block), data=metadata.E)
summary(Mod.CRM)
anova(Mod.CRM)
########## Nlevel
#                    Sum Sq Mean Sq NumDF  DenDF F.value    Pr(>F)    
#Nlevel              65.873  65.873     1 75.290  74.601 7.092e-13 ***
#Sterility            0.533   0.533     1 75.405   0.603    0.4398    
#PresenceofMoldAug31  0.247   0.247     1 77.334   0.280    0.5982    
#InoculumType         2.536   2.536     1 76.168   2.872    0.0942 .  
#Nlevel:Sterility     0.078   0.078     1 75.271   0.089    0.7667     








#########################################
###Re-run with mold term removed
##########################################


################################################################################################
###Shoot mass-log transformed data looks good (sort of...)
####################################################################################################
#the samples sizes are slightly uneven-type 3 appropriate?
#Autoclaved, LowN, Liquid (9) and solid (10)
#Autoclaved, HighN, Liquid (11) and solid (11)
#Not autoclaved, LowN, Liquid (10) and solid (12)
#Not autoclaved, HighN, Liquid (11) and solid (11)

#Full Model shoot mass

#log transformation improves model fit
Mod.shootmass<-lmer (log( TotalShootMass) ~ Nlevel*Sterility  + InoculumType +(1|G5Block), data=metadata.E)
summary(Mod.shootmass)
anova(Mod.shootmass)
#Analysis of Variance Table of type III  with  Satterthwaite 
#approximation for degrees of freedom

#                   Sum Sq Mean Sq NumDF  DenDF F.value Pr(>F)    
#Nlevel           116.779 116.779     1 76.065 170.671 <2e-16 ***
#  Sterility          0.320   0.320     1 76.186   0.467 0.4964    
#InoculumType       1.864   1.864     1 76.246   2.724 0.1029    
#Nlevel:Sterility   0.277   0.277     1 76.192   0.405 0.5265     
#remove inoculum type



################################################################################################
###R.S 
####################################################################################################
#Full Model R.S 

Mod.R.S <-lmer (log(R.S ) ~ Nlevel*Sterility  + InoculumType +(1|G5Block), data=metadata.E)
summary(Mod.R.S )
anova(Mod.R.S )


#Analysis of Variance Table of type III  with  Satterthwaite 
#approximation for degrees of freedom
#                      Sum Sq Mean Sq NumDF DenDF F.value    Pr(>F)    
#Nlevel           6.9484  6.9484     1    79 14.7361 0.0002479 ***
#  Sterility        0.0388  0.0388     1    79  0.0823 0.7750011    
#InoculumType     0.0434  0.0434     1    79  0.0921 0.7622751    
#Nlevel:Sterility 0.0289  0.0289     1    79  0.0612 0.8052506   





################################################################################################
###Belowground Traits
####################################################################################################






################################################################################################
###Core root length
####################################################################################################

#####Full Model
Mod.CRL<-lmer (log(CoreRootLength)~ Nlevel*Sterility + InoculumType +(1|G5Block), data=metadata.E)
summary(Mod.CRL)
anova(Mod.CRL)

########## N level is the important determinant of CRL 
#Analysis of Variance Table of type III  with  Satterthwaite 
#approximation for degrees of freedom
#Sum Sq Mean Sq NumDF  DenDF F.value    Pr(>F)    
#Nlevel           43.450  43.450     1 76.074  57.432 6.996e-11 ***
#  Sterility         0.000   0.000     1 76.311   0.000    0.9920    
#InoculumType      0.481   0.481     1 76.427   0.635    0.4278    
#Nlevel:Sterility  0.514   0.514     1 76.322   0.679    0.4123    


################################################################################################
###Core_Mean_Root_Diameter
####################################################################################################

#####Full Model
Mod.CRD<-lmer (CoreMeanRootDiameter ~ Nlevel*Sterility  + InoculumType +(1|G5Block), data=metadata.E)
summary(Mod.CRD)
anova(Mod.CRD)


##############################################

#Analysis of Variance Table of type III  with  Satterthwaite 
#approximation for degrees of freedom
#                    Sum Sq    Mean Sq NumDF  DenDF F.value   Pr(>F)   
#Nlevel           1.0533e-04 1.0533e-04     1 76.170 10.1023 0.002141 **
#  Sterility        1.2229e-05 1.2229e-05     1 76.624  1.1728 0.282222   
#InoculumType     3.6950e-06 3.6950e-06     1 76.843  0.3544 0.553370   
#Nlevel:Sterility 1.4911e-05 1.4911e-05     1 76.642  1.4300 0.235445   


################################################################################################
###Core_Specific_Root_Length
####################################################################################################

#####Full Model
Mod.SRL<-lmer (CoreSpecificRootLength ~ Nlevel*Sterility  + InoculumType +(1|G5Block), data=metadata.E)
summary(Mod.SRL)
anova(Mod.SRL)



#                       Sum Sq Mean Sq NumDF  DenDF F.value   Pr(>F)   
#Nlevel           1660.69 1660.69     1 75.687  7.8229 0.006537 **
#  Sterility         303.63  303.63     1 76.281  1.4303 0.235425   
#InoculumType      323.36  323.36     1 76.565  1.5232 0.220904   
#Nlevel:Sterility   19.77   19.77     1 76.304  0.0931 0.761091  

################################################################################################
##Root Tissue Density
####################################################################################################

#####Full Model
Mod.RTD<-lmer (Root_Tissue_Density  ~ Nlevel*Sterility + InoculumType +(1|G5Block), data=metadata.E)
summary(Mod.RTD)
anova(Mod.RTD)
########## Nlevel
#                       Sum Sq    Mean Sq NumDF  DenDF F.value Pr(>F)
#Nlevel           0.00010775 0.00010775     1 75.928 0.30735 0.5809
#Sterility        0.00067787 0.00067787     1 76.237 1.93360 0.1684
#InoculumType     0.00058941 0.00058941     1 76.389 1.68127 0.1987
#Nlevel:Sterility 0.00001088 0.00001088     1 76.251 0.03102 0.8607



################################################################################################
###CoreRootMass 
####################################################################################################


#####Full Model
Mod.CRM<-lmer (log(CoreRootMass) ~ Nlevel*Sterility  + InoculumType +(1|G5Block), data=metadata.E)
summary(Mod.CRM)
anova(Mod.CRM)
########## Nlevel
#                    Sum Sq Mean Sq NumDF  DenDF F.value    Pr(>F)    
#Nlevel           66.653  66.653     1 76.138  76.139 4.414e-13 ***
#  Sterility         0.547   0.547     1 76.372   0.625    0.4318    
#InoculumType      2.288   2.288     1 76.487   2.614    0.1100    
#Nlevel:Sterility  0.111   0.111     1 76.383   0.126    0.7233

###Re-run with inoculum type term removed
##########################################


################################################################################################

#Now we arrive at the analysis in the paper

################################################################################################

# ==== CJ is frowning now that she has finally found where the real analysis starts! ====

################################################################################################
###Shoot mass-log transformed data looks good (sort of...)
####################################################################################################
#the samples sizes are slightly uneven-type 3 appropriate?
#Autoclaved, LowN, Liquid (9) and solid (10)=19
#Autoclaved, HighN, Liquid (11) and solid (11)=22
#Not autoclaved, LowN, Liquid (10) and solid (12)=22
#Not autoclaved, HighN, Liquid (11) and solid (11)=22

#Full Model shoot mass

#log transformation improves model fit
Mod.shootmass<-lmer (log( TotalShootMass) ~ Nlevel*Sterility  +(1|G5Block), data=metadata.E)
summary(Mod.shootmass)
anova(Mod.shootmass)
#Analysis of Variance Table of type III  with  Satterthwaite 
#approximation for degrees of freedom

#                   Sum Sq Mean Sq NumDF  DenDF F.value Pr(>F)    
#Nlevel           116.160 116.160     1 77.068 165.855 <2e-16 ***
#  Sterility          0.358   0.358     1 77.194   0.511 0.4770    
#Nlevel:Sterility   0.304   0.304     1 77.209   0.435 0.5117    

#Mod.shootmass<-lmer (log( TotalShootMass) ~ Nlevel*Sterility + PresenceofMoldAug31 + InoculumType + (1|G5Block), data=data)
IC.emm <- emmeans(Mod.shootmass, ~ Nlevel*Sterility)

contrast(IC.emm, "consec", simple = "each", combine = TRUE, adjust = "mvt")
# Sterility   Nlevel contrast                estimate    SE  df t.ratio p.value
# Unperturbed .      Low - High               -2.1540 0.164 259 -13.133 <.0001 
# Perturbed   .      Low - High               -2.1591 0.164 259 -13.131 <.0001 
# .           High   Perturbed - Unperturbed  -0.0866 0.163 259  -0.532 0.9732 
# .           Low    Perturbed - Unperturbed  -0.0917 0.163 259  -0.564 0.9668 
# 
# Results are averaged over some or all of the levels of: InoculumType 
# Results are given on the log (not the response) scale. 
# P value adjustment: sidak method for 4 tests 

################################################################################################
###R.S 
####################################################################################################
#Full Model R.S 

Mod.R.S <-lmer (log(R.S ) ~ Nlevel*Sterility +(1|G5Block), data=metadata.E)
summary(Mod.R.S )
anova(Mod.R.S )


#Analysis of Variance Table of type III  with  Satterthwaite 
#approximation for degrees of freedom
#                      Sum Sq Mean Sq NumDF DenDF F.value    Pr(>F)    
#Nlevel           6.9769  6.9769     1    80 14.9664 0.000222 ***
#  Sterility        0.0407  0.0407     1    80  0.0873 0.768412    
#Nlevel:Sterility 0.0273  0.0273     1    80  0.0585 0.809440  

IC.emm <- emmeans(Mod.R.S, ~ Nlevel*Sterility)
contrast(IC.emm, "consec", simple = "each", combine = TRUE, adjust = "mvt")

# Sterility   Nlevel contrast                estimate    SE   df t.ratio p.value
# Unperturbed .      Low - High                0.5424 0.206 77.3  2.631  0.0357 
# Perturbed   .      Low - High                0.6148 0.217 77.4  2.829  0.0214 
# .           High   Perturbed - Unperturbed  -0.0804 0.207 77.9 -0.389  0.9737 
# .           Low    Perturbed - Unperturbed  -0.0080 0.217 77.4 -0.037  1.0000 
# 
# Results are given on the log (not the response) scale. 
# P value adjustment: mvt method for 4 tests

################################################################################################
###Belowground Traits
####################################################################################################






################################################################################################
###Core root length
####################################################################################################

#####Full Model
Mod.CRL<-lmer (log(CoreRootLength)~ Nlevel*Sterility  +(1|G5Block), data=metadata.E)
summary(Mod.CRL)
anova(Mod.CRL)

########## N level is the important determinant of CRL 
#Analysis of Variance Table of type III  with  Satterthwaite 
#approximation for degrees of freedom
#Sum Sq Mean Sq NumDF  DenDF F.value    Pr(>F)    
#Nlevel           43.264  43.264     1 77.063  57.463 6.481e-11 ***
 # Sterility         0.001   0.001     1 77.296   0.001    0.9764    
#Nlevel:Sterility  0.534   0.534     1 77.324   0.709    0.4025  

IC.emm <- emmeans(Mod.CRL, ~ Nlevel*Sterility)
# contrast(IC.emm, "consec", simple = "each", combine = TRUE, adjust = "mvt")
# 
# Sterility   Nlevel contrast                estimate    SE   df t.ratio p.value
# Unperturbed .      Low - High                -1.281 0.262 77.1 -4.888  <.0001 
# Perturbed   .      Low - High                -1.601 0.276 77.1 -5.797  <.0001 
# .           High   Perturbed - Unperturbed    0.155 0.263 77.4  0.589  0.9170 
# .           Low    Perturbed - Unperturbed   -0.166 0.276 77.1 -0.601  0.9124 
# 
# Results are given on the log (not the response) scale. 
# P value adjustment: mvt method for 4 tests

################################################################################################
###Core_Mean_Root_Diameter
####################################################################################################

#####Full Model
Mod.CRD<-lmer (CoreMeanRootDiameter ~ Nlevel*Sterility + (1|G5Block), data=metadata.E)
summary(Mod.CRD)
anova(Mod.CRD)


##############################################

#Analysis of Variance Table of type III  with  Satterthwaite 
#approximation for degrees of freedom
#                    Sum Sq    Mean Sq NumDF  DenDF F.value   Pr(>F)   
#Nlevel           1.0448e-04 1.0448e-04     1 77.193 10.0874 0.002146 **
#  Sterility        1.2572e-05 1.2572e-05     1 77.656  1.2139 0.273966   
#Nlevel:Sterility 1.4660e-05 1.4660e-05     1 77.707  1.4155 0.237774   

IC.emm <- emmeans(Mod.CRD, ~ Nlevel*Sterility)
contrast(IC.emm, "consec", simple = "each", combine = TRUE, adjust = "mvt")

# Sterility   Nlevel contrast                 estimate       SE   df t.ratio p.value
# Unperturbed .      Low - High              -3.08e-03 0.000972 77.3 -3.168  0.0081 
# Perturbed   .      Low - High              -1.40e-03 0.001024 77.3 -1.366  0.4601 
# .           High   Perturbed - Unperturbed -1.62e-03 0.000974 77.7 -1.659  0.2945 
# .           Low    Perturbed - Unperturbed  6.21e-05 0.001024 77.3  0.061  0.9999 
# 
# P value adjustment: mvt method for 4 tests 

################################################################################################
###Core_Specific_Root_Length
####################################################################################################

#####Full Model
Mod.SRL<-lmer (CoreSpecificRootLength ~ Nlevel*Sterility +(1|G5Block), data=metadata.E)
summary(Mod.SRL)
anova(Mod.SRL)



#                       Sum Sq Mean Sq NumDF  DenDF F.value   Pr(>F)   
# Nlevel           1627.44 1627.44     1 76.665  7.5809 0.007362 **
#  Sterility         322.72  322.72     1 77.332  1.5033 0.223884   
#  Nlevel:Sterility   17.85   17.85     1 77.402  0.0831 0.773870  


IC.emm <- emmeans(Mod.SRL, ~ Nlevel*Sterility)
contrast(IC.emm, "consec", simple = "each", combine = TRUE, adjust = "mvt")
# Sterility   Nlevel contrast                estimate   SE   df t.ratio p.value
# Unperturbed .      Low - High                  9.76 4.42 77.3 2.207   0.1008 
# Perturbed   .      Low - High                  7.91 4.66 77.4 1.696   0.2764 
# .           High   Perturbed - Unperturbed     4.86 4.44 77.9 1.096   0.6367 
# .           Low    Perturbed - Unperturbed     3.01 4.66 77.4 0.645   0.8944 
# 
# P value adjustment: mvt method for 4 tests 

################################################################################################
##Root Tissue Density
####################################################################################################

#####Full Model
Mod.RTD<-lmer (Root_Tissue_Density  ~ Nlevel*Sterility + (1|G5Block), data=metadata.E)
summary(Mod.RTD)
anova(Mod.RTD)
########## Nlevel
#                       Sum Sq    Mean Sq NumDF  DenDF F.value Pr(>F)
#Nlevel           0.00011974 0.00011974     1 76.934 0.33785 0.5628
#Sterility        0.00071104 0.00071104     1 77.261 2.00612 0.1607
#Nlevel:Sterility 0.00001360 0.00001360     1 77.298 0.03838 0.8452


IC.emm <- emmeans(Mod.RTD, ~ Nlevel*Sterility)
contrast(IC.emm, "consec", simple = "each", combine = TRUE, adjust = "mvt")

# Sterility   Nlevel contrast                estimate      SE   df t.ratio p.value
# Unperturbed .      Low - High               0.00321 0.00568 77.2  0.564  0.9260 
# Perturbed   .      Low - High               0.00159 0.00599 77.2  0.265  0.9913 
# .           High   Perturbed - Unperturbed -0.00504 0.00570 77.5 -0.884  0.7709 
# .           Low    Perturbed - Unperturbed -0.00666 0.00599 77.2 -1.111  0.6267 
# 
# P value adjustment: mvt method for 4 tests

################################################################################################
###CoreRootMass 
####################################################################################################


#####Full Model
Mod.CRM<-lmer (log(CoreRootMass) ~ Nlevel*Sterility +(1|G5Block), data=metadata.E)
summary(Mod.CRM)
anova(Mod.CRM)
########## Nlevel
#                    Sum Sq Mean Sq NumDF  DenDF F.value    Pr(>F)    
#Nlevel           66.120  66.120     1 77.139  73.894 7.112e-13 ***
#  Sterility         0.600   0.600     1 77.383   0.671    0.4153    
#Nlevel:Sterility  0.131   0.131     1 77.411   0.147    0.7026 

IC.emm <- emmeans(Mod.CRM, ~ Nlevel*Sterility)
contrast(IC.emm, "consec", simple = "each", combine = TRUE, adjust = "mvt")

# Sterility   Nlevel contrast                estimate    SE   df t.ratio p.value
# Unperturbed .      Low - High               -1.7016 0.286 77.1 -5.958  <.0001 
# Perturbed   .      Low - High               -1.8607 0.301 77.2 -6.180  <.0001 
# .           High   Perturbed - Unperturbed  -0.0904 0.286 77.4 -0.316  0.9855 
# .           Low    Perturbed - Unperturbed  -0.2494 0.301 77.2 -0.829  0.8031 

##################################################
#adjusting for mc
###############################################

########################################
## adjustment for mc
########################################




#Pvalues for OTUs =14 comparisons
Ntrait<-c(0.0000000000000002,0.000222, 0.00000000006481, 0.002146, 0.007362, 0.5628, 0.0000000000007112 )


BHNtrait = p.adjust(Ntrait, method = "BH")
print(BHNtrait)

# 1.400000e-15 3.885000e-04 1.512233e-10 3.004400e-03 8.589000e-03 5.628000e-01 2.489200e-12
#all still sig after correcting for multiple comparisons

# ==== CJ: remaking figures ====
#########

metadata.E <- read.csv("~/Dropbox/PanicumMicrobiome_CleanFolder/data/Traits/metadata.E2.csv")
myColors <- c("grey36", "gray74")
names(myColors) <- levels(metadata.E$Nlevel)
colScale <- scale_fill_manual(name = "Nitrogen \nLevel", values = myColors)

table(metadata.E$Sterility)
levels(metadata.E$Sterility) <- c("Unperturbed", "Perturbed")
#theme_set(theme_cowplot())
##########

# ==== Total Shoot mass fig ====
TSMfig <- ggplot(metadata.E, aes(Sterility, TotalShootMass, fill = Nlevel)) + geom_boxplot() + colScale  + theme(legend.position = "none", axis.title.x = element_blank()) + labs(y = "Shoot mass (g)") + theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 13)) + geom_signif(xmin = c(0.8, 1.8), xmax = c(1.2,2.2), y_position = 2.5, annotation = c("***", "***"), tip_length = 0.015, textsize = 5)



# ==== R:S fig ====
RSfig <-ggplot(metadata.E, aes(Sterility, R.S, fill = Nlevel)) + geom_boxplot() + colScale  + labs(y = "Root : Shoot") + theme(legend.position = "none", axis.title.x = element_blank()) + theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 13))+ geom_signif(xmin = c(0.8, 1.8), xmax = c(1.2,2.2), y_position =0.2, annotation = c("*", "*"), tip_length = 0.015, textsize = 5) + coord_cartesian(ylim = c(0,.25))

# ==== CRL fig ====
CRLfig <- ggplot(metadata.E, aes(Sterility, CoreRootLength, fill = Nlevel)) + geom_boxplot() + colScale  + labs(y = "Root length (cm)") + theme(legend.position = "none", axis.title.x = element_blank()) + theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 13))+ geom_signif(xmin = c(0.8, 1.8), xmax = c(1.2,2.2), y_position =2000, annotation = c("***", "***"), tip_length = 0.015, textsize = 5)

# ==== CRD fig ====
CRDfig <- ggplot(metadata.E, aes(Sterility, CoreMeanRootDiameter, fill = Nlevel)) + geom_boxplot() + colScale  + labs(y = "Root Diameter (cm)") + theme(legend.position = "none", axis.title.x = element_blank())+ theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 13))+ geom_signif(xmin = 0.8, xmax = 1.2, y_position =0.035, annotation = "**", tip_length = 0.015, textsize = 5)+ coord_cartesian(ylim=c(0, .04))

# ==== CSRL ====
CSRLfig <- ggplot(metadata.E, aes(Sterility, CoreSpecificRootLength, fill = Nlevel)) + geom_boxplot()  + colScale + labs(y =  expression(paste("Specific root length (g ", cm^{-3}, ")"))) + theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 13)) + theme(legend.position = "none", axis.title.x = element_blank())# + geom_signif(xmin = c(0.8, 1.8), xmax = c(1.2,2.2), y_position =75, annotation = c("NS", "NS"), tip_length = 0.015, textsize = 5)

# ==== CRM ====

CRMfig <- ggplot(metadata.E, aes(Sterility, CoreRootMass, fill = Nlevel)) + geom_boxplot() + colScale + labs(y = "Root mass (g)") + theme(legend.position = "none", axis.title.x = element_blank()) + theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 13))+ geom_signif(xmin = c(0.8, 1.8), xmax = c(1.2,2.2), y_position =0.08, annotation = c("***", "***"), tip_length = 0.015, textsize = 5)

# ==== RTD ====

RTDfig <-ggplot(metadata.E, aes(Sterility, Root_Tissue_Density, fill = Nlevel)) + geom_boxplot()  + colScale  + labs(y = expression(paste("Root tissue density (g ", cm^{-3}, ")"))) + theme(axis.title.y = element_text(size = 18), axis.text = element_text(size = 13))+ theme(axis.title.x = element_blank())# + geom_signif(xmin = c(0.8, 1.8), xmax = c(1.2,2.2), y_position =.13, annotation = c("NS", "NS"), tip_length = 0.015, textsize = 3)

# put all together
Mylegend <- get_legend(RTDfig)
p1 <-plot_grid(TSMfig, RSfig, CRLfig, CRDfig, CSRLfig, CRMfig, RTDfig + theme(legend.position = "none"), Mylegend, ncol = 4, labels = c("A", "B", "C", "D", "E", "F", "G"), align = "hv", axis = "l")

title <- ggdraw() + draw_label("Microbial community", fontface = "bold", fontfamily = "")

p2 <-plot_grid(p1,title, ncol=1, rel_heights=c(1, .07)) # rel_heights values control title margin
#ggdraw(add_sub(p2, "Label", vpadding=grid::unit(0,"lines"),y=-0.5, x=0.5, vjust=0

save_plot("Traits5Sept2019.tiff", dpi = 1000, p2, ncol = 1, nrow = 2, base_width = 14)
