#title: "OTU effects on Panicum phenotype using top relative abundances"
#author: "Renee H. Petipas"
#date: August 5, 2019
#To clear the R environment
#Re-running top abundance by OTU analysis with relative abundance and transformed values (consistent with transformations in trait analysis) and also adding random effects
rm(list=ls(all=TRUE))
#Load libraries
library(phyloseq)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(utils)
library(stargazer)
library(car)
library(kableExtra)
library(gridExtra)
library(vegan)
library(microbiome)
library(lme4)
library(lmerTest)
#I couldn't read tree before I added this package (WTF?)
library(ape)

###################################################################

#Belowground traits
#
######################################################################
# 
setwd("~/Dropbox/PanicumMicrobiome_CleanFolder/data/MicrobiomexTrait")
df<-read.csv("MetaTableTotal_RelAbun_FINAL_05Aug2019.csv", row.names = 1)
str(df)
#Trying this analysis again using a filtered data set (ask Emily for specifications)
df$G5Block<-as.factor(df$G5Block)
df$PresenceofMoldAug31 <-as.factor(df$PresenceofMoldAug31 )
df$InoculumVolume<-as.factor(df$InoculumVolume)

#Calculate R.S and RTD
df$RootTissueDensity  <- df$CoreRootMass /df$CoreRootVolume
df$R.S <- df$CoreRootMass /df$TotalShootMass


########################################
## CoreRootLength
#log transform to be consistent with trait analysis
########################################

#OTU 1712
#  X726a6245b76a9e74caea7993598fd969
#Wow goes from being significant to NS
#I played around with graphing and it lookslike block four has shortest roots and lowest abundance of this OTU so 
#probably what is driving the relationship :0 OMG!!


options(contrasts=c('contr.sum','contr.poly'))
modotu1<-lmer(log(CoreRootLength)~X726a6245b76a9e74caea7993598fd969*Nlevel*Sterility+ (1|G5Block), data=df)
summary(modotu1)
anova(modotu1)


#Type III Analysis of Variance Table with Satterthwaite's method
#                                                    Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#X726a6245b76a9e74caea7993598fd969                   0.0035  0.0035     1 29.481  0.0065 0.9361371    
#Nlevel                                             10.4105 10.4105     1 29.261 19.5906 0.0001225 ***
#Sterility                                           0.0150  0.0150     1 29.534  0.0282 0.8678686    
#X726a6245b76a9e74caea7993598fd969:Nlevel            0.0302  0.0302     1 29.528  0.0569 0.8131362    
#X726a6245b76a9e74caea7993598fd969:Sterility         0.0000  0.0000     1 30.363  0.0000 0.9990972    
#Nlevel:Sterility                                    0.7173  0.7173     1 29.670  1.3497 0.2545857    
#X726a6245b76a9e74caea7993598fd969:Nlevel:Sterility  0.1096  0.1096     1 30.671  0.2062 0.6529449                       


#Nlevel effects



               
#OTU 1301
#  X7d3587bbb834c6ac86656e1d2715c331

modotu2<-lmer (log(CoreRootLength)~ X7d3587bbb834c6ac86656e1d2715c331*Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu2)
anova(modotu2)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#X7d3587bbb834c6ac86656e1d2715c331                  0.1169  0.1169     1 30.535  0.2526 0.6188258    
#Nlevel                                             8.9190  8.9190     1 31.033 19.2675 0.0001222 ***
#Sterility                                          0.9813  0.9813     1 29.721  2.1199 0.1558785    
#X7d3587bbb834c6ac86656e1d2715c331:Nlevel           0.4180  0.4180     1 29.906  0.9030 0.3495872    
#X7d3587bbb834c6ac86656e1d2715c331:Sterility        0.3459  0.3459     1 30.547  0.7473 0.3940850    
#Nlevel:Sterility                                   3.6365  3.6365     1 31.226  7.8560 0.0086281 ** 
#X7d3587bbb834c6ac86656e1d2715c331:Nlevel:Sterility 0.8086  0.8086     1 29.940  1.7468 0.1962857 
         
#Nlevel effects and Nlevel*Sterility NS after correction

#OTU 1857
#a90d5b4368338dea15569269e8d33fb8 

modotu3<-lmer(log(CoreRootLength)~a90d5b4368338dea15569269e8d33fb8 *Nlevel*Sterility +  (1|G5Block), data=df)
summary(modotu3)
anova(modotu3)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
#a90d5b4368338dea15569269e8d33fb8                  2.29462 2.29462     1 30.772  4.9416 0.03370 *
#Nlevel                                            1.58460 1.58460     1 28.710  3.4125 0.07503 .
#Sterility                                         0.02458 0.02458     1 29.599  0.0529 0.81962  
#a90d5b4368338dea15569269e8d33fb8:Nlevel           2.68379 2.68379     1 13.501  5.7797 0.03121 *
#a90d5b4368338dea15569269e8d33fb8:Sterility        0.02227 0.02227     1 31.041  0.0480 0.82807  
#Nlevel:Sterility                                  1.26311 1.26311     1 31.969  2.7202 0.10888  
#a90d5b4368338dea15569269e8d33fb8:Nlevel:Sterility 0.50530 0.50530     1 26.539  1.0882 0.30628       

#OTU 1857 significant and OTUxNlevel but NS after correcting for MC


#####################################
#This OTU was no longer in the most abundant data set
#################################
#OTU 416
#  bea3f263c468bf3d3fb6645562737c70 

#modotu4<-lm (log(CoreRootLength)~ bea3f263c468bf3d3fb6645562737c70 *Nlevel*Sterility, data=df)
#summary(modotu4)
#anova(modotu4)
   
#bea3f263c468bf3d3fb6645562737c70                   1   75467   75467  0.3782 0.542906    
#Nlevel                                             1 2869012 2869012 14.3787 0.000626 ***
#  Sterility                                          1  347361  347361  1.7409 0.196390    
#bea3f263c468bf3d3fb6645562737c70:Nlevel            1    8609    8609  0.0431 0.836769    
#bea3f263c468bf3d3fb6645562737c70:Sterility         1    1437    1437  0.0072 0.932894    
#Nlevel:Sterility                                   1  550347  550347  2.7582 0.106526    
#bea3f263c468bf3d3fb6645562737c70:Nlevel:Sterility  1   29324   29324  0.1470 0.703987 

#Nlevel effects


##################################
#only OTU that is significant after correcting for MC
#####################################
#OTU 3749
# X071e2a4c88faf23787b614210e527ce3
#presence of mold is NS
modotu5<-lmer (log(CoreRootLength)~ X071e2a4c88faf23787b614210e527ce3*Nlevel*Sterility +  (1|G5Block), data=df)
summary(modotu5)
anova(modotu5)

#Type III Analysis of Variance Table with Satterthwaite's method
                                                       #Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#X071e2a4c88faf23787b614210e527ce3                  4.2761  4.2761     1 31.456  9.7208 0.0038780 ** 
#Nlevel                                             7.6159  7.6159     1 29.773 17.3130 0.0002478 ***
#Sterility                                          0.0026  0.0026     1 31.938  0.0059 0.9390521    
#X071e2a4c88faf23787b614210e527ce3:Nlevel           1.9039  1.9039     1 28.657  4.3281 0.0465390 *  
#X071e2a4c88faf23787b614210e527ce3:Sterility        0.0052  0.0052     1 31.969  0.0118 0.9142609    
#Nlevel:Sterility                                   1.5117  1.5117     1 29.844  3.4365 0.0736795 .  
#X071e2a4c88faf23787b614210e527ce3:Nlevel:Sterility 0.1335  0.1335     1 29.720  0.3034 0.5858810      

#OTU 3749 is significant even after correcting for MC
#Also N level effects
#OTU by N interaction but NS after correcting for MC


modN<-lm (X071e2a4c88faf23787b614210e527ce3~Nlevel*Sterility, data=df)
summary(modN)
anova(modN)

#Relative abundance of Xo7 is determined by N level 
head(df)
p_rl2 <- ggplot(df, aes(x= X071e2a4c88faf23787b614210e527ce3, y=CoreRootLength)) + geom_point(aes(x=X071e2a4c88faf23787b614210e527ce3, y=CoreRootLength, shape=Sterility, color=Nlevel)) + theme_classic() + xlab("Micromonospora sp.") + ylab("Root Length")+ stat_smooth(method= "lm", col = "black")
p_rl2


plot(log(CoreRootLength)~ X071e2a4c88faf23787b614210e527ce3, data=df)

#Micromonospora still negatively correlated with root length but in both sterile and non-sterile treatments

#abundance of X07 is affected by N level...



#OTU 5891
#X2bb5ef7086fb6a049aa217591531d386 

modotu6<-lmer (log(CoreRootLength)~ X2bb5ef7086fb6a049aa217591531d386 *Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu6)
anova(modotu6)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                    Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#X2bb5ef7086fb6a049aa217591531d386                   0.7022  0.7022     1 29.350  1.3868    0.2484    
#Nlevel                                             10.2276 10.2276     1 31.813 20.1999 8.671e-05 ***
#Sterility                                           0.1303  0.1303     1 29.682  0.2574    0.6157    
#X2bb5ef7086fb6a049aa217591531d386:Nlevel            1.1126  1.1126     1 31.895  2.1974    0.1481    
#X2bb5ef7086fb6a049aa217591531d386:Sterility         0.0964  0.0964     1 30.177  0.1903    0.6658    
#Nlevel:Sterility                                    0.6552  0.6552     1 31.245  1.2941    0.2639    
#X2bb5ef7086fb6a049aa217591531d386:Nlevel:Sterility  0.4161  0.4161     1 30.482  0.8217    0.3718  

#Nlevel effects


#OTU 718
#X3899d68a62bd9b58f6d28de79e6f5d1d

modotu7<-lmer (log(CoreRootLength)~ X3899d68a62bd9b58f6d28de79e6f5d1d*Nlevel*Sterility +  (1|G5Block), data=df)
summary(modotu7)
anova(modotu7)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
#X3899d68a62bd9b58f6d28de79e6f5d1d                  0.2810  0.2810     1 31.417  0.5549 0.461868   
#Nlevel                                             5.5308  5.5308     1 29.592 10.9200 0.002497 **
#Sterility                                          0.0003  0.0003     1 31.160  0.0007 0.979522   
#X3899d68a62bd9b58f6d28de79e6f5d1d:Nlevel           0.0881  0.0881     1 29.794  0.1740 0.679574   
#X3899d68a62bd9b58f6d28de79e6f5d1d:Sterility        0.1254  0.1254     1 31.827  0.2476 0.622206   
#Nlevel:Sterility                                   0.0139  0.0139     1 30.204  0.0274 0.869625   
#X3899d68a62bd9b58f6d28de79e6f5d1d:Nlevel:Sterility 0.7729  0.7729     1 30.582  1.5260 0.226119

#Nlevel effects


#OTU 741
#X453d47fbe68f91cd6b3aec11fc3b964b

modotu8<-lmer (log(CoreRootLength)~ X453d47fbe68f91cd6b3aec11fc3b964b *Nlevel*Sterility+  (1|G5Block), data=df)
summary(modotu8)
anova(modotu8)
#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
#X453d47fbe68f91cd6b3aec11fc3b964b                  0.0425  0.0425     1 31.308  0.0778 0.782147   
#Nlevel                                             5.8507  5.8507     1 29.855 10.7151 0.002688 **
#Sterility                                          0.7028  0.7028     1 30.145  1.2872 0.265507   
#X453d47fbe68f91cd6b3aec11fc3b964b:Nlevel           0.3509  0.3509     1 30.237  0.6426 0.429020   
#X453d47fbe68f91cd6b3aec11fc3b964b:Sterility        0.4646  0.4646     1 30.918  0.8510 0.363426   
#Nlevel:Sterility                                   0.6922  0.6922     1 31.770  1.2677 0.268628   
#X453d47fbe68f91cd6b3aec11fc3b964b:Nlevel:Sterility 0.0048  0.0048     1 31.861  0.0089 0.925519   


#Nlevel effects



#OTU 6838
#X69d3d7cc762732aaf09fd2959046648a

modotu9<-lmer (log(CoreRootLength)~ X69d3d7cc762732aaf09fd2959046648a*Nlevel*Sterility +  (1|G5Block), data=df)
summary(modotu9)
anova(modotu9)
#Type III Analysis of Variance Table with Satterthwaite's method
#                                                    Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
#X69d3d7cc762732aaf09fd2959046648a                  0.40766 0.40766     1 30.482  0.7281 0.40015  
#Nlevel                                             2.40201 2.40201     1 30.695  4.2901 0.04683 *
#Sterility                                          0.83484 0.83484     1 31.686  1.4911 0.23106  
#X69d3d7cc762732aaf09fd2959046648a:Nlevel           0.78447 0.78447     1 31.288  1.4011 0.24545  
#X69d3d7cc762732aaf09fd2959046648a:Sterility        0.79963 0.79963     1 31.941  1.4282 0.24086  
#Nlevel:Sterility                                   0.44790 0.44790     1 30.549  0.8000 0.37810  
#X69d3d7cc762732aaf09fd2959046648a:Nlevel:Sterility 0.41830 0.41830     1 31.972  0.7471 0.39383  

#Nlevel effects

#OTU 6371
#X7fb63308d1a53c4af530608e9ccee61d

modotu10<-lmer (log(CoreRootLength)~ X7fb63308d1a53c4af530608e9ccee61d*Nlevel*Sterility  + (1|G5Block), data=df)
summary(modotu10)
anova(modotu10)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                    Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
#X7fb63308d1a53c4af530608e9ccee61d                  0.07770 0.07770     1 30.694  0.1721 0.68113  
#Nlevel                                             1.93449 1.93449     1 30.733  4.2849 0.04695 *
#Sterility                                          1.91869 1.91869     1 30.809  4.2498 0.04778 *
#X7fb63308d1a53c4af530608e9ccee61d:Nlevel           0.34437 0.34437     1 31.726  0.7628 0.38902  
#X7fb63308d1a53c4af530608e9ccee61d:Sterility        0.91421 0.91421     1 31.589  2.0249 0.16454  
#Nlevel:Sterility                                   2.30552 2.30552     1 31.590  5.1067 0.03087 *
#X7fb63308d1a53c4af530608e9ccee61d:Nlevel:Sterility 0.15781 0.15781     1 31.941  0.3495 0.55854                  

#Nlevel effects and N level by sterility interaction NS after correcting for MC


#OTU 1192
#X8864a30ab9e33ba131278de4a262ee0a

modotu11<-lmer(log(CoreRootLength)~ X8864a30ab9e33ba131278de4a262ee0a*Nlevel*Sterility  +  (1|G5Block), data=df)
summary(modotu11)
anova(modotu11)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                    Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#X8864a30ab9e33ba131278de4a262ee0a                   0.7971  0.7971     1 31.861  1.4526   0.23699    
#Nlevel                                             16.1459 16.1459     1 29.201 29.4224 7.641e-06 ***
#Sterility                                           0.4680  0.4680     1 30.217  0.8528   0.36308    
#X8864a30ab9e33ba131278de4a262ee0a:Nlevel            0.5258  0.5258     1 30.419  0.9581   0.33538    
#X8864a30ab9e33ba131278de4a262ee0a:Sterility         0.0027  0.0027     1 29.722  0.0050   0.94436    
#Nlevel:Sterility                                    1.7655  1.7655     1 29.369  3.2172   0.08317 .  
#X8864a30ab9e33ba131278de4a262ee0a:Nlevel:Sterility  0.0590  0.0590     1 30.690  0.1075   0.74521                

#Nlevel effects



#OTU 6710
#d68451a2b18e9fd8e6738547074179fd
modotu12<-lmer(log(CoreRootLength)~ d68451a2b18e9fd8e6738547074179fd*Nlevel*Sterility  + (1|G5Block), data=df)
summary(modotu12)
anova(modotu12)


#Type III Analysis of Variance Table with Satterthwaite's method
#                                                  Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)   
#d68451a2b18e9fd8e6738547074179fd                  0.1512  0.1512     1 29.765  0.2858 0.59687  
#Nlevel                                            5.0829  5.0829     1 29.518  9.6094 0.00423 **
#Sterility                                         0.3778  0.3778     1 30.746  0.7143 0.40456   
#d68451a2b18e9fd8e6738547074179fd:Nlevel           0.0075  0.0075     1 29.312  0.0141 0.90616   
#d68451a2b18e9fd8e6738547074179fd:Sterility        0.1077  0.1077     1 31.139  0.2036 0.65492   
#Nlevel:Sterility                                  1.0392  1.0392     1 29.640  1.9647 0.17140   
#d68451a2b18e9fd8e6738547074179fd:Nlevel:Sterility 0.0241  0.0241     1 30.753  0.0455 0.83253    

#Nlevel effects

#OTU 2558
#d84b98e91b181d4c2879e8ed860ef94e 
modotu13<-lmer(log(CoreRootLength)~ d84b98e91b181d4c2879e8ed860ef94e *Nlevel*Sterility  + (1|G5Block), data=df)
summary(modotu13)
anova(modotu13)

#Nlevel effects and OTU effects but NS after correcting for MC




#Type III Analysis of Variance Table with Satterthwaite's method
#                                                  Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#d84b98e91b181d4c2879e8ed860ef94e                  3.0162  3.0162     1 30.582  7.1002   0.01219 *  
#Nlevel                                            8.9781  8.9781     1 31.085 21.1342 6.754e-05 ***
#Sterility                                         0.2562  0.2562     1 30.451  0.6032   0.44337    
#d84b98e91b181d4c2879e8ed860ef94e:Nlevel           0.2864  0.2864     1 31.576  0.6742   0.41776    
#d84b98e91b181d4c2879e8ed860ef94e:Sterility        1.2879  1.2879     1 30.819  3.0316   0.09163 .  
#Nlevel:Sterility                                  0.4960  0.4960     1 31.088  1.1675   0.28823    
#d84b98e91b181d4c2879e8ed860ef94e:Nlevel:Sterility 0.0117  0.0117     1 31.998  0.0275   0.86931               



#OTU 6708
#ec4beab3963241c35c22f57155a96813
modotu14<-lmer(log(CoreRootLength)~ ec4beab3963241c35c22f57155a96813 *Nlevel * Sterility  + (1|G5Block), data=df)
summary(modotu14)
anova(modotu14)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#ec4beab3963241c35c22f57155a96813                   0.3056  0.3056     1 31.998  0.5562    0.4613    
#Nlevel                                            11.2016 11.2016     1 30.631 20.3840 8.777e-05 ***
#Sterility                                          0.4613  0.4613     1 30.841  0.8394    0.3667    
#ec4beab3963241c35c22f57155a96813:Nlevel            0.3302  0.3302     1 31.830  0.6009    0.4439    
#ec4beab3963241c35c22f57155a96813:Sterility         0.1998  0.1998     1 31.996  0.3635    0.5508    
#Nlevel:Sterility                                   1.3988  1.3988     1 30.364  2.5455    0.1210    
#ec4beab3963241c35c22f57155a96813:Nlevel:Sterility  0.0832  0.0832     1 31.611  0.1514    0.6998  

#Nlevel effects




#######Checked for mold effects in all these models and none



########################################
## adjustment for mc -double checked all p-values =correct
########################################

#Root length
RLOTUP<-c(0.9361371,0.6188258, 0.03370, 0.0038780,0.2484, 0.461868, 0.782147, 0.40015, 0.68113, 0.23699 , 0.59687, 0.01219, 0.4613)

RLOTUN<-c( 0.8131362,0.3495872, 0.03121, 0.0465390, 0.1481, 0.679574, 0.429020, 0.24545, 0.38902, 0.33538, 0.90616, 0.41776, 0.4439)

RLOTUS<-c( 0.9990972, 0.3940850, 0.82807, 0.9142609, 0.6658, 0.622206, 0.363426, 0.24086, 0.16454,  0.94436, 0.65492, 0.09163, 0.5508)

RL3Way<-c(0.6529449, 0.1962857, 0.30628, 0.5858810 , 0.3718 , 0.226119, 0.925519, 0.39383 , 0.55854, 0.74521, 0.83253, 0.86931 ,0.6998) 

#Fourth OTU is sig
BHRLOTUP = p.adjust(RLOTUP, method = "BH")
print(BHRLOTUP)
#0.9361371 0.8044735 0.1460333 0.0504140 0.6458400 0.7505355 0.8473259 0.7505355 0.8049718 0.6458400 0.8044735 0.0792350 0.7505355

BonRLOTUP = p.adjust(RLOTUP, method = "bonferroni")
print(BonRLOTUP)
#also fourth is still sig by bonferroni

#All NS
BHRLOTUN = p.adjust(RLOTUN, method = "BH")
print(BHRLOTUN)
#0.8808975 0.5770700 0.3025035 0.3025035 0.5770700 0.8031329 0.5770700 0.5770700 0.5770700 0.5770700 0.9061600 0.5770700 0.5770700


#All NS
BHRLOTUS = p.adjust(RLOTUS, method = "BH")
print(BHRLOTUS)
# 0.9990972 0.9617111 0.9990972 0.9990972 0.9617111 0.9617111 0.9617111 0.9617111 0.9617111 0.9990972 0.9617111 0.9617111 0.9617111





#all 3 way interactions NS
BHRL3Way = p.adjust(RL3Way, method = "BH")
print(BHRL3Way)

#0.925519 0.925519 0.925519 0.925519 0.925519 0.925519 0.925519 0.925519 0.925519 0.925519 0.925519 0.925519 0.925519

########################################
## CoreMeanRootDiameter 

########################################

########################################

#OTU 1
#  X726a6245b76a9e74caea7993598fd969

modotu1<-lmer(CoreMeanRootDiameter~  X726a6245b76a9e74caea7993598fd969*Nlevel*Sterility  + (1|G5Block), data=df)
summary(modotu1)
anova(modotu1)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                       Sum Sq    Mean Sq NumDF DenDF F value  Pr(>F)  
#X726a6245b76a9e74caea7993598fd969                  4.1685e-06 4.1685e-06     1    32  0.4642 0.50058  
#Nlevel                                             1.8610e-07 1.8610e-07     1    32  0.0207 0.88644  
#Sterility                                          3.1432e-05 3.1432e-05     1    32  3.5000 0.07054 .
#X726a6245b76a9e74caea7993598fd969:Nlevel           1.9975e-05 1.9975e-05     1    32  2.2243 0.14565  
#X726a6245b76a9e74caea7993598fd969:Sterility        9.4549e-06 9.4549e-06     1    32  1.0528 0.31255  
#Nlevel:Sterility                                   4.1400e-08 4.1400e-08     1    32  0.0046 0.94627  
#X726a6245b76a9e74caea7993598fd969:Nlevel:Sterility 1.7738e-05 1.7738e-05     1    32  1.9752 0.16953               

#All NS




#OTU 2
#  X7d3587bbb834c6ac86656e1d2715c331

modotu2<-lmer (CoreMeanRootDiameter~ X7d3587bbb834c6ac86656e1d2715c331*Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu2)
anova(modotu2)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                       Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#X7d3587bbb834c6ac86656e1d2715c331                  2.7850e-07 2.7850e-07     1    32  0.0257 0.8738
#Nlevel                                             2.8358e-06 2.8358e-06     1    32  0.2612 0.6128
#Sterility                                          9.1023e-06 9.1023e-06     1    32  0.8384 0.3667
#X7d3587bbb834c6ac86656e1d2715c331:Nlevel           1.3810e-07 1.3810e-07     1    32  0.0127 0.9109
#X7d3587bbb834c6ac86656e1d2715c331:Sterility        8.0900e-08 8.0900e-08     1    32  0.0075 0.9317
#Nlevel:Sterility                                   1.4565e-06 1.4565e-06     1    32  0.1342 0.7166
#X7d3587bbb834c6ac86656e1d2715c331:Nlevel:Sterility 1.3530e-07 1.3530e-07     1    32  0.0125 0.9118

#All NS

#OTU 3
#a90d5b4368338dea15569269e8d33fb8 

modotu3<-lmer(CoreMeanRootDiameter~a90d5b4368338dea15569269e8d33fb8 *Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu3)
anova(modotu3)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                      Sum Sq    Mean Sq NumDF DenDF F value  Pr(>F)  
#a90d5b4368338dea15569269e8d33fb8                  7.2500e-07 7.2500e-07     1    32  0.0701 0.79283  
#Nlevel                                            9.4010e-06 9.4010e-06     1    32  0.9093 0.34744  
#Sterility                                         3.5186e-05 3.5186e-05     1    32  3.4036 0.07433 .
#a90d5b4368338dea15569269e8d33fb8:Nlevel           4.3970e-06 4.3970e-06     1    32  0.4253 0.51897  
#a90d5b4368338dea15569269e8d33fb8:Sterility        7.2750e-06 7.2750e-06     1    32  0.7037 0.40777  
#Nlevel:Sterility                                  1.2805e-05 1.2805e-05     1    32  1.2386 0.27403  
#a90d5b4368338dea15569269e8d33fb8:Nlevel:Sterility 1.3998e-05 1.3998e-05     1    32  1.3540 0.25318                                    

#All NS


#OTU 4
#  bea3f263c468bf3d3fb6645562737c70 
#################################
#Not in the top abundances after relative fitness and rarefaction
#########################################
#modotu4<-lm (CoreMeanRootDiameter~ bea3f263c468bf3d3fb6645562737c70 *Nlevel*Sterility, data=df)
#summary(modotu4)
#anova(modotu4)

#Response: CoreMeanRootDiameter

#bea3f263c468bf3d3fb6645562737c70                   1 0.00000014 1.4300e-07  0.0131 0.90958  

#bea3f263c468bf3d3fb6645562737c70:Nlevel            1 0.00000249 2.4880e-06  0.2280 0.63628  

#bea3f263c468bf3d3fb6645562737c70:Sterility         1 0.00000070 7.0100e-07  0.0643 0.80152  

#bea3f263c468bf3d3fb6645562737c70:Nlevel:Sterility  1 0.00000119 1.1910e-06  0.1091 0.74329  


#OTU 5
# X071e2a4c88faf23787b614210e527ce3

modotu5<-lmer (CoreMeanRootDiameter~ X071e2a4c88faf23787b614210e527ce3*Nlevel*Sterility  + (1|G5Block), data=df)
summary(modotu5)
anova(modotu5)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                       Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#X071e2a4c88faf23787b614210e527ce3                  1.1654e-06 1.1654e-06     1    32  0.1108 0.7414
#Nlevel                                             2.8641e-05 2.8641e-05     1    32  2.7223 0.1087
#Sterility                                          1.2888e-05 1.2888e-05     1    32  1.2250 0.2766
#X071e2a4c88faf23787b614210e527ce3:Nlevel           1.0286e-05 1.0286e-05     1    32  0.9777 0.3302
#X071e2a4c88faf23787b614210e527ce3:Sterility        6.6000e-09 6.6000e-09     1    32  0.0006 0.9801
#Nlevel:Sterility                                   9.3000e-09 9.3000e-09     1    32  0.0009 0.9765
#X071e2a4c88faf23787b614210e527ce3:Nlevel:Sterility 1.7980e-07 1.7980e-07     1    32  0.0171 0.8968

#All NS





#OTU 6
#X2bb5ef7086fb6a049aa217591531d386 

modotu6<-lmer (CoreMeanRootDiameter~ X2bb5ef7086fb6a049aa217591531d386 *Nlevel*Sterility +(1|G5Block), data=df)
summary(modotu6)
anova(modotu6)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                       Sum Sq    Mean Sq NumDF  DenDF F value Pr(>F)
#X2bb5ef7086fb6a049aa217591531d386                  1.1943e-05 1.1943e-05     1 31.986  1.1450 0.2926
#Nlevel                                             1.4002e-05 1.4002e-05     1 31.986  1.3425 0.2552
#Sterility                                          2.3173e-05 2.3173e-05     1 31.992  2.2217 0.1459
#X2bb5ef7086fb6a049aa217591531d386:Nlevel           4.2000e-09 4.2000e-09     1 31.980  0.0004 0.9841
#X2bb5ef7086fb6a049aa217591531d386:Sterility        3.2680e-07 3.2680e-07     1 31.998  0.0313 0.8606
#Nlevel:Sterility                                   1.4583e-06 1.4583e-06     1 31.999  0.1398 0.7109
#X2bb5ef7086fb6a049aa217591531d386:Nlevel:Sterility 5.3380e-06 5.3380e-06     1 31.999  0.5118 0.4796
#All NS



#OTU 7
#OTU 718 with revised OTU names 
#X3899d68a62bd9b58f6d28de79e6f5d1d

modotu7<-lmer (CoreMeanRootDiameter~ X3899d68a62bd9b58f6d28de79e6f5d1d*Nlevel*Sterility  + (1|G5Block), data=df)
summary(modotu7)
anova(modotu7)
#Type III Analysis of Variance Table with Satterthwaite's method
#                                                       Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#X3899d68a62bd9b58f6d28de79e6f5d1d                  3.1515e-06 3.1515e-06     1    32  0.3085 0.5825
#Nlevel                                             7.5830e-07 7.5830e-07     1    32  0.0742 0.7870
#Sterility                                          5.5111e-06 5.5111e-06     1    32  0.5394 0.4680
#X3899d68a62bd9b58f6d28de79e6f5d1d:Nlevel           8.2408e-06 8.2408e-06     1    32  0.8066 0.3758
#X3899d68a62bd9b58f6d28de79e6f5d1d:Sterility        1.3300e-08 1.3300e-08     1    32  0.0013 0.9715
#Nlevel:Sterility                                   1.5300e-08 1.5300e-08     1    32  0.0015 0.9694
#X3899d68a62bd9b58f6d28de79e6f5d1d:Nlevel:Sterility 1.3170e-07 1.3170e-07     1    32  0.0129 0.9103


#All NS   


#OTU 8
#X453d47fbe68f91cd6b3aec11fc3b964b

modotu8<-lmer (CoreMeanRootDiameter~ X453d47fbe68f91cd6b3aec11fc3b964b *Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu8)
anova(modotu8)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                       Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#X453d47fbe68f91cd6b3aec11fc3b964b                  3.4208e-06 3.4208e-06     1    32  0.3329 0.5680
#Nlevel                                             9.0840e-07 9.0840e-07     1    32  0.0884 0.7682
#Sterility                                          5.5080e-07 5.5080e-07     1    32  0.0536 0.8184
#X453d47fbe68f91cd6b3aec11fc3b964b:Nlevel           3.0320e-07 3.0320e-07     1    32  0.0295 0.8647
#X453d47fbe68f91cd6b3aec11fc3b964b:Sterility        3.6281e-06 3.6281e-06     1    32  0.3530 0.5566
#Nlevel:Sterility                                   5.9500e-07 5.9500e-07     1    32  0.0579 0.8114
#X453d47fbe68f91cd6b3aec11fc3b964b:Nlevel:Sterility 1.1438e-06 1.1438e-06     1    32  0.1113 0.7409
   
#All NS




#OTU 9
#X69d3d7cc762732aaf09fd2959046648a

modotu9<-lmer (CoreMeanRootDiameter~ X69d3d7cc762732aaf09fd2959046648a*Nlevel*Sterility  + (1|G5Block), data=df)
summary(modotu9)
anova(modotu9)
 
#Type III Analysis of Variance Table with Satterthwaite's method
#                                                       Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#X69d3d7cc762732aaf09fd2959046648a                  1.0568e-06 1.0568e-06     1    32  0.1025 0.7510
#Nlevel                                             7.0950e-06 7.0950e-06     1    32  0.6879 0.4130
#Sterility                                          2.9281e-05 2.9281e-05     1    32  2.8388 0.1017
#X69d3d7cc762732aaf09fd2959046648a:Nlevel           1.5230e-07 1.5230e-07     1    32  0.0148 0.9040
#X69d3d7cc762732aaf09fd2959046648a:Sterility        7.5131e-06 7.5131e-06     1    32  0.7284 0.3997
#Nlevel:Sterility                                   6.5030e-07 6.5030e-07     1    32  0.0630 0.8034
#X69d3d7cc762732aaf09fd2959046648a:Nlevel:Sterility 5.6361e-06 5.6361e-06     1    32  0.5464 0.4652

#All NS

#OTU 10
#X7fb63308d1a53c4af530608e9ccee61d

modotu10<-lmer (CoreMeanRootDiameter~ X7fb63308d1a53c4af530608e9ccee61d*Nlevel*Sterility  + (1|G5Block), data=df)
summary(modotu10)
anova(modotu10)

#Type III Analysis of Variance Table with Satterthwaite's method
 #                                                      Sum Sq    Mean Sq NumDF DenDF F value  Pr(>F)  
#X7fb63308d1a53c4af530608e9ccee61d                  2.0260e-06 2.0260e-06     1    32  0.2006 0.65725  
#Nlevel                                             4.6485e-05 4.6485e-05     1    32  4.6025 0.03962 *
#Sterility                                          1.7141e-05 1.7141e-05     1    32  1.6971 0.20196  
#X7fb63308d1a53c4af530608e9ccee61d:Nlevel           3.0591e-05 3.0591e-05     1    32  3.0288 0.09141 .
#X7fb63308d1a53c4af530608e9ccee61d:Sterility        1.5090e-06 1.5090e-06     1    32  0.1494 0.70171  
#Nlevel:Sterility                                   7.9400e-07 7.9400e-07     1    32  0.0786 0.78098  
#X7fb63308d1a53c4af530608e9ccee61d:Nlevel:Sterility 8.8000e-06 8.8000e-06     1    32  0.8713 0.35759  

#All NS




#OTU 11
#X8864a30ab9e33ba131278de4a262ee0a

modotu11<-lmer (CoreMeanRootDiameter~ X8864a30ab9e33ba131278de4a262ee0a*Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu11)
anova(modotu11)
#Type III Analysis of Variance Table with Satterthwaite's method
 #                                                      Sum Sq    Mean Sq NumDF DenDF F value  Pr(>F)  
#X8864a30ab9e33ba131278de4a262ee0a                  1.7663e-05 1.7663e-05     1    32  1.8903 0.17872  
#Nlevel                                             7.5790e-06 7.5790e-06     1    32  0.8111 0.37452  
#Sterility                                          6.3024e-05 6.3024e-05     1    32  6.7449 0.01409 *
#X8864a30ab9e33ba131278de4a262ee0a:Nlevel           4.2710e-06 4.2710e-06     1    32  0.4571 0.50383  
#X8864a30ab9e33ba131278de4a262ee0a:Sterility        3.6400e-06 3.6400e-06     1    32  0.3896 0.53693  
#Nlevel:Sterility                                   7.6880e-06 7.6880e-06     1    32  0.8228 0.37114  
#X8864a30ab9e33ba131278de4a262ee0a:Nlevel:Sterility 8.9870e-06 8.9870e-06     1    32  0.9618 0.33409                                        

#All NS




#OTU 12
#d68451a2b18e9fd8e6738547074179fd
modotu12<-lmer (CoreMeanRootDiameter~ d68451a2b18e9fd8e6738547074179fd*Nlevel*Sterility  + (1|G5Block), data=df)
summary(modotu12)
anova(modotu12)


#Type III Analysis of Variance Table with Satterthwaite's method
 #                                                     Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#d68451a2b18e9fd8e6738547074179fd                  3.1405e-06 3.1405e-06     1    32  0.3066 0.5836
#Nlevel                                            1.6874e-06 1.6874e-06     1    32  0.1647 0.6875
#Sterility                                         2.2228e-05 2.2228e-05     1    32  2.1701 0.1505
#d68451a2b18e9fd8e6738547074179fd:Nlevel           8.2700e-08 8.2700e-08     1    32  0.0081 0.9290
#d68451a2b18e9fd8e6738547074179fd:Sterility        2.3657e-06 2.3657e-06     1    32  0.2310 0.6341
#Nlevel:Sterility                                  4.2000e-07 4.2000e-07     1    32  0.0410 0.8408
#d68451a2b18e9fd8e6738547074179fd:Nlevel:Sterility 3.3904e-06 3.3904e-06     1    32  0.3310 0.5691

#All NS


#OTU 13
#d84b98e91b181d4c2879e8ed860ef94e 
modotu13<-lmer (CoreMeanRootDiameter~ d84b98e91b181d4c2879e8ed860ef94e *Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu13)
anova(modotu13)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                      Sum Sq    Mean Sq NumDF  DenDF F value Pr(>F)
#d84b98e91b181d4c2879e8ed860ef94e                  6.4216e-06 6.4216e-06     1 23.108  0.6078 0.4435
#Nlevel                                            3.2633e-06 3.2633e-06     1 31.995  0.3089 0.5822
#Sterility                                         1.1195e-05 1.1195e-05     1 30.758  1.0596 0.3113
#d84b98e91b181d4c2879e8ed860ef94e:Nlevel           5.1940e-07 5.1940e-07     1 29.230  0.0492 0.8261
#d84b98e91b181d4c2879e8ed860ef94e:Sterility        8.9900e-08 8.9900e-08     1 31.543  0.0085 0.9271
#Nlevel:Sterility                                  1.5830e-06 1.5830e-06     1 31.969  0.1498 0.7013
#d84b98e91b181d4c2879e8ed860ef94e:Nlevel:Sterility 6.2900e-06 6.2900e-06     1 27.098  0.5954 0.4470
#All NS



#OTU 14
#ec4beab3963241c35c22f57155a96813
modotu14<-lmer (CoreMeanRootDiameter~ ec4beab3963241c35c22f57155a96813 *Nlevel * Sterility + (1|G5Block), data=df)
summary(modotu14)
anova(modotu14)


#Type III Analysis of Variance Table with Satterthwaite's method
#                                                      Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#ec4beab3963241c35c22f57155a96813                  2.0019e-06 2.0019e-06     1    32  0.2017 0.6563
#Nlevel                                            4.4500e-08 4.4500e-08     1    32  0.0045 0.9470
#Sterility                                         1.3335e-05 1.3335e-05     1    32  1.3439 0.2549
#ec4beab3963241c35c22f57155a96813:Nlevel           2.6156e-06 2.6156e-06     1    32  0.2636 0.6112
#ec4beab3963241c35c22f57155a96813:Sterility        1.5050e-07 1.5050e-07     1    32  0.0152 0.9027
#Nlevel:Sterility                                  1.3550e-07 1.3550e-07     1    32  0.0137 0.9077
#ec4beab3963241c35c22f57155a96813:Nlevel:Sterility 2.0406e-06 2.0406e-06     1    32  0.2056 0.6533

#All NS


#######Checked for mold effects in all these models and none

########################################
##
#All NS so need to check for MC

########################################
RDOTU<-c( )   
RDOTUN<-c( ) 
RDOTUS<-c()
RD3Way<-c()
#OTU effects are NS after correction
BHRDOTU = p.adjust(RDOTU, method = "BH")
print(BHRDOTU)

#OUT x Nitrogen effects
BHRDOTUN = p.adjust(RDOTUN, method = "BH")
print(BHRDOTUN)

#OTUx Sterility effects
BHRDOTUS = p.adjust(RDOTUS, method = "BH")
print(BHRDOTUS)

BHRD3Way = p.adjust(RD3Way, method = "BH")
print(BHRD3Way)



########################################
##CoreSpecificRootLength

########################################
#Column called P-check is where I am double checking I used the right p-values in MC correction-ignore for interpreting model


#OTU 1
#  X726a6245b76a9e74caea7993598fd969

modotu1<-lmer(CoreSpecificRootLength~  X726a6245b76a9e74caea7993598fd969*Nlevel*Sterility +(1|G5Block) , data=df)
summary(modotu1)
anova(modotu1)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  P-check
#X726a6245b76a9e74caea7993598fd969                  407.49  407.49     1    32  1.7319 0.19752  0.19752
#Nlevel                                               0.01    0.01     1    32  0.0000 0.99445  
#Sterility                                          386.43  386.43     1    32  1.6424 0.20921  
#X726a6245b76a9e74caea7993598fd969:Nlevel           895.26  895.26     1    32  3.8049 0.05991 .0.05991
#X726a6245b76a9e74caea7993598fd969:Sterility        521.97  521.97     1    32  2.2184 0.14616  0.14616
#Nlevel:Sterility                                    93.13   93.13     1    32  0.3958 0.53373  
#X726a6245b76a9e74caea7993598fd969:Nlevel:Sterility 643.14  643.14     1    32  2.7334 0.10805 0.10805







#OTU 2
#  X7d3587bbb834c6ac86656e1d2715c331

modotu2<-lmer (CoreSpecificRootLength~ X7d3587bbb834c6ac86656e1d2715c331*Nlevel*Sterility  + (1|G5Block), data=df)
summary(modotu2)
anova(modotu2)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF DenDF F value Pr(>F) P-check
#X7d3587bbb834c6ac86656e1d2715c331                  46.058  46.058     1    32  0.1694 0.6834 0.6834
#Nlevel                                              2.447   2.447     1    32  0.0090 0.9250
#Sterility                                           5.660   5.660     1    32  0.0208 0.8862
#X7d3587bbb834c6ac86656e1d2715c331:Nlevel           10.985  10.985     1    32  0.0404 0.8420   0.8420
#X7d3587bbb834c6ac86656e1d2715c331:Sterility        30.867  30.867     1    32  0.1135 0.7383  0.7383
#Nlevel:Sterility                                   58.436  58.436     1    32  0.2150 0.6460
#X7d3587bbb834c6ac86656e1d2715c331:Nlevel:Sterility 23.111  23.111     1    32  0.0850 0.7725   0.7725



#OTU 3
#a90d5b4368338dea15569269e8d33fb8 

modotu3<-lmer(CoreSpecificRootLength~a90d5b4368338dea15569269e8d33fb8 *Nlevel*Sterility  + (1|G5Block), data=df)
summary(modotu3)
anova(modotu3)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                  Sum Sq Mean Sq NumDF DenDF F value Pr(>F) p-check
#a90d5b4368338dea15569269e8d33fb8                  107.19  107.19     1    32  0.4026 0.5303 0.5303
#Nlevel                                            374.17  374.17     1    32  1.4052 0.2446
#Sterility                                         589.82  589.82     1    32  2.2151 0.1465
#a90d5b4368338dea15569269e8d33fb8:Nlevel           143.99  143.99     1    32  0.5408 0.4675  0.4675
#a90d5b4368338dea15569269e8d33fb8:Sterility        316.71  316.71     1    32  1.1894 0.2836  0.2836
#Nlevel:Sterility                                   12.88   12.88     1    32  0.0484 0.8273
#a90d5b4368338dea15569269e8d33fb8:Nlevel:Sterility 176.86  176.86     1    32  0.6642 0.4211  0.4211


##########################################
# OTU not found in most abundant list
##################################

#OTU 4
#  bea3f263c468bf3d3fb6645562737c70 

#modotu4<-lm (CoreSpecificRootLength~ bea3f263c468bf3d3fb6645562737c70 *Nlevel*Sterility, data=df)
#summary(modotu4)
#anova(modotu4)

#Response: CoreSpecificRootLength

#bea3f263c468bf3d3fb6645562737c70                   1   16.3   16.28  0.0591 0.8094

#bea3f263c468bf3d3fb6645562737c70:Nlevel            1   90.2   90.16  0.3274 0.5712
#bea3f263c468bf3d3fb6645562737c70:Sterility         1    1.3    1.26  0.0046 0.9465

#bea3f263c468bf3d3fb6645562737c70:Nlevel:Sterility  1  180.9  180.90  0.6570 0.4236


#OTU 5
# X071e2a4c88faf23787b614210e527ce3

modotu5<-lmer (CoreSpecificRootLength~ X071e2a4c88faf23787b614210e527ce3*Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu5)
anova(modotu5)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF DenDF F value Pr(>F) p-check
#X071e2a4c88faf23787b614210e527ce3                   94.18   94.18     1    32  0.3553 0.5553  0.5553
#Nlevel                                             599.64  599.64     1    32  2.2622 0.1424
#Sterility                                          317.47  317.47     1    32  1.1977 0.2820
#X071e2a4c88faf23787b614210e527ce3:Nlevel            61.77   61.77     1    32  0.2330 0.6326  0.6326
#X071e2a4c88faf23787b614210e527ce3:Sterility         38.11   38.11     1    32  0.1438 0.7071 0.7071
#Nlevel:Sterility                                    96.60   96.60     1    32  0.3644 0.5503
#X071e2a4c88faf23787b614210e527ce3:Nlevel:Sterility   5.60    5.60     1    32  0.0211 0.8854 0.8854



#OTU 6
#X2bb5ef7086fb6a049aa217591531d386 

modotu6<-lmer (CoreSpecificRootLength~ X2bb5ef7086fb6a049aa217591531d386 *Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu6)
anova(modotu6)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF DenDF F value Pr(>F) p-check
#X2bb5ef7086fb6a049aa217591531d386                  602.82  602.82     1    32  2.3597 0.1343 0.1343
#Nlevel                                             457.00  457.00     1    32  1.7889 0.1905
#Sterility                                          171.42  171.42     1    32  0.6710 0.4188
#X2bb5ef7086fb6a049aa217591531d386:Nlevel             1.01    1.01     1    32  0.0040 0.9502  0.9502
#X2bb5ef7086fb6a049aa217591531d386:Sterility         72.36   72.36     1    32  0.2833 0.5983 0.5983
#Nlevel:Sterility                                   356.60  356.60     1    32  1.3959 0.2461
#X2bb5ef7086fb6a049aa217591531d386:Nlevel:Sterility 282.07  282.07     1    32  1.1042 0.3012  0.3012

#OTU 7
#X3899d68a62bd9b58f6d28de79e6f5d1d

modotu7<-lmer (CoreSpecificRootLength~ X3899d68a62bd9b58f6d28de79e6f5d1d*Nlevel*Sterility  + (1|G5Block), data=df)
summary(modotu7)
anova(modotu7)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                    Sum Sq Mean Sq NumDF DenDF F value Pr(>F) p-check
#X3899d68a62bd9b58f6d28de79e6f5d1d                  172.872 172.872     1    32  0.6499 0.4261 0.4261
#Nlevel                                               6.671   6.671     1    32  0.0251 0.8752
#Sterility                                          183.257 183.257     1    32  0.6889 0.4127
#X3899d68a62bd9b58f6d28de79e6f5d1d:Nlevel           141.299 141.299     1    32  0.5312 0.4714 0.4714
#X3899d68a62bd9b58f6d28de79e6f5d1d:Sterility         67.495  67.495     1    32  0.2537 0.6179 0.6179
#Nlevel:Sterility                                   106.008 106.008     1    32  0.3985 0.5323
#X3899d68a62bd9b58f6d28de79e6f5d1d:Nlevel:Sterility  31.639  31.639     1    32  0.1189 0.7324 0.7324

#OTU 8
#X453d47fbe68f91cd6b3aec11fc3b964b

modotu8<-lmer (CoreSpecificRootLength~ X453d47fbe68f91cd6b3aec11fc3b964b *Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu8)
anova(modotu8)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF DenDF F value Pr(>F) p-check
#X453d47fbe68f91cd6b3aec11fc3b964b                  104.52  104.52     1    32  0.4232 0.5200 0.5200
#Nlevel                                              94.98   94.98     1    32  0.3846 0.5396
#Sterility                                           40.72   40.72     1    32  0.1648 0.6874
#X453d47fbe68f91cd6b3aec11fc3b964b:Nlevel            16.21   16.21     1    32  0.0656 0.7995 0.7995
#X453d47fbe68f91cd6b3aec11fc3b964b:Sterility        174.21  174.21     1    32  0.7054 0.4072 0.4072
#Nlevel:Sterility                                   369.01  369.01     1    32  1.4941 0.2305
#X453d47fbe68f91cd6b3aec11fc3b964b:Nlevel:Sterility 227.22  227.22     1    32  0.9200 0.3447 0.3447
 
#OTU 9
#X69d3d7cc762732aaf09fd2959046648a

modotu9<-lmer (CoreSpecificRootLength~ X69d3d7cc762732aaf09fd2959046648a*Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu9)
anova(modotu9)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                    Sum Sq Mean Sq NumDF DenDF F value   Pr(>F) p-check
#X69d3d7cc762732aaf09fd2959046648a                    15.90   15.90     1    32  0.0733 0.788332  0.788332  
#Nlevel                                              838.30  838.30     1    32  3.8654 0.058018 . 
#Sterility                                          1643.53 1643.53     1    32  7.5783 0.009653 **
#X69d3d7cc762732aaf09fd2959046648a:Nlevel            303.35  303.35     1    32  1.3987 0.245647  0.245647  
#X69d3d7cc762732aaf09fd2959046648a:Sterility        1161.92 1161.92     1    32  5.3576 0.027204 * 0.027204
#Nlevel:Sterility                                     29.44   29.44     1    32  0.1357 0.714990   
#X69d3d7cc762732aaf09fd2959046648a:Nlevel:Sterility  331.80  331.80     1    32  1.5299 0.225124  0.225124




#OTU 10
#X7fb63308d1a53c4af530608e9ccee61d

modotu10<-lmer (CoreSpecificRootLength~ X7fb63308d1a53c4af530608e9ccee61d*Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu10)
anova(modotu10)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                    Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  p-check
#X7fb63308d1a53c4af530608e9ccee61d                   107.31  107.31     1 31.803  0.4615 0.50184 0.50184
#Nlevel                                             1610.34 1610.34     1 31.757  6.9255 0.01300 *
#Sterility                                           320.24  320.24     1 31.998  1.3772 0.24923  
#X7fb63308d1a53c4af530608e9ccee61d:Nlevel           1475.64 1475.64     1 30.844  6.3462 0.01716 *0.01716
#X7fb63308d1a53c4af530608e9ccee61d:Sterility         190.89  190.89     1 30.512  0.8210 0.37200  0.37200
#Nlevel:Sterility                                    660.69  660.69     1 31.395  2.8414 0.10178  
#X7fb63308d1a53c4af530608e9ccee61d:Nlevel:Sterility 1079.17 1079.17     1 28.571  4.6412 0.03979 *0.03979

#OTU 11
#X8864a30ab9e33ba131278de4a262ee0a

modotu11<-lmer (CoreSpecificRootLength~ X8864a30ab9e33ba131278de4a262ee0a*Nlevel*Sterility +(1|G5Block), data=df)
summary(modotu11)
anova(modotu11)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  p-check
#X8864a30ab9e33ba131278de4a262ee0a                  300.30  300.30     1    32  1.1822 0.28502 0.28502 
#Nlevel                                             106.77  106.77     1    32  0.4203 0.52139  
#Sterility                                          775.68  775.68     1    32  3.0538 0.09014 .
#X8864a30ab9e33ba131278de4a262ee0a:Nlevel           106.52  106.52     1    32  0.4194 0.52187  0.52187
#X8864a30ab9e33ba131278de4a262ee0a:Sterility        156.94  156.94     1    32  0.6179 0.43762  0.43762
#Nlevel:Sterility                                     8.43    8.43     1    32  0.0332 0.85659  
#X8864a30ab9e33ba131278de4a262ee0a:Nlevel:Sterility  31.49   31.49     1    32  0.1240 0.72709  0.72709




#OTU 12
#d68451a2b18e9fd8e6738547074179fd
modotu12<-lmer (CoreSpecificRootLength~ d68451a2b18e9fd8e6738547074179fd*Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu12)
anova(modotu12)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                  Sum Sq Mean Sq NumDF DenDF F value Pr(>F) p-check
#d68451a2b18e9fd8e6738547074179fd                  112.81  112.81     1    32  0.4822 0.4924  0.4924
#Nlevel                                             40.12   40.12     1    32  0.1715 0.6816
#Sterility                                         336.00  336.00     1    32  1.4362 0.2396
#d68451a2b18e9fd8e6738547074179fd:Nlevel             0.79    0.79     1    32  0.0034 0.9540 0.9540
#d68451a2b18e9fd8e6738547074179fd:Sterility        259.68  259.68     1    32  1.1099 0.3000 0.3000
#Nlevel:Sterility                                  491.13  491.13     1    32  2.0992 0.1571
#d68451a2b18e9fd8e6738547074179fd:Nlevel:Sterility 470.00  470.00     1    32  2.0089 0.1660 0.1660



#OTU 13
#d84b98e91b181d4c2879e8ed860ef94e 
modotu13<-lmer (CoreSpecificRootLength~ d84b98e91b181d4c2879e8ed860ef94e *Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu13)
anova(modotu13)

#Analysis of Variance Table

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                  Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  p-check
#d84b98e91b181d4c2879e8ed860ef94e                  956.35  956.35     1 25.433  4.0225 0.05565 .0.05565
#Nlevel                                            526.67  526.67     1 30.875  2.2153 0.14680  
#Sterility                                         352.90  352.90     1 29.138  1.4844 0.23287  
#d84b98e91b181d4c2879e8ed860ef94e:Nlevel           193.33  193.33     1 31.912  0.8132 0.37394  0.37394
#d84b98e91b181d4c2879e8ed860ef94e:Sterility        141.51  141.51     1 30.065  0.5952 0.44643  0.44643
#Nlevel:Sterility                                  642.07  642.07     1 30.920  2.7006 0.11044  
#d84b98e91b181d4c2879e8ed860ef94e:Nlevel:Sterility 774.20  774.20     1 31.277  3.2564 0.08078 .0.08078



#OTU 14
#ec4beab3963241c35c22f57155a96813
modotu14<-lmer (CoreSpecificRootLength~ ec4beab3963241c35c22f57155a96813 *Nlevel * Sterility  + (1|G5Block), data=df)
summary(modotu14)
anova(modotu14)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF DenDF F value Pr(>F) p-check
#ec4beab3963241c35c22f57155a96813                  102.283 102.283     1    32  0.3917 0.5358 0.5358
#Nlevel                                             23.124  23.124     1    32  0.0886 0.7679
#Sterility                                         115.063 115.063     1    32  0.4407 0.5116
#ec4beab3963241c35c22f57155a96813:Nlevel            10.390  10.390     1    32  0.0398 0.8432 0.8432  
#ec4beab3963241c35c22f57155a96813:Sterility          0.062   0.062     1    32  0.0002 0.9878  0.9878   
#Nlevel:Sterility                                  139.887 139.887     1    32  0.5357 0.4695
#ec4beab3963241c35c22f57155a96813:Nlevel:Sterility  58.251  58.251     1    32  0.2231 0.6399 0.6399 





#######Checked for mold effects in all these models and none

########################################
##Double checked all p-values are correct 

#Specific Root Length
SRLOTU<-c(0.19752, 0.6834, 0.5303, 0.5553, 0.1343, 0.4261, 0.5200, 0.788332, 0.50184, 0.28502, 0.4924, 0.05565, 0.5358 )   
SRLOTUN<-c( 0.05991, 0.8420, 0.4675, 0.6326, 0.9502, 0.4714, 0.7995, 0.245647, 0.01716, 0.52187, 0.9540, 0.37394, 0.8432   ) 
SRLOTUS<-c(0.14616, 0.7383, 0.2836, 0.7071, 0.5983, 0.6179, 0.4072, 0.027204, 0.37200, 0.43762, 0.3000, 0.44643, 0.9878    )
SRL3Way<-c(0.10805, 0.7725, 0.4211, 0.8854, 0.3012, 0.7324, 0.3447, 0.225124, 0.03979, 0.72709 , 0.1660, 0.08078, 0.6399 )


#OTU effects are NS after correction
BHSRLOTU = p.adjust(SRLOTU, method = "BH")
print(BHSRLOTU)
#0.6562636 0.7403500 0.6562636 0.6562636 0.6562636 0.6562636 0.6562636 0.7883320 0.6562636 0.6562636 0.6562636 0.6562636 0.6562636


#N x OTU effects are NS after correction
BHSRLOTUN = p.adjust(SRLOTUN, method = "BH")
print(BHSRLOTUN)
#0.389415 0.954000 0.954000 0.954000 0.954000 0.954000 0.954000 0.954000 0.223080 0.954000 0.954000 0.954000 0.954000


#OTU effects are NS after correction
BHSRLOTUS = p.adjust(SRLOTUS, method = "BH")
print(BHSRLOTUS)
#0.7254487 0.7998250 0.7254487 0.7998250 0.7998250 0.7998250 0.7254487 0.3536520 0.7254487 0.7254487 0.7254487 0.7254487 0.9878000


#3 way effects are NS after correction
BHSRL3Way = p.adjust(SRL3Way, method = "BH")
print(BHSRL3Way)

#0.4682167 0.8368750 0.6842875 0.8854000 0.6401571 0.8368750 0.6401571 0.5853224 0.4682167 0.8368750 0.5395000 0.4682167 0.8368750
########################################
# all are NS after correcting for multiple comparisons 


########################################
## R.S 

########################################
str(df)


#OTU 1
#  X726a6245b76a9e74caea7993598fd969

modotu1<-lmer( log(R.S) ~  X726a6245b76a9e74caea7993598fd969*Nlevel*Sterility + (1|G5Block) , data=df)
summary(modotu1)
anova(modotu1)

#Response: R.S


#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)    
#X726a6245b76a9e74caea7993598fd969                  0.6661  0.6661     1 31.195  2.5479  0.12052    
#Nlevel                                             6.5624  6.5624     1 30.456 25.1015 2.18e-05 ***
#Sterility                                          1.1771  1.1771     1 30.844  4.5026  0.04198 *  
#X726a6245b76a9e74caea7993598fd969:Nlevel           1.0480  1.0480     1 31.525  4.0088  0.05393 .  
#X726a6245b76a9e74caea7993598fd969:Sterility        0.7926  0.7926     1 31.970  3.0317  0.09127 .  
#Nlevel:Sterility                                   0.1391  0.1391     1 30.749  0.5319  0.47133    
#X726a6245b76a9e74caea7993598fd969:Nlevel:Sterility 0.4582  0.4582     1 32.430  1.7525  0.19482  
#NS

#OTU 2
#  X7d3587bbb834c6ac86656e1d2715c331

modotu2<-lmer (log( R.S) ~ X7d3587bbb834c6ac86656e1d2715c331*Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu2)
anova(modotu2)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
#X7d3587bbb834c6ac86656e1d2715c331                  1.0735  1.0735     1    33  3.8093 0.0594984 .  
#Nlevel                                             4.5705  4.5705     1    33 16.2188 0.0003114 ***
#Sterility                                          0.0025  0.0025     1    33  0.0087 0.9260397    
#X7d3587bbb834c6ac86656e1d2715c331:Nlevel           1.0950  1.0950     1    33  3.8859 0.0571277 .  
#X7d3587bbb834c6ac86656e1d2715c331:Sterility        1.1520  1.1520     1    33  4.0880 0.0513587 .  
#Nlevel:Sterility                                   0.1363  0.1363     1    33  0.4838 0.4915643    
#X7d3587bbb834c6ac86656e1d2715c331:Nlevel:Sterility 1.0088  1.0088     1    33  3.5800 0.0672821 .  



#OTU 3
#a90d5b4368338dea15569269e8d33fb8 

modotu3<-lmer(log(R.S) ~a90d5b4368338dea15569269e8d33fb8 *Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu3)
anova(modotu3)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
#a90d5b4368338dea15569269e8d33fb8                  0.22746 0.22746     1 31.711  0.8348 0.36778  
#Nlevel                                            1.22275 1.22275     1 31.733  4.4877 0.04207 *
#Sterility                                         0.96521 0.96521     1 31.130  3.5425 0.06919 .
#a90d5b4368338dea15569269e8d33fb8:Nlevel           0.31268 0.31268     1 21.214  1.1476 0.29609  
#a90d5b4368338dea15569269e8d33fb8:Sterility        0.06205 0.06205     1 32.025  0.2277 0.63645  
#Nlevel:Sterility                                  0.04475 0.04475     1 32.931  0.1643 0.68789  
#a90d5b4368338dea15569269e8d33fb8:Nlevel:Sterility 0.84893 0.84893     1 30.674  3.1157 0.08750 .

#####################################################
#Not in the most abundant file
#############################################

#OTU 4
#  bea3f263c468bf3d3fb6645562737c70 

#modotu4<-lm ( log(R.S) ~ bea3f263c468bf3d3fb6645562737c70 *Nlevel*Sterility, data=df)
#summary(modotu4)
#anova(modotu4)

#Response: R.S

#bea3f263c468bf3d3fb6645562737c70                   1 0.002259 0.0022590  2.1424 0.1527369    

#bea3f263c468bf3d3fb6645562737c70:Nlevel            1 0.000534 0.0005341  0.5066 0.4816332    

#bea3f263c468bf3d3fb6645562737c70:Sterility         1 0.001028 0.0010277  0.9747 0.3307048    
    
#bea3f263c468bf3d3fb6645562737c70:Nlevel:Sterility  1 0.000006 0.0000058  0.0055 0.9413773  

#OTU 5
# X071e2a4c88faf23787b614210e527ce3

modotu5<-lmer(log(R.S) ~ X071e2a4c88faf23787b614210e527ce3*Nlevel*Sterility ++ (1|G5Block), data=df)
summary(modotu5)
anova(modotu5)
 

#Type III Analysis of Variance Table with Satterthwaite's method
 #                                                   Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
#X071e2a4c88faf23787b614210e527ce3                  0.00367 0.00367     1 32.681  0.0128 0.9106
#Nlevel                                             0.11047 0.11047     1 31.692  0.3848 0.5395
#Sterility                                          0.03116 0.03116     1 32.849  0.1085 0.7439
#X071e2a4c88faf23787b614210e527ce3:Nlevel           0.66624 0.66624     1 30.243  2.3205 0.1381
#X071e2a4c88faf23787b614210e527ce3:Sterility        0.10255 0.10255     1 32.986  0.3572 0.5542
#Nlevel:Sterility                                   0.25121 0.25121     1 31.606  0.8750 0.3567
#X071e2a4c88faf23787b614210e527ce3:Nlevel:Sterility 0.06064 0.06064     1 31.172  0.2112 0.6490

#OTU 6
#X2bb5ef7086fb6a049aa217591531d386 

modotu6<-lmer (log(R.S) ~ X2bb5ef7086fb6a049aa217591531d386 *Nlevel*Sterility +(1|G5Block), data=df)
summary(modotu6)
anova(modotu6)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                    Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
#X2bb5ef7086fb6a049aa217591531d386                  0.28687 0.28687     1    33  1.0443 0.31426  
#Nlevel                                             0.47992 0.47992     1    33  1.7471 0.19533  
#Sterility                                          0.12651 0.12651     1    33  0.4605 0.50210  
#X2bb5ef7086fb6a049aa217591531d386:Nlevel           0.59467 0.59467     1    33  2.1648 0.15067  
#X2bb5ef7086fb6a049aa217591531d386:Sterility        0.17595 0.17595     1    33  0.6405 0.42925  
#Nlevel:Sterility                                   1.02968 1.02968     1    33  3.7484 0.06146 .
#X2bb5ef7086fb6a049aa217591531d386:Nlevel:Sterility 0.15783 0.15783     1    33  0.5745 0.45384 



#OTU 7
#X3899d68a62bd9b58f6d28de79e6f5d1d

modotu7<-lmer (log(R.S) ~ X3899d68a62bd9b58f6d28de79e6f5d1d*Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu7)
anova(modotu7)

#Type III Analysis of Variance Table with Satterthwaite's method
 #                                                   Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
#X3899d68a62bd9b58f6d28de79e6f5d1d                  0.19986 0.19986     1 32.991  0.6755 0.417042   
#Nlevel                                             2.73542 2.73542     1 30.497  9.2454 0.004817 **
#Sterility                                          0.19402 0.19402     1 32.997  0.6558 0.423857   
#X3899d68a62bd9b58f6d28de79e6f5d1d:Nlevel           0.41961 0.41961     1 30.814  1.4182 0.242785   
#X3899d68a62bd9b58f6d28de79e6f5d1d:Sterility        0.03619 0.03619     1 31.964  0.1223 0.728820   
#Nlevel:Sterility                                   0.00850 0.00850     1 31.957  0.0287 0.866489   
#X3899d68a62bd9b58f6d28de79e6f5d1d:Nlevel:Sterility 0.24515 0.24515     1 32.354  0.8286 0.369415  

#OTU 8
#X453d47fbe68f91cd6b3aec11fc3b964b

modotu8<-lmer (log(R.S) ~ X453d47fbe68f91cd6b3aec11fc3b964b *Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu8)
anova(modotu8)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                    Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
#X453d47fbe68f91cd6b3aec11fc3b964b                  0.27164 0.27164     1    33  0.9797 0.32949  
#Nlevel                                             1.79249 1.79249     1    33  6.4646 0.01588 *
#Sterility                                          0.00092 0.00092     1    33  0.0033 0.95446  
#X453d47fbe68f91cd6b3aec11fc3b964b:Nlevel           0.03937 0.03937     1    33  0.1420 0.70871  
#X453d47fbe68f91cd6b3aec11fc3b964b:Sterility        0.05412 0.05412     1    33  0.1952 0.66152  
#Nlevel:Sterility                                   1.84469 1.84469     1    33  6.6529 0.01454 *
#X453d47fbe68f91cd6b3aec11fc3b964b:Nlevel:Sterility 1.23601 1.23601     1    33  4.4577 0.04240 *


#OTU 9
#X69d3d7cc762732aaf09fd2959046648a

modotu9<-lmer (log(R.S) ~ X69d3d7cc762732aaf09fd2959046648a*Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu9)
anova(modotu9)


#Type III Analysis of Variance Table with Satterthwaite's method
 #                                                   Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
#X69d3d7cc762732aaf09fd2959046648a                  0.57583 0.57583     1 32.045  2.2275 0.145353   
#Nlevel                                             2.47097 2.47097     1 31.361  9.5586 0.004156 **
#Sterility                                          1.54366 1.54366     1 32.088  5.9714 0.020219 * 
#X69d3d7cc762732aaf09fd2959046648a:Nlevel           0.28444 0.28444     1 30.948  1.1003 0.302322   
##X69d3d7cc762732aaf09fd2959046648a:Sterility        1.01004 1.01004     1 32.379  3.9072 0.056648 . 
#Nlevel:Sterility                                   0.12717 0.12717     1 30.278  0.4920 0.488415   
#X69d3d7cc762732aaf09fd2959046648a:Nlevel:Sterility 0.12434 0.12434     1 30.923  0.4810 0.493152 







#OTU 10
#X7fb63308d1a53c4af530608e9ccee61d

modotu10<-lmer (log(R.S) ~ X7fb63308d1a53c4af530608e9ccee61d*Nlevel*Sterility +(1|G5Block), data=df)
summary(modotu10)
anova(modotu10)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                    Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
#X7fb63308d1a53c4af530608e9ccee61d                  0.01349 0.01349     1    33  0.0480 0.82795  
#Nlevel                                             0.81032 0.81032     1    33  2.8822 0.09898 .
#Sterility                                          0.33781 0.33781     1    33  1.2016 0.28095  
#X7fb63308d1a53c4af530608e9ccee61d:Nlevel           0.44102 0.44102     1    33  1.5687 0.21921  
#X7fb63308d1a53c4af530608e9ccee61d:Sterility        0.12358 0.12358     1    33  0.4396 0.51194  
#Nlevel:Sterility                                   1.85634 1.85634     1    33  6.6028 0.01489 *
#X7fb63308d1a53c4af530608e9ccee61d:Nlevel:Sterility 1.49782 1.49782     1    33  5.3275 0.02740 *

#OTU 11
#X8864a30ab9e33ba131278de4a262ee0a

modotu11<-lmer (log(R.S) ~ X8864a30ab9e33ba131278de4a262ee0a*Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu11)
anova(modotu11)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#X8864a30ab9e33ba131278de4a262ee0a                  0.4834  0.4834     1 32.968  1.8749 0.1801678    
#Nlevel                                             4.7220  4.7220     1 30.705 18.3157 0.0001696 ***
#Sterility                                          1.7998  1.7998     1 32.403  6.9809 0.0125831 *  
#X8864a30ab9e33ba131278de4a262ee0a:Nlevel           0.0080  0.0080     1 32.258  0.0309 0.8615528    
#X8864a30ab9e33ba131278de4a262ee0a:Sterility        0.5981  0.5981     1 31.773  2.3201 0.1376066    
#Nlevel:Sterility                                   0.3313  0.3313     1 31.170  1.2849 0.2656245    
#X8864a30ab9e33ba131278de4a262ee0a:Nlevel:Sterility 0.0006  0.0006     1 32.291  0.0023 0.9624527  

#OTU 12
#d68451a2b18e9fd8e6738547074179fd
modotu12<-lmer (log(R.S) ~ d68451a2b18e9fd8e6738547074179fd*Nlevel*Sterility +(1|G5Block), data=df)
summary(modotu12)
anova(modotu12)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                  Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#d68451a2b18e9fd8e6738547074179fd                  0.0034  0.0034     1 30.923  0.0131 0.9097411    
#Nlevel                                            4.4119  4.4119     1 30.601 16.9934 0.0002647 ***
#Sterility                                         0.0852  0.0852     1 32.193  0.3280 0.5708207    
#d68451a2b18e9fd8e6738547074179fd:Nlevel           0.2826  0.2826     1 30.331  1.0886 0.3050208    
#d68451a2b18e9fd8e6738547074179fd:Sterility        0.0027  0.0027     1 32.570  0.0105 0.9191609    
#Nlevel:Sterility                                  0.6274  0.6274     1 30.792  2.4165 0.1302816    
#d68451a2b18e9fd8e6738547074179fd:Nlevel:Sterility 0.1392  0.1392     1 32.156  0.5362 0.4693047



#OTU 13
#d84b98e91b181d4c2879e8ed860ef94e 
modotu13<-lmer (log(R.S) ~ d84b98e91b181d4c2879e8ed860ef94e *Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu13)
anova(modotu13)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                   Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)  
#d84b98e91b181d4c2879e8ed860ef94e                  0.10761 0.10761     1 18.465  0.3666 0.5522  
#Nlevel                                            1.47450 1.47450     1 31.478  5.0231 0.0322 *
#Sterility                                         0.65630 0.65630     1 28.628  2.2358 0.1458  
#d84b98e91b181d4c2879e8ed860ef94e:Nlevel           0.04464 0.04464     1 32.843  0.1521 0.6991  
#d84b98e91b181d4c2879e8ed860ef94e:Sterility        0.13500 0.13500     1 29.128  0.4599 0.5030  
#Nlevel:Sterility                                  0.71998 0.71998     1 31.495  2.4527 0.1273  
#d84b98e91b181d4c2879e8ed860ef94e:Nlevel:Sterility 0.29733 0.29733     1 30.450  1.0129 0.3221    


#OTU 14
#ec4beab3963241c35c22f57155a96813
modotu14<-lmer (log(R.S) ~ ec4beab3963241c35c22f57155a96813 *Nlevel * Sterility + (1|G5Block), data=df)
summary(modotu14)
anova(modotu14)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                  Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#ec4beab3963241c35c22f57155a96813                  0.3565  0.3565     1 31.032  1.2660 0.2691389    
#Nlevel                                            5.4001  5.4001     1 32.157 19.1762 0.0001185 ***
#Sterility                                         0.3916  0.3916     1 32.839  1.3906 0.2467782    
#ec4beab3963241c35c22f57155a96813:Nlevel           0.8513  0.8513     1 32.470  3.0232 0.0915559 .  
#ec4beab3963241c35c22f57155a96813:Sterility        0.3273  0.3273     1 30.191  1.1624 0.2895079    
#Nlevel:Sterility                                  1.2264  1.2264     1 32.213  4.3550 0.0448971 *  
#ec4beab3963241c35c22f57155a96813:Nlevel:Sterility 0.8512  0.8512     1 32.954  3.0226 0.0914456 .  





#######Checked for mold effects in all these models and none

########################################
##
RSOTU<-c( 0.12052, 0.0594984,0.36778, 0.9106, 0.31426, 0.417042, 0.32949, 0.145353, 0.82795, 0.1801678 , 0.9097411,  0.5522, 0.2691389)   

RSOTUN<-c(0.05393 ,0.0571277, 0.29609, 0.1381, 0.15067, 0.242785, 0.70871, 0.302322, 0.21921 , 0.8615528, 0.3050208, 0.6991, 0.0915559  ) 

RSOTUS<-c(0.09127, 0.0513587, 0.63645, 0.5542, 0.42925, 0.728820, 0.66152, 0.056648, 0.51194 , 0.1376066 , 0.9191609, 0.5030 ,0.2895079) 

RS3way  <-c(0.19482, 0.0672821, 0.08750, 0.6490, 0.45384, 0.369415, 0.04240, 0.493152 , 0.02740, 0.9624527, 0.4693047, 0.3221 ,0.0914456 ) 

#OTU only
#All NS after MC
BHRSOTU = p.adjust(RSOTU, method = "BH")
print(BHRSOTU)
# [1] 0.5855453 0.5855453 0.5976425 0.9106000 0.5976425 0.6023940 0.5976425 0.5855453 0.9106000 0.5855453 0.9106000 0.7178600
#[13] 0.5976425

#OTU by N
#all NA after MC

BHRSOTUN = p.adjust(RSOTUN, method = "BH")
print(BHRSOTUN)
#0.3713300 0.3713300 0.3965270 0.3917420 0.3917420 0.3965270 0.7677692 0.3965270 0.3965270 0.8615528 0.3965270 0.7677692
#0.3917420

#OTU by sterility interactions 
#NS after correcting for multiple comparisons
BHRSOTUS = p.adjust(RSOTUS, method = "BH")
print(BHRSOTUS)
# 0.3955033 0.3682120 0.7817964 0.7817964 0.7817964 0.7895550 0.7817964 0.3682120 0.7817964 0.4472214 0.9191609 0.7817964
# 0.7527205


#3-Way
#NS after adjusting for multiple comparisons
BHRS3way = p.adjust(RS3way, method = "BH")
print(BHRS3way)
#0.4221100 0.2377586 0.2377586 0.7030833 0.5828160 0.5828160 0.2377586 0.5828160 0.2377586 0.9624527 0.5828160 0.5828160
#0.2377586
########################################



########################################
##  RootTissueDensity 

########################################



#OTU 1
#  X726a6245b76a9e74caea7993598fd969

modotu1<-lmer(RootTissueDensity  ~  X726a6245b76a9e74caea7993598fd969*Nlevel*Sterility  + (1|G5Block), data=df)
summary(modotu1)
anova(modotu1)

#Type III Analysis of Variance Table with Satterthwaite's method
   #                                                    Sum Sq    Mean Sq NumDF DenDF F value  Pr(>F)  
#X726a6245b76a9e74caea7993598fd969                  0.00075027 0.00075027     1    32  3.0864 0.08851 .
#Nlevel                                             0.00088109 0.00088109     1    32  3.6246 0.06596 .
#Sterility                                          0.00013878 0.00013878     1    32  0.5709 0.45543  
#X726a6245b76a9e74caea7993598fd969:Nlevel           0.00138817 0.00138817     1    32  5.7106 0.02292 *
#X726a6245b76a9e74caea7993598fd969:Sterility        0.00002724 0.00002724     1    32  0.1121 0.73998  
#Nlevel:Sterility                                   0.00020178 0.00020178     1    32  0.8301 0.36907  
#X726a6245b76a9e74caea7993598fd969:Nlevel:Sterility 0.00009046 0.00009046     1    32  0.3721 0.54616  




#OTU 2
#  X7d3587bbb834c6ac86656e1d2715c331

modotu2<-lmer (RootTissueDensity  ~ X7d3587bbb834c6ac86656e1d2715c331*Nlevel*Sterility +(1|G5Block), data=df)
summary(modotu2)
anova(modotu2)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                       Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#X7d3587bbb834c6ac86656e1d2715c331                  0.00005868 0.00005868     1    32  0.2055 0.6534
#Nlevel                                             0.00032824 0.00032824     1    32  1.1496 0.2916
#Sterility                                          0.00005485 0.00005485     1    32  0.1921 0.6641
#X7d3587bbb834c6ac86656e1d2715c331:Nlevel           0.00002316 0.00002316     1    32  0.0811 0.7776
#X7d3587bbb834c6ac86656e1d2715c331:Sterility        0.00005340 0.00005340     1    32  0.1870 0.6683
##Nlevel:Sterility                                   0.00056021 0.00056021     1    32  1.9621 0.1709
#X7d3587bbb834c6ac86656e1d2715c331:Nlevel:Sterility 0.00005181 0.00005181     1    32  0.1815 0.6730


#OTU 3
#a90d5b4368338dea15569269e8d33fb8 

modotu3<-lmer(  RootTissueDensity  ~a90d5b4368338dea15569269e8d33fb8 *Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu3)
anova(modotu3)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                      Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#a90d5b4368338dea15569269e8d33fb8                  0.00064579 0.00064579     1    32  2.4778 0.1253
#Nlevel                                            0.00001517 0.00001517     1    32  0.0582 0.8109
#Sterility                                         0.00061805 0.00061805     1    32  2.3713 0.1334
#a90d5b4368338dea15569269e8d33fb8:Nlevel           0.00007800 0.00007800     1    32  0.2993 0.5881
#a90d5b4368338dea15569269e8d33fb8:Sterility        0.00064661 0.00064661     1    32  2.4809 0.1251
#Nlevel:Sterility                                  0.00051199 0.00051199     1    32  1.9644 0.1707
#a90d5b4368338dea15569269e8d33fb8:Nlevel:Sterility 0.00026220 0.00026220     1    32  1.0060 0.3234


#OTU 4
#  bea3f263c468bf3d3fb6645562737c70 

#modotu4<-lm (  RootTissueDensity  ~ bea3f263c468bf3d3fb6645562737c70 *Nlevel*Sterility, data=df)
#summary(modotu4)
#anova(modotu4)

#Response: RootTissueDensity

#bea3f263c468bf3d3fb6645562737c70                   1 0.0001472 0.00014722  0.5120 0.4795

#bea3f263c468bf3d3fb6645562737c70:Nlevel            1 0.0000009 0.00000089  0.0031 0.9560

#bea3f263c468bf3d3fb6645562737c70:Sterility         1 0.0000026 0.00000259  0.0090 0.9250

#bea3f263c468bf3d3fb6645562737c70:Nlevel:Sterility  1 0.0004290 0.00042904  1.4920 0.2308

#OTU 5
# X071e2a4c88faf23787b614210e527ce3

modotu5<-lmer( RootTissueDensity  ~ X071e2a4c88faf23787b614210e527ce3*Nlevel*Sterility +(1|G5Block), data=df)
summary(modotu5)
anova(modotu5)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                       Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#X071e2a4c88faf23787b614210e527ce3                  0.00023431 0.00023431     1    32  0.9446 0.3384
#Nlevel                                             0.00017621 0.00017621     1    32  0.7104 0.4056
#Sterility                                          0.00049254 0.00049254     1    32  1.9856 0.1684
#X071e2a4c88faf23787b614210e527ce3:Nlevel           0.00005613 0.00005613     1    32  0.2263 0.6375
#X071e2a4c88faf23787b614210e527ce3:Sterility        0.00000522 0.00000522     1    32  0.0210 0.8856
#Nlevel:Sterility                                   0.00018521 0.00018521     1    32  0.7467 0.3940
#X071e2a4c88faf23787b614210e527ce3:Nlevel:Sterility 0.00024208 0.00024208     1    32  0.9759 0.3306

#OTU 6
#X2bb5ef7086fb6a049aa217591531d386 

modotu6<-lmer( RootTissueDensity  ~ X2bb5ef7086fb6a049aa217591531d386 *Nlevel*Sterility +(1|G5Block), data=df)
summary(modotu6)
anova(modotu6)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                       Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#X2bb5ef7086fb6a049aa217591531d386                  0.00014466 0.00014466     1    32  0.5388 0.4683
#Nlevel                                             0.00010651 0.00010651     1    32  0.3967 0.5333
#Sterility                                          0.00033655 0.00033655     1    32  1.2535 0.2712
#X2bb5ef7086fb6a049aa217591531d386:Nlevel           0.00006467 0.00006467     1    32  0.2409 0.6269
#X2bb5ef7086fb6a049aa217591531d386:Sterility        0.00000276 0.00000276     1    32  0.0103 0.9199
#Nlevel:Sterility                                   0.00016252 0.00016252     1    32  0.6053 0.4423
#X2bb5ef7086fb6a049aa217591531d386:Nlevel:Sterility 0.00044571 0.00044571     1    32  1.6601 0.2068


#OTU 7
#X3899d68a62bd9b58f6d28de79e6f5d1d

modotu7<-lmer ( RootTissueDensity  ~ X3899d68a62bd9b58f6d28de79e6f5d1d*Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu7)
anova(modotu7)

#Type III Analysis of Variance Table with Satterthwaite's method
 #                                                      Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#X3899d68a62bd9b58f6d28de79e6f5d1d                  0.00021283 0.00021283     1    32  0.7718 0.3862
#Nlevel                                             0.00000008 0.00000008     1    32  0.0003 0.9869
#Sterility                                          0.00058146 0.00058146     1    32  2.1087 0.1562
#X3899d68a62bd9b58f6d28de79e6f5d1d:Nlevel           0.00002042 0.00002042     1    32  0.0741 0.7873
#X3899d68a62bd9b58f6d28de79e6f5d1d:Sterility        0.00026741 0.00026741     1    32  0.9698 0.3321
#Nlevel:Sterility                                   0.00011845 0.00011845     1    32  0.4296 0.5169
#X3899d68a62bd9b58f6d28de79e6f5d1d:Nlevel:Sterility 0.00029619 0.00029619     1    32  1.0742 0.3078 


#OTU 8
#X453d47fbe68f91cd6b3aec11fc3b964b

modotu8<-lmer( RootTissueDensity  ~ X453d47fbe68f91cd6b3aec11fc3b964b *Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu8)
anova(modotu8)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                       Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#X453d47fbe68f91cd6b3aec11fc3b964b                  0.00000006 0.00000006     1    32  0.0002 0.9880
#Nlevel                                             0.00002411 0.00002411     1    32  0.0872 0.7697
#Sterility                                          0.00020037 0.00020037     1    32  0.7245 0.4010
#X453d47fbe68f91cd6b3aec11fc3b964b:Nlevel           0.00015457 0.00015457     1    32  0.5589 0.4602
#X453d47fbe68f91cd6b3aec11fc3b964b:Sterility        0.00063099 0.00063099     1    32  2.2816 0.1407
#Nlevel:Sterility                                   0.00000989 0.00000989     1    32  0.0358 0.8512
#X453d47fbe68f91cd6b3aec11fc3b964b:Nlevel:Sterility 0.00007094 0.00007094     1    32  0.2565 0.6160
 

#OTU 9
#X69d3d7cc762732aaf09fd2959046648a

modotu9<-lmer ( RootTissueDensity  ~ X69d3d7cc762732aaf09fd2959046648a*Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu9)
anova(modotu9)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                       Sum Sq    Mean Sq NumDF DenDF F value   Pr(>F)   
#X69d3d7cc762732aaf09fd2959046648a                  0.00016911 0.00016911     1    32  0.6899 0.412359   
#Nlevel                                             0.00016958 0.00016958     1    32  0.6918 0.411724   
#Sterility                                          0.00199362 0.00199362     1    32  8.1328 0.007553 **
#X69d3d7cc762732aaf09fd2959046648a:Nlevel           0.00059227 0.00059227     1    32  2.4161 0.129926   
#X69d3d7cc762732aaf09fd2959046648a:Sterility        0.00160361 0.00160361     1    32  6.5418 0.015476 * 
#Nlevel:Sterility                                   0.00001634 0.00001634     1    32  0.0667 0.797930   
#X69d3d7cc762732aaf09fd2959046648a:Nlevel:Sterility 0.00054097 0.00054097     1    32  2.2069 0.147184  


#OTU 10
#X7fb63308d1a53c4af530608e9ccee61d

modotu10<-lmer ( RootTissueDensity  ~ X7fb63308d1a53c4af530608e9ccee61d*Nlevel*Sterility +(1|G5Block), data=df)
summary(modotu10)
anova(modotu10)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                       Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#X7fb63308d1a53c4af530608e9ccee61d                  0.00002838 0.00002838     1    32  0.1034 0.7498
#Nlevel                                             0.00024279 0.00024279     1    32  0.8850 0.3539
#Sterility                                          0.00049775 0.00049775     1    32  1.8144 0.1874
#X7fb63308d1a53c4af530608e9ccee61d:Nlevel           0.00063965 0.00063965     1    32  2.3317 0.1366
#X7fb63308d1a53c4af530608e9ccee61d:Sterility        0.00021353 0.00021353     1    32  0.7784 0.3842
#Nlevel:Sterility                                   0.00007435 0.00007435     1    32  0.2710 0.6062
#X7fb63308d1a53c4af530608e9ccee61d:Nlevel:Sterility 0.00057391 0.00057391     1    32  2.0920 0.1578

#OTU 11
#X8864a30ab9e33ba131278de4a262ee0a

modotu11<-lmer (  RootTissueDensity  ~ X8864a30ab9e33ba131278de4a262ee0a*Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu11)
anova(modotu11)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                       Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#X8864a30ab9e33ba131278de4a262ee0a                  0.00069108 0.00069108     1    32  2.5564 0.1197
#Nlevel                                             0.00025439 0.00025439     1    32  0.9410 0.3393
#Sterility                                          0.00064968 0.00064968     1    32  2.4033 0.1309
#X8864a30ab9e33ba131278de4a262ee0a:Nlevel           0.00035300 0.00035300     1    32  1.3058 0.2616
#X8864a30ab9e33ba131278de4a262ee0a:Sterility        0.00000040 0.00000040     1    32  0.0015 0.9695
#Nlevel:Sterility                                   0.00015271 0.00015271     1    32  0.5649 0.4578
#X8864a30ab9e33ba131278de4a262ee0a:Nlevel:Sterility 0.00007948 0.00007948     1    32  0.2940 0.5914


#OTU 12
#d68451a2b18e9fd8e6738547074179fd
modotu12<-lmer( RootTissueDensity  ~ d68451a2b18e9fd8e6738547074179fd*Nlevel*Sterility +(1|G5Block), data=df)
summary(modotu12)
anova(modotu12)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                      Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#d68451a2b18e9fd8e6738547074179fd                  0.00005183 0.00005183     1    32  0.1957 0.6612
#Nlevel                                            0.00000809 0.00000809     1    32  0.0306 0.8623
#Sterility                                         0.00042303 0.00042303     1    32  1.5972 0.2154
#d68451a2b18e9fd8e6738547074179fd:Nlevel           0.00001304 0.00001304     1    32  0.0493 0.8258
#d68451a2b18e9fd8e6738547074179fd:Sterility        0.00025085 0.00025085     1    32  0.9471 0.3378
#Nlevel:Sterility                                  0.00009860 0.00009860     1    32  0.3723 0.5461
#d68451a2b18e9fd8e6738547074179fd:Nlevel:Sterility 0.00046245 0.00046245     1    32  1.7460 0.1958




#OTU 13
#d84b98e91b181d4c2879e8ed860ef94e 
modotu13<-lmer (  RootTissueDensity ~ d84b98e91b181d4c2879e8ed860ef94e *Nlevel*Sterility + (1|G5Block), data=df)
summary(modotu13)
anova(modotu13)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                      Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#d84b98e91b181d4c2879e8ed860ef94e                  0.00019793 0.00019793     1    32  0.6866 0.4134
#Nlevel                                            0.00001767 0.00001767     1    32  0.0613 0.8060
#Sterility                                         0.00037943 0.00037943     1    32  1.3163 0.2598
#d84b98e91b181d4c2879e8ed860ef94e:Nlevel           0.00011576 0.00011576     1    32  0.4016 0.5308
#d84b98e91b181d4c2879e8ed860ef94e:Sterility        0.00008408 0.00008408     1    32  0.2917 0.5929
#Nlevel:Sterility                                  0.00002343 0.00002343     1    32  0.0813 0.7774
#d84b98e91b181d4c2879e8ed860ef94e:Nlevel:Sterility 0.00025020 0.00025020     1    32  0.8680 0.3585

#OTU 14
#ec4beab3963241c35c22f57155a96813
modotu14<-lmer ( RootTissueDensity  ~ ec4beab3963241c35c22f57155a96813 *Nlevel * Sterility +(1|G5Block), data=df)
summary(modotu14)
anova(modotu14)

#Type III Analysis of Variance Table with Satterthwaite's method
#                                                      Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#ec4beab3963241c35c22f57155a96813                  1.4425e-05 1.4425e-05     1    32  0.0483 0.8275
#Nlevel                                            3.1808e-05 3.1808e-05     1    32  0.1064 0.7464
#Sterility                                         6.0794e-05 6.0794e-05     1    32  0.2034 0.6550
#ec4beab3963241c35c22f57155a96813:Nlevel           2.4520e-06 2.4520e-06     1    32  0.0082 0.9284
#ec4beab3963241c35c22f57155a96813:Sterility        3.5628e-05 3.5628e-05     1    32  0.1192 0.7321
#Nlevel:Sterility                                  5.4787e-05 5.4787e-05     1    32  0.1833 0.6714
#ec4beab3963241c35c22f57155a96813:Nlevel:Sterility 1.3560e-06 1.3560e-06     1    32  0.0045 0.9467


########################################
##NS after multicomparison correction
RTDOTU<-c( 0.08851, 0.6534, 0.1253, 0.3384, 0.4683, 0.3862, 0.9880,  0.412359, 0.7498, 0.1197, 0.6612, 0.4134,  0.8275)
RTDOTUN<-c(0.02292, 0.7776, 0.5881, 0.6375, 0.6269, 0.7873, 0.4602, 0.129926, 0.1366, 0.2616, 0.8258, 0.5308, 0.9284)
RTDOTUS<-c(0.73998, 0.6683, 0.1251, 0.8856, 0.9199, 0.3321, 0.1407, 0.015476, 0.3842, 0.9695, 0.3378, 0.5929,0.7321)
RTD3way  <-c(0.54616, 0.6730, 0.3234, 0.3306, 0.2068, 0.3078, 0.6160, 0.147184, 0.1578, 0.5914, 0.1958, 0.3585, 0.9467 )



#OTU 
BHRTDOTU = p.adjust(RTDOTU, method = "BH")
print(BHRTDOTU)

#NS after adjusting for MC
#0.5429667 0.8595600 0.5429667 0.7609875 0.7609875 0.7609875 0.9880000 0.7609875 0.8861273 0.5429667 0.8595600 0.7609875
#0.8964583



#N xOTU

BHRTDOTUN = p.adjust(RTDOTUN, method = "BH")
print(BHRTDOTUN)
#0.2979600 0.8946167 0.8946167 0.8946167 0.8946167 0.8946167 0.8946167 0.5919333 0.5919333 0.8502000 0.8946167 0.8946167
# 0.9284000
#NS after correcting for MC

#OTU x Sterility

BHRTDOTUS = p.adjust(RTDOTUS, method = "BH")
print(BHRTDOTUS)
#0.9619740 0.9619740 0.6097000 0.9695000 0.9695000 0.8324333 0.6097000 0.2011880 0.8324333 0.9695000 0.8324333 0.9619740
#0.9619740
#NS after correcting for MC





BHRTD3way = p.adjust(RTD3way, method = "BH")
print(BHRTD3way)
#0.7280000 0.7290833 0.5825625 0.5825625 0.5825625 0.5825625 0.7280000 0.5825625 0.5825625 0.7280000 0.5825625 0.5825625
#0.9467000
########################################



