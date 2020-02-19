###########################
#July 22, 2019
#Renee H. Petipas
#AVSelectionAnalysis_Revisions1
#Re-run selection analysis using all the samples including liquid inoculum that was 
#not included in the first submission. Made a new clean folder in my Dropbox so I can
# get the code and everything ready for submission
# Selection analysis for IJPS manuscript
# the final files used in the publication are version 8 (I tested for an effect of liquid inoculum and didn't find it so removed from models)
rm(list=ls(all=TRUE)) # clear R environment

#set working directory
setwd("~/Dropbox/PanicumMicrobiome_CleanFolder/data/Traits")

##########################

library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(MASS)
library(matrixStats)
library(lm.beta)



trait1<-read.csv("metadata.E2.csv")
#this file includes OTUs
str(trait1)
#All this code was to try and figure out samples different between 16s trait data and full trait data set
#missing from trait1 file are samples 136, 52, and 59
#manually added those back from the _Expt_Subset file
#traitfull<-read.csv("AV_Metadata_Expt_Subset.csv")

trait1 <- dplyr:: filter(trait1, CoreRootLength != "NA") # Remove NAs

str(trait1)

trait1$G5Block<-as.factor(trait1$G5Block)
trait1$InoculumVolume<-as.factor(trait1$InoculumVolume)
trait1$PresenceofMoldAug31<-as.factor(trait1$PresenceofMoldAug31) 

#Calculate R.S and RTD
trait1$Root_Tissue_Density <- trait1$CoreRootMass /trait1$CoreRootVolume
trait1$R.S <- trait1$CoreRootMass /trait1$TotalShootMass
head(trait1)
#trait1 <- trait1[,c(1:8,22,27:29,32,34:35)]

simple <- dplyr::filter(trait1, Sterility == "Sterile") 
complex <- dplyr::filter(trait1, Sterility == "NonSter")
simpleHigh <- dplyr::filter(simple, Nlevel == "High") # 22 obs. of  33 variables:
simpleLow <- dplyr::filter(simple, Nlevel == "Low") # 18 obs. of  33 variables:
complexHigh <- dplyr::filter(complex, Nlevel == "High") # 22 obs. of  33 variables:
complexLow <- dplyr::filter(complex, Nlevel == "Low") # 22 obs. of  33 variables:


#I don't think we need even sample sizes for lmer
#set.seed(777)

#simpleHigh <- dplyr::sample_n(simpleHigh, 9, replace = FALSE)
#complexHigh <- dplyr::sample_n(complexHigh, 9, replace = FALSE)
#complexLow <- dplyr::sample_n(complexLow, 9, replace = FALSE)

# Recombine for the full dataset

simples <- dplyr::union(simpleHigh, simpleLow)
complexes <- dplyr::union(complexHigh, complexLow)
trait2 <- dplyr::union(simples, complexes)

# Calculate relative fitness and scale to mean 0 SD 1 

trait2$Relative_Fitness <- trait2$TotalShootMass/mean(trait2$TotalShootMass)
head(trait2)

colnames(trait2 )[28] 
colnames(trait2 )[35] 
#this was previously 11 oops did I scale total shoot mass when I shouldn't have
#actually I had already calculated relative fitness so maybe it is fine 
trait2[,c(28:36)] <- scale(trait2[,c(28:36)])
head(trait2)







###########################################################################################

# Selection Differentials - Full


###########################################################################################
#removed SLA not contributing to the story and not a belowground trait

#singular fit
full_mod_diff <- lmer(Relative_Fitness ~ CoreRootLength  +(1|G5Block), data = trait2)
summary(full_mod_diff)
diff_cor <- summary(full_mod_diff)$coefficients

full_mod_diff <- lmer(Relative_Fitness ~ CoreMeanRootDiameter  + (1|G5Block)  , data = trait2)
diff_cor <- rbind(diff_cor, summary(full_mod_diff)$coefficients)

full_mod_diff <- lmer(Relative_Fitness ~ CoreSpecificRootLength  + (1|G5Block), data = trait2)
diff_cor <- rbind(diff_cor, summary(full_mod_diff)$coefficients)

full_mod_diff <- lmer(Relative_Fitness ~  Root_Tissue_Density   + (1|G5Block), data = trait2)
diff_cor <- rbind(diff_cor, summary(full_mod_diff)$coefficients)

full_mod_diff <- lmer(Relative_Fitness ~ R.S + (1|G5Block), data = trait2)
summary(full_mod_diff)
diff_cor <- rbind(diff_cor, summary(full_mod_diff)$coefficients)



write.table(diff_cor, file = "Sel_Diff_Full_V8.csv", sep = ",")

# Scale and relative fitness within groups
head(simpleHigh)
colnames(simpleHigh )[27] 
colnames(simpleHigh )[35] 

simpleHigh$Relative_Fitness <- simpleHigh$TotalShootMass /mean(simpleHigh$TotalShootMass )

simpleHigh[,c(28:35)] <- scale(simpleHigh[,c(28:35)])


simpleLow$Relative_Fitness <- simpleLow$TotalShootMass /mean(simpleLow$TotalShootMass )
simpleLow[,c(28:35)]<- scale(simpleLow[,c(28:35)])
head(simpleLow)

complexHigh$Relative_Fitness <- complexHigh$TotalShootMass /mean(complexHigh$TotalShootMass )
complexHigh[,c(28:35)] <- scale(complexHigh[,c(28:35)])


complexLow$Relative_Fitness <- complexLow$TotalShootMass /mean(complexLow$TotalShootMass )
complexLow[,c(28:35)]<- scale(complexLow[,c(28:35)])

# Selection Differentials - Separate

# simpleHigh (AKA sterile high N)
#removed SLA
head(simpleHigh)
sep_mod_diff <- lmer(Relative_Fitness ~  CoreRootLength  + (1|G5Block), data = simpleHigh)
diff_cor <- summary(sep_mod_diff)$coefficients

sep_mod_diff <- lmer(Relative_Fitness ~ CoreMeanRootDiameter + (1|G5Block), data = simpleHigh)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)


sep_mod_diff <- lmer(Relative_Fitness ~ CoreSpecificRootLength + (1|G5Block), data = simpleHigh)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)


sep_mod_diff <- lmer(Relative_Fitness ~ Root_Tissue_Density  + (1|G5Block), data = simpleHigh)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)


sep_mod_diff <- lmer(Relative_Fitness ~ R.S  + (1|G5Block), data = simpleHigh)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)



write.table(diff_cor, file = "Sel_Diff_simpleHigh_V8.csv", sep = ",")

# simpleLow
head(simpleLow)
#Singular fit
sep_mod_diff <- lmer(Relative_Fitness ~  CoreRootLength+ (1|G5Block), data = simpleLow)
diff_cor <- summary(sep_mod_diff)$coefficients

sep_mod_diff <- lmer(Relative_Fitness ~ CoreMeanRootDiameter  + (1|G5Block), data = simpleLow)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)

#Singular fit
sep_mod_diff <- lmer(Relative_Fitness ~ CoreSpecificRootLength + (1|G5Block), data = simpleLow)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)

sep_mod_diff <- lmer(Relative_Fitness ~ Root_Tissue_Density    + (1|G5Block), data = simpleLow)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)

#Singular fit
sep_mod_diff <- lmer(Relative_Fitness ~ R.S   + (1|G5Block), data = simpleLow)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)
summary(sep_mod_diff)

write.table(diff_cor, file = "Sel_Diff_simpleLow_V8.csv", sep = ",")


# complexHigh



sep_mod_diff <- lmer(Relative_Fitness ~  CoreRootLength + (1|G5Block), data = complexHigh)
diff_cor <- summary(sep_mod_diff)$coefficients

sep_mod_diff <- lmer(Relative_Fitness ~ CoreMeanRootDiameter + (1|G5Block), data = complexHigh)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)

sep_mod_diff <- lmer(Relative_Fitness ~ CoreSpecificRootLength + (1|G5Block), data = complexHigh)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)


sep_mod_diff <- lmer(Relative_Fitness ~ Root_Tissue_Density  + (1|G5Block), data = complexHigh)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)


sep_mod_diff <- lmer(Relative_Fitness ~ R.S + (1|G5Block), data = complexHigh)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)


write.table(diff_cor, file = "Sel_Diff_complexHigh_V8.csv", sep = ",")


# complexLow



#Singular fit
sep_mod_diff <- lmer(Relative_Fitness ~  CoreRootLength + (1|G5Block), data =  complexLow)
diff_cor <- summary(sep_mod_diff)$coefficients

sep_mod_diff <- lmer(Relative_Fitness ~ CoreMeanRootDiameter  + (1|G5Block), data =  complexLow)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)

sep_mod_diff <- lmer(Relative_Fitness ~ CoreSpecificRootLength  + (1|G5Block), data =  complexLow)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)


sep_mod_diff <- lmer(Relative_Fitness ~ Root_Tissue_Density  + (1|G5Block), data =  complexLow)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)


sep_mod_diff <- lmer(Relative_Fitness ~ R.S  + (1|G5Block), data =  complexLow)
diff_cor <- rbind(diff_cor, summary(sep_mod_diff)$coefficients)
summary(sep_mod_diff)
write.table(diff_cor, file = "Sel_Diff_complexLow_V8.csv", sep = ",")

# Selection gradients - Both full and separate
# Try removing date planted
#singular fit

full_grad_mod <- lmer(Relative_Fitness ~ CoreRootLength  + 
                        CoreMeanRootDiameter + CoreSpecificRootLength + 
                        Root_Tissue_Density  + R.S + (1|G5Block),data = trait2, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000000)))
summary(full_grad_mod)
anova(full_grad_mod)

# # 
isSingular(full_grad_mod, tol = 1e-05)
# # Removed selection gradient analysis for the separate groups-because small sample sizes

# # 
# # 
# # 



grad_coeff <- summary(full_grad_mod)$coefficients


write.table(grad_coeff, file = "Sel_Grads_V7.csv", sep = ",")

# ANCOVAs 

# Differentials

# Removed SLA
# Core_Root_Length         
# Core_Mean_Root_Diameter  
# Core_Specific_Root_Length
# Root_Tissue_Density
# Root_Shoot_Ratio
#The ancova is asking if the relationship between relative fitness and the trait changes under the different nutrient/microbe treatments
#
str(trait2)
diff_mod2 <- lmer(Relative_Fitness ~ CoreRootLength * Sterility * Nlevel   +  (1|G5Block), data = trait2)
anova(diff_mod2)
diff_mod3 <- lmer(Relative_Fitness ~ CoreMeanRootDiameter * Sterility * Nlevel  + (1|G5Block), data = trait2)
anova(diff_mod3)
diff_mod4 <- lmer(Relative_Fitness ~ CoreSpecificRootLength * Sterility * Nlevel + (1|G5Block), data = trait2)
anova(diff_mod4)
diff_mod5 <- lmer(Relative_Fitness ~  Root_Tissue_Density  * Sterility * Nlevel  + (1|G5Block), data = trait2)
anova(diff_mod5)
diff_mod6 <- lmer(Relative_Fitness ~ R.S * Sterility *Nlevel   + (1|G5Block), data = trait2)
anova(diff_mod6)



ancova_diff <- anova(diff_mod2)
ancova_diff <- rbind(ancova_diff, anova(diff_mod3))
ancova_diff <- rbind(ancova_diff, anova(diff_mod4))
ancova_diff <- rbind(ancova_diff, anova(diff_mod5))
ancova_diff <- rbind(ancova_diff, anova(diff_mod6))


write.table(ancova_diff, file = "ancova_diff_V8.csv", sep = ",")
library(emmeans)

emmeans(diff_mod5, list(pairwise~Sterility*Nlevel), adjust = "tukey")

# Gradients
str(trait2)
# SLA
# Core_Root_Length         
# Core_Mean_Root_Diameter  
# Core_Specific_Root_Length
# Root_Tissue_Density
# Root_Shoot_Ratio

#Singular fit
grad_mod1 <- lmer(Relative_Fitness ~ (CoreRootLength  +
                                        CoreMeanRootDiameter + CoreSpecificRootLength +
                                        Root_Tissue_Density  + R.S ) * Sterility *
                    Nlevel  + (1|G5Block), data = trait2)

ancova_grad <- anova(grad_mod1)




write.table(ancova_grad, file = "ancova_grad_V8.csv", sep = ",")





##################Graphs for Publication##########################################################
#load libraries
library(ggplot2)
library(cowplot)
#########
myColors <- c("grey36", "gray74")
names(myColors) <- levels(trait2$Nlevel)
str(trait2)
colScale <- scale_fill_manual(name = "Nitrogen \nLevel", values = myColors)

table(trait2$Sterility)
levels(trait2$Sterility) <- c("Unperturbed", "Perturbed")
theme_set(theme_cowplot())
##########

#Plot1
# Make ggplot of relative fitness x root length
p_rl3 <- ggplot(trait2, aes(x=CoreRootLength, y=Relative_Fitness)) + geom_point(aes(x=CoreRootLength, y=Relative_Fitness, shape=Sterility, color=Nlevel)) + theme_classic() + xlab("Root Length") + ylab("Relative Fitness") + stat_smooth(method= "lm", col = "black") 
p_rl3



#Plot2
# Make ggplot of relative fitness x root length facetted by nitrogen

p_rs1 <- ggplot(trait2, aes(x=R.S, y=Relative_Fitness)) + geom_point(aes(x=R.S, y=Relative_Fitness, shape=Sterility, color=Nlevel)) + theme_classic() + xlab("Root to Shoot Ratio") + ylab("Relative Fitness") + stat_smooth(method= "lm", col = "black") 

# Divide by day, going horizontally and wrapping with 2 columns
p_rs1 <- p_rs1 + facet_wrap( ~ Nlevel, ncol=2)
p_rs1


#Plot3
# Make ggplot of relative fitness x RTD
p_rtd1 <- ggplot(trait2, aes(x=Root_Tissue_Density, y=Relative_Fitness)) + geom_point(aes(x=Root_Tissue_Density, y=Relative_Fitness, shape=Sterility, color=Nlevel)) + theme_classic() + xlab("Root Tissue Density (g cm-3)") + ylab("Relative Fitness") + stat_smooth(method= "lm", col = "black") 

# Divide by day, going horizontally and wrapping with 2 columns
p_rtd1 <- p_rtd1 + facet_wrap( ~ Sterility*Nlevel, ncol=2)
p_rtd1


# ==== CJ The real real plots for publication ====
# Figure 5
#  Micro vs root length

modotu5<-lmer (log(CoreRootLength)~ X071e2a4c88faf23787b614210e527ce3*Nlevel*Sterility +  (1|G5Block), data=df)

df<-read.csv("~/Dropbox/PanicumMicrobiome_CleanFolder/data/MicrobiomexTrait/MetaTableTotal_RelAbun_FINAL_05Aug2019.csv")
str(df)

df$G5Block<-as.factor(df$G5Block)
df$PresenceofMoldAug31 <-as.factor(df$PresenceofMoldAug31 )
df$InoculumVolume<-as.factor(df$InoculumVolume)

myColors2 <- c("grey36", "gray74")
names(myColors2) <- levels(df$Nlevel)
colScale2 <- scale_colour_manual(name = "Nitrogen \nLevel", values = myColors2)

table(df$Sterility)
levels(df$Sterility) <- c("Unperturbed", "Perturbed")

Fig5A <- ggplot(df, aes(x = X071e2a4c88faf23787b614210e527ce3, y = CoreRootLength))+ colScale2  + geom_point(aes(colour = Nlevel, shape = Sterility), size = 2.5) + labs(x = expression(paste(italic("Micromonosporo "), sp.)), y = "Root Length (cm)") + stat_smooth(method = "lm", colour = "black", fill = "grey80")+ guides(colour = guide_legend(order = 1, override.aes = list(shape = 15)), shape = guide_legend(order = 2)) + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 13))

# Root lenght relative fitness
trait1<-read.csv("~/Dropbox/PanicumMicrobiome_CleanFolder/data/Traits/metadata.E2.csv")
trait1 <- dplyr:: filter(trait1, CoreRootLength != "NA") # Remove NAs
trait1$G5Block<-as.factor(trait1$G5Block)
trait1$InoculumVolume<-as.factor(trait1$InoculumVolume)
trait1$PresenceofMoldAug31<-as.factor(trait1$PresenceofMoldAug31)
trait1$Root_Tissue_Density <- trait1$CoreRootMass /trait1$CoreRootVolume
trait1$R.S <- trait1$CoreRootMass /trait1$TotalShootMass

simple <- dplyr::filter(trait1, Sterility == "Sterile")
complex <- dplyr::filter(trait1, Sterility == "NonSter")
simpleHigh <- dplyr::filter(simple, Nlevel == "High") # 22 obs. of  33 variables:
simpleLow <- dplyr::filter(simple, Nlevel == "Low") # 18 obs. of  33 variables:
complexHigh <- dplyr::filter(complex, Nlevel == "High") # 22 obs. of  33 variables:
complexLow <- dplyr::filter(complex, Nlevel == "Low") # 22 obs. of  33 variables:
simples <- dplyr::union(simpleHigh, simpleLow)
complexes <- dplyr::union(complexHigh, complexLow)
trait2 <- dplyr::union(simples, complexes)
trait2$Relative_Fitness <- trait2$TotalShootMass/mean(trait2$TotalShootMass)
trait2[,c(28:36)] <- scale(trait2[,c(28:36)])

Fig5B <-ggplot(trait2, aes(x = CoreRootLength, y = Relative_Fitness))+ colScale2  + geom_point(aes(colour = Nlevel, shape = Sterility), size = 2.5) + stat_smooth(method = "lm", colour = "black", fill = "grey80")+ labs(x = "Root Length (cm)", y = "Relative Fitness") + guides(colour = guide_legend(order = 1, override.aes = list(shape = 15)), shape = guide_legend(order = 2)) + theme(legend.position = "none") + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 13))

Mylegend <- get_legend(Fig5A + theme(legend.box.margin = margin(0, 0, 0, 5))) 

Fig5.1 <-plot_grid(Fig5A + theme(legend.position = "none"), Fig5B, nrow = 1, labels = c("A", "B"),  align = "hv", axis = "l")

Fig5<- plot_grid(Fig5.1, Mylegend, rel_widths = c(2, .25))

save_plot("Fig5_12Aug2019.pdf", Fig5, nrow = 1, base_width = 14)

# Fig 6- Rel fitne vs r:s

levels(trait2$Sterility) <- c("Unperturbed", "Perturbed")

ggplot(trait2, aes(x = R.S, y = Relative_Fitness, colour = Nlevel))  + labs(x = "Root to Shoot Ratio", y = "Relative Fitness") + geom_point(aes(shape = Sterility), size = 2.5) + geom_smooth(aes(group = Nlevel),method = "lm", fill = "grey80", colour = "black") + colScale2 + guides(colour = guide_legend(order = 1, override.aes = list(shape = 15)), shape = guide_legend(order = 2))+ theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14))

ggsave("Fig6_12Aug2019.pdf")



ggplot(trait2, aes(x = Root_Tissue_Density, y = Relative_Fitness, colour = Nlevel, shape = Sterility))  + labs(x = "Root Length (cm)", y = "Relative Fitness") + geom_point(size = 2.5) + geom_smooth(method = "lm", fill = "grey80", colour = "black") + facet_wrap(~Sterility)# + colScale2 + guides(colour = guide_legend(order = 1, override.aes = list(shape = 15)), shape = guide_legend(order = 2))


# Figure 7

trait2$Nster <- interaction(trait2$Nlevel, trait2$Sterility, sep = ".")

Fig7A <- ggplot(subset(trait2, Sterility == "Unperturbed"), aes(x = Root_Tissue_Density, y = Relative_Fitness, colour = Nlevel, shape = Sterility))  + labs(x = "Root tissue density", y = "Relative Fitness") + geom_point(size = 2.5) + geom_smooth(aes(group = Nster),method = "lm", fill = "grey80", colour = "black") + coord_cartesian(ylim = c(0,6), xlim = c(-2,6))+ colScale2 + guides(colour = guide_legend(order = 1, override.aes = list(shape = 15)), shape = guide_legend(order = 2)) + theme(legend.position = "none") + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14))

Fig7B <- ggplot(subset(trait2, Sterility == "Perturbed"), aes(x = Root_Tissue_Density, y = Relative_Fitness, colour = Nlevel, shape = Sterility))  + labs(x = "Root tissue density", y = "Relative Fitness") + geom_point(shape = 17,size = 2.5) + geom_smooth(aes(group = Nster),method = "lm", fill = "grey80", colour = "black") + colScale2 + guides(colour = guide_legend(order = 1, override.aes = list(shape = 15)), shape = guide_legend(order = 2)) + coord_cartesian(ylim = c(0,6), xlim = c(-2,6)) + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14))

Mylegend <- get_legend(Fig5A + theme(legend.box.margin = margin(0, 0, 0, 5))) 

Fig7.1 <-plot_grid(Fig7A , Fig7B+ theme(legend.position = "none"), nrow = 1, labels = c("A", "B"),  align = "hv", axis = "l")

Fig7<- plot_grid(Fig7.1, Mylegend, rel_widths = c(2, .25))

save_plot("Fig7_12Aug2019.pdf", Fig7, nrow = 1, base_width = 14)

# #lm model without random effects
# grad_mod1 <- lm(Relative_Fitness ~ (CoreRootLength  + 
#                                         CoreMeanRootDiameter + CoreSpecificRootLength + 
#                                         RootTissueDensity  + R.S + X071e2a4c88faf23787b614210e527ce3 +
#                                         X69d3d7cc762732aaf09fd2959046648a) * Sterility * 
#                     Nlevel, data = trait2)
# 
# ancova_grad <- anova(grad_mod1)
# 
# 
# 
# 
# write.table(ancova_grad, file = "ancova_grad_V5.csv", sep = ",")
#lm models=no random effects
# Selection gradients - Both full and separate
# 
# full_grad_mod <- lm(Relative_Fitness ~ CoreRootLength  +
#                          CoreMeanRootDiameter + CoreSpecificRootLength +
#                          RootTissueDensity  + R.S + X071e2a4c88faf23787b614210e527ce3 +
#                          X69d3d7cc762732aaf09fd2959046648a, data = trait2)
# 
# 
# sH_selgrad <- lm(Relative_Fitness ~ CoreRootLength  +
#                       CoreMeanRootDiameter + CoreSpecificRootLength +
#                      RootTissueDensity  + R.S + X071e2a4c88faf23787b614210e527ce3 +
#                       X69d3d7cc762732aaf09fd2959046648a, data = simpleHigh)
# 
# 
# sL_selgrad <- lm(Relative_Fitness ~ CoreRootLength  +
#                      CoreMeanRootDiameter + CoreSpecificRootLength +
#                       RootTissueDensity  + R.S + X071e2a4c88faf23787b614210e527ce3 +
#                      X69d3d7cc762732aaf09fd2959046648a, data = simpleLow)
# 
# 
# #
# cH_selgrad <- lm(Relative_Fitness ~ CoreRootLength  +
#                       CoreMeanRootDiameter + CoreSpecificRootLength +
#                       RootTissueDensity  + R.S + X071e2a4c88faf23787b614210e527ce3 +
#                       X69d3d7cc762732aaf09fd2959046648a, data = complexHigh)
# 
# 
# cL_selgrad <- lm(Relative_Fitness ~ CoreRootLength  +
#                       CoreMeanRootDiameter + CoreSpecificRootLength +
#                      RootTissueDensity  + R.S + X071e2a4c88faf23787b614210e527ce3 +
#                      X69d3d7cc762732aaf09fd2959046648a, data = complexLow)
# 

#