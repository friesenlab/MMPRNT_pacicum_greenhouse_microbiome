
#title: "AV_Richness_Ordination_July2019"
#author: "Renee H. Petipas"
#final pre-publication code used to generate NMDS plots, run PERMANOV, and do richness analysis
#date: August 1, 2019
#Emily re-ran all the qiime2 qc steps 


#To clear the R environment
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
library(scatterplot3d)
library(microbiome)
library(dplyr)
#I couldn't read tree before I added this package (WTF?)
library(ape)
library(vegan)
library(lme4)
library(lmerTest)


setwd("~/Dropbox/PanicumMicrobiome_CleanFolder/data/Microbiome/PhyloseqFiles")
#### Creating Trim0 Phyloseq ####
#Import samples

metadata.E <-read.csv("AV_TraitDataE_Aug2019.csv", row.names = 1)
head(metadata.E)
metadata.E[metadata.E == "na"] <- NA
str(metadata.E)
write.csv(metadata.E, "metadata.E2.csv")

metadata.E <-read.csv("metadata.E2.csv", row.names = 1)
str(metadata.E)

#Calling things as.numeric changed the values :0 :0 also it is just easier to read in a new csv file

"""
metadata.E$NumberofExpandedLeavesAug31<-as.integer(metadata.E$NumberofExpandedLeavesAug31)
metadata.E$PresenceofMoldAug31<-as.factor(metadata.E$PresenceofMoldAug31) 
metadata.E$HeightAug31<-as.numeric(metadata.E$HeightAug31) 
metadata.E$LengthLargestLeafAug31  <-as.numeric(metadata.E$LengthLargestLeafAug31  ) 
metadata.E$ HeightSep23<-as.numeric(metadata.E$ HeightSep23) 
metadata.E$ HeightSep29<-as.numeric(metadata.E$ HeightSep29) 
metadata.E$ NumberofStemsSep29<-as.integer(metadata.E$ NumberofStemsSep29) 
metadata.E$ TotalHeightSep29<-as.numeric(metadata.E$ TotalHeightSep29) 
metadata.E$ SoilNO3 <-as.numeric(metadata.E$SoilNO3 ) 
metadata.E$ SoilNH4 <-as.numeric(metadata.E$SoilNH4 ) 
metadata.E$ Soildrymattercontent <-as.numeric(metadata.E$Soildrymattercontent) 
metadata.E$ TopLeafMass <-as.numeric(metadata.E$TopLeafMass) 
metadata.E$ TopLeafArea <-as.numeric(metadata.E$TopLeafArea) 
metadata.E$ SpecificLeafArea  <-as.numeric(metadata.E$ SpecificLeafArea )
metadata.E$ RemainingLeafMass  <-as.numeric(metadata.E$ RemainingLeafMass)
metadata.E$ TotalLeafMass   <-as.numeric(metadata.E$ TotalLeafMass )
metadata.E$ Stems   <-as.integer(metadata.E$ Stems )
metadata.E$ Flowers    <-as.numeric(metadata.E$ Flowers )
metadata.E$ TotalShootMass    <-as.numeric(metadata.E$ TotalShootMass)
metadata.E$ CoreRootLength    <-as.numeric(metadata.E$ CoreRootLength )
metadata.E$  CoreMeanRootDiameter    <-as.numeric(metadata.E$  CoreMeanRootDiameter )
metadata.E$  CoreRootSurfaceArea<-as.numeric(metadata.E$ CoreRootSurfaceArea)
metadata.E$  CoreRootVolume     <-as.numeric(metadata.E$ CoreRootVolume  )

metadata.E$  CoreRootVolume    <-as.numeric(metadata.E$  CoreRootVolume )
metadata.E$  CoreSpecificRootLength    <-as.numeric(metadata.E$ CoreSpecificRootLength )
metadata.E$  CoreRootMass    <-as.numeric(metadata.E$ CoreRootMass )
"""
#change so R recognizes factors as factors
metadata.E$G5Block<-as.factor(metadata.E$G5Block)
metadata.E$InoculumVolume<-as.factor(metadata.E$InoculumVolume)
metadata.E$PresenceofMoldAug31<-as.factor(metadata.E$PresenceofMoldAug31) 



#recheck structure of the data
str(metadata.E)



#OTU table
otu.table.E <- read.csv("AV_OtuTable0E_Aug2019.csv", row.names = 1)






#because PERMANOVA requires equal sample sizes I randomly removed samples 128, 135, 137, 63,65 (picked by R :)) 
#It doesn't make much a differences for the analysis so I am using the full data file


#This OTU table was with the samples removed
#otu.table.E <- read.csv("20190219_AV001_OtuTable0E_001_Subset.csv", row.names = 1)


tax.table <- read.csv("AV_Taxonomy0_Aug2019.csv", row.names = 1)

tree.data <- read.tree("AV_Rooted0_Aug2019.nwk")

#Convert necessary files to Matrices

otu.mat.E <- as.matrix(otu.table.E)
tax.mat <- as.matrix(tax.table)

#Convert files to Phyloseq Objects

METAE <- sample_data(metadata.E)
TAX   <- tax_table(tax.mat)

OTUE  <- otu_table(otu.mat.E, taxa_are_rows = TRUE)
#make sure OTU names are consistent across objects
taxa_names(TAX)
taxa_names(OTUE)
taxa_names(tree.data)
#make sure sample names are consistent across objects

sample_names(OTUE)
sample_names(METAE)

#Make Phyloseq-Class Object

Ephylo0 <- phyloseq(OTUE, TAX, METAE, tree.data)

#Check that making phyloseq object worked

Ephylo0

#before removing chloroplasts and mitochondria there are 7091 taxa


#Filtering out Chloroplasts & Mitochondria in Phyloseq
#Remove chloroplast from taxonomy

Ephylo0 <-subset_taxa(Ephylo0, !(class %in% c("Chloroplast")))
#Remove mitochondria from taxonomy
Ephylo0 <-subset_taxa(Ephylo0, !(family %in% c("mitochondria")))
#Results in 163 fewer OTUs
Ephylo0


#6928 taxa after removing

#see what sampling depth to rarefy to
sample_sums(otu_table(Ephylo0))
set.seed(394679238)
#rarefy to even depth
#where did this 12,940 come from???
#Ephylo0 <- rarefy_even_depth(Ephylo0, sample.size = 12940)
Ephylo0 <- rarefy_even_depth(Ephylo0, sample.size =  15558, rngseed=F)
#I think we were rarefying to the wrong value-not sure where the 12,000 came from but when I checked it -the lowest value was 15,000 
#so I changed it-fortunately it doesn't really change the results-minor changes to F and R2
Ephylo0
#4234 taxa after rarefaction

#Write richness table for analysis
rich_table <- estimate_richness(Ephylo0, measures = c("Observed","Simpson", "InvSimpson", "Shannon", "Chao1", "ACE","se.ACE"))

write.csv(rich_table, file = "richness_table.csv")

#also tested not rarefied richness (as suggested by authors of package) and it is not qualitatively different than rarefied analysis


############################
#Introduce a step to calculate relative abundance????
#I don't think this appropriate when using bray curtis...
############################

# Calculate compositional version of the full data set
# (relative abundances)
#Ephylo0 <- microbiome::transform(Ephylo0, "compositional")

######################################################################
#Create df for permanova using phyloseq object
######################################################################

bacplantdat <- data.frame(sample_data(Ephylo0))
bacotudat <- data.frame(t(otu_table(Ephylo0)))
head(bacplantdat)
#based on comparisons below use count data and even sample sizes (although it doesn't matter too much)
#Pretty sure that PERMANOVA (if multilevel requires equal sample sizes)
#PERMANOVA for bray
#strata need to be balanced
#G5block and date planted are overlapping for example all G5 #3 were planted on the third...well that explains some weirdness...
adonisbray<-adonis(bacotudat ~ Nlevel * Sterility + G5Block + PresenceofMoldAug31, bacplantdat, perm=999, distance="bray", contr.unordered = "contr.sum")
adonisbray



##Final model included in the manuscript
#Using rarefied count data (but not removed rare samples), uneven sample sizes, ignoring date planted (because it overlaps with G5 Block)-and including mold (Although very few non-mold samples)
#Permutation: free
#Number of permutations: 999




#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#                  Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Nlevel            1    0.5320 0.53198  2.2262 0.04959  0.001 ***
#  Sterility       1    0.7768 0.77681  3.2508 0.07241  0.001 ***
##  G5Block        3    1.0243 0.34143  1.4288 0.09548  0.001 ***
# PresenceofMold   1    0.2388 0.23882  0.9994 0.02226  0.439    
#Nlevel:Sterilit   1    0.2701 0.27009  1.1303 0.02518  0.198    
#Residuals         33    7.8858 0.23896         0.73508           
#Total             40   10.7278                 1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Using count data, even sample sizes, ignoring date planted (because it overlaps with G5 Block)-and including mold (Although very few non-mold samples)
#Terms added sequentially (first to last)

#                       Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#  Nlevel               1    0.4529 0.45289  2.1777 0.05350  0.002 ** 
#  Sterility            1    0.8174 0.81737  3.9302 0.09656  0.001 ***
#  G5Block              3    0.9150 0.30500  1.4666 0.10810  0.003 ** 
#  PresenceofMoldAug31  1    0.2046 0.20455  0.9836 0.02417  0.484    
#  Nlevel:Sterility     1    0.2519 0.25190  1.2112 0.02976  0.157    
#Residuals           28    5.8231 0.20797         0.68792           
#Total               35    8.4648                 1.00000           







#Using relative abundance with even sample sizes
#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Nlevel            1    0.4668 0.46675  2.4278 0.05906  0.001 ***
#  Sterility         1    0.7722 0.77224  4.0168 0.09772  0.001 ***
#  G5Block           3    0.8676 0.28920  1.5043 0.10979  0.001 ***
#  Nlevel:Sterility  1    0.2208 0.22080  1.1485 0.02794  0.191    
#Residuals        29    5.5754 0.19225         0.70550           
#Total            35    7.9028                 1.00000          


#Using count data with even sample sizes
#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Nlevel            1    0.4529 0.45289  2.1843 0.05350  0.002 ** 
#  Sterility         1    0.8174 0.81737  3.9423 0.09656  0.001 ***
#  G5Block           3    0.9150 0.30500  1.4711 0.10810  0.002 ** 
#  Nlevel:Sterility  1    0.2669 0.26691  1.2874 0.03153  0.089 .  
#Residuals        29    6.0127 0.20733         0.71031           
#Total            35    8.4648                 1.00000     

#Using count data with uneven sample sizes
#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Nlevel            1    0.5111 0.51114  2.4342 0.05266  0.001 ***
#  Sterility         1    0.8267 0.82669  3.9370 0.08517  0.001 ***
#  G5Block           3    0.9728 0.32427  1.5443 0.10022  0.003 ** 
#  Nlevel:Sterility  1    0.2565 0.25648  1.2214 0.02642  0.120    
#Residuals        34    7.1394 0.20998         0.73553           
#Total            40    9.7065                 1.00000        



#Using relative abundance data with uneven sample sizes
#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Nlevel            1    0.5287 0.52872  2.7158 0.05826  0.001 ***
#  Sterility         1    0.7742 0.77423  3.9768 0.08531  0.001 ***
#  G5Block           3    0.9269 0.30898  1.5871 0.10214  0.001 ***
#  Nlevel:Sterility  1    0.2263 0.22626  1.1622 0.02493  0.193    
#Residuals        34    6.6193 0.19468         0.72936           
#Total            40    9.0754                 1.00000     


#Using count data with even sample sizes but not trimming the data to remove rare OTUs
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Nlevel            1    0.4669 0.46688  1.9099 0.04832  0.001 ***
#  Sterility         1    0.8046 0.80457  3.2913 0.08327  0.001 ***
#  G5Block           3    1.0078 0.33593  1.3742 0.10430  0.003 ** 
#  Nlevel:Sterility  1    0.2936 0.29361  1.2011 0.03039  0.121    
#Residuals        29    7.0891 0.24445         0.73371           
#Total            35    9.6619                 1.00000      




#Removed rare samples to clean-up ordination-this doesn't affect the PERMANOVA but I moved it to after just to clean up Mds
#Remove OTUs that do not show appear more than 2 times in more than 10% of the samples
#
wh0 = genefilter_sample(Ephylo0, filterfun_sample(function(x) x > 2), A=0.1*nsamples(Ephylo0))
Ephylo0 = prune_taxa(wh0, Ephylo0)
#Not conservative enough-although went from 4234 taxa to 692 but my thresholds are not as low as in tutorial
Ephylo0
#now only 692 samples left

#Plotting using NMDS  ******************************* cj
vare.mds<-metaMDS(bacotudat, distance="bray", trace=FALSE)
vare.mds
plot(vare.mds)
vare.mds<-metaMDS(bacotudat, trace=FALSE)
vare.mds
plot(vare.mds)


#only used ordihull in the plot 'labeledAMFdiversity_NMDSall'
with(bacplantdat, ordihull(vare.mds, Sterility, col="blue", lty=6, label=TRUE))
#read data-using 12 that actually have soils data associated


#colvec<-c("red", "gold2")
colvec<- c("grey36", "gray74")
#type="n" clears away all points for rebuilding
#plot(vare.mds, type="n", xaxt='n', yaxt='n', xlab="NMDS axis 1", ylab="NMDS axis 2")
plot(vare.mds, type="n")
str(bacplantdat)
with(bacplantdat, points(vare.mds, display="sites", col=colvec[Nlevel],pch=21, bg=colvec[Nlevel], cex=1.0))  

#to double check colors are actually aligning with treatments
head(with(bacplantdat, colvec[Nlevel]))
#to surround the points-further delineates
with(bacplantdat, ordihull(vare.mds, Sterility, col="black", lty=1, label=TRUE))
#lty means line type 1=solid increasing numbers=increasinly dashed
with(bacplantdat, legend ("bottomleft", legend=levels(Nlevel), bty="n", col=colvec, pch=21, pt.bg=colvec, cex=0.8))


# CJ: remaking figure **********
library(data.table)
df <- as.data.frame(vare.mds[["points"]])
setDT(df,keep.rownames = "SampleID.1")

df$Sterility <- bacplantdat[match(df$SampleID.1, bacplantdat$SampleID.1), "Sterility"]
df$Nitrogen <- bacplantdat[match(df$SampleID.1, bacplantdat$SampleID.1), "Nlevel"]

table(df$Sterility)
levels(df$Sterility) <- c("Unperturbed", "Perturbed")
myColors2 <- c("grey36", "gray74")
names(myColors2) <- levels(df$Nitrogen)
colScale <- scale_colour_manual(name = "Nitrogen \nLevel", values = myColors2)
ggplot(df, aes(x = MDS1, y = MDS2, colour = Nitrogen, shape = Sterility)) + theme_classic() + geom_point(size = 3) + colScale + scale_shape(name = "Microbial Community") + labs(x = "NMDS1", y = "NMDS2") + stat_ellipse( aes(x = MDS1, y = MDS2,  shape = Sterility, linetype = Sterility),type = "t", inherit.aes = F, show.legend = F)+ theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14))+ guides(colour = guide_legend(order = 1, override.aes = list(shape = 15)), shape = guide_legend(order = 2))
ggsave(file = "NMDSFigure5Sep2019.tiff", dpi = 1000)


#+ theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14)








######################################################################
#Alphadiversity
######################################################################
#Found a mistake in the previous draft OMG! We accidentally included mitochondria and chloroplasts in richness analysis
#the richness file below was written from the code above and should be correct (done after proper upstream steps)

#Need to remember to add sample ID name to the column with sample IDs or it won't work to join by SampleID
richness <- read.csv("richness_table.csv")

colnames(richness)[colnames(richness)=="X"] <- "SampleID.1"

str(richness)


str(metadata.E)
#to join dataframes I had to change row names to a column with a column header =ID



#Unrarefied is qualitatively the same so I will stick with rarefied for consistency

metarichness<-inner_join(metadata.E,richness, by = "SampleID.1")
metarichness$SampleID.1<-as.factor(metarichness$SampleID.1)
str(metarichness)






############################################################################################
##Assumptions of linear regression Observed
############################################################################################



M<-lmer (Observed ~ Nlevel*Sterility 
         + PresenceofMoldAug31 + (1|G5Block), data=metarichness)
plot(M, add.smooth = FALSE)
E<-resid(M)
E <- resid(M)
## Histogram of residual frequency-seems to evidence of non-normality
hist(E, xlab = "Residuals", main = "")


#### log might fit a little better-tradeoff skewness vs. homoegeneity
Mlog<-lmer (log(Observed) ~ Nlevel*Sterility 
            + PresenceofMoldAug31 + (1|G5Block), data=metarichness)
plot(Mlog, add.smooth = FALSE)
E<-resid(Mlog)
E <- resid(Mlog)
## Histogram of residual frequency-seems to evidence of non-normality
hist(E, xlab = "Residuals", main = "")

## QQplots of sample vs. theoretical quantiles also indicate non-normality
qqnorm(richness$Observed)
qqline(richness$Observed)


#####Full Model Observed richness-core root length
Mod.obs<-lmer (Observed ~ Nlevel*Sterility + PresenceofMoldAug31 + (1|G5Block), data=metarichness)
summary(Mod.obs)
anova(Mod.obs)
########## Doesn't seem to be much going on with observed richeness and the traits
#Type III Analysis of Variance Table with Satterthwaite's method
 #                    Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
#Nlevel               4116.3  4116.3     1 33.711  0.5789 0.45202  
#Sterility            7674.4  7674.4     1 33.881  1.0793 0.30620  
#PresenceofMoldAug31 24347.8 24347.8     1 35.994  3.4243 0.07247 .
#Nlevel:Sterility     1035.0  1035.0     1 33.871  0.1456 0.70519  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


############################################################################################
##Assumptions of linear regression Shannon
############################################################################################




M<-lmer (Shannon ~ Nlevel*Sterility 
         + PresenceofMoldAug31 + (1|G5Block), data=metarichness)
plot(M, add.smooth = FALSE)
E<-resid(M)
E <- resid(M)
## Histogram of residual frequency-seems to evidence of non-normality
hist(E, xlab = "Residuals", main = "")


#### log might fit a little better-tradeoff skewness vs. homoegeneity
Mlog<-lmer (log(Shannon) ~ Nlevel*Sterility 
            + PresenceofMoldAug31 + (1|G5Block), data=metarichness)
plot(Mlog, add.smooth = FALSE)
E<-resid(Mlog)
E <- resid(Mlog)
## Histogram of residual frequency-seems to evidence of non-normality
hist(E, xlab = "Residuals", main = "")

## QQplots of sample vs. theoretical quantiles also indicate non-normality
qqnorm(richness$Shannon)
qqline(richness$Shannon)



#####
Mod.shan<-lmer (Shannon ~ Nlevel*Sterility + PresenceofMoldAug31 + (1|G5Block), data=metarichness)
summary(Mod.shan)
anova(Mod.shan)
########## Doesn't seem to be much going on with shannon and traits but are N level effects on shannon
#Type III Analysis of Variance Table with Satterthwaite's method
#                      Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)  
#Nlevel              0.044663 0.044663     1 33.332  0.9149 0.34570  
#Sterility           0.142929 0.142929     1 33.490  2.9279 0.09631 .
#PresenceofMoldAug31 0.175737 0.175737     1 35.976  3.6000 0.06583 .
#Nlevel:Sterility    0.133835 0.133835     1 33.486  2.7416 0.10711  
#
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

################################
#Richness Figures
################################
#########
myColors <- c("grey36", "gray74")
names(myColors) <- levels(metarichness$Nlevel)
str(metarichness)
colScale <- scale_fill_manual(name = "Nitrogen \nLevel", values = myColors)
library(cowplot)
table(metarichness$Sterility)
levels(metarichness$Sterility) <- c("Unperturbed", "Perturbed")
theme_set(theme_cowplot())
##########

Obsfig <- ggplot(metarichness, aes(Sterility, Observed, fill = Nlevel)) + geom_boxplot() + colScale  + theme(legend.position = "none", axis.title.x = element_blank()) + labs(y = "Observed Richness") + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 15))
Shannfig <-ggplot(metarichness, aes(Sterility, Shannon, fill = Nlevel)) + geom_boxplot()  + colScale  + labs(y = "Shannon Diversity") + theme(axis.title.x = element_blank()) + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 15))
# put all together
Mylegend <- get_legend(Shannfig)
p1 <-plot_grid(Obsfig, Shannfig + theme(legend.position = "none"), ncol = 2, labels = "AUTO", align = "hv")
p2 <- plot_grid(p1, Mylegend, rel_widths = c(6,1))
title <- ggdraw() + draw_label("Microbial community", fontface = "bold", fontfamily = "", size = 18)
p3 <-plot_grid(p2,title, ncol=1, rel_heights=c(1, .07)) # rel_heights values control title margin
#ggdraw(add_sub(p2, "Label", vpadding=grid::unit(0,"lines"),y=-0.5, x=0.5, vjust=0))

save_plot("Rich12Aug2019.pdf", p3, ncol = 1, nrow = 2, base_width = 14)
