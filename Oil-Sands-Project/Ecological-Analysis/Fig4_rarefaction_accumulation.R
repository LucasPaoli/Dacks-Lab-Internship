############################################
# OTU CLUSTERING DATA DESCRIPTION ANALYSIS #
############################################

library(vegan)

################
# LOADING DATA #
################

# Loading BML data
workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Data/BML2/usearch/output_BML2_95"
setwd(workingdirectory)

Class_BML=read.delim(file="Class_filtered.txt",row.names=1)
Species_BML=read.delim(file="Species_filtered.txt",row.names=1)
OTU_BML=read.delim(file="otutab.txt",row.names=1)

# Reordering the samples
Class_BML=Class_BML[,c("SMPL0","SMPL1","SMPL2","SMPL3","SMPL4","SMPL5","SMPL6",
                       "SMPL7","SMPL8","SMPL9","SMPL10","SMPL11","SMPL12","SMPL13",
                       "SMPL14","SMPL15","SMPL16","SMPL17","SMPL18","SMPL19","SMPL20",
                       "SMPL21","SMPL22","SMPL23","SMPL24","SMPL25","SMPL26","SMPL27",
                       "SMPL28","SMPL29","SMPL30","SMPL31","SMPL32","SMPL33","SMPL34",
                       "SMPL35","SMPL36","SMPL37","SMPL38","SMPL39","SMPL40","SMPL41",
                       "SMPL42","SMPL43","SMPL44","SMPL45")]
Species_BML=Species_BML[,c("SMPL0","SMPL1","SMPL2","SMPL3","SMPL4","SMPL5","SMPL6",
                           "SMPL7","SMPL8","SMPL9","SMPL10","SMPL11","SMPL12","SMPL13",
                           "SMPL14","SMPL15","SMPL16","SMPL17","SMPL18","SMPL19","SMPL20",
                           "SMPL21","SMPL22","SMPL23","SMPL24","SMPL25","SMPL26","SMPL27",
                           "SMPL28","SMPL29","SMPL30","SMPL31","SMPL32","SMPL33","SMPL34",
                           "SMPL35","SMPL36","SMPL37","SMPL38","SMPL39","SMPL40","SMPL41",
                           "SMPL42","SMPL43","SMPL44","SMPL45")]
OTU_BML=OTU_BML[,c("SMPL0","SMPL1","SMPL2","SMPL3","SMPL4","SMPL5","SMPL6",
                   "SMPL7","SMPL8","SMPL9","SMPL10","SMPL11","SMPL12","SMPL13",
                   "SMPL14","SMPL15","SMPL16","SMPL17","SMPL18","SMPL19","SMPL20",
                   "SMPL21","SMPL22","SMPL23","SMPL24","SMPL25","SMPL26","SMPL27",
                   "SMPL28","SMPL29","SMPL30","SMPL31","SMPL32","SMPL33","SMPL34",
                   "SMPL35","SMPL36","SMPL37","SMPL38","SMPL39","SMPL40","SMPL41",
                   "SMPL42","SMPL43","SMPL44","SMPL45")]


# Taking relevant samples (Excluding point that are not replicated)
# Ex : June depths 1,2,3,5,6,7,9,10m

OTU_sub=OTU_BML[,c("SMPL0","SMPL1","SMPL2",
                   "SMPL3","SMPL4","SMPL5",
                   "SMPL6","SMPL7","SMPL8",
                   "SMPL9","SMPL10","SMPL11",
                   "SMPL12","SMPL13","SMPL14",
                   "SMPL15","SMPL16","SMPL17",
                   "SMPL18","SMPL22","SMPL26",
                   "SMPL27","SMPL31","SMPL35",
                   "SMPL37","SMPL41","SMPL45")]

Class_sub=Class_BML[,c("SMPL0","SMPL1","SMPL2",
                       "SMPL3","SMPL4","SMPL5",
                       "SMPL6","SMPL7","SMPL8",
                       "SMPL9","SMPL10","SMPL11",
                       "SMPL12","SMPL13","SMPL14",
                       "SMPL15","SMPL16","SMPL17",
                       "SMPL18","SMPL22","SMPL26",
                       "SMPL27","SMPL31","SMPL35",
                       "SMPL37","SMPL41","SMPL45")]

Species_sub=Species_BML[,c("SMPL0","SMPL1","SMPL2",
                           "SMPL3","SMPL4","SMPL5",
                           "SMPL6","SMPL7","SMPL8",
                           "SMPL9","SMPL10","SMPL11",
                           "SMPL12","SMPL13","SMPL14",
                           "SMPL15","SMPL16","SMPL17",
                           "SMPL18","SMPL22","SMPL26",
                           "SMPL27","SMPL31","SMPL35",
                           "SMPL37","SMPL41","SMPL45")]

# Getting rid of the unknown sequences
Species_sub=Species_sub[-nrow(Species_sub),]
Class_sub=Class_sub[-nrow(Class_sub),]
OTU_sub=OTU_sub[-nrow(OTU_sub),]

Species.season=rbind(colSums(t(Species_sub)[1:9,]),
                     colSums(t(Species_sub)[10:18,]),
                     colSums(t(Species_sub)[19:27,]))

Class.season=rbind(colSums(t(Class_sub)[1:9,]),
                   colSums(t(Class_sub)[10:18,]),
                   colSums(t(Class_sub)[19:27,]))

OTU.season=rbind(colSums(t(OTU_sub)[1:9,]),
                 colSums(t(OTU_sub)[10:18,]),
                 colSums(t(OTU_sub)[19:27,]))


# Loading Soil data
workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Data/Soil/D1_G2/usearch_output_subsample_0.95/"
setwd(workingdirectory)

Class_soil=read.delim(file="Class_filtered.txt",row.names=1)
Species_soil=read.delim(file="Species_filtered.txt",row.names=1)
OTU_soil=read.delim(file="otutab.txt",row.names=1)

# Get rid of Unknown
Species_soil=Species_soil[-nrow(Species_soil),]
Class_soil=Class_soil[-nrow(Class_soil),]
OTU_soil=OTU_soil[-nrow(OTU_soil),]

Species.reclamation=rbind(colSums(t(Species_soil)[1:6,]),
                          colSums(t(Species_soil)[7:12,]),
                          colSums(t(Species_soil)[13:18,]))

Class.reclamation=rbind(colSums(t(Class_soil)[1:6,]),
                        colSums(t(Class_soil)[7:12,]),
                        colSums(t(Class_soil)[13:18,]))

OTU.reclamation=rbind(colSums(t(OTU_soil)[1:6,]),
                      colSums(t(OTU_soil)[7:12,]),
                      colSums(t(OTU_soil)[13:18,]))



#######################
# ACCUMULATION CURVES #
#######################

workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Report/Figures"
setwd(workingdirectory)

pdf(file = "Figure4.pdf", width=10, height=5)
par(mfrow=c(1,2),oma=c(0,0,2.5,0))

#######
# BML #
#######

raremax.bml <- min(rowSums(t(Species_sub)))
row.names(Species.season)=c("May","September","June")

rarecurve(Species.season,sample=9*raremax.bml,xlab="Sampling size (Number of sequences)")
title("A - Rarefaction curve and effect of rarefaction\non the lake for each season")

########
# SOIL #
########

raremax.soil <- min(rowSums(t(Species_soil)))
row.names(Species.reclamation)=c("Undisturbed","Tailings","Overburden")

rarecurve(Species.reclamation,sample=6*raremax.soil,xlab="Sampling size (Number of sequences)")
title("B - Rarefaction curve and effect of rarefaction\non the soil dataset for each treatment")

mtext("Rarefaction curves",outer=T,cex=1.5)

dev.off()



