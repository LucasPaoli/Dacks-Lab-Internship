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

pdf(file = "FigureS2.pdf", width=16, height=10)
par(mfrow=c(2,3),mar=c(5,5,5,3),oma=c(0,4,4,0))

#######
# BML #
#######

raremax.bml <- min(rowSums(t(Species_sub)))
row.names(Species.season)=c("May","September","June")

plot(specaccum(t(Species_sub),method="random",permutations=100),
     ci.type="poly", col="steelblue", lwd=2, ci.lty=0, ci.col="lightblue",
     ylab="Number of Species",
     xlab="Number of sampling sites",cex.lab=1.7,cex.axis=1.5)

title("A - Accumulation curve for the whole lake\nover the three time points",cex.main=2)

mtext(side=2,padj=-5,"BASE MINE LAKE",cex=1.5)

rarecurve(t(Species_sub), step = 20, sample = raremax.bml, col = "steelblue",xlab="Sampling size (Number of sequences)",cex=1.1,cex.lab=1.7,cex.axis=1.5)
title("B - Rarefaction curve and effect of\nrarefaction on each sample",cex.main=2)

M = Species_sub
M = sweep(M,2,colSums(M),"/")
M[M==0]=NA

M1 = t(rrarefy(t(Species_sub),raremax.bml))
M1 = sweep(M1,2,colSums(M1),"/")
M1[M1==0]=NA

plot(sort(M[,1],decreasing=T),type="n",ylim=c(0,0.75),xlim=c(1,50),ylab="Relative abundance",
     main="C - Decreasing relative abundances\nof species per sample",
     xlab="Species index",log="x",
     cex.lab=1.7,cex.axis=1.5,cex.main=2)
for (i in 1:ncol(M)){lines(sort(M[,i],decreasing=T),col="black",lwd=2)}
for (i in 1:ncol(M1)){lines(sort(M1[,i],decreasing=T),col="steelblue",lwd=2)}

legend("topright", # places a legend at the appropriate place 
       c("Non rarefied","Rarefied"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),
       col=c("black","steelblue"),
       cex=1.2) # gives the legend lines the correct color and width


########
# SOIL #
########

raremax.soil <- min(rowSums(t(Species_soil)))
row.names(Species.reclamation)=c("Undisturbed","Tailings","Overburden")

plot(specaccum(t(Species_soil),method="random",permutations=100),
     ci.type="poly", col="steelblue", lwd=2, ci.lty=0, ci.col="lightblue",
     ylab="Number of Species",
     xlab="Number of sampling sites",cex.lab=1.7,cex.axis=1.5)

title("D - Accumulation curve for\nthe soils dataset",cex.main=2)

mtext(side=2,padj=-5,"SOIL DATASET",cex=1.5)

rarecurve(t(Species_soil), step = 20, sample = raremax.soil, col = "steelblue",xlab="Sampling size (Number of sequences)",cex=1.1,cex.lab=1.7,cex.axis=1.5)
title("E - Rarefaction curve and effect of\nrarefaction on each sample",cex.main=2)

M = Species_soil
M = sweep(M,2,colSums(M),"/")
M[M==0]=NA

M1 = t(rrarefy(t(Species_soil),raremax.soil))
M1 = sweep(M1,2,colSums(M1),"/")
M1[M1==0]=NA

plot(sort(M[,1],decreasing=T),type="n",ylim=c(0,0.5),xlim=c(1,350),ylab="Relative abundance",
     main="F - Decreasing relative abundances\nof species per sample",
     xlab="Species index",
     log="x",
     cex.lab=1.7,cex.axis=1.5,cex.main=2)
for (i in 1:ncol(M)){lines(sort(M[,i],decreasing=T),col="black",lwd=2)}
for (i in 1:ncol(M1)){lines(sort(M1[,i],decreasing=T),col="steelblue",lwd=2)}

legend("topright", # places a legend at the appropriate place 
       c("Non rarefied","Rarefied"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),
       col=c("black","steelblue"),
       cex=1.2) # gives the legend lines the correct color and width


mtext("Rarefaction effect and accumulation curves",outer=T,cex=2,padj=-0.5)

dev.off()



