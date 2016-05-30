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

Species_BML=read.delim(file="Species_filtered.txt",row.names=1)

# Reordering the samples
Species_BML=Species_BML[,c("SMPL0","SMPL1","SMPL2","SMPL3","SMPL4","SMPL5","SMPL6",
                       "SMPL7","SMPL8","SMPL9","SMPL10","SMPL11","SMPL12","SMPL13",
                       "SMPL14","SMPL15","SMPL16","SMPL17","SMPL18","SMPL19","SMPL20",
                       "SMPL21","SMPL22","SMPL23","SMPL24","SMPL25","SMPL26","SMPL27",
                       "SMPL28","SMPL29","SMPL30","SMPL31","SMPL32","SMPL33","SMPL34",
                       "SMPL35","SMPL36","SMPL37","SMPL38","SMPL39","SMPL40","SMPL41",
                       "SMPL42","SMPL43","SMPL44","SMPL45")]

# Taking relevant samples (Excluding point that are not replicated)
# Ex : June depths 1,2,3,5,6,7,9,1m0

Species_sub=Species_BML[,c("SMPL0","SMPL1","SMPL2",
                       "SMPL3","SMPL4","SMPL5",
                       "SMPL6","SMPL7","SMPL8",
                       "SMPL9","SMPL10","SMPL11",
                       "SMPL12","SMPL13","SMPL14",
                       "SMPL15","SMPL16","SMPL17",
                       "SMPL18","SMPL22","SMPL26",
                       "SMPL27","SMPL31","SMPL35",
                       "SMPL37","SMPL41","SMPL45")]


Species_sub=Species_sub[-nrow(Species_sub),]

######################
# DIVERSITY BARPLOTS #
######################

workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Report/Figures"
setwd(workingdirectory)

raremax.bml <- min(rowSums(t(Species_sub)))
Species.r=rrarefy(t(Species_sub),raremax.bml)
shan = diversity (Species.r, index="shannon")
simp = diversity (Species.r, index="simpson")

shan
simp

m0May=c(shan["SMPL6"],2*simp["SMPL6"],shan["SMPL0"],2*simp["SMPL0"],shan["SMPL3"],2*simp["SMPL3"])
m4May=c(shan["SMPL37"],2*simp["SMPL37"],shan["SMPL18"],2*simp["SMPL18"],shan["SMPL27"],2*simp["SMPL27"])
m8May=c(shan["SMPL15"],2*simp["SMPL15"],shan["SMPL9"],2*simp["SMPL9"],shan["SMPL12"],2*simp["SMPL12"])

m0June=c(shan["SMPL7"],2*simp["SMPL7"],shan["SMPL1"],2*simp["SMPL1"],shan["SMPL4"],2*simp["SMPL4"])
m4June=c(shan["SMPL41"],2*simp["SMPL41"],shan["SMPL22"],2*simp["SMPL22"],shan["SMPL31"],2*simp["SMPL31"])
m8June=c(shan["SMPL16"],2*simp["SMPL16"],shan["SMPL10"],2*simp["SMPL10"],shan["SMPL13"],2*simp["SMPL13"])

m0Sept=c(shan["SMPL8"],2*simp["SMPL8"],shan["SMPL2"],2*simp["SMPL2"],shan["SMPL5"],2*simp["SMPL5"])
m4Sept=c(shan["SMPL45"],2*simp["SMPL45"],shan["SMPL26"],2*simp["SMPL26"],shan["SMPL35"],2*simp["SMPL35"])
m8Sept=c(shan["SMPL17"],2*simp["SMPL17"],shan["SMPL11"],2*simp["SMPL11"],shan["SMPL3"],2*simp["SMPL14"])

colors=rep(c("steelblue","grey40"),3)

pdf(file = "Figure7A.pdf", width=8, height=7)
par(mfrow=c(3,3),oma=c(5,6,6,5),mar=c(1,0.25,2,0.25))

barplot(m0May,border=NA,names.arg=c(3,3,1,1,2,2),space=0.05,ylim=c(0,2),axes=F,xpd=T,col=colors)
axis(2, col="steelblue",lwd=2)
mtext(side = 2, text = "Shannon Index",padj=-3,col="steelblue",cex=0.8)
mtext(side=3,"MAY",cex=1,padj=-1)
mtext(side = 2, text = "Depth : 0m",padj=-5)
barplot(m4May,border=NA,names.arg=c(3,3,1,1,2,2),space=0.05,ylim=c(0,2),axes=F,xpd=T,col=colors)
mtext(side=3,"JUNE",cex=1,padj=-1)
barplot(m8May,border=NA,names.arg=c(3,3,1,1,2,2),space=0.05,ylim=c(0,2),axes=F,xpd=T,col=colors)
mtext(side=3,"SEPTEMBER",cex=1,padj=-1)
axis(4, col="grey40",lwd=2,at=c(0,0.5,1,1.5,2),labels=c(0,0.25,0.5,0.75,1))
mtext(side = 4, text = "Simpson Index",padj=3,col="grey40",cex=0.8)

barplot(m0June,border=NA,names.arg=c(3,3,1,1,2,2),space=0.05,ylim=c(0,2),axes=F,xpd=T,col=colors)
axis(2, col="steelblue",lwd=2)
mtext(side = 2, text = "Shannon Index",padj=-3,col="steelblue",cex=0.8)
mtext(side = 2, text = "Depth : 4m",padj=-5)
barplot(m4June,border=NA,names.arg=c(3,3,1,1,2,2),space=0.05,ylim=c(0,2),axes=F,xpd=T,col=colors)
barplot(m8June,border=NA,names.arg=c(3,3,1,1,2,2),space=0.05,ylim=c(0,2),axes=F,xpd=T,col=colors)
axis(4, col="grey40",lwd=2,at=c(0,0.5,1,1.5,2),labels=c(0,0.25,0.5,0.75,1))
mtext(side = 4, text = "Simpson Index",padj=3,col="grey40",cex=0.8)

barplot(m0Sept,border=NA,names.arg=c(3,3,1,1,2,2),space=0.05,ylim=c(0,2),axes=F,xpd=T,col=colors)
axis(2, col="steelblue",lwd=2)
mtext(side = 2, text = "Shannon Index",padj=-3,col="steelblue",cex=0.8)
mtext(side = 2, text = "Depth : 8m",padj=-5)
mtext(side = 1, text = "Platforms",padj=3)
barplot(m4Sept,border=NA,names.arg=c(3,3,1,1,2,2),space=0.05,ylim=c(0,2),axes=F,xpd=T,col=colors)
mtext(side = 1, text = "Platforms",padj=3)
barplot(m8Sept,border=NA,names.arg=c(3,3,1,1,2,2),space=0.05,ylim=c(0,2),axes=F,xpd=T,col=colors)
axis(4, col="grey40",lwd=2,at=c(0,0.5,1,1.5,2),labels=c(0,0.25,0.5,0.75,1))
mtext(side = 4, text = "Simpson Index",padj=3,col="grey40",cex=0.8)
mtext(side = 1, text = "Platforms",padj=3)

mtext("A - Diversity barplot, Base Mine Lake",outer=T,cex=1.5,padj=-2)

dev.off()

