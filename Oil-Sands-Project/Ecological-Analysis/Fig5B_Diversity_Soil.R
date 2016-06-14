############################################
# OTU CLUSTERING DATA DESCRIPTION ANALYSIS #
############################################

library(vegan)

################
# LOADING DATA #
################

# Loading BML data
workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Data/Soil/D1_G2/usearch_output_subsample_0.95/"
setwd(workingdirectory)

Species_soil=read.delim(file="Species_filtered.txt",row.names=1)
Species_protists=read.delim(file="Species_protists.txt",row.names=1)

Species_soil=Species_soil[-nrow(Species_soil),]
Species_protists=Species_protists[-nrow(Species_protists),]

# Loading Fungi
Species_fungi=read.delim(file="Species_fungi.txt")

######################
# DIVERSITY BARPLOTS #
######################

workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Report/Figures"
setwd(workingdirectory)

raremax.soil <- min(rowSums(t(Species_soil)))
Species.r=rrarefy(t(Species_soil),raremax.soil)
shan = diversity (Species.r, index="shannon")
simp = diversity (Species.r, index="simpson")

shan
simp

allundis=c(shan["SMPL0"],4*simp["SMPL0"],shan["SMPL1"],4*simp["SMPL1"],shan["SMPL2"],4*simp["SMPL2"],
shan["SMPL3"],4*simp["SMPL3"],shan["SMPL4"],4*simp["SMPL4"],shan["SMPL5"],4*simp["SMPL5"])
allover=c(shan["SMPL12"],4*simp["SMPL12"],shan["SMPL13"],4*simp["SMPL13"],shan["SMPL14"],4*simp["SMPL14"],
shan["SMPL15"],4*simp["SMPL15"],shan["SMPL16"],4*simp["SMPL16"],shan["SMPL17"],4*simp["SMPL17"])
alltail=c(shan["SMPL6"],4*simp["SMPL6"],shan["SMPL7"],4*simp["SMPL7"],shan["SMPL8"],4*simp["SMPL8"],
shan["SMPL9"],4*simp["SMPL9"],shan["SMPL10"],4*simp["SMPL10"],shan["SMPL11"],4*simp["SMPL11"])

raremax.soil <- min(rowSums(t(Species_protists)))
Species.r=rrarefy(t(Species_protists),raremax.soil)
shan = diversity (Species.r, index="shannon")
simp = diversity (Species.r, index="simpson")

shan
simp

protundis=c(shan["SMPL0"],4*simp["SMPL0"],shan["SMPL1"],4*simp["SMPL1"],shan["SMPL2"],4*simp["SMPL2"],
shan["SMPL3"],4*simp["SMPL3"],shan["SMPL4"],4*simp["SMPL4"],shan["SMPL5"],4*simp["SMPL5"])
protover=c(shan["SMPL12"],4*simp["SMPL12"],shan["SMPL13"],4*simp["SMPL13"],shan["SMPL14"],4*simp["SMPL14"],
shan["SMPL15"],4*simp["SMPL15"],shan["SMPL16"],4*simp["SMPL16"],shan["SMPL17"],4*simp["SMPL17"])
prottail=c(shan["SMPL6"],4*simp["SMPL6"],shan["SMPL7"],4*simp["SMPL7"],shan["SMPL8"],4*simp["SMPL8"],
shan["SMPL9"],4*simp["SMPL9"],shan["SMPL10"],4*simp["SMPL10"],shan["SMPL11"],4*simp["SMPL11"])

raremax.soil <- min(rowSums(t(Species_fungi)))
Species.r=rrarefy(t(Species_fungi),raremax.soil)
shan = diversity (Species.r, index="shannon")
simp = diversity (Species.r, index="simpson")

shan
simp

fungundis=c(shan["SMPL0"],4*simp["SMPL0"],shan["SMPL1"],4*simp["SMPL1"],shan["SMPL2"],4*simp["SMPL2"],
shan["SMPL3"],4*simp["SMPL3"],shan["SMPL4"],4*simp["SMPL4"],shan["SMPL5"],4*simp["SMPL5"])
fungover=c(shan["SMPL12"],4*simp["SMPL12"],shan["SMPL13"],4*simp["SMPL13"],shan["SMPL14"],4*simp["SMPL14"],
shan["SMPL15"],4*simp["SMPL15"],shan["SMPL16"],4*simp["SMPL16"],shan["SMPL17"],4*simp["SMPL17"])
fungtail=c(shan["SMPL6"],4*simp["SMPL6"],shan["SMPL7"],4*simp["SMPL7"],shan["SMPL8"],4*simp["SMPL8"],
shan["SMPL9"],4*simp["SMPL9"],shan["SMPL10"],4*simp["SMPL10"],shan["SMPL11"],4*simp["SMPL11"])

colors=rep(c("steelblue","grey40"),3)

pdf(file = "Figure5B.pdf", width=8, height=7)
par(mfrow=c(3,3),oma=c(5,6,6,5),mar=c(1,0.25,2,0.25))

barplot(allundis,border=NA,names.arg=c(1,1,2,2,3,3,1,1,2,2,3,3),space=0.05,ylim=c(0,4),axes=F,xpd=T,col=colors)
axis(2, col="steelblue",lwd=2)
mtext(side = 2, text = "Shannon Index",padj=-3,col="steelblue",cex=0.8)
mtext(side=3,"UNDISTURBED",cex=1,padj=-1)
mtext(side = 2, text = "All taxa",padj=-5)
barplot(allover,border=NA,names.arg=c(1,1,2,2,3,3,1,1,2,2,3,3),space=0.05,ylim=c(0,4),axes=F,xpd=T,col=colors)
mtext(side=3,"OVERBURDEN",cex=1,padj=-1)
barplot(alltail,border=NA,names.arg=c(1,1,2,2,3,3,1,1,2,2,3,3),space=0.05,ylim=c(0,4),axes=F,xpd=T,col=colors)
mtext(side=3,"TAILINGS",cex=1,padj=-1)
axis(4, col="grey40",lwd=2,at=c(0,1,2,3,4),labels=c(0,0.25,0.5,0.75,1))
mtext(side = 4, text = "Simpson Index",padj=3,col="grey40",cex=0.8)

barplot(protundis,border=NA,names.arg=c(1,1,2,2,3,3,1,1,2,2,3,3),space=0.05,ylim=c(0,4),axes=F,xpd=T,col=colors)
axis(2, col="steelblue",lwd=2)
mtext(side = 2, text = "Shannon Index",padj=-3,col="steelblue",cex=0.8)
mtext(side = 2, text = "Non fungal supergroups",padj=-5)
barplot(protover,border=NA,names.arg=c(1,1,2,2,3,3,1,1,2,2,3,3),space=0.05,ylim=c(0,4),axes=F,xpd=T,col=colors)
barplot(prottail,border=NA,names.arg=c(1,1,2,2,3,3,1,1,2,2,3,3),space=0.05,ylim=c(0,4),axes=F,xpd=T,col=colors)
axis(4, col="grey40",lwd=2,at=c(0,1,2,3,4),labels=c(0,0.25,0.5,0.75,1))
mtext(side = 4, text = "Simpson Index",padj=3,col="grey40",cex=0.8)

barplot(fungundis,border=NA,names.arg=c(1,1,2,2,3,3,1,1,2,2,3,3),space=0.05,ylim=c(0,4),axes=F,xpd=T,col=colors)
axis(2, col="steelblue",lwd=2)
mtext(side = 2, text = "Shannon Index",padj=-3,col="steelblue",cex=0.8)
mtext(side = 2, text = "Fungal classes",padj=-5)
mtext(side = 1, text = "Mineral     /     Organic",padj=3)
barplot(fungover,border=NA,names.arg=c(1,1,2,2,3,3,1,1,2,2,3,3),space=0.05,ylim=c(0,4),axes=F,xpd=T,col=colors)
mtext(side = 1, text = "Mineral     /     Organic",padj=3)
barplot(fungtail,border=NA,names.arg=c(1,1,2,2,3,3,1,1,2,2,3,3),space=0.05,ylim=c(0,4),axes=F,xpd=T,col=colors)
axis(4, col="grey40",lwd=2,at=c(0,1,2,3,4),labels=c(0,0.25,0.5,0.75,1))
mtext(side = 4, text = "Simpson Index",padj=3,col="grey40",cex=0.8)
mtext(side = 1, text = "Mineral     /     Organic",padj=3)

mtext("B - Diversity barplot, Soil Dataset",outer=T,cex=1.5,padj=-2)

dev.off()

