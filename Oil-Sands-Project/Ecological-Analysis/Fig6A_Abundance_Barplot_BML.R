############################################
# OTU CLUSTERING DATA DESCRIPTION ANALYSIS #
############################################

################
# LOADING DATA #
################

# Loading BML data
workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Data/BML2/usearch/output_BML2_95"
setwd(workingdirectory)

Class_BML=read.delim(file="Class_filtered.txt",row.names=1)

# Reordering the samples
Class_BML=Class_BML[,c("SMPL0","SMPL1","SMPL2","SMPL3","SMPL4","SMPL5","SMPL6",
                       "SMPL7","SMPL8","SMPL9","SMPL10","SMPL11","SMPL12","SMPL13",
                       "SMPL14","SMPL15","SMPL16","SMPL17","SMPL18","SMPL19","SMPL20",
                       "SMPL21","SMPL22","SMPL23","SMPL24","SMPL25","SMPL26","SMPL27",
                       "SMPL28","SMPL29","SMPL30","SMPL31","SMPL32","SMPL33","SMPL34",
                       "SMPL35","SMPL36","SMPL37","SMPL38","SMPL39","SMPL40","SMPL41",
                       "SMPL42","SMPL43","SMPL44","SMPL45")]

# Taking relevant samples (Excluding point that are not replicated)
# Ex : June depths 1,2,3,5,6,7,9,10m

Class_sub=Class_BML[,c("SMPL0","SMPL1","SMPL2",
                       "SMPL3","SMPL4","SMPL5",
                       "SMPL6","SMPL7","SMPL8",
                       "SMPL9","SMPL10","SMPL11",
                       "SMPL12","SMPL13","SMPL14",
                       "SMPL15","SMPL16","SMPL17",
                       "SMPL18","SMPL22","SMPL26",
                       "SMPL27","SMPL31","SMPL35",
                       "SMPL37","SMPL41","SMPL45")]

Class_sub=Class_sub[c("Choanoflagellida",
                      "Mesomycetozoa",
                      "Fungi",
                      "Hilomonadea",
                      "Apusomonadidae",
                      "Amoebozoa_X",
                      "Lobosa",
                      "Conosa",
                      "Chlorophyta",
                      "Discoba",
                      "Centroheliozoa",
                      "Telonemia",
                      "Cercozoa",
                      "Stramenopiles_X",
                      "Ochrophyta",
                      "Dinophyta",
                      "Apicomplexa",
                      "Perkinsea",
                      "Ciliophora"),]


#######################
# ABUNDANCES BARPLOTS #
#######################

workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Report/Figures"
setwd(workingdirectory)

Class_sub_N=sweep(Class_sub,2,colSums(Class_sub),"/")

pdf(file = "Figure6A.pdf", width=8, height=6)
par(mfrow=c(3,4),oma=c(5,2,6,0),mar=c(0.5,0.25,0.5,0.25))

barplot(border=NA,names.arg=c(3,1,2),space=0.05,axes=F,as.matrix(Class_sub_N[,c("SMPL6","SMPL0","SMPL3")]),col=rainbow(length(row.names(Class_sub_N))))
mtext(side=3,"MAY",cex=1,padj=-0.5)
mtext(side = 2, text = "Depth : 0m",padj=-1)
barplot(border=NA,names.arg=c(3,1,2),space=0.05,axes=F,as.matrix(Class_sub_N[,c("SMPL37","SMPL18","SMPL27")]),col=rainbow(length(row.names(Class_sub_N))))
mtext(side=3,"JUNE",cex=1,padj=-0.5)
barplot(border=NA,names.arg=c(3,1,2),space=0.05,axes=F,as.matrix(Class_sub_N[,c("SMPL15","SMPL9","SMPL12")]),col=rainbow(length(row.names(Class_sub_N))))
mtext(side=3,"SEPTEMBER",cex=1,padj=-0.5)

frame()

barplot(border=NA,names.arg=c(3,1,2),space=0.05,axes=F,as.matrix(Class_sub_N[,c("SMPL7","SMPL1","SMPL4")]),col=rainbow(length(row.names(Class_sub_N))))
mtext(side = 2, text = "Depth : 4m",padj=-1)
barplot(border=NA,names.arg=c(3,1,2),space=0.05,axes=F,as.matrix(Class_sub_N[,c("SMPL41","SMPL22","SMPL31")]),col=rainbow(length(row.names(Class_sub_N))))
barplot(border=NA,names.arg=c(3,1,2),space=0.05,axes=F,as.matrix(Class_sub_N[,c("SMPL16","SMPL10","SMPL13")]),col=rainbow(length(row.names(Class_sub_N))))

frame()


barplot(border=NA,names.arg=c(3,1,2),space=0.05,axes=F,as.matrix(Class_sub_N[,c("SMPL8","SMPL2","SMPL5")]),col=rainbow(length(row.names(Class_sub_N))))
mtext(side = 2, text = "Depth : 8m",padj=-1)
mtext(side = 1, text = "Platforms",padj=3)
barplot(border=NA,names.arg=c(3,1,2),space=0.05,axes=F,as.matrix(Class_sub_N[,c("SMPL45","SMPL26","SMPL35")]),col=rainbow(length(row.names(Class_sub_N))))
mtext(side = 1, text = "Platforms",padj=3)
barplot(border=NA,names.arg=c(3,1,2),space=0.05,axes=F,as.matrix(Class_sub_N[,c("SMPL17","SMPL11","SMPL14")]),col=rainbow(length(row.names(Class_sub_N))))
mtext(side = 1, text = "Platforms",padj=3)

frame()

par(mfrow=c(1,4),oma=c(5,2,6,0),mar=c(0.5,0.25,0.5,0.25))
legend("left",legend = rev(row.names(Class_sub_N)),
fill = rev(rainbow(length(row.names(Class_sub_N))))
, ncol=1,
cex = 1.5,
border=NA,
bty="n")

mtext("A - Abundance barplot, Base Mine Lake",outer=T,cex=1.5,padj=-2)

dev.off()

workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Data/BML2/usearch/output_BML2_95"
setwd(workingdirectory)

# Displaying the taxonomic assignment by relative abundance
Class_BML=read.delim(file="Species_filtered.txt",row.names=1)

# Taking relevant samples (Excluding point that are not replicated)
# Ex : June depths 1,2,3,5,6,7,9,10m

Class_sub=Class_BML[,c("SMPL0","SMPL1","SMPL2",
"SMPL3","SMPL4","SMPL5",
"SMPL6","SMPL7","SMPL8",
"SMPL9","SMPL10","SMPL11",
"SMPL12","SMPL13","SMPL14",
"SMPL15","SMPL16","SMPL17",
"SMPL18","SMPL22","SMPL26",
"SMPL27","SMPL31","SMPL35",
"SMPL37","SMPL41","SMPL45")]


Class_sub_N=rowSums(Class_sub)/sum(rowSums(Class_sub))

sort(Class_sub_N)

