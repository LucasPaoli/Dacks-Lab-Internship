############################################
# OTU CLUSTERING DATA DESCRIPTION ANALYSIS #
############################################

################
# LOADING DATA #
################

# Loading data
workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Data/Soil/D1_G2/usearch_output_subsample_0.95/"
setwd(workingdirectory)

Class_soil=read.delim(file="Class_filtered.txt",row.names=1)
Class_protists=read.delim(file="Class_protists.txt",row.names=1)

# Reordering the Classes
Class_soil=Class_soil[c("Eukaryota_XX",
                        "Opisthokonta_X",
                        "Choanoflagellida",
                        "Fungi",
                        "Mesomycetozoa",
                        "Hilomonadea",
                        "Apusomonadidae",
                        "Amoebozoa_X",
                        "Lobosa",
                        "Conosa",
                        "Chlorophyta",
                        "Rhodophyta",
                        "Discoba",
                        "Malawimonadidae",
                        "Centroheliozoa",
                        "Telonemia",
                        "Haptophyta",
                        "Cryptophyta",
                        "Picobiliphyta",
                        "Cercozoa",
                        "Foraminifera",
                        "Radiolaria",
                        "Stramenopiles_X",
                        "Ochrophyta",
                        "Dinophyta",
                        "Apicomplexa",
                        "Perkinsea",
                        "Ciliophora"),]

# Reordering the Classes
Class_protists=Class_protists[c("Eukaryota_XX",
                        "Opisthokonta_X",
                        "Choanoflagellida",
                        "Mesomycetozoa",
                        "Hilomonadea",
                        "Apusomonadidae",
                        "Amoebozoa_X",
                        "Lobosa",
                        "Conosa",
                        "Chlorophyta",
                        "Rhodophyta",
                        "Discoba",
                        "Malawimonadidae",
                        "Centroheliozoa",
                        "Telonemia",
                        "Haptophyta",
                        "Cryptophyta",
                        "Picobiliphyta",
                        "Cercozoa",
                        "Foraminifera",
                        "Radiolaria",
                        "Stramenopiles_X",
                        "Ochrophyta",
                        "Dinophyta",
                        "Apicomplexa",
                        "Perkinsea",
                        "Ciliophora"),]

# Loading Fungi
Species_fungi=read.delim(file="Species_fungi_.txt")

Species_fungi[,1]=gsub("Eukaryota_Opisthokonta_Fungi_","",Species_fungi[,1])
Species_fungi[,1]=gsub("_.*","",Species_fungi[,1])

Fungi=as.data.frame(matrix(0,length(unique(Species_fungi[,1])),(length(Species_fungi)-1)))
row.names(Fungi)=unique(Species_fungi[,1])
names(Fungi)=names(Species_fungi[,2:length(Species_fungi)])

for (i in row.names(Fungi)){
  A=subset(Species_fungi,names==i)[,2:length(Species_fungi)]
  Fungi[i,]=colSums(A)
}

#######################
# ABUNDANCES BARPLOTS #
#######################

workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Report/Figures"
setwd(workingdirectory)

Class_soil_N=sweep(Class_soil,2,colSums(Class_soil),"/")
Fungi_N = sweep(Fungi,2,colSums(Fungi),"/")
Class_protists_N=sweep(Class_protists,2,colSums(Class_protists),"/")

pdf(file = "Figure6B.pdf", width=9, height=9)
par(mfrow=c(3,4),oma=c(5,2,6,0),mar=c(0.5,0.25,0.5,0.25))

barplot(border=NA,names.arg=rep(c(1,2,3),2),space=0.05,axes=F,
        as.matrix(Class_soil_N[,c("SMPL0","SMPL1","SMPL2","SMPL3","SMPL4","SMPL5")]),
        col=rainbow(length(row.names(Class_soil_N))))
mtext(side=3,"UNDISTURBED",cex=1,padj=-0.5)
mtext(side = 2, text = "All Taxa",padj=-1)
barplot(border=NA,names.arg=rep(c(1,2,3),2),space=0.05,axes=F,
        as.matrix(Class_soil_N[,c("SMPL12","SMPL13","SMPL15","SMPL16","SMPL14","SMPL17")]),
        col=rainbow(length(row.names(Class_soil_N))))
mtext(side=3,"OVERBURDEN",cex=1,padj=-0.5)
barplot(border=NA,names.arg=rep(c(1,2,3),2),space=0.05,axes=F,
        as.matrix(Class_soil_N[,c("SMPL6","SMPL7","SMPL9","SMPL10","SMPL8","SMPL11")]),
        col=rainbow(length(row.names(Class_soil_N))))
mtext(side=3,"TAILINGS",cex=1,padj=-0.5)

frame()

barplot(border=NA,names.arg=rep(c(1,2,3),2),space=0.05,axes=F,
        as.matrix(Class_protists_N[,c("SMPL0","SMPL1","SMPL2","SMPL3","SMPL4","SMPL5")]),
        col=rainbow(length(row.names(Class_protists_N))))
mtext(side = 2, text = "Non fungal supergroups",padj=-1)
barplot(border=NA,names.arg=rep(c(1,2,3),2),space=0.05,axes=F,
        as.matrix(Class_protists_N[,c("SMPL12","SMPL13","SMPL15","SMPL16","SMPL14","SMPL17")]),
        col=rainbow(length(row.names(Class_protists_N))))
barplot(border=NA,names.arg=rep(c(1,2,3),2),space=0.05,axes=F,
        as.matrix(Class_protists_N[,c("SMPL6","SMPL7","SMPL9","SMPL10","SMPL8","SMPL11")]),
        col=rainbow(length(row.names(Class_protists_N))))


frame()

barplot(border=NA,names.arg=rep(c(1,2,3),2),space=0.05,axes=F,
        as.matrix(Fungi_N[,c("SMPL0","SMPL1","SMPL2","SMPL3","SMPL4","SMPL5")]),
        col=rainbow(length(row.names(Fungi_N))))
mtext(side = 2, text = "Fungal classes",padj=-1)
mtext(side = 1, text = "Mineral     /     Organic",padj=3)
barplot(border=NA,names.arg=rep(c(1,2,3),2),space=0.05,axes=F,
        as.matrix(Fungi_N[,c("SMPL12","SMPL13","SMPL15","SMPL16","SMPL14","SMPL17")]),
        col=rainbow(length(row.names(Fungi_N))))
mtext(side = 1, text = "Mineral     /     Organic",padj=3)
barplot(border=NA,names.arg=rep(c(1,2,3),2),space=0.05,axes=F,
        as.matrix(Fungi_N[,c("SMPL6","SMPL7","SMPL9","SMPL10","SMPL8","SMPL11")]),
        col=rainbow(length(row.names(Fungi_N))))
mtext(side = 1, text = "Mineral     /     Organic",padj=3)

frame()

par(mfrow=c(1,4),oma=c(5,2,6,0),mar=c(0.5,0.25,0.5,0.25))
legend("topleft",legend = rev(row.names(Class_soil_N)),
       fill = rev(rainbow(length(row.names(Class_soil_N))))
       , ncol=1,
       cex = 1.25,
       border=NA,
       bty="n")
legend("bottomleft",legend = rev(row.names(Fungi_N)),
       fill = rev(rainbow(length(row.names(Fungi_N))))
       , ncol=1,
       cex = 1.4,
       border=NA,
       bty="n")

mtext("B - Abundance barplot, Soil Dataset",outer=T,cex=1.5,padj=-2)

dev.off()


workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Data/Soil/D1_G2/usearch_output_subsample_0.95/"
setwd(workingdirectory)

Class_soil=read.delim(file="Species_filtered.txt",row.names=1)


Class_sub_N=rowSums(Class_soil)/sum(rowSums(Class_soil))

sort(Class_sub_N)