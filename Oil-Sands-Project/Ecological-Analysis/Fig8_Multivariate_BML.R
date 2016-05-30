##################################
# BML DATA MULTIVARIATE ANALYSIS #
##################################


library(vegan)
library(ggbiplot)


################
# LOADING DATA #
################

# Loading BML data
workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Data/BML2/usearch/output_BML2_95"
setwd(workingdirectory)

Species_BML=read.delim(file="Species_filtered.txt",row.names=1)
Class_BML=read.delim(file="Class_filtered.txt",row.names=1)

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

Class_sub=Class_BML[,c("SMPL0","SMPL1","SMPL2",
"SMPL3","SMPL4","SMPL5",
"SMPL6","SMPL7","SMPL8",
"SMPL9","SMPL10","SMPL11",
"SMPL12","SMPL13","SMPL14",
"SMPL15","SMPL16","SMPL17",
"SMPL18","SMPL22","SMPL26",
"SMPL27","SMPL31","SMPL35",
"SMPL37","SMPL41","SMPL45")]


Species_sub=Species_sub[-nrow(Species_sub),]
Class_sub=Class_sub[-nrow(Class_sub),]

# Environmental Matrix
env.matrix=data.frame(matrix(0,ncol=length(Species_sub),nrow=3),row.names=c("Platform","Depth","Season"))
names(env.matrix)=names(Species_sub)
env.matrix[1,]=rep(c(rep("Platform#1",3),rep("Platform#2",3),rep("Platform#3",3)),3)
env.matrix[2,]=rep(c(0,4,8),9)
env.matrix[3,]=c(rep("May",9),rep("September",9),rep("June",9))
env.matrix=t(env.matrix)

# Environmental Matrix with dummy variables
env.response=as.data.frame(cbind(env.matrix[,1],env.matrix[,1],env.matrix[,1],
                               env.matrix[,3],env.matrix[,3],env.matrix[,3],
                               env.matrix[,2],env.matrix[,1]))
names(env.response)=c("Platform#1","Platform#2","Platform#3","May","September","June","Depth","Platform")

env.response[,"Platform#1"]=gsub("Platform#1",1,env.response[,"Platform#1"])
env.response[,"Platform#1"]=gsub("Platform#.",0,env.response[,"Platform#1"])
env.response[,"Platform#2"]=gsub("Platform#2",1,env.response[,"Platform#2"])
env.response[,"Platform#2"]=gsub("Platform#.",0,env.response[,"Platform#2"])
env.response[,"Platform#3"]=gsub("Platform#3",1,env.response[,"Platform#3"])
env.response[,"Platform#3"]=gsub("Platform#.",0,env.response[,"Platform#3"])

env.response[,"May"]=gsub("May",1,env.response[,"May"])
env.response[,"May"]=gsub("September",0,env.response[,"May"])
env.response[,"May"]=gsub("June",0,env.response[,"May"])
env.response[,"June"]=gsub("May",0,env.response[,"June"])
env.response[,"June"]=gsub("September",0,env.response[,"June"])
env.response[,"June"]=gsub("June",1,env.response[,"June"])
env.response[,"September"]=gsub("May",0,env.response[,"September"])
env.response[,"September"]=gsub("September",1,env.response[,"September"])
env.response[,"September"]=gsub("June",0,env.response[,"September"])

env.response[,"Platform"]=rep(c(rep(1,3),rep(2,3),rep(3,3)),3)

env=c(rep("May",9),rep("September",9),rep("June",9))
names(env)=names(Species_sub)

#########################
# MULTIVARIATE ANALYSIS #
#########################

workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Report/Figures"
setwd(workingdirectory)

raremax.bml <- min(rowSums(t(Species_sub)))

# Class Matrix transformation
Class.r=rrarefy(t(Class_sub),raremax.bml)
# Species Hellinger transformation
Class.H=decostand(Class.r,method="hellinger")

# Matrix transformation
Species.r=rrarefy(t(Species_sub),raremax.bml)
# Hellinger transformation
Species.H=decostand(Species.r,method="hellinger")

########
# NMDS #
########

nmds=metaMDS(Species.H, dist="bray")

#env.matrixal Fit
ef <- envfit(nmds, as.data.frame(env.matrix), permu = 999)
ef

pdf(file = "FigureS4A.pdf", width=15, height=5)
par(mfcol=c(1,3),oma=c(0,0,4,0))
ordiplot(nmds,display = "sites")
ordihull(nmds,groups=as.character(env.matrix[,3]),draw="polygon",col="grey90",label=T,cex=1)
title("Season effect : p-value < 0.001")
ordiplot(nmds,display = "sites")
ordihull(nmds,groups=as.character(env.matrix[,1]),draw="polygon",col="grey90",label=T,cex=1)
title("Platform effect : p-value = 0.902")
ordiplot(nmds,display = "sites")
ordihull(nmds,groups=as.character(env.matrix[,2]),draw="polygon",col="grey90",label=T,cex=1)
title("Depth effect : p-value = 0.838")
mtext("A - Exploratory Analysis of BML : NMDS",outer=T,cex=1.5,padj=-1)
dev.off()

#######
# CCA #
#######

#PLOT
pdf(file = "Figure8A.pdf",width=6,height=6)
par(mfcol=c(1,1))
Platform3=env.response[,3]
colvec <- c("red2", "green4", "mediumblue")
cca=cca(formula = Species.H ~ May+June+Platform3, env.response)
anova(cca)
anova(cca, by = "margin")
plot(cca, display = "sites", type = "n")
ordisurf(cca, diversity(Species.H), add = TRUE)
ordiellipse(cca,groups=as.character(env),draw="polygon",label=F,cex=1)
with(env.response, points(cca, pch=16,col=colvec[as.numeric(env.response[,"Platform"])]))
text(cca, display="bp",col="blue",cex=1,labels=c("May","June","Platform3"))
with(env.response, legend("topleft", c("Platform#1","Platform#2","Platform#3"), pch = 16,col=colvec,cex=1))
title("CCA, Species abuncance ~ Month + Platform3")
mtext("Month : p-value < 0.001 | Platform3 : p-value = 0.03",cex=1,padj=-0.75)
dev.off()

#######
# PCA #
#######

#PCA
if (sum(Class.H[,ncol(Class.H)])==0){
Class.H=Class.H[,-ncol(Class.H)]}

Class.H

pdf(file = "Figure8B.pdf",width=6,height=6)
par(mfcol=c(1,1))
# PCA using prcomp function
pca.class <- prcomp(Class.H)

g <- ggbiplot(pca.class,
              groups = env.matrix[,"Season"], 
              ellipse = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
g <- g + labs(title="PCA of supergroups abundance matrix")
print(g)
dev.off()




