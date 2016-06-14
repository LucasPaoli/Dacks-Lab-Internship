###################################
# SOIL DATA MULTIVARIATE ANALYSIS #
###################################


library(vegan)
library(ggbiplot)


################
# LOADING DATA #
################

# Loading data
workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Data/Soil/D1_G2/usearch_output_subsample_0.95/"
setwd(workingdirectory)

Species_soil=read.delim(file="Species_filtered.txt",row.names=1)
Class_soil=read.delim(file="Class_filtered.txt",row.names=1)

Species_soil=Species_soil[-nrow(Species_soil),]
Class_soil=Class_soil[-nrow(Class_soil),]

# Environmental Matrix
env.matrix=data.frame(matrix(0,ncol=length(Species_soil),nrow=3),row.names=c("layer","Condition","Site"))
names(env.matrix)=names(Species_soil)
env.matrix[1,]=rep(c(rep("Mineral",3),rep("Organic",3)),3)
env.matrix[2,]=c(rep("Overburden",6),rep("Undisturbed",6),rep("Tailings",6))
env.matrix[3,]=c(rep("Mineral_Overburden",3),rep("Organic_Overburden",3),
rep("Mineral_Undisturbed",3),rep("Organic_Undisturbed",3),
rep("Mineral_Tailings",3),rep("Organic_Tailings",3))
env.matrix=t(env.matrix)

# Environmental Matrix with dummy variables
env.response=as.data.frame(cbind(env.matrix[,1],env.matrix[,1],env.matrix[,2],
env.matrix[,2]))

names(env.response)=c("Organic","Mineral","Tailings","Overburden")

env.response[,"Organic"]=rep(c(rep(0,3),rep(1,3)),3)
env.response[,"Mineral"]=rep(c(rep(1,3),rep(0,3)),3)
env.response[,"Tailings"]=c(rep(0,6),rep(0,6),rep(1,6))
env.response[,"Overburden"]=c(rep(1,6),rep(0,6),rep(0,6))

env=c(rep("Overburden",6),rep("Undisturbed",6),rep("Tailings",6))
names(env)=names(Species_soil)

#########################
# MULTIVARIATE ANALYSIS #
#########################

workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Report/Figures"
setwd(workingdirectory)

raremax.soil <- min(min(rowSums(t(Species_soil))),min(rowSums(t(Class_soil))))

# Class Matrix transformation
Class.r=rrarefy(t(Class_soil),raremax.soil)
# Species Hellinger transformation
Class.H=decostand(Class.r,method="hellinger")

# Matrix transformation
Species.r=rrarefy(t(Species_soil),raremax.soil)
# Hellinger transformation
Species.H=decostand(Species.r,method="hellinger")

########
# NMDS #
########

nmds=metaMDS(Species.H, dist="bray")

#env.matrixal Fit
ef <- envfit(nmds, as.data.frame(env.matrix), permu = 999)
ef

pdf(file = "FigureS4B.pdf", width=15, height=5)
par(mfcol=c(1,3),oma=c(0,0,4,0),mar=c(4,5,4,3))
ordiplot(nmds,display = "sites",cex=1.5,cex.lab=1.5,cex.axis=1.5)
ordihull(nmds,groups=as.character(env.matrix[,1]),draw="polygon",col="grey90",label=T,cex=1.5)
title("Layer effect : p-value = 0.938",cex.main=2)
ordiplot(nmds,display = "sites",cex=1.5,cex.lab=1.5,cex.axis=1.5)
ordihull(nmds,groups=as.character(env.matrix[,2]),draw="polygon",col="grey90",label=T,cex=1.5)
title("Treatment effect : p-value < 0.001",cex.main=2)
ordiplot(nmds,display = "sites",cex=1.5,cex.lab=1.5,cex.axis=1.5)
ordihull(nmds,groups=as.character(env.matrix[,3]),draw="polygon",col="grey90",label=T,cex=1.5)
title("Treatment*Layer effect : p-value < 0.001",cex.main=2)
mtext("B - Exploratory Analysis of Soil Dataset : NMDS",outer=T,cex=2,padj=-0.5)
dev.off()

#######
# CCA #
#######

#PLOT
pdf(file = "Figure8C.pdf",width=6,height=6)
par(mfcol=c(1,1))
colvec <- c("steelblue", "green4")
cca=cca(formula = Species.H ~ Tailings+Overburden, env.response)
anova(cca)
anova(cca, by = "margin")
plot(cca, display = "sites", type = "n")
ordisurf(cca, diversity(Species.H), add = TRUE)
ordiellipse(cca,groups=as.character(env),draw="polygon",label=F,cex=1)
with(env.response, points(cca, pch=16,col=colvec[as.numeric(env.response[,"Organic"])+1]))
with(env.response, legend("topleft", c("Mineral","Organic"), pch = 16,col=colvec,cex=1))
text(cca, display="bp",col="blue",cex=1,labels=c("Tailings","Overburden"))
title("CCA, Species abuncance ~ Reclamation")
mtext("p-value < 0.001",cex=1)
dev.off()

#######
# PCA #
#######

#PCA
if (sum(Class.H[,ncol(Class.H)])==0){
Class.H=Class.H[,-ncol(Class.H)]}

Class.H

pdf(file = "Figure8D.pdf",width=6,height=6)
par(mfcol=c(1,1))
# PCA using prcomp function
pca.class <- prcomp(Class.H)

g <- ggbiplot(pca.class,
              groups = env.matrix[,"Condition"],
              ellipse = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
g <- g + labs(title="PCA of supergroups abundance matrix")
print(g)
dev.off()




