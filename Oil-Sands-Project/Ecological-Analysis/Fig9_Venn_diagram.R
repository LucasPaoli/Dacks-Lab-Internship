#########################
# VENN DIAGRAM ANALYSIS #
#########################


library(vegan)
library(ggbiplot)
library(VennDiagram)

################
# LOADING DATA #
################

#Soil
workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Data/Soil/D1_G2/usearch_output_subsample_0.95/"
setwd(workingdirectory)

Species_soil_init=read.delim(file="Species_filtered.txt",row.names=1)
Species_soil_init=Species_soil_init[-nrow(Species_soil_init),]
Species.soil=as.data.frame(cbind(rowSums((Species_soil_init)[,1:6]),
                                 rowSums((Species_soil_init)[,13:18]),
                                 rowSums((Species_soil_init)[,7:12])))

names(Species.soil)=c("Undisturbed.soil","Overburden.soil","Tailings.soil")

Species.soil=as.data.frame(Species.soil>0)

#Base Mine Lake
workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Data/BML2/usearch/output_BML2_95"
setwd(workingdirectory)

Species_BML=read.delim(file="Species_filtered.txt",row.names=1)
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
Species_BML=as.data.frame(rowSums(Species_sub))
names(Species_BML)="BML"
Species_BML=as.data.frame(Species_BML>0)

workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Data/MLSB_metagenome"
setwd(workingdirectory)

Species_tailings=read.delim(file="species_jeu_blast_PR2.csv",header=T,sep=",")
Species_tailings_init=as.data.frame(cbind(as.character(unique(Species_tailings[,1])),rep(TRUE,length(unique(Species_tailings[,1])))))

Species_tailings=as.data.frame(Species_tailings_init[,2])
row.names(Species_tailings)=Species_tailings_init[,1]
names(Species_tailings)="Tailings"


N=unique(c(row.names(Species.soil),row.names(Species_tailings),row.names(Species_BML)))
Venn=data.frame(matrix(data=FALSE,nrow=length(N),ncol=5),row.names=N)
names(Venn)=c("Undisturbed.soil","Overburden.soil","Tailings.soil","Tailings","BML")
for (i in 1:nrow(Species.soil)){
  Venn[row.names(Species.soil)[i],"Undisturbed.soil"]=Species.soil[i,"Undisturbed.soil"]
  Venn[row.names(Species.soil)[i],"Overburden.soil"]=Species.soil[i,"Overburden.soil"]
  Venn[row.names(Species.soil)[i],"Tailings.soil"]=Species.soil[i,"Tailings.soil"]
}
for (i in 1:nrow(Species_tailings)){
  Venn[row.names(Species_tailings)[i],"Tailings"]=Species_tailings[i,"Tailings"]
}
for (i in 1:nrow(Species_BML)){
  Venn[row.names(Species_BML)[i],"BML"]=Species_BML[i,"BML"]
}

setwd("/Users/Lucas/Documents/ENS/BIO_M1_stage/Report/Figures")

pdf(file = "Figure9.pdf", width=7, height=7)

grid.newpage()
draw.quintuple.venn(area1 = nrow(subset(Venn, Undisturbed.soil == 1)), 
                 area2 = nrow(subset(Venn, Overburden.soil == 1)),
                 area3 = nrow(subset(Venn, Tailings.soil == 1)),
                 area4 = nrow(subset(Venn, Tailings == 1)),
                 area5 = nrow(subset(Venn, BML == 1)),
                 n12 = nrow(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1)),
                 n13 = nrow(subset(Venn, Undisturbed.soil == 1 & Tailings.soil == 1)),
                 n14 = nrow(subset(Venn, Undisturbed.soil == 1 & Tailings == 1)),
                 n15 = nrow(subset(Venn, Undisturbed.soil == 1 & BML == 1)),
                 n23 = nrow(subset(Venn, Overburden.soil == 1 & Tailings.soil == 1)),
                 n24 = nrow(subset(Venn, Overburden.soil == 1 & Tailings == 1)),
                 n25 = nrow(subset(Venn, Overburden.soil == 1 & BML == 1)),
                 n34 = nrow(subset(Venn, Tailings.soil == 1 & Tailings == 1)),
                 n35 = nrow(subset(Venn, Tailings.soil == 1 & BML == 1)),
                 n45 = nrow(subset(Venn, Tailings == 1 & BML == 1)), 
                 n123 = nrow(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings.soil == 1)),
                 n124 = nrow(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings == 1)),
                 n125 = nrow(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & BML == 1)),
                 n134 = nrow(subset(Venn, Undisturbed.soil == 1 & Tailings.soil == 1 & Tailings == 1)),
                 n135 = nrow(subset(Venn, Undisturbed.soil == 1 & Tailings.soil == 1 & BML == 1)),
                 n145 = nrow(subset(Venn, Undisturbed.soil == 1 & Tailings == 1 & BML == 1)),
                 n234 = nrow(subset(Venn, Overburden.soil == 1 & Tailings.soil == 1 & Tailings == 1)),
                 n235 = nrow(subset(Venn, Overburden.soil == 1 & Tailings.soil == 1 & BML == 1)),
                 n245 = nrow(subset(Venn, Overburden.soil == 1 & Tailings == 1 & BML == 1)),
                 n345 = nrow(subset(Venn, Tailings.soil == 1 & Tailings == 1 & BML == 1)),
                 n1234 = nrow(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings.soil == 1 & Tailings == 1)),
                 n1235 = nrow(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings.soil == 1 & BML == 1)),
                 n1245 = nrow(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings == 1 & BML == 1)),
                 n1345 = nrow(subset(Venn, Undisturbed.soil == 1 & Tailings.soil == 1 & Tailings == 1 & BML == 1)),
                 n2345 = nrow(subset(Venn, Overburden.soil == 1 & Tailings.soil == 1 & Tailings == 1 & BML == 1)),
                 n12345 = nrow(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings.soil == 1 & Tailings == 1 & BML == 1)),
                 category = c("Undisturbed.soil","Overburden.soil","Tailings.soil","Tailings","BML"), 
                 lty = "blank", 
                 fill = c("sandybrown","maroon", "mediumorchid","pink1","skyblue"),
                 margin=0.1)
dev.off()

workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Data/Soil/D1_G2/usearch_output_subsample_0.95/"
setwd(workingdirectory)

Species_soil_init=read.delim(file="Species_filtered.txt",row.names=1)
Species_soil_init=Species_soil_init[-nrow(Species_soil_init),]
Species.soil=as.data.frame(cbind(rowSums((Species_soil_init)[,1:6]),
                                 rowSums((Species_soil_init)[,13:18]),
                                 rowSums((Species_soil_init)[,7:12])))

names(Species.soil)=c("Undisturbed.soil","Overburden.soil","Tailings.soil")
Species.soil=sweep(Species.soil,2,colSums(Species.soil),"/")


workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Data/BML2/usearch/output_BML2_95"
setwd(workingdirectory)

Species_BML=read.delim(file="Species_filtered.txt",row.names=1)
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
Species_BML=as.data.frame(rowSums(Species_sub))
names(Species_BML)="BML"

Species_BML=sweep(Species_BML,2,colSums(Species_BML),"/")

workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Data/MLSB_metagenome"
setwd(workingdirectory)

Species_tailings_init=read.delim(file="Species.tailings.csv",header=T,sep=";")
Species_tailings=as.data.frame(rep(0,length(unique(Species_tailings_init[,1]))))
row.names(Species_tailings)=unique(Species_tailings_init[,1])
names(Species_tailings)="Tailings"

for (i in 1:nrow(Species_tailings_init)){
  Species_tailings[as.character(Species_tailings_init[i,1]),1]=Species_tailings[Species_tailings_init[i,1],1]+Species_tailings_init[i,2]
}


Species_tailings=sweep(Species_tailings,2,colSums(Species_tailings),"/")

print("Attribution shared by all the samples")

row.names(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings.soil == 1 & Tailings == 1 & BML == 1))
colSums(100*Species.soil[row.names(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings.soil == 1 & Tailings == 1 & BML == 1)),])
sum(100*Species_BML[row.names(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings.soil == 1 & Tailings == 1 & BML == 1)),])
sum(100*Species_tailings[row.names(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings.soil == 1 & Tailings == 1 & BML == 1)),])

print("Attribution shared by Soils")
row.names(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings.soil == 1 & Tailings == 0 & BML == 0))
colSums(100*Species.soil[row.names(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings.soil == 1 & Tailings == 0 & BML == 0)),])

print("Attribution shared by Soils&Tailings")
row.names(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings.soil == 1 & Tailings == 1 & BML == 0))
colSums(100*Species.soil[row.names(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings.soil == 1 & Tailings == 1 & BML == 0)),])
sum(100*Species_tailings[row.names(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings.soil == 1 & Tailings == 1 & BML == 0)),])

print("Attribution shared by Lake&Tailings")
row.names(subset(Venn, Undisturbed.soil == 0 & Overburden.soil == 0 & Tailings.soil == 0 & Tailings == 1 & BML == 1))
print("lake")
sum(100*Species_BML[row.names(subset(Venn, Undisturbed.soil == 0 & Overburden.soil == 0 & Tailings.soil == 0 & Tailings == 1 & BML == 1)),])
print("tailings")
sum(100*Species_tailings[row.names(subset(Venn, Undisturbed.soil == 0 & Overburden.soil == 0 & Tailings.soil == 0 & Tailings == 1 & BML == 1)),])

print("Soils+BML")
row.names(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings.soil == 1 & Tailings == 0 & BML == 1))
print("Abundance in BML")
sum(100*Species_BML[row.names(subset(Venn, Undisturbed.soil == 1 & Overburden.soil == 1 & Tailings.soil == 1 & Tailings == 0 & BML == 1)),])

print("BML only")
sum(100*Species_BML[row.names(subset(Venn, Undisturbed.soil == 0 & Overburden.soil == 0 & Tailings.soil == 0 & Tailings == 0 & BML == 1)),])
print("tailings only")
sum(100*Species_tailings[row.names(subset(Venn, Undisturbed.soil == 0 & Overburden.soil == 0 & Tailings.soil == 0 & Tailings == 1 & BML == 0)),])


