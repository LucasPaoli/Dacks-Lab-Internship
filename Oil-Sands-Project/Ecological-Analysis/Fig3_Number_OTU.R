#############################################################
# OTU CLUSTERING DATA DESCRIPTION ANALYSIS : NUMBER OF OTUS #
#############################################################


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

Species.reclamation=rbind(colSums(t(Species_soil)[1:6,]),
                          colSums(t(Species_soil)[7:12,]),
                          colSums(t(Species_soil)[13:18,]))

Class.reclamation=rbind(colSums(t(Class_soil)[1:6,]),
                          colSums(t(Class_soil)[7:12,]),
                          colSums(t(Class_soil)[13:18,]))

OTU.reclamation=rbind(colSums(t(OTU_soil)[1:6,]),
                          colSums(t(OTU_soil)[7:12,]),
                          colSums(t(OTU_soil)[13:18,]))

#######################################
# NUMBER OF OTU, TAXONOMIC ASSIGNMENT #
#######################################

workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Report/Figures"
setwd(workingdirectory)

sum(rowSums(OTU_sub)>0)
sum(rowSums(OTU_soil)>0)

pdf(file = "Figure3.pdf", width=11, height=5)

par(mfrow=c(1,2),oma=c(0,0,4,0))

assignment <- function(x)  # From a Sample (vector), get the proportion of attributed and non attributed reads
{
  N=sum(x)
  assigned=sum(x[1:length(x)-1])
  non.assigned=x["Unknown"]
  y=c(assigned/N,non.assigned/N)
  return(y)
}

assignment.matrix <- function(x,y)  # From a Sample set (matrix), get the numbers of attributed and non attributed reads
{
  z=matrix(0,nrow=2,ncol=ncol(x))
  for (i in 1:ncol(x)){
    N=sum(x[,i])
    z[1,i]=sum(x[1:nrow(x)-1,i])/N
    z[2,i]=x["Unknown",i]/N
    z[,i]=z[,i]*sum(y[,i]>=1)
  }
  return(z)
}

#######
# BML #
#######

legend_BML=c("Total",
         "May","June","September")

Total=sum(rowSums(OTU_sub)>=1)
May=sum(t(OTU.season)[,1]>=1)
June=sum(t(OTU.season)[,3]>=1)
September=sum(t(OTU.season)[,2]>=1)

Barplot = barplot(
  cbind(
  assignment(rowSums(Class_sub))*Total,
  assignment(t(Class.season)[,1])*May,
  assignment(t(Class.season)[,3])*June,
  assignment(t(Class.season)[,2])*September
  ),
space=0.05,
  ylab="Number of OTUs",
  ylim=c(0,350),
  col=c("steelblue","grey40"),
  main="A - Base Mine Lake")
legend("topright", # places a legend at the appropriate place
       c("Assigned","Unassigned"), # puts text in the legend
       pch=c(15,15), # gives the legend appropriate symbol
       col=c("steelblue","grey40"),
       cex=1,border=NA,bty="n",title="% of sequences")
text(cex=0.8, x=Barplot, y = par("usr")[3]-15, legend_BML, xpd=TRUE, srt=20,adj=1)

########
# SOIL #
########

legend_soil=c("Total",
              "Undisturbed",
         "Overburden",
         "Tailings")

Total=sum(rowSums(OTU_soil)>=1)
Undisturbed=sum(t(OTU.reclamation)[,1]>=1)
Tailings=sum(t(OTU.reclamation)[,2]>=1)
Overburden=sum(t(OTU.reclamation)[,3]>=1)

Barplot = barplot(
  cbind(
    assignment(rowSums(Class_soil))*Total,
    assignment(t(Class.reclamation)[,1])*Undisturbed,
    assignment(t(Class.reclamation)[,3])*Overburden,
    assignment(t(Class.reclamation)[,2])*Tailings
  ),
  ylab="Number of OTUs",
  ylim=c(0,2500),space=0.05,
  col=c("steelblue","grey40"),
  main="B - Soil dataset")
legend("topright", # places a legend at the appropriate place 
       c("Assigned","Unassigned"), # puts text in the legend
       pch=c(15,15), # gives the legend appropriate symbol
       col=c("steelblue","grey40"),
       cex=1,border=NA,bty="n",title="% of sequences")

text(cex=0.8, x=Barplot, y = par("usr")[3]-107, legend_soil, xpd=TRUE, srt=20,adj=1)

mtext("Number of OTUs and proportion of OTUs\nwith a resolved taxonomic assignment",outer=T,cex=1.5)


dev.off()
