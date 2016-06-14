##################################################
# Statistical Analysis of Dinos ARN editing Data #
##################################################

## We focus on 4 organisms : kv, km, sm, pl
## And 4 genes family, atp, psa, psb, pet

library(lattice)

##############################
# Loading, transforming Data #
##############################

# /!\
# Optional if already done.
# /!\

###########
# ENTROPY #
###########

## LOADING DATA

setwd("/Users/Lucas/dropbox/Chris/Very\ Final\ Stuff/Entropy/")

names=c("alignment.position","genomic.position","edited.or.not",
        "positional.entropy","plastid",
        "organism")

kv=read.delim(file="Kveneficum/Spreadsheets/kv_entropy.csv",sep=",")
kv=cbind(kv,rep("kv",length(kv$alignment.position)))
kv=cbind(kv,rep("kv",length(kv$alignment.position)))
kv=na.omit(kv)
names(kv)=names

km=read.delim(file="Kmikimotoi/Spreadsheets/km_entropy.csv",sep=",")
km=cbind(km,rep("km",length(km$alignment.position)))
km=cbind(km,rep("km",length(km$alignment.position)))
names(km)=names

pl=read.delim(file="Plunula/Spreadsheets/pl_entropy.csv",sep=",")
pl=cbind(pl,rep("pl",length(pl$alignment.position)))
pl=cbind(pl,rep("pl",length(pl$alignment.position)))
names(pl)=names

sm=read.delim(file="Sminutum/Spreadsheets/sm_entropy.csv",sep=",")
sm=cbind(sm,rep("sm",length(sm$alignment.position)))
sm=cbind(sm,rep("sm",length(sm$alignment.position)))
names(sm)=names

table=rbind(kv,km,pl,sm)

table.0=subset(table,table$edited.or.not==0)
table.1=subset(table,table$edited.or.not==1)


## DISTRIBUTION ANALYSIS

hist(table$edited.or.not,breaks=seq(0,1,0.1))
hist(table$positional.entropy,breaks=seq(0,1,0.05))

hist(table.0$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in non edited codons",
     xlab="Postitional Entropy",
     col="lightgrey")
hist(table.1$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons",
     xlab="Postitional Entropy",
     col="lightgrey")

ks.test(table.0$positional.entropy,table.1$positional.entropy)

pdf(file = "Entropy distribtion.pdf", width=14, height=8)
par(mfcol=c(1,2))
hist(table.0$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in non edited codons",
     xlab="Postitional Entropy",
     col="lightgrey")
hist(table.1$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons",
     xlab="Postitional Entropy",
     col="lightgrey")
dev.off()

boxplot(table.0$positional.entropy,table.1$positional.entropy,names=c("Non edited positions", "Edited positions"))

pdf(file = "Entropy Boxplot.pdf", width=8, height=8)
boxplot(table.0$positional.entropy,table.1$positional.entropy,names=c("Non edited positions", "Edited positions"),
        main="Boxplot of positional entropy",
        ylab="positional entropy",
        col="lightgrey")
dev.off()

## PLASTID EFFECT ?
table.1.P=subset(table,table$plastid=="kv" | table$plastid=="km")
table.1.F=subset(table,table$plastid=="pl" | table$plastid=="sm")

boxplot(table.1.F$positional.entropy,table.1.P$positional.entropy,names=c("Fucoxanthine", "Peredinine"))
hist(table.1.F$positional.entropy,breaks=seq(0,1,0.05))
hist(table.1.P$positional.entropy,breaks=seq(0,1,0.05))

pdf(file = "Entropy Boxplot Plastid split.pdf", width=8, height=8)
boxplot(table.1.F$positional.entropy,table.1.P$positional.entropy,names=c("Fucoxanthine", "Peredinine"),
        main="Boxplot of positional entropy",
        ylab="positional entropy",
        col="lightgrey")
dev.off()
