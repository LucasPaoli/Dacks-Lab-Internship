##################################################
# Statistical Analysis of Dinos ARN editing Data #
##################################################

## We focus on 4 organisms : kv, km, sm, pl
## And 4 genes family, atp, psa, psb, pet

##########################
# EDITING SCORE ANALYSIS #
##########################

library(lattice)
library(ggbiplot)

setwd("/Users/Lucas/dropbox/Chris/Very\ Final\ Stuff/Amino_Acid_Analysis_V3")
table=read.delim(file="table_analysis_global_lucas.csv",sep=",")


#####################################
# NORMALITY TESTS AND MISCELLEANOUS #
#####################################

# Testing the normal distribution
hist(table$average.edit.score.diff,breaks=seq(-7,7,0.5))
shapiro.test(table$average.edit.score.diff)
table.F=subset(table,plastid=="fuco")
table.P=subset(table,plastid=="pere")
shapiro.test(table.F$average.edit.score.diff)
hist(table.F$average.edit.score.diff,breaks=seq(-7,7,0.5))
shapiro.test(table.P$average.edit.score.diff)
hist(table.P$average.edit.score.diff,breaks=seq(-7,7,0.5))

# testing for organism effect
hist(table.F$average.edit.score.diff,breaks=seq(-7,7,0.5))
shapiro.test(table.F$average.edit.score.diff)
table.kv=subset(table,organism=="kv")
table.km=subset(table,organism=="km")
shapiro.test(table.kv$average.edit.score.diff)
shapiro.test(table.km$average.edit.score.diff)
hist(table.kv$average.edit.score.diff,breaks=seq(-7,7,0.5))
hist(table.km$average.edit.score.diff,breaks=seq(-7,7,0.5))
hist(table.P$average.edit.score.diff,breaks=seq(-7,7,0.5))
shapiro.test(table.P$average.edit.score.diff)
table.sm=subset(table,organism=="sm")
table.pl=subset(table,organism=="pl")
shapiro.test(table.sm$average.edit.score.diff)
shapiro.test(table.pl$average.edit.score.diff)
hist(table.sm$average.edit.score.diff,breaks=seq(-7,7,0.5))
hist(table.pl$average.edit.score.diff,breaks=seq(-7,7,0.5))

# Testing reference effect

boxplot(table$average.edit.score.diff~table$reference)
boxplot(table.F$average.edit.score.diff~table.F$reference)
boxplot(table.P$average.edit.score.diff~table.P$reference)

kruskal.test(table$average.edit.score.diff~table$reference)
kruskal.test(table.F$average.edit.score.diff~table.F$reference)
kruskal.test(table.P$average.edit.score.diff~table.P$reference)

# As their is no ref effect : get the average for the 5 ref

########################
# AVERAGING REFERENCES #
########################
table=read.delim(file="table_analysis_global_lucas.csv",sep=",")
table=cbind(table,rep(0,nrow(table)))
names(table)=c(names(table)[1:length(names(table))-1],"Number")

table.ac=subset(table,table$reference=="Ac")
table.ct=subset(table,table$reference=="Ct")
table.pt=subset(table,table$reference=="Pt")
table.eh=subset(table,table$reference=="Eh")
table.vb=subset(table,table$reference=="Vb")

flag=c(table.ac, table.ct, table.pt, table.eh, table.vb)
Number=data.frame(0)
table_sub=cbind(table[0,])


averaging.reference <- function(table_sub, j){
  for (i in 1:nrow(j)){
    line=subset(table_sub, gene==j[i,]$gene & organism==j[i,]$organism)
    if (nrow(line)==0){
      table_sub=rbind(table_sub,j[i,])
    }else{
      table_sub[row.names(line),"num.AA.changes"]=
        table_sub[row.names(line),"num.AA.changes"]+j[i,"num.AA.changes"]
      table_sub[row.names(line),"average.edit.score.diff"]=
        table_sub[row.names(line),"average.edit.score.diff"]+j[i,"average.edit.score.diff"]
      table_sub[row.names(line),"Number"]=
        table_sub[row.names(line),"Number"]+1
    }
  }
  return(table_sub)
}

j=table.ac
table_sub=averaging.reference(table_sub,j)
j=table.ct
table_sub=averaging.reference(table_sub,j)
j=table.vb
table_sub=averaging.reference(table_sub,j)
j=table.eh
table_sub=averaging.reference(table_sub,j)
j=table.pt
table_sub=averaging.reference(table_sub,j)

for (i in 1:nrow(table_sub)){
  table_sub[i,"num.AA.changes"]=table_sub[i,"num.AA.changes"]/(table_sub[i,"Number"]+1)
  table_sub[i,"average.edit.score.diff"]=table_sub[i,"average.edit.score.diff"]/(table_sub[i,"Number"]+1)
}

summary(table.ac)
summary(table.vb)
summary(table.pt)
summary(table.eh)
summary(table.ct)

summary(table_sub)

# redifining variables

table_sub.F=subset(table_sub,plastid=="fuco")
table_sub.P=subset(table_sub,plastid=="pere")
table_sub.sm=subset(table_sub,organism=="sm")
table_sub.pl=subset(table_sub,organism=="pl")
table_sub.kv=subset(table_sub,organism=="kv")
table_sub.km=subset(table_sub,organism=="km")

##############################################
# CORRELATION EDITING SCORE NBER OF AA EDITS #
##############################################

## NB : DONE FOR THE AVERAGE OF REFERENCES

# Overall
cor.test(table_sub[,3],table_sub[,4])
xyplot(table_sub[,3]~table_sub[,4])

# Fucoxanthin
cor.test(table_sub.F[,3],table_sub.F[,4])
xyplot(table_sub.F[,3]~table_sub.F[,4])
# Peridinin
cor.test(table_sub.P[,3],table_sub.P[,4])
xyplot(table_sub.P[,3]~table_sub.P[,4])

###########
# Figures #
###########

pdf(file= "Editing_score_summary.pdf", width=28, height=12)
par(mfrow=c(3,7),mar=c(4,7,7,4))

############ Reference ############

boxplot(table$average.edit.score.diff~table$reference,
        col="lightgrey",
        main="Editing scores depending on references organisms",
        ylab="Editing score",
        ylim=c(-5,5),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:5,
     labels=c(paste0("Ac\n n = ",nrow(subset(table,reference=="Ac"))),
              paste0("Ct\n n = ",nrow(subset(table,reference=="Ct"))),
              paste0("Eh\n n = ",nrow(subset(table,reference=="Eh"))),
              paste0("Pt\n n = ",nrow(subset(table,reference=="Pt"))),
              paste0("Vb\n n = ",nrow(subset(table,reference=="Vb")))
              )
     )
kruskal.test(table$average.edit.score.diff~table$reference)
shapiro.test(table$average.edit.score.diff)
object=aov(table$average.edit.score.diff~table$reference)
summary(object)
TukeyHSD(object)
mtext("OVERALL",side=3,padj=-5.5)
mtext("REFERENCE EFFECT",side=2,padj=-5)

boxplot(table.F$average.edit.score.diff~table.F$reference,
        col="lightgrey",
        main="Editing scores depending on references organisms \nwithin Peridinin plastids",
        ylab="Editing score",
        ylim=c(-5,5),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:5,
     labels=c(paste0("Ac\n n = ",nrow(subset(table.F,reference=="Ac"))),
              paste0("Ct\n n = ",nrow(subset(table.F,reference=="Ct"))),
              paste0("Eh\n n = ",nrow(subset(table.F,reference=="Eh"))),
              paste0("Pt\n n = ",nrow(subset(table.F,reference=="Pt"))),
              paste0("Vb\n n = ",nrow(subset(table.F,reference=="Vb")))
              )
     )
kruskal.test(table.F$average.edit.score.diff~table.F$reference)
shapiro.test(table.F$average.edit.score.diff)
object=aov(table.F$average.edit.score.diff~table.F$reference)
summary(object)
TukeyHSD(object)
mtext("FUCOXANTHIN",side=3,padj=-5.5)

boxplot(table.P$average.edit.score.diff~table.P$reference,
        col="lightgrey",
        main="Editing scores depending on references organisms \nwithin Peridinin plastids",
        ylab="Editing score",
        ylim=c(-5,5),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:5,
     labels=c(paste0("Ac\n n = ",nrow(subset(table.P,reference=="Ac"))),
                     paste0("Ct\n n = ",nrow(subset(table.P,reference=="Ct"))),
                     paste0("Eh\n n = ",nrow(subset(table.P,reference=="Eh"))),
                     paste0("Pt\n n = ",nrow(subset(table.P,reference=="Pt"))),
                     paste0("Vb\n n = ",nrow(subset(table.P,reference=="Vb")))
              )
     )
kruskal.test(table.P$average.edit.score.diff~table.P$reference)
shapiro.test(table.P$average.edit.score.diff)
object=aov(table.P$average.edit.score.diff~table.P$reference)
summary(object)
TukeyHSD(object)
mtext("PERIDININ",side=3,padj=-5.5)

formule=subset(table, organism=="kv")$average.edit.score.diff~
  subset(table, organism=="kv")$reference
boxplot(formule,
        col="lightgrey",
        main="Editing scores depending on references organisms \nfor KV",
        ylab="Editing score",
        ylim=c(-5,5),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:5,
     labels=c(paste0("Ac\n n = ",nrow(subset(subset(table, organism=="kv"),reference=="Ac"))),
              paste0("Ct\n n = ",nrow(subset(subset(table, organism=="kv"),reference=="Ct"))),
              paste0("Eh\n n = ",nrow(subset(subset(table, organism=="kv"),reference=="Eh"))),
              paste0("Pt\n n = ",nrow(subset(subset(table, organism=="kv"),reference=="Pt"))),
              paste0("Vb\n n = ",nrow(subset(subset(table, organism=="kv"),reference=="Vb")))
     )
)
kruskal.test(formule)
shapiro.test(subset(table, organism=="kv")$average.edit.score.diff)
object=aov(formule)
summary(object)
TukeyHSD(object)
mtext("KV",side=3,padj=-5.5)

formule=subset(table, organism=="km")$average.edit.score.diff~
  subset(table, organism=="km")$reference
boxplot(formule,
        col="lightgrey",
        main="Editing scores depending on references organisms \nfor KM",
        ylab="Editing score",
        ylim=c(-5,5),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:5,
     labels=c(paste0("Ac\n n = ",nrow(subset(subset(table, organism=="km"),reference=="Ac"))),
              paste0("Ct\n n = ",nrow(subset(subset(table, organism=="km"),reference=="Ct"))),
              paste0("Eh\n n = ",nrow(subset(subset(table, organism=="km"),reference=="Eh"))),
              paste0("Pt\n n = ",nrow(subset(subset(table, organism=="km"),reference=="Pt"))),
              paste0("Vb\n n = ",nrow(subset(subset(table, organism=="km"),reference=="Vb")))
     )
)
kruskal.test(formule)
shapiro.test(subset(table, organism=="km")$average.edit.score.diff)
object=aov(formule)
summary(object)
TukeyHSD(object)
mtext("KM",side=3,padj=-5.5)

formule=subset(table, organism=="sm")$average.edit.score.diff~
  subset(table, organism=="sm")$reference
boxplot(formule,
        col="lightgrey",
        main="Editing scores depending on references organisms \nfor SM",
        ylab="Editing score",
        ylim=c(-5,5),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:5,
     labels=c(paste0("Ac\n n = ",nrow(subset(subset(table, organism=="sm"),reference=="Ac"))),
              paste0("Ct\n n = ",nrow(subset(subset(table, organism=="sm"),reference=="Ct"))),
              paste0("Eh\n n = ",nrow(subset(subset(table, organism=="sm"),reference=="Eh"))),
              paste0("Pt\n n = ",nrow(subset(subset(table, organism=="sm"),reference=="Pt"))),
              paste0("Vb\n n = ",nrow(subset(subset(table, organism=="sm"),reference=="Vb")))
     )
)
kruskal.test(formule)
shapiro.test(subset(table, organism=="sm")$average.edit.score.diff)
object=aov(formule)
summary(object)
TukeyHSD(object)
mtext("SM",side=3,padj=-5.5)

formule=subset(table, organism=="pl")$average.edit.score.diff~
  subset(table, organism=="pl")$reference
boxplot(formule,
        col="lightgrey",
        main="Editing scores depending on references organisms \nfor PL",
        ylab="Editing score",
        ylim=c(-5,5),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:5,
     labels=c(paste0("Ac\n n = ",nrow(subset(subset(table, organism=="pl"),reference=="Ac"))),
              paste0("Ct\n n = ",nrow(subset(subset(table, organism=="pl"),reference=="Ct"))),
              paste0("Eh\n n = ",nrow(subset(subset(table, organism=="pl"),reference=="Eh"))),
              paste0("Pt\n n = ",nrow(subset(subset(table, organism=="pl"),reference=="Pt"))),
              paste0("Vb\n n = ",nrow(subset(subset(table, organism=="pl"),reference=="Vb")))
     )
)
kruskal.test(formule)
shapiro.test(subset(table, organism=="pl")$average.edit.score.diff)
#hist(subset(table, organism=="kv")$average.edit.score.diff)
object=aov(formule)
summary(object)
TukeyHSD(object)
mtext("PL",side=3,padj=-5.5)

############ Lineage ############

boxplot(table_sub.F$average.edit.score.diff,table_sub.P$average.edit.score.diff,
        main="Editing scores for each plastid type",
        ylab="Editing score",
        col="lightgrey",
        ylim=c(-6,7),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:2,
     labels=c(paste0("Fucoxanthin\n n = ",nrow(table_sub.F)),
              paste0("Peridinin\n n = ",nrow(table_sub.P))
     )
)
t.test(table_sub.F$average.edit.score.diff,table_sub.P$average.edit.score.diff)
lines(x=c(1,2),y=c(6,6))
text(x=1.5, y=6, "**", pos=3, cex=1)
mtext("LINEAGE EFFECT",side=2,padj=-5)

frame()
frame()

boxplot(table_sub$average.edit.score.diff~table_sub$organism,
        col="lightgrey",
        main="Editing scores for each organism",
        ylab="Editing score",
        ylim=c(-6,9),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("Km\n n = ",nrow(subset(table_sub, organism=="km"))),
              paste0("Kv\n n = ",nrow(subset(table_sub, organism=="kv"))),
              paste0("Pl\n n = ",nrow(subset(table_sub, organism=="pl"))),
              paste0("Sm\n n = ",nrow(subset(table_sub, organism=="sm")))
     )
)
shapiro.test(table_sub.kv$average.edit.score.diff)
shapiro.test(table_sub.km$average.edit.score.diff)
shapiro.test(table_sub.sm$average.edit.score.diff)
shapiro.test(table_sub.pl$average.edit.score.diff)
#hist(table_sub.pl$average.edit.score.diff)
#hist(table_sub.sm$average.edit.score.diff)
t.test(table_sub.kv$average.edit.score.diff,table_sub.km$average.edit.score.diff)
t.test(table_sub.sm$average.edit.score.diff,table_sub.pl$average.edit.score.diff)
object=aov(table_sub$average.edit.score.diff~table_sub$organism)
summary(object)
TukeyHSD(object)
t.test(table_sub.kv$average.edit.score.diff,table_sub.pl$average.edit.score.diff)
wilcox.test(table_sub.kv$average.edit.score.diff,table_sub.pl$average.edit.score.diff)
t.test(table_sub.kv$average.edit.score.diff,table_sub.sm$average.edit.score.diff)
wilcox.test(table_sub.kv$average.edit.score.diff,table_sub.sm$average.edit.score.diff)
t.test(table_sub.km$average.edit.score.diff,table_sub.pl$average.edit.score.diff)
wilcox.test(table_sub.km$average.edit.score.diff,table_sub.pl$average.edit.score.diff)
t.test(table_sub.km$average.edit.score.diff,table_sub.sm$average.edit.score.diff)
wilcox.test(table_sub.km$average.edit.score.diff,table_sub.sm$average.edit.score.diff)
lines(x=c(2,3),y=c(6,6))
text(x=2.5, y=5.5, "**", pos=3, cex=1)
lines(x=c(2,4),y=c(7,7))
text(x=3, y=6.5, "*", pos=3, cex=1)
lines(x=c(1,3),y=c(8,8))
text(x=2, y=7.5, "*", pos=3, cex=1)


frame()
frame()
frame()

############ Gene family ############

boxplot(table_sub$average.edit.score.diff~table_sub$gene_family,
        col="lightgrey",
        main="Editing scores depending on gene's family",
        ylab="Editing score",
        ylim=c(-6,6),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("atp\n n = ",nrow(subset(table_sub, gene_family=="atp"))),
              paste0("pet\n n = ",nrow(subset(table_sub, gene_family=="pet"))),
              paste0("psa\n n = ",nrow(subset(table_sub, gene_family=="psa"))),
              paste0("psb\n n = ",nrow(subset(table_sub, gene_family=="psb")))
     )
)
object=aov(table_sub$average.edit.score.diff~table_sub$gene_family)
summary(object)
TukeyHSD(object)
mtext("GENE FAMILY EFFECT",side=2,padj=-5)


boxplot(table_sub.F$average.edit.score.diff~table_sub.F$gene_family,        
        col="lightgrey",
        main="Editing scores depending on gene's family \nwithin Fucoxanthin plastids",
        ylab="Editing score",
        ylim=c(-6,6),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("atp\n n = ",nrow(subset(table_sub.F, gene_family=="atp"))),
              paste0("pet\n n = ",nrow(subset(table_sub.F, gene_family=="pet"))),
              paste0("psa\n n = ",nrow(subset(table_sub.F, gene_family=="psa"))),
              paste0("psb\n n = ",nrow(subset(table_sub.F, gene_family=="psb")))
     )
)
object=aov(table_sub.P$average.edit.score.diff~table_sub.P$gene_family)
summary(object)
TukeyHSD(object)

boxplot(table_sub.P$average.edit.score.diff~table_sub.P$gene_family,
        col="lightgrey",
        main="Editing scores depending on gene's family \nwithin Peridinin plastids",
        ylab="Editing score",
        ylim=c(-6,6),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("atp\n n = ",nrow(subset(table_sub.P, gene_family=="atp"))),
              paste0("pet\n n = ",nrow(subset(table_sub.P, gene_family=="pet"))),
              paste0("psa\n n = ",nrow(subset(table_sub.P, gene_family=="psa"))),
              paste0("psb\n n = ",nrow(subset(table_sub.P, gene_family=="psb")))
     )
)
object=aov(table_sub.F$average.edit.score.diff~table_sub.F$gene_family)
summary(object)
TukeyHSD(object)

formule=subset(table_sub, organism=="kv")$average.edit.score.diff~
  subset(table_sub, organism=="kv")$gene_family
boxplot(formule,
        col="lightgrey",
        main="Editing scores depending on gene's family \nfor KV",
        ylab="Editing score",
        ylim=c(-6,6),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("atp\n n = ",nrow(subset(subset(table_sub, organism=="kv"), gene_family=="atp"))),
              paste0("pet\n n = ",nrow(subset(subset(table_sub, organism=="kv"), gene_family=="pet"))),
              paste0("psa\n n = ",nrow(subset(subset(table_sub, organism=="kv"), gene_family=="psa"))),
              paste0("psb\n n = ",nrow(subset(subset(table_sub, organism=="kv"), gene_family=="psb")))
     )
)
object=aov(formule)
summary(object)
TukeyHSD(object)

formule=subset(table_sub, organism=="km")$average.edit.score.diff~
  subset(table_sub, organism=="km")$gene_family
boxplot(formule,
        col="lightgrey",
        main="Editing scores depending on gene's family \nfor KM",
        ylab="Editing score",
        ylim=c(-6,6),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("atp\n n = ",nrow(subset(subset(table_sub, organism=="km"), gene_family=="atp"))),
              paste0("pet\n n = ",nrow(subset(subset(table_sub, organism=="km"), gene_family=="pet"))),
              paste0("psa\n n = ",nrow(subset(subset(table_sub, organism=="km"), gene_family=="psa"))),
              paste0("psb\n n = ",nrow(subset(subset(table_sub, organism=="km"), gene_family=="psb")))
     )
)
object=aov(formule)
summary(object)
TukeyHSD(object)
lines(x=c(1,4),y=c(4,4))
text(x=2.5, y=4, "*", pos=3, cex=1)

formule=subset(table_sub, organism=="sm")$average.edit.score.diff~
  subset(table_sub, organism=="sm")$gene_family
boxplot(formule,
        col="lightgrey",
        main="Editing scores depending on gene's family \nfor SM",
        ylab="Editing score",
        ylim=c(-6,6),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("atp\n n = ",nrow(subset(subset(table_sub, organism=="sm"), gene_family=="atp"))),
              paste0("pet\n n = ",nrow(subset(subset(table_sub, organism=="sm"), gene_family=="pet"))),
              paste0("psa\n n = ",nrow(subset(subset(table_sub, organism=="sm"), gene_family=="psa"))),
              paste0("psb\n n = ",nrow(subset(subset(table_sub, organism=="sm"), gene_family=="psb")))
     )
)
object=aov(formule)
summary(object)
TukeyHSD(object)

formule=subset(table_sub, organism=="pl")$average.edit.score.diff~
  subset(table_sub, organism=="pl")$gene_family
boxplot(formule,
        col="lightgrey",
        main="Editing scores depending on gene's family \nfor PL",
        ylab="Editing score",
        ylim=c(-6,6),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("atp\n n = ",nrow(subset(subset(table_sub, organism=="pl"), gene_family=="atp"))),
              paste0("pet\n n = ",nrow(subset(subset(table_sub, organism=="pl"), gene_family=="pet"))),
              paste0("psa\n n = ",nrow(subset(subset(table_sub, organism=="pl"), gene_family=="psa"))),
              paste0("psb\n n = ",nrow(subset(subset(table_sub, organism=="pl"), gene_family=="psb")))
     )
)
object=aov(formule)
summary(object)
TukeyHSD(object)

############ Gene Group ############

dev.off()



#Correlation with number of AA edit
pdf(file="Edit_score_number_edit.pdf", width=18, height=6)
# Overall
cor.test(table_sub[,3],table_sub[,4])
plot1=xyplot(table_sub[,3]~table_sub[,4],
       groups=table_sub$plastid,
       pch=16,
       col=c("green3","blue2"),
       main="Correlation between edit score and number of edits",
       ylab="Number of Edits",
       xlab="Average edit score differences",
       sub="Pearson corr. coeff. = 0.08261125 , p-value = 0.5268",
       auto.key=list(text=c("Fucoxanthin","Peridinin"),title="Plastid :",column=2,cex=0.8,points=F,col=c("green3","blue2"))
)

# Fuco
cor.test(table_sub.F[,3],table_sub.F[,4])
plot2=xyplot(table_sub.F[,3]~table_sub.F[,4],
       groups=table_sub$organism,
       col=c("green1","green4"),
       pch=16,
       main="Correlation between edit score and number of edits\nwithin Fucoxanthin plastids",
       ylab="Number of Edits",
       xlab="Average edit score differences",
       sub="Pearson corr. coeff. = 0.02741075 , p-value = 0.8702",
       auto.key=list(text=c("kv","km"),title="Organism :",column=2,cex=0.8,points=F,col=c("green1","green4"))
)

# Per
cor.test(table_sub.P[,3],table_sub.P[,4])
plot3=xyplot(table_sub.P[,3]~table_sub.P[,4],
       groups=table_sub$organism,
       col=c("lightblue","blue3"),
       pch=16,
       main="Correlation between edit score and number of edits\nwithin Fucoxanthin plastids",
       ylab="Number of Edits",
       xlab="Average edit score differences",
       sub="Pearson corr. coeff. = -0.1754453 , p-value = 0.4233",
       auto.key=list(text=c("pl","sm"),title="Organism :",column=2,cex=0.8,points=F,col=c("lightblue","blue3"))
)

print(plot1, position = c(0, 0, 0.33, 1), more = TRUE)
print(plot2, position = c(0.33, 0, 0.66, 1), more = TRUE)
print(plot3, position = c(0.66, 0, 1, 1))

dev.off()


#################################
## LINEAR MODELS (DNU FOR NOW) ##
#################################

model1=lm(table_sub$average.edit.score.diff ~ table_sub$plastid)
summary(model1)
aov(model1)
anova(model1)
plot(model1)
plot(resid(model1))

model1=lm(table_sub$average.edit.score.diff ~ table_sub$gene_family)
summary(model1)
aov(model1)
anova(model1)
plot(model1)
plot(resid(model1))

model1=lm(table_sub$average.edit.score.diff ~ table_sub$plastid * table_sub$gene_family)
summary(model1)
aov(model1)
anova(model1)
plot(model1)
plot(resid(model1))
