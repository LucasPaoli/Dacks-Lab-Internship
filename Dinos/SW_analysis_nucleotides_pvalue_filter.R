library(lattice)

###########################
# SLIDING WINDOW ANALYSIS #
###########################

setwd("/Users/Lucas/dropbox/Chris/Very\ Final\ Stuff/Very_Final_Sliding_Window/Lucas")
SW=read.delim(file="SW_analysis_global_lucas.csv",sep=",")

SW=subset(SW,dt(SW$nucleotide_t_value,df=SW$DF_.N.2.)<=0.05)

# Testing the normal distribution
hist(SW$nucleotide_pearson_correlation,breaks=seq(-2,1,0.2))

shapiro.test(SW$nucleotide_pearson_correlation)
SW.F=subset(SW, organism== "kv" | organism== "km")
SW.P=subset(SW, organism== "sm" | organism== "pl")
shapiro.test(SW.F$nucleotide_pearson_correlation)
hist(SW.F$nucleotide_pearson_correlation,breaks=seq(-2,1,0.2))
shapiro.test(SW.P$nucleotide_pearson_correlation)
hist(SW.P$nucleotide_pearson_correlation,breaks=seq(-2,1,0.2))

################################
# Testing for reference effect #
################################

# Overall
boxplot(SW$nucleotide_pearson_correlation~SW$reference)
kruskal.test(SW$nucleotide_pearson_correlation~SW$reference)
summary(aov(SW$nucleotide_pearson_correlation~SW$reference))
TukeyHSD(aov(SW$nucleotide_pearson_correlation~SW$reference))

# Each organism

#KV
SW.organism=subset(SW, organism=="kv")
boxplot(SW.organism$nucleotide_pearson_correlation~SW.organism$reference)
kruskal.test(SW.organism$nucleotide_pearson_correlation~SW.organism$reference)
summary(aov(SW.organism$nucleotide_pearson_correlation~SW.organism$reference))
TukeyHSD(aov(SW.organism$nucleotide_pearson_correlation~SW.organism$reference))

#KM
SW.organism=subset(SW, organism=="km")
boxplot(SW.organism$nucleotide_pearson_correlation~SW.organism$reference)
kruskal.test(SW.organism$nucleotide_pearson_correlation~SW.organism$reference)
summary(aov(SW.organism$nucleotide_pearson_correlation~SW.organism$reference))
TukeyHSD(aov(SW.organism$nucleotide_pearson_correlation~SW.organism$reference))

#PL
SW.organism=subset(SW, organism=="pl")
boxplot(SW.organism$nucleotide_pearson_correlation~SW.organism$reference)
kruskal.test(SW.organism$nucleotide_pearson_correlation~SW.organism$reference)
summary(aov(SW.organism$nucleotide_pearson_correlation~SW.organism$reference))
TukeyHSD(aov(SW.organism$nucleotide_pearson_correlation~SW.organism$reference))

#SM
SW.organism=subset(SW, organism=="sm")
boxplot(SW.organism$nucleotide_pearson_correlation~SW.organism$reference)
kruskal.test(SW.organism$nucleotide_pearson_correlation~SW.organism$reference)
summary(aov(SW.organism$nucleotide_pearson_correlation~SW.organism$reference))
TukeyHSD(aov(SW.organism$nucleotide_pearson_correlation~SW.organism$reference))

########################
# AVERAGING REFERENCES #
########################

SW=read.delim(file="SW_analysis_global_lucas.csv",sep=",")
SW=subset(SW,dt(SW$nucleotide_t_value,df=SW$DF_.N.2.)<=0.05)
SW=cbind(SW,rep(0,nrow(SW)))
names(SW)=c(names(SW)[1:length(names(SW))-1],"Number")

SW.ac=subset(SW,SW$reference=="Ac")
SW.ct=subset(SW,SW$reference=="Ct")
SW.pt=subset(SW,SW$reference=="Pt")
SW.eh=subset(SW,SW$reference=="Eh")
SW.vb=subset(SW,SW$reference=="Vb")

flag=c(SW.ac, SW.ct, SW.pt, SW.eh, SW.vb)
Number=data.frame(0)
SW_sub=cbind(SW[0,])


averaging.reference <- function(SW_sub, j){
  for (i in 1:nrow(j)){
    line=subset(SW_sub, gene==j[i,]$gene & organism==j[i,]$organism)
    if (nrow(line)==0){
      SW_sub=rbind(SW_sub,j[i,])
    }else{
      SW_sub[row.names(line),"nucleotide_pearson_correlation"]=
        SW_sub[row.names(line),"nucleotide_pearson_correlation"]+j[i,"nucleotide_pearson_correlation"]
      SW_sub[row.names(line),"Number"]=
        SW_sub[row.names(line),"Number"]+1
    }
  }
  return(SW_sub)
}

j=SW.ac
SW_sub=averaging.reference(SW_sub,j)
j=SW.ct
SW_sub=averaging.reference(SW_sub,j)
j=SW.vb
SW_sub=averaging.reference(SW_sub,j)
j=SW.eh
SW_sub=averaging.reference(SW_sub,j)
j=SW.pt
SW_sub=averaging.reference(SW_sub,j)

for (i in 1:nrow(SW_sub)){
  SW_sub[i,"nucleotide_pearson_correlation"]=SW_sub[i,"nucleotide_pearson_correlation"]/(SW_sub[i,"Number"]+1)
}

summary(SW.ac)
summary(SW.vb)
summary(SW.pt)
summary(SW.eh)
summary(SW.ct)

summary(SW_sub)

# redifining variables

SW_sub.F=subset(SW_sub,plastid=="fuco")
SW_sub.P=subset(SW_sub,plastid=="pere")
SW_sub.sm=subset(SW_sub,organism=="sm")
SW_sub.pl=subset(SW_sub,organism=="pl")
SW_sub.kv=subset(SW_sub,organism=="kv")
SW_sub.km=subset(SW_sub,organism=="km")

##################
# Lineage effect #
##################

#PLastid effect
hist(SW_sub$nucleotide_pearson_correlation)
t.test(SW_sub.F$nucleotide_pearson_correlation,SW_sub.P$nucleotide_pearson_correlation)
wilcox.test(SW_sub.F$nucleotide_pearson_correlation,SW_sub.P$nucleotide_pearson_correlation)
boxplot(SW_sub.F$nucleotide_pearson_correlation,SW_sub.P$nucleotide_pearson_correlation)

# testing for organism effect

# FUCOXANTHIN

hist(SW_sub.F$nucleotide_pearson_correlation,breaks=seq(-2,1,0.2))
shapiro.test(SW_sub.F$nucleotide_pearson_correlation)

shapiro.test(SW_sub.kv$nucleotide_pearson_correlation)
shapiro.test(SW_sub.km$nucleotide_pearson_correlation)
hist(SW_sub.kv$nucleotide_pearson_correlation,breaks=seq(-2,1,0.2))
hist(SW_sub.km$nucleotide_pearson_correlation,breaks=seq(-2,1,0.2))

t.test(SW_sub.kv$nucleotide_pearson_correlation,SW_sub.km$nucleotide_pearson_correlation)
boxplot(SW_sub.kv$nucleotide_pearson_correlation,SW_sub.km$nucleotide_pearson_correlation)
summary(aov(SW_sub.F$nucleotide_pearson_correlation~SW_sub.F$organism))
# In case of non normality :
wilcox.test(SW_sub.kv$nucleotide_pearson_correlation,SW_sub.km$nucleotide_pearson_correlation)

## PERIDININ

hist(SW_sub.P$nucleotide_pearson_correlation,breaks=seq(-2,1,0.2))
shapiro.test(SW_sub.P$nucleotide_pearson_correlation)

shapiro.test(SW_sub.pl$nucleotide_pearson_correlation)
shapiro.test(SW_sub.sm$nucleotide_pearson_correlation)
hist(SW_sub.pl$nucleotide_pearson_correlation,breaks=seq(-2,1,0.2))
hist(SW_sub.sm$nucleotide_pearson_correlation,breaks=seq(-2,1,0.2))

t.test(SW_sub.pl$nucleotide_pearson_correlation,SW_sub.sm$nucleotide_pearson_correlation)
boxplot(SW_sub.pl$nucleotide_pearson_correlation,SW_sub.sm$nucleotide_pearson_correlation)
summary(aov(SW_sub.P$nucleotide_pearson_correlation~SW_sub.P$organism))
# In case of non normality :
wilcox.test(SW_sub.pl$nucleotide_pearson_correlation,SW_sub.sm$nucleotide_pearson_correlation)

######################
# Gene family effect #
######################

boxplot(SW_sub$nucleotide_pearson_correlation~SW_sub$gene_family)
boxplot(SW_sub.F$nucleotide_pearson_correlation~SW_sub.F$gene_family)
boxplot(SW_sub.P$nucleotide_pearson_correlation~SW_sub.P$gene_family)

kruskal.test(SW_sub$nucleotide_pearson_correlation~SW_sub$gene_family)
kruskal.test(SW_sub.F$nucleotide_pearson_correlation~SW_sub.F$gene_family)
kruskal.test(SW_sub.P$nucleotide_pearson_correlation~SW_sub.P$gene_family)

oneway.test(SW_sub$nucleotide_pearson_correlation~SW_sub$gene_family)
object=aov(SW_sub$nucleotide_pearson_correlation~SW_sub$gene_family)
summary(object)
TukeyHSD(object)

oneway.test(SW_sub.F$nucleotide_pearson_correlation~SW_sub.F$gene_family)
object=aov(SW_sub.F$nucleotide_pearson_correlation~SW_sub.F$gene_family)
summary(object)
TukeyHSD(object)

oneway.test(SW_sub.P$nucleotide_pearson_correlation~SW_sub.P$gene_family)
object=aov(SW_sub.P$nucleotide_pearson_correlation~SW_sub.P$gene_family)
summary(object)
TukeyHSD(object)

###########
# Figures #
###########


pdf(file= "SW_summary_nucleotide_pvalue.pdf", width=28, height=12)
par(mfrow=c(3,7),mar=c(4,7,7,4))

############ Reference ############

boxplot(SW$nucleotide_pearson_correlation~SW$reference,
        col="lightgrey",
        main="Sliding Window Pearson correlation \ndepending on references organisms",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:5,
     labels=c(paste0("Ac\n n = ",nrow(subset(SW,reference=="Ac"))),
              paste0("Ct\n n = ",nrow(subset(SW,reference=="Ct"))),
              paste0("Eh\n n = ",nrow(subset(SW,reference=="Eh"))),
              paste0("Pt\n n = ",nrow(subset(SW,reference=="Pt"))),
              paste0("Vb\n n = ",nrow(subset(SW,reference=="Vb")))
     )
)
kruskal.test(SW$nucleotide_pearson_correlation~SW$reference)
shapiro.test(SW$nucleotide_pearson_correlation)
object=aov(SW$nucleotide_pearson_correlation~SW$reference)
summary(object)
TukeyHSD(object)
mtext("OVERALL",side=3,padj=-5.5)
mtext("REFERENCE EFFECT",side=2,padj=-5)

boxplot(SW.F$nucleotide_pearson_correlation~SW.F$reference,
        col="lightgrey",
        main="Sliding Window Pearson correlation \ndepending on references organisms within Peridinin plastids",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:5,
     labels=c(paste0("Ac\n n = ",nrow(subset(SW.F,reference=="Ac"))),
              paste0("Ct\n n = ",nrow(subset(SW.F,reference=="Ct"))),
              paste0("Eh\n n = ",nrow(subset(SW.F,reference=="Eh"))),
              paste0("Pt\n n = ",nrow(subset(SW.F,reference=="Pt"))),
              paste0("Vb\n n = ",nrow(subset(SW.F,reference=="Vb")))
     )
)
kruskal.test(SW.F$nucleotide_pearson_correlation~SW.F$reference)
shapiro.test(SW.F$nucleotide_pearson_correlation)
object=aov(SW.F$nucleotide_pearson_correlation~SW.F$reference)
summary(object)
TukeyHSD(object)
mtext("FUCOXANTHIN",side=3,padj=-5.5)

boxplot(SW.P$nucleotide_pearson_correlation~SW.P$reference,
        col="lightgrey",
        main="Sliding Window Pearson correlation \ndepending on references organisms within Peridinin plastids",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:5,
     labels=c(paste0("Ac\n n = ",nrow(subset(SW.P,reference=="Ac"))),
              paste0("Ct\n n = ",nrow(subset(SW.P,reference=="Ct"))),
              paste0("Eh\n n = ",nrow(subset(SW.P,reference=="Eh"))),
              paste0("Pt\n n = ",nrow(subset(SW.P,reference=="Pt"))),
              paste0("Vb\n n = ",nrow(subset(SW.P,reference=="Vb")))
     )
)
kruskal.test(SW.P$nucleotide_pearson_correlation~SW.P$reference)
shapiro.test(SW.P$nucleotide_pearson_correlation)
object=aov(SW.P$nucleotide_pearson_correlation~SW.P$reference)
summary(object)
TukeyHSD(object)
mtext("PERIDININ",side=3,padj=-5.5)

formule=subset(SW, organism=="kv")$nucleotide_pearson_correlation~
  subset(SW, organism=="kv")$reference
boxplot(formule,
        col="lightgrey",
        main="Sliding Window Pearson correlation \ndepending on references organisms for KV",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:5,
     labels=c(paste0("Ac\n n = ",nrow(subset(subset(SW, organism=="kv"),reference=="Ac"))),
              paste0("Ct\n n = ",nrow(subset(subset(SW, organism=="kv"),reference=="Ct"))),
              paste0("Eh\n n = ",nrow(subset(subset(SW, organism=="kv"),reference=="Eh"))),
              paste0("Pt\n n = ",nrow(subset(subset(SW, organism=="kv"),reference=="Pt"))),
              paste0("Vb\n n = ",nrow(subset(subset(SW, organism=="kv"),reference=="Vb")))
     )
)
kruskal.test(formule)
shapiro.test(subset(SW, organism=="kv")$nucleotide_pearson_correlation)
object=aov(formule)
summary(object)
TukeyHSD(object)
mtext("KV",side=3,padj=-5.5)

formule=subset(SW, organism=="km")$nucleotide_pearson_correlation~
  subset(SW, organism=="km")$reference
boxplot(formule,
        col="lightgrey",
        main="Sliding Window Pearson correlation \ndepending on references organisms for KM",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:5,
     labels=c(paste0("Ac\n n = ",nrow(subset(subset(SW, organism=="km"),reference=="Ac"))),
              paste0("Ct\n n = ",nrow(subset(subset(SW, organism=="km"),reference=="Ct"))),
              paste0("Eh\n n = ",nrow(subset(subset(SW, organism=="km"),reference=="Eh"))),
              paste0("Pt\n n = ",nrow(subset(subset(SW, organism=="km"),reference=="Pt"))),
              paste0("Vb\n n = ",nrow(subset(subset(SW, organism=="km"),reference=="Vb")))
     )
)
kruskal.test(formule)
shapiro.test(subset(SW, organism=="km")$nucleotide_pearson_correlation)
object=aov(formule)
summary(object)
TukeyHSD(object)
mtext("KM",side=3,padj=-5.5)

formule=subset(SW, organism=="sm")$nucleotide_pearson_correlation~
  subset(SW, organism=="sm")$reference
boxplot(formule,
        col="lightgrey",
        main="Sliding Window Pearson correlation \ndepending on references organisms for SM",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:5,
     labels=c(paste0("Ac\n n = ",nrow(subset(subset(SW, organism=="sm"),reference=="Ac"))),
              paste0("Ct\n n = ",nrow(subset(subset(SW, organism=="sm"),reference=="Ct"))),
              paste0("Eh\n n = ",nrow(subset(subset(SW, organism=="sm"),reference=="Eh"))),
              paste0("Pt\n n = ",nrow(subset(subset(SW, organism=="sm"),reference=="Pt"))),
              paste0("Vb\n n = ",nrow(subset(subset(SW, organism=="sm"),reference=="Vb")))
     )
)
kruskal.test(formule)
shapiro.test(subset(SW, organism=="sm")$nucleotide_pearson_correlation)
object=aov(formule)
summary(object)
TukeyHSD(object)
mtext("SM",side=3,padj=-5.5)

formule=subset(SW, organism=="pl")$nucleotide_pearson_correlation~
  subset(SW, organism=="pl")$reference
boxplot(formule,
        col="lightgrey",
        main="Sliding Window Pearson correlation \ndepending on references organisms for PL",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:5,
     labels=c(paste0("Ac\n n = ",nrow(subset(subset(SW, organism=="pl"),reference=="Ac"))),
              paste0("Ct\n n = ",nrow(subset(subset(SW, organism=="pl"),reference=="Ct"))),
              paste0("Eh\n n = ",nrow(subset(subset(SW, organism=="pl"),reference=="Eh"))),
              paste0("Pt\n n = ",nrow(subset(subset(SW, organism=="pl"),reference=="Pt"))),
              paste0("Vb\n n = ",nrow(subset(subset(SW, organism=="pl"),reference=="Vb")))
     )
)
kruskal.test(formule)
shapiro.test(subset(SW, organism=="pl")$nucleotide_pearson_correlation)
#hist(subset(SW, organism=="kv")$nucleotide_pearson_correlation)
object=aov(formule)
summary(object)
TukeyHSD(object)
mtext("PL",side=3,padj=-5.5)

############ Lineage ############

boxplot(SW_sub.F$nucleotide_pearson_correlation,SW_sub.P$nucleotide_pearson_correlation,
        col="lightgrey",
        main="Sliding Window Pearson correlation \nfor each plastid type",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:2,
     labels=c(paste0("Fucoxanthin\n n = ",nrow(SW_sub.F)),
              paste0("Peridinin\n n = ",nrow(SW_sub.P))
     )
)
t.test(SW_sub.F$nucleotide_pearson_correlation,SW_sub.P$nucleotide_pearson_correlation)
mtext("LINEAGE EFFECT",side=2,padj=-5)

frame()
frame()

boxplot(SW_sub$nucleotide_pearson_correlation~SW_sub$organism,
        col="lightgrey",
        main="Sliding Window Pearson correlation \nfor each organism",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("Km\n n = ",nrow(subset(SW_sub, organism=="km"))),
              paste0("Kv\n n = ",nrow(subset(SW_sub, organism=="kv"))),
              paste0("Pl\n n = ",nrow(subset(SW_sub, organism=="pl"))),
              paste0("Sm\n n = ",nrow(subset(SW_sub, organism=="sm")))
     )
)
shapiro.test(SW_sub.kv$nucleotide_pearson_correlation)
shapiro.test(SW_sub.km$nucleotide_pearson_correlation)
shapiro.test(SW_sub.sm$nucleotide_pearson_correlation)
shapiro.test(SW_sub.pl$nucleotide_pearson_correlation)
#hist(SW_sub.pl$nucleotide_pearson_correlation)
#hist(SW_sub.sm$nucleotide_pearson_correlation)
#hist(SW_sub.kv$nucleotide_pearson_correlation)
#hist(SW_sub.km$nucleotide_pearson_correlation)
t.test(SW_sub.kv$nucleotide_pearson_correlation,SW_sub.km$nucleotide_pearson_correlation)
t.test(SW_sub.sm$nucleotide_pearson_correlation,SW_sub.pl$nucleotide_pearson_correlation)
wilcox.test(SW_sub.sm$nucleotide_pearson_correlation,SW_sub.pl$nucleotide_pearson_correlation)
object=aov(SW_sub$nucleotide_pearson_correlation~SW_sub$organism)
summary(object)
TukeyHSD(object)
t.test(SW_sub.kv$nucleotide_pearson_correlation,SW_sub.pl$nucleotide_pearson_correlation)
wilcox.test(SW_sub.kv$nucleotide_pearson_correlation,SW_sub.pl$nucleotide_pearson_correlation)
t.test(SW_sub.kv$nucleotide_pearson_correlation,SW_sub.sm$nucleotide_pearson_correlation)
wilcox.test(SW_sub.kv$nucleotide_pearson_correlation,SW_sub.sm$nucleotide_pearson_correlation)
t.test(SW_sub.km$nucleotide_pearson_correlation,SW_sub.pl$nucleotide_pearson_correlation)
wilcox.test(SW_sub.km$nucleotide_pearson_correlation,SW_sub.pl$nucleotide_pearson_correlation)
t.test(SW_sub.km$nucleotide_pearson_correlation,SW_sub.sm$nucleotide_pearson_correlation)
wilcox.test(SW_sub.km$nucleotide_pearson_correlation,SW_sub.sm$nucleotide_pearson_correlation)


frame()
frame()
frame()

############ Gene family ############

boxplot(SW_sub$nucleotide_pearson_correlation~SW_sub$gene_family,
        col="lightgrey",
        main="Sliding Window Pearson correlation depending on \ngene's family",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("atp\n n = ",nrow(subset(SW_sub, gene_family=="atp"))),
              paste0("pet\n n = ",nrow(subset(SW_sub, gene_family=="pet"))),
              paste0("psa\n n = ",nrow(subset(SW_sub, gene_family=="psa"))),
              paste0("psb\n n = ",nrow(subset(SW_sub, gene_family=="psb")))
     )
)
object=aov(SW_sub$nucleotide_pearson_correlation~SW_sub$gene_family)
summary(object)
TukeyHSD(object)
mtext("GENE FAMILY EFFECT",side=2,padj=-5)


boxplot(SW_sub.F$nucleotide_pearson_correlation~SW_sub.F$gene_family,        
        col="lightgrey",
        main="Sliding Window Pearson correlation depending on \ngene's family within Fucoxanthin plastids",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("atp\n n = ",nrow(subset(SW_sub.F, gene_family=="atp"))),
              paste0("pet\n n = ",nrow(subset(SW_sub.F, gene_family=="pet"))),
              paste0("psa\n n = ",nrow(subset(SW_sub.F, gene_family=="psa"))),
              paste0("psb\n n = ",nrow(subset(SW_sub.F, gene_family=="psb")))
     )
)
object=aov(SW_sub.P$nucleotide_pearson_correlation~SW_sub.P$gene_family)
summary(object)
TukeyHSD(object)

boxplot(SW_sub.P$nucleotide_pearson_correlation~SW_sub.P$gene_family,
        col="lightgrey",
        main="Sliding Window Pearson correlation depending on \ngene's family within Peridinin plastids",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("atp\n n = ",nrow(subset(SW_sub.P, gene_family=="atp"))),
              paste0("pet\n n = ",nrow(subset(SW_sub.P, gene_family=="pet"))),
              paste0("psa\n n = ",nrow(subset(SW_sub.P, gene_family=="psa"))),
              paste0("psb\n n = ",nrow(subset(SW_sub.P, gene_family=="psb")))
     )
)
object=aov(SW_sub.F$nucleotide_pearson_correlation~SW_sub.F$gene_family)
summary(object)
TukeyHSD(object)

formule=subset(SW_sub, organism=="kv")$nucleotide_pearson_correlation~
  subset(SW_sub, organism=="kv")$gene_family
boxplot(formule,
        col="lightgrey",
        main="Sliding Window Pearson correlation depending on \ngene's family for KV",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("atp\n n = ",nrow(subset(subset(SW_sub, organism=="kv"), gene_family=="atp"))),
              paste0("pet\n n = ",nrow(subset(subset(SW_sub, organism=="kv"), gene_family=="pet"))),
              paste0("psa\n n = ",nrow(subset(subset(SW_sub, organism=="kv"), gene_family=="psa"))),
              paste0("psb\n n = ",nrow(subset(subset(SW_sub, organism=="kv"), gene_family=="psb")))
     )
)
object=aov(formule)
summary(object)
TukeyHSD(object)

formule=subset(SW_sub, organism=="km")$nucleotide_pearson_correlation~
  subset(SW_sub, organism=="km")$gene_family
boxplot(formule,
        col="lightgrey",
        main="Sliding Window Pearson correlation depending on \ngene's family for KM",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("atp\n n = ",nrow(subset(subset(SW_sub, organism=="km"), gene_family=="atp"))),
              paste0("pet\n n = ",nrow(subset(subset(SW_sub, organism=="km"), gene_family=="pet"))),
              paste0("psa\n n = ",nrow(subset(subset(SW_sub, organism=="km"), gene_family=="psa"))),
              paste0("psb\n n = ",nrow(subset(subset(SW_sub, organism=="km"), gene_family=="psb")))
     )
)
object=aov(formule)
summary(object)
TukeyHSD(object)

lines(x=c(1,2),y=c(0.4,0.4))
text(x=1.5, y=0.4, "*", pos=3, cex=1)

formule=subset(SW_sub, organism=="sm")$nucleotide_pearson_correlation~
  subset(SW_sub, organism=="sm")$gene_family
boxplot(formule,
        col="lightgrey",
        main="Sliding Window Pearson correlation depending on \ngene's family for SM",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("atp\n n = ",nrow(subset(subset(SW_sub, organism=="sm"), gene_family=="atp"))),
              paste0("pet\n n = ",nrow(subset(subset(SW_sub, organism=="sm"), gene_family=="pet"))),
              paste0("psa\n n = ",nrow(subset(subset(SW_sub, organism=="sm"), gene_family=="psa"))),
              paste0("psb\n n = ",nrow(subset(subset(SW_sub, organism=="sm"), gene_family=="psb")))
     )
)
object=aov(formule)
summary(object)
TukeyHSD(object)

formule=subset(SW_sub, organism=="pl")$nucleotide_pearson_correlation~
  subset(SW_sub, organism=="pl")$gene_family
boxplot(formule,
        col="lightgrey",
        main="Sliding Window Pearson correlation depending on \ngene's family for PL",
        ylab="Sliding Window Pearson correlation",
        ylim=c(-1,1),
        xaxt="n"
)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(paste0("atp\n n = ",nrow(subset(subset(SW_sub, organism=="pl"), gene_family=="atp"))),
              paste0("pet\n n = ",nrow(subset(subset(SW_sub, organism=="pl"), gene_family=="pet"))),
              paste0("psa\n n = ",nrow(subset(subset(SW_sub, organism=="pl"), gene_family=="psa"))),
              paste0("psb\n n = ",nrow(subset(subset(SW_sub, organism=="pl"), gene_family=="psb")))
     )
)
object=aov(formule)
summary(object)
TukeyHSD(object)

############ Gene Group ############

dev.off()

