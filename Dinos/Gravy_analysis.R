##################################################
# Statistical Analysis of Dinos ARN editing Data #
##################################################

## We focus on 4 organisms : kv, km, sm, pl
## And 4 genes family, atp, psa, psb, pet

library(lattice)


########################
# GRAVY SCORE ANALYSIS #
########################


names=c("gene","weight","gravy",
        "rna.mol.weight","rna.gravy","plastid",
        "organism")
  
  setwd("/Users/Lucas/dropbox/Chris/Very\ Final\ Stuff/Gravy_Scores")
  
  kv=read.delim(file="Kveneficum/Kv_gravy.csv",sep=",")
  kv=cbind(kv,rep("kv",length(kv$gene)))
  kv=cbind(kv,rep("kv",length(kv$gene)))
  names(kv)=names

  km=read.delim(file="Kmikimotoi/Km_gravy.csv",sep=",")
  km=cbind(km,rep("km",length(km$gene)))
  km=cbind(km,rep("km",length(km$gene)))
  names(km)=names
  
  pl=read.delim(file="Plunula/Pl_gravy.csv",sep=",")
  pl=cbind(pl,rep("pl",length(pl$gene)))
  pl=cbind(pl,rep("pl",length(pl$gene)))
  names(pl)=names
  
  sm=read.delim(file="Sminutum/Sm_gravy.csv",sep=",")
  sm=cbind(sm,rep("sm",length(sm$gene)))
  sm=cbind(sm,rep("sm",length(sm$gene)))
  names(sm)=names
  
  table_trans=rbind(kv,km,pl,sm)
  
  table=cbind(table_trans[,1],table_trans[,1:length(table_trans)])
  names(table)=c("gene","gene_family","weight","gravy",
                 "rna.mol.weight","rna.gravy","plastid",
                 "organism")
  
  table[,2]=gsub(".*atp.*","atp",table[,2])
  table[,2]=gsub(".*psa.*","psa",table[,2])
  table[,2]=gsub(".*psb.*","psb",table[,2])
  table[,2]=gsub(".*pet.*","pet",table[,2])
  table[,2]=gsub(".*rpl.*","rpl",table[,2])
  table[,2]=gsub(".*rps.*","rps",table[,2])
  table[,2]=gsub(".*ycf.*","ycf",table[,2])
  table[,2]=gsub(".*rbc.*","rbc",table[,2])
  table[,2]=gsub(".*rpo.*","rpo",table[,2])
  table[,2]=gsub(".*sec.*","sec",table[,2])
  
  table[,length(table)-1]=gsub("kv","fucoxanthin",table[,length(table)-1])
  table[,length(table)-1]=gsub("km","fucoxanthin",table[,length(table)-1])
  table[,length(table)-1]=gsub("sm","peridinin",table[,length(table)-1])
  table[,length(table)-1]=gsub("pl","peridinin",table[,length(table)-1])
  
  table[,length(table)+1]=table[,"rna.mol.weight"]-table[,"weight"]
  table[,length(table)+1]=table[,"rna.gravy"]-table[,"gravy"]
  
 names(table)=c("gene","gene_family","weight","gravy",
               "rna.mol.weight","rna.gravy","plastid",
               "organism","weight","gravy")

  table_sub=subset(table,gene_family== "psa" | gene_family== "psb" | gene_family== "atp" | gene_family== "pet")

## CORRELATION RNA

cor.test(table[,"gen.mol.weight"],table[,"rna.mol.weight"])
cor.test(table[,"gen.gravy"],table[,"rna.gravy"])

xyplot(table[,"gen.gravy"]~table[,"rna.gravy"],groups=table$plastid,pch=16,col=c("green","blue"))
xyplot(table[,"gen.mol.weight"]~table[,"rna.mol.weight"],groups=table$plastid,pch=16,col=c("green","blue"))


# ALL

cor.test(table[,"weight"],table[,"gravy"])
xyplot(table[,"weight"]~table[,"gravy"],groups=table$plastid,pch=16,col=c("green","blue"))

## PHOTO GENES ONLY

table_sub=subset(table,gene_family== "psa" | gene_family== "psb" | gene_family== "atp" | gene_family== "pet")
cor.test(table_sub[,"weight"],table_sub[,"gravy"])
xyplot(table_sub[,"weight"]~table_sub[,"gravy"],groups=table_sub$plastid,pch=16,col=c("green","blue"))

pdf(file = "correlation_gravy_score_general.pdf", width=16, height=8)
par(mfcol=c(1,2))
plot1=xyplot(table[,"weight"]~table[,"gravy"],
       groups=table$plastid,
       pch=16,col=c("green3","blue"),
       main="Correlation for all the genes",
       ylab="Molecular weight",
       xlab="Gravy score",
       sub="Pearson corr. coeff. = 0.06120272 , p-value = 0.5253",
       auto.key=list(title="Plastid :",column=2,cex=0.8,points=F,col=c("green3","blue"))
)

plot2=xyplot(table_sub[,"weight"]~table_sub[,"gravy"],
       groups=table_sub$plastid,
       pch=16,
       col=c("green3","blue"),
       main="Correlation for conserved genes (photosynthetic genes)",
       ylab="Molecular weight",
       xlab="Gravy score",
       sub="Pearson corr. coeff. = 0.0496924 , p-value = 0.7037",
       auto.key=list(title="Plastid :",column=2,cex=0.8,points=F,col=c("green3","blue"))
)

print(plot1, position = c(0, 0, 0.5, 1), more = TRUE)
print(plot2, position = c(0.5, 0, 1, 1))

dev.off()

# PERIDININ

table.P=subset(table, organism== "sm" | organism== "pl")
cor.test(table.P[,"weight"],table.P[,"gravy"])
xyplot(table.P[,"weight"]~table.P[,"gravy"],groups=table.P$gene_family,pch=16)

pdf(file = "correlation_gravy_score_peridinin.pdf", width=8, height=8)
xyplot(table.P[,"weight"]~table.P[,"gravy"],groups=table.P$gene_family,
       main="Correlation for Peridinin plastids",
       ylab="Molecular weight",
       xlab="Gravy score",
       sub="Pearson corr. coeff. = 0.2275462, p-value = 0.2964",
       par.settings = list(superpose.symbol = list(pch = 16)),
       auto.key=list(title="Gene family :",column=4,cex=0.8)
)
dev.off()

# FUCOXANTHIN

table.F=subset(table, organism== "kv" | organism== "km")
cor.test(table.F[,"weight"],table.F[,"gravy"])
xyplot(table.F[,"weight"]~table.F[,"gravy"],groups=table.F$gene_family,pch=16)

table_sub.F=subset(table_sub, organism== "kv" | organism== "km")
cor.test(table_sub.F[,"weight"],table_sub.F[,"gravy"])
xyplot(table_sub.F[,"weight"]~table_sub.F[,"gravy"],groups=table_sub.F$gene_family,pch=16)

pdf(file = "correlation_gravy_score_photo_fucoxanthin.pdf", width=8, height=8)
xyplot(table_sub.F[,"weight"]~table_sub.F[,"gravy"],groups=table_sub.F$gene_family,
       main="Correlation for Fucoxanthin plastids, photosynthetic genes only",
       ylab="Molecular weight",
       xlab="Gravy score",
       sub="Pearson corr. coeff. = 0.03263602, p-value = 0.7641",
       par.settings = list(superpose.symbol = list(pch = 16)),
       auto.key=list(title="Gene family :",column=4,cex=0.8)
)
dev.off()

table.sub=subset(table,gene_family!= "psa" & gene_family!= "psb" & gene_family!= "atp" & gene_family!= "pet")
cor.test(table.sub[,"weight"],table.sub[,"gravy"])
xyplot(table.sub[,"weight"]~table.sub[,"gravy"],groups=table.sub$gene_family,pch=16)

pdf(file = "correlation_gravy_score_noncvgene_fucoxanthin.pdf", width=8, height=8)
xyplot(table.sub[,"weight"]~table.sub[,"gravy"],groups=table.sub$gene_family,
       main=c("Correlation for Fucoxanthin plastids, housekeeping genes"),
       ylab="Molecular weight",
       xlab="Gravy score",
       sub="Pearson corr. coeff. = -0.03024132, p-value = 0.857",
       par.settings = list(superpose.symbol = list(pch = 16)),
       auto.key=list(title="Gene family :",column=7,cex=0.8)
)
dev.off()

## PLASTID FIGURE

pdf(file = "correlation_gravy_score_plastids.pdf", width=24, height=8)
peridinin = xyplot(table.P[,"weight"]~table.P[,"gravy"],groups=table.P$gene_family,
                   main="Correlation for Peridinin plastids",
                   ylab="Molecular weight",
                   xlab="Gravy score",
                   sub=" Pearson corr. coeff. = 0.2275462, p-value = 0.2964",
                   par.settings = list(superpose.symbol = list(pch = 16)),
                   auto.key=list(title="Gene family :",column=4,cex=0.8)
)
fuco1=xyplot(table_sub.F[,"weight"]~table_sub.F[,"gravy"],groups=table_sub.F$gene_family,
             main="Correlation for Fucoxanthin plastids, photosynthetic genes only",
             ylab="Molecular weight",
             xlab="Gravy score",
             sub="Pearson corr. coeff. = 0.03263602, p-value = 0.7641",
             par.settings = list(superpose.symbol = list(pch = 16)),
             auto.key=list(title="Gene family :",column=4,cex=0.8)
)
fuco2=xyplot(table.sub[,"weight"]~table.sub[,"gravy"],groups=table.sub$gene_family,
             main=c("Correlation for Fucoxanthin plastids, housekeeping genes only"),
             ylab="Molecular weight",
             xlab="Gravy score",
             sub="Pearson corr. coeff. = -0.03024132, p-value = 0.857",
             par.settings = list(superpose.symbol = list(pch = 16)),
             auto.key=list(title="Gene family :",column=7,cex=0.8)
)
print(peridinin, position = c(0, 0, 0.33, 1), more = TRUE)
print(fuco1, position = c(0.33, 0, 0.66, 1), more = TRUE)
print(fuco2, position = c(0.66, 0, 1, 1))
dev.off()


## ORGANISM

pdf(file = "correlation_gravy_score_organism.pdf", width=24, height=16)

cor.test(subset(table.sub,organism=="kv")$weight,
         subset(table.sub,organism=="kv")$gravy)
kv1=xyplot(subset(table.sub,organism=="kv")$weight~
         subset(table.sub,organism=="kv")$gravy,
       groups=subset(table.sub,organism=="kv")$gene_family,
       main="Correlation for KV, housekeeping genes",
       ylab="Molecular weight",
       xlab="Gravy score",
       sub="Pearson corr. coeff. = 0.2200501, p-value = 0.1783",
       par.settings = list(superpose.symbol = list(pch = 16)),
       auto.key=list(title="Gene family :",column=4,cex=0.8),
       cex=1.2
)
cor.test(subset(table_sub,organism=="kv")$weight,
         subset(table_sub,organism=="kv")$gravy)
kv2=xyplot(subset(table_sub,organism=="kv")$weight~
             subset(table_sub,organism=="kv")$gravy,
           groups=subset(table_sub,organism=="kv")$gene_family,
           main="Correlation for KV, photosynthetic genes",
           ylab="Molecular weight",
           xlab="Gravy score",
           sub="Pearson corr. coeff. = -0.1065756, p-value = 0.6284",
           par.settings = list(superpose.symbol = list(pch = 16)),
           auto.key=list(title="Gene family :",column=4,cex=0.8),
           cex=1.2
)

cor.test(subset(table.sub,organism=="km")$weight,
         subset(table.sub,organism=="km")$gravy)
km1=xyplot(subset(table.sub,organism=="km")$weight~
         subset(table.sub,organism=="km")$gravy,
       groups=subset(table.sub,organism=="km")$gene_family,
       main="Correlation for KM, housekeeping genes",
       ylab="Molecular weight",
       xlab="Gravy score",
       sub="Pearson corr. coeff. = -0.3116398, p-value = 0.3807",
       par.settings = list(superpose.symbol = list(pch = 16)),
       auto.key=list(title="Gene family :",column=4,cex=0.8),
       cex=1.2
)

cor.test(subset(table_sub,organism=="km")$weight,
         subset(table_sub,organism=="km")$gravy)
km2=xyplot(subset(table_sub,organism=="km")$weight~
             subset(table_sub,organism=="km")$gravy,
           groups=subset(table_sub,organism=="km")$gene_family,
           main="Correlation for KM, photosynthetic genes",
           ylab="Molecular weight",
           xlab="Gravy score",
           sub="Pearson corr. coeff. = 0.007388237, p-value = 0.9792",
           par.settings = list(superpose.symbol = list(pch = 16)),
           auto.key=list(title="Gene family :",column=4,cex=0.8),
           cex=1.2
)

cor.test(subset(table,organism=="sm")$weight,
         subset(table,organism=="sm")$gravy)
sm1=xyplot(subset(table,organism=="sm")$weight~
         subset(table,organism=="sm")$gravy,
       groups=subset(table,organism=="sm")$gene_family,
       main="Correlation for SM, photosynthetic genes",
       ylab="Molecular weight",
       xlab="Gravy score",
       sub="Pearson corr. coeff. = -0.6849792, p-value = 0.01397",
       par.settings = list(superpose.symbol = list(pch = 16)),
       auto.key=list(title="Gene family :",column=4,cex=0.8),
       cex=1.2
)

cor.test(subset(table,organism=="pl")$weight,
         subset(table,organism=="pl")$gravy)
pl1=xyplot(subset(table,organism=="pl")$weight~
         subset(table,organism=="pl")$gravy,
       groups=subset(table,organism=="pl")$gene_family,
       main="Correlation for PL, photosynthetic genes",
       ylab="Molecular weight",
       xlab="Gravy score",
       sub="Pearson corr. coeff. = 0.6499008, p-value = 0.03042",
       par.settings = list(superpose.symbol = list(pch = 16)),
       auto.key=list(title="Gene family :",column=4,cex=0.8),
       cex=1.2
)
print(kv1, position = c(0, 0.52, 0.33, 1), more = TRUE)
print(km1, position = c(0, 0, 0.33, 0.48), more = TRUE)
print(kv2, position = c(0.33, 0.52, 0.66, 1), more = TRUE)
print(km2, position = c(0.33, 0, 0.66, 0.48), more = TRUE)
print(sm1, position = c(0.66, 0.52, 1, 1), more = TRUE)
print(pl1, position = c(0.66, 0, 1, 0.48))
dev.off()




