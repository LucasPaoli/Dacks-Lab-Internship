##################################################
# Statistical Analysis of Dinos ARN editing Data #
##################################################

###########################
# SECOND ENTROPY ANALYSIS #
###########################

## We focus on 4 organisms : kv, km, sm, pl
## And 4 genes family, atp, psa, psb, pet

library(lattice)


###########
# ENTROPY #
###########

################
# LOADING DATA #
################

setwd("/Users/Lucas/dropbox/Chris/Very\ Final\ Stuff/Entropy/")

names=c("alignment.position","genomic.position","edited.or.not",
        "positional.entropy","score","plastid",
        "organism")

kv=read.delim(file="Kveneficum/Spreadsheets_score/kv_entropy.csv",sep=",")
kv=cbind(kv,rep("kv",length(kv$alignment.position)))
kv=cbind(kv,rep("kv",length(kv$alignment.position)))
kv=na.omit(kv)
names(kv)=names

km=read.delim(file="Kmikimotoi/Spreadsheets_score/km_entropy.csv",sep=",")
km=cbind(km,rep("km",length(km$alignment.position)))
km=cbind(km,rep("km",length(km$alignment.position)))
names(km)=names

pl=read.delim(file="Plunula/Spreadsheets_score/pl_entropy.csv",sep=",")
pl=cbind(pl,rep("pl",length(pl$alignment.position)))
pl=cbind(pl,rep("pl",length(pl$alignment.position)))
names(pl)=names

sm=read.delim(file="Sminutum/Spreadsheets_score/sm_entropy.csv",sep=",")
sm=cbind(sm,rep("sm",length(sm$alignment.position)))
sm=cbind(sm,rep("sm",length(sm$alignment.position)))
names(sm)=names

table=rbind(kv,km,pl,sm)

table.0=subset(table,table$edited.or.not==0)
table.1=subset(table,table$edited.or.not==1)

######################
# EDITS DISTRIBUTION #
######################

hist(table.0$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in non edited codons",
     xlab="Postitional Entropy",
     col="lightgrey")
hist(table.1$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons",
     xlab="Postitional Entropy",
     col="lightgrey")

ks.test(table.0$positional.entropy/sum(table.0$positional.entropy),table.1$positional.entropy/sum(table.1$positional.entropy))

## PLASTID EFFECT

table.1.P=subset(table,table$plastid=="pl" | table$plastid=="sm")
table.1.F=subset(table,table$plastid=="kv" | table$plastid=="km")
table.0.P=subset(table,table$plastid=="pl" | table$plastid=="sm")
table.0.F=subset(table,table$plastid=="kv" | table$plastid=="km")

hist(table.1.F$positional.entropy,breaks=seq(0,1,0.05))
hist(table.1.P$positional.entropy,breaks=seq(0,1,0.05))
hist(table.0.F$positional.entropy,breaks=seq(0,1,0.05))
hist(table.0.P$positional.entropy,breaks=seq(0,1,0.05))

## ORGANISM EFFECT

table.1.kv=subset(table,table$plastid=="kv")
table.1.km=subset(table,table$plastid=="km")
table.1.pl=subset(table,table$plastid=="pl")
table.1.sm=subset(table,table$plastid=="sm")
table.0.kv=subset(table,table$plastid=="kv")
table.0.km=subset(table,table$plastid=="km")
table.0.pl=subset(table,table$plastid=="pl")
table.0.sm=subset(table,table$plastid=="sm")

hist(table.1.kv$positional.entropy,breaks=seq(0,1,0.05))
hist(table.1.km$positional.entropy,breaks=seq(0,1,0.05))
hist(table.1.pl$positional.entropy,breaks=seq(0,1,0.05))
hist(table.1.sm$positional.entropy,breaks=seq(0,1,0.05))
hist(table.0.kv$positional.entropy,breaks=seq(0,1,0.05))
hist(table.0.km$positional.entropy,breaks=seq(0,1,0.05))
hist(table.0.pl$positional.entropy,breaks=seq(0,1,0.05))
hist(table.0.sm$positional.entropy,breaks=seq(0,1,0.05))

ks.test(table.1.kv$positional.entropy,table.1.km$positional.entropy)
ks.test(table.1.sm$positional.entropy,table.1.km$positional.entropy)
ks.test(table.1.pl$positional.entropy,table.1.km$positional.entropy)
ks.test(table.1.sm$positional.entropy,table.1.pl$positional.entropy)


# Figure

pdf(file = "Entropy_distribution_organisms.pdf", width=10, height=8)
par(mfcol=c(4,3))
hist(table.0$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy distribution in non edited codons",
     xlab="Postitional Entropy",
     col="lightgrey")
hist(table.1$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy distribution in edited codons",
     xlab="Postitional Entropy",
     col="lightgrey")

hist(table.0.kv$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy distribution in non edited codons\nfor KV",
     xlab="Postitional Entropy",
     col="lightgrey")
hist(table.1.kv$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy distribution in edited codons\nfor KV",
     xlab="Postitional Entropy",
     col="lightgrey")

hist(table.0.km$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy distribution in non edited codons\nfor KM",
     xlab="Postitional Entropy",
     col="lightgrey")
hist(table.1.km$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy distribution in edited codons\nfor KM",
     xlab="Postitional Entropy",
     col="lightgrey")

hist(table.0.pl$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy distribution in non edited codons\nfor PL",
     xlab="Postitional Entropy",
     col="lightgrey")
hist(table.1.pl$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy distribution in edited codons\nfor PL",
     xlab="Postitional Entropy",
     col="lightgrey")

hist(table.0.sm$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy distribution in non edited codons\nfor SM",
     xlab="Postitional Entropy",
     col="lightgrey")
hist(table.1.sm$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy distribution in edited codons\nfor SM",
     xlab="Postitional Entropy",
     col="lightgrey")

dev.off()


## INFLUENCE OF SCORE ON DISTRIBUTION

# BINARY

table.1.conserv=subset(table.1,table.1$score>0)
table.1.unconserv=subset(table.1,table.1$score<0)

hist(table.1.conserv$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency Editing Score diff > 0",
     xlab="Postitional Entropy",
     col="lightgrey")
hist(table.1.unconserv$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency Editing Score diff < 0",
     xlab="Postitional Entropy",
     col="lightgrey")

ks.test(table.1.conserv$positional.entropy,table.1.unconserv$positional.entropy)

summary(table.1$score)

# MORE PRECISE

table.1.max=subset(table.1,table.1$score>=5)
table.1.rest=subset(table.1,table.1$score<5)
table.1.5=subset(table.1,table.1$score<5 & table.1$score>0)
table.1_5=subset(table.1,table.1$score<=0 & table.1$score>-5)
table.1.min=subset(table.1,table.1$score<=-5)

hist(table.1.min$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons",
     xlab="Postitional Entropy",
     col="lightgrey")
hist(table.1_5$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons",
     xlab="Postitional Entropy",
     col="lightgrey")
hist(table.1.5$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons",
     xlab="Postitional Entropy",
     col="lightgrey")
hist(table.1.max$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons",
     xlab="Postitional Entropy",
     col="lightgrey")
hist(table.1.rest$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons",
     xlab="Postitional Entropy",
     col="lightgrey")

ks.test(table.1.max$positional.entropy,table.1.rest$positional.entropy)

boxplot(table.1.max$positional.entropy,
        table.1.5$positional.entropy,
        table.1_5$positional.entropy,
        table.1.min$positional.entropy,
        names=c("Max scores", "slight positive","slight negative","max divergence"),
        main="Boxplot of positional entropy",
        ylab="positional entropy",
        col="lightgrey",
        las=2)

# Figure

pdf(file = "Entropy_distribution_score.pdf", width=14, height=7)
par(mfcol=c(2,4))

hist(table.1.conserv$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons for\nEditing Score difference > 0",
     xlab="Postitional Entropy",
     col="lightgrey")
hist(table.1.unconserv$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons for\nEditing Score difference < 0",
     xlab="Postitional Entropy",
     col="lightgrey")

hist(table.1.min$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons with",
     xlab="Postitional Entropy",
     col="lightgrey")
mtext(expression("Diff"<="-5"))
hist(table.1_5$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons with",
     xlab="Postitional Entropy",
     col="lightgrey")
mtext(expression("-5 < Diff" <= "0"))
hist(table.1.5$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons with",
     xlab="Postitional Entropy",
     col="lightgrey")
mtext(expression("0 < Diff < 5"))
hist(table.1.max$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons with",
     xlab="Postitional Entropy",
     col="lightgrey")
mtext(expression("Diff">="5"))

boxplot(table.1.max$positional.entropy,
        table.1.5$positional.entropy,
        table.1_5$positional.entropy,
        table.1.min$positional.entropy,
        ylim=c(0.2,1.1),
        main="Boxplot of positional entropy depending on\nEditing Score differences",
        ylab="Positional Entropy",
        col="lightgrey",
        las=2,
        xaxt="n")
text(x=1, y=1, "***", pos=3, cex=1)
text(x=4, y=1, "***", pos=3, cex=1)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(expression("Diff">="5"),
              expression("0 < Diff" < "5"),
              expression("-5 < Diff" <= "0"),
              expression("Diff"<="-5"))
)
axis(side=1,
     padj=2,
     at=1:4,
     lwd=0,
     lwd.ticks=0,
     labels=c(paste0("n = ",nrow(table.1.max)),
              paste0("n = ",nrow(table.1.5)),
              paste0("n = ",nrow(table.1_5)),
              paste0("n = ",nrow(table.1.min))
     )
)
dev.off()

shapiro.test(table.1.max$positional.entropy)
shapiro.test(table.1.5$positional.entropy)
shapiro.test(table.1.min$positional.entropy)
shapiro.test(table.1_5$positional.entropy)

t.test(table.1.max$positional.entropy,
       table.1.5$positional.entropy)
t.test(table.1.min$positional.entropy,
       table.1_5$positional.entropy)
t.test(table.1.5$positional.entropy,
       table.1_5$positional.entropy)
t.test(table.1.min$positional.entropy,
       table.1.max$positional.entropy)

# PLASTID EFFECT
table.1.P=subset(table.1,table.1$plastid=="pl" 
                 | table.1$plastid=="sm")
table.1.F=subset(table.1,table.1$plastid=="kv" 
                 | table.1$plastid=="km")

table.1.P.max=subset(table.1.P,table.1.P$score>=5)
table.1.P.rest=subset(table.1.P,table.1.P$score<5)
table.1.P.5=subset(table.1.P,table.1.P$score<5 & table.1.P$score>0)
table.1.P_5=subset(table.1.P,table.1.P$score<=0 & table.1.P$score>-5)
table.1.P.min=subset(table.1.P,table.1.P$score<=-5)

table.1.F.max=subset(table.1.F,table.1.F$score>=5)
table.1.F.rest=subset(table.1.F,table.1.F$score<5)
table.1.F.5=subset(table.1.F,table.1.F$score<5 & table.1.F$score>0)
table.1.F_5=subset(table.1.F,table.1.F$score<=0 & table.1.F$score>-5)
table.1.F.min=subset(table.1.F,table.1.F$score<=-5)

pdf(file = "Entropy_distribution_score_plastid.pdf", width=23, height=8)
par(mfrow=c(2,5),mar=c(6,6,6,6))

hist(table.1.P.min$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons with",
     xlab="Postitional Entropy",
     col="lightgrey")
mtext(expression("Diff"<="-5"))
hist(table.1.P_5$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons with",
     xlab="Postitional Entropy",
     col="lightgrey")
mtext(expression("-5 < Diff" <= "0"))
hist(table.1.P.5$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons with",
     xlab="Postitional Entropy",
     col="lightgrey")
mtext(expression("0 < Diff < 5"))
hist(table.1.P.max$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons with",
     xlab="Postitional Entropy",
     col="lightgrey")
mtext(expression("Diff">="5"))

boxplot(table.1.P.max$positional.entropy,
        table.1.P.5$positional.entropy,
        table.1.P_5$positional.entropy,
        table.1.P.min$positional.entropy,
        ylim=c(0.2,1.1),
        main="Boxplot of positional entropy depending on\nEditing Score differences for Peridinin plastids",
        ylab="Positional Entropy",
        col="lightgrey",
        xaxt="n")
text(x=1, y=1, "***", pos=3, cex=1)
text(x=4, y=1, "**", pos=3, cex=1)
mtext("PERIDININ PLASTIDS",side=4,padj=2)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(expression("Diff">="5"),
              expression("0 < Diff" < "5"),
              expression("-5 < Diff" <= "0"),
              expression("Diff"<="-5"))
)
axis(side=1,
     padj=2,
     at=1:4,
     lwd=0,
     lwd.ticks=0,
     labels=c(paste0("n = ",nrow(table.1.P.max)),
              paste0("n = ",nrow(table.1.P.5)),
              paste0("n = ",nrow(table.1.P_5)),
              paste0("n = ",nrow(table.1.P.min))
     )
)

hist(table.1.F.min$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons with",
     xlab="Postitional Entropy",
     col="lightgrey")
mtext(expression("Diff"<="-5"))
hist(table.1.F_5$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons with",
     xlab="Postitional Entropy",
     col="lightgrey")
mtext(expression("-5 < Diff" <= "0"))
hist(table.1.F.5$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons with",
     xlab="Postitional Entropy",
     col="lightgrey")
mtext(expression("0 < Diff < 5"))
hist(table.1.F.max$positional.entropy,breaks=seq(0,1,0.05),
     main="Entropy frequency in edited codons with",
     xlab="Postitional Entropy",
     col="lightgrey")
mtext(expression("Diff">="5"))
boxplot(table.1.F.max$positional.entropy,
        table.1.F.5$positional.entropy,
        table.1.F_5$positional.entropy,
        table.1.F.min$positional.entropy,
        ylim=c(0.2,1.1),
        main="Boxplot of positional entropy depending on\nEditing Score differences for Fucoxathin plastids",
        ylab="Positional Entropy",
        col="lightgrey",
        xaxt="n")
text(x=1, y=1, "***", pos=3, cex=1)
text(x=4, y=1, "***", pos=3, cex=1)
mtext("FUCOXANTHIN PLASTIDS",side=4,padj=2)
axis(side=1,
     padj=0.5,
     at=1:4,
     labels=c(expression("Diff">="5"),
              expression("0 < Diff" < "5"),
              expression("-5 < Diff" <= "0"),
              expression("Diff"<="-5"))
)
axis(side=1,
     padj=2,
     at=1:4,
     lwd=0,
     lwd.ticks=0,
     labels=c(paste0("n = ",nrow(table.1.F.max)),
              paste0("n = ",nrow(table.1.F.5)),
              paste0("n = ",nrow(table.1.F_5)),
              paste0("n = ",nrow(table.1.F.min))
     )
)


dev.off()

shapiro.test(table.1.P.max$positional.entropy)
shapiro.test(table.1.P.5$positional.entropy)
shapiro.test(table.1.P.min$positional.entropy)
shapiro.test(table.1.P_5$positional.entropy)

shapiro.test(table.1.F.max$positional.entropy)
shapiro.test(table.1.F.5$positional.entropy)
shapiro.test(table.1.F.min$positional.entropy)
shapiro.test(table.1.F_5$positional.entropy)


t.test(table.1.P.max$positional.entropy,
       table.1.P.5$positional.entropy)
t.test(table.1.P.min$positional.entropy,
       table.1.P_5$positional.entropy)
t.test(table.1.P.5$positional.entropy,
       table.1.P_5$positional.entropy)
t.test(table.1.P.min$positional.entropy,
       table.1.P.max$positional.entropy)

t.test(table.1.F.max$positional.entropy,
       table.1.F.5$positional.entropy)
t.test(table.1.F.min$positional.entropy,
       table.1.F_5$positional.entropy)
t.test(table.1.F.5$positional.entropy,
       table.1.F_5$positional.entropy)
t.test(table.1.F.min$positional.entropy,
       table.1.F.max$positional.entropy)

##############################################
# INFLUENCE OF POS. ENTROPY ON EDITING SCORE #
##############################################

## OVERALL

table.1.pos.entropy=subset(table.1,table.1$positional.entropy > 0.95)
table.1.neg.entropy=subset(table.1,table.1$positional.entropy <= 0.95)

boxplot(table.1.pos.entropy$score,
        table.1.neg.entropy$score,
        main="Boxplot of positional entropy",
        names=c("Entropy > 0.95","Entropy > 0.95"),
        ylab="Editing Score difference",
        col="lightgrey",
        las=2)

shapiro.test(table.1.neg.entropy$score)
shapiro.test(table.1.pos.entropy$score)

t.test(table.1.pos.entropy$score,table.1.neg.entropy$score)

## PLASTID EFFECT

table.1.pos.entropy.P=subset(table.1.pos.entropy,table.1.pos.entropy$plastid=="pl" 
                             | table.1.pos.entropy$plastid=="sm")
table.1.pos.entropy.F=subset(table.1.pos.entropy,table.1.pos.entropy$plastid=="kv" 
                             | table.1.pos.entropy$plastid=="km")
table.1.neg.entropy.P=subset(table.1.neg.entropy,table.1.neg.entropy$plastid=="pl" 
                             | table.1.neg.entropy$plastid=="sm")
table.1.neg.entropy.F=subset(table.1.neg.entropy,table.1.neg.entropy$plastid=="kv" 
                             | table.1.neg.entropy$plastid=="km")

boxplot(table.1.pos.entropy.P$score,
        table.1.pos.entropy.F$score,
        table.1.neg.entropy.P$score,
        table.1.neg.entropy.F$score,
        main="Boxplot of editing score according to the positional entropy",
        names=c("peridinin\nentropy > 0.95","Fucoxanthin\nentropy > 0.95","peridinin\nentropy ≤ 0.95","Fucoxanthin\nentropy ≤ 0.95"),
        ylab="Editing Score difference",
        col="lightgrey",
        ylim=c(-15,20))
lines(x=c(1,3),y=c(13,13))
text(x=2, y=10, "***", pos=3, cex=1)
lines(x=c(1,2),y=c(17,17))
text(x=1.5, y=14, "***", pos=3, cex=1)

shapiro.test(table.1.pos.entropy.P$score)
shapiro.test(table.1.pos.entropy.F$score)
shapiro.test(table.1.neg.entropy.P$score)
shapiro.test(table.1.neg.entropy.F$score)

t.test(table.1.pos.entropy.F$score,table.1.neg.entropy.F$score)
t.test(table.1.pos.entropy.P$score,table.1.neg.entropy.P$score)
t.test(table.1.pos.entropy.P$score,table.1.pos.entropy.F$score)
t.test(table.1.neg.entropy.P$score,table.1.neg.entropy.F$score)

# Figure

pdf(file = "Entropy_Editing_score_plastids.pdf", width=10, height=8)
par(mar=c(9,4,6,4))
boxplot(table.1.pos.entropy.P$score,
        table.1.pos.entropy.F$score,
        table.1.neg.entropy.P$score,
        table.1.neg.entropy.F$score,
        main="Boxplot of editing score difference according to the positional entropy",
        ylab="Editing Score difference",
        col="lightgrey",
        ylim=c(-15,20),
        xaxt="n")
lines(x=c(1,3),y=c(13,13))
text(x=2, y=12.5, "***", pos=3, cex=1)
lines(x=c(1,2),y=c(17,17))
text(x=1.5, y=16.5, "***", pos=3, cex=1)

axis(side=1,
     padj=1,
     at=1:4,
     labels=c(expression("Peridinin\nentropy">"0.95"),
              expression("Fucoxanthin\nentropy">"0.95"),
              expression("Peridinin\nentropy"<="0.95"),
              expression("Fucoxanthin\nentropy"<="0.95"))
)
axis(side=1,
     padj=3,
     at=1:4,
     lwd=0,
     lwd.ticks=0,
     labels=c(paste0("n = ",nrow(table.1.pos.entropy.P)),
              paste0("n = ",nrow(table.1.pos.entropy.F)),
              paste0("n = ",nrow(table.1.neg.entropy.P)),
              paste0("n = ",nrow(table.1.neg.entropy.F))
     )
)

legend("topright", # places a legend at the appropriate place
       c("***  t-test p-value < 0.001"), # puts text in the legend
)
text(1, summary(table.1.pos.entropy.P$score)["Mean"]+1,labels=summary(table.1.pos.entropy.P$score)["Mean"])
text(2, summary(table.1.pos.entropy.F$score)["Mean"]+2,labels=summary(table.1.pos.entropy.F$score)["Mean"])
text(3, summary(table.1.neg.entropy.P$score)["Mean"]+0.7,labels=summary(table.1.neg.entropy.P$score)["Mean"])
text(4, summary(table.1.neg.entropy.F$score)["Mean"]+0.8,labels=summary(table.1.neg.entropy.F$score)["Mean"])

dev.off()

