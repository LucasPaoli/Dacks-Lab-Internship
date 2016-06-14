##############
# GC CONTENT #
##############

setwd("/Users/Lucas/Dropbox/Chris/Very Final Stuff/Pipeline_CSV_files_for_Lucas/Position_specific_GC_outfiles")
table.GC=read.delim(file="Global_GC_position.csv")
setwd("/Users/Lucas/Dropbox/Chris/Very Final Stuff/Pipeline_CSV_files_for_Lucas/Position_specific_edit_type_outfiles")
table.types=read.delim(file="Global_post_editing_types.csv")



pdf(file = "GC_barplot.pdf", width=24, height=12)
par(mfrow=c(2,4),mar=c(8,6,4,4)+1)

increase.first=mean(subset(table.types, position=="first")$A.to.C)+
  mean(subset(table.types, position=="first")$A.to.G)+
  mean(subset(table.types, position=="first")$T.to.C)+
  mean(subset(table.types, position=="first")$T.to.G)-
  mean(subset(table.types, position=="first")$C.to.A)-
  mean(subset(table.types, position=="first")$C.to.T)-
  mean(subset(table.types, position=="first")$G.to.A)-
  mean(subset(table.types, position=="first")$G.to.T)
increase.second=mean(subset(table.types, position=="second")$A.to.C)+
  mean(subset(table.types, position=="second")$A.to.G)+
  mean(subset(table.types, position=="second")$T.to.C)+
  mean(subset(table.types, position=="second")$T.to.G)-
  mean(subset(table.types, position=="second")$C.to.A)-
  mean(subset(table.types, position=="second")$C.to.T)-
  mean(subset(table.types, position=="second")$G.to.A)-
  mean(subset(table.types, position=="second")$G.to.T)
increase.third=mean(subset(table.types, position=="third")$A.to.C)+
  mean(subset(table.types, position=="third")$A.to.G)+
  mean(subset(table.types, position=="third")$T.to.C)+
  mean(subset(table.types, position=="third")$T.to.G)-
  mean(subset(table.types, position=="third")$C.to.A)-
  mean(subset(table.types, position=="third")$C.to.T)-
  mean(subset(table.types, position=="third")$G.to.A)-
  mean(subset(table.types, position=="third")$G.to.T)

barplot.pos=barplot(c(mean(table.GC$X1st.GC.before),
                      mean(table.GC$X1st.GC.after),
                      increase.first,
                      mean(table.GC$X2nd.GC.before),
                      mean(table.GC$X2nd.GC.after),
                      increase.second,
                      mean(table.GC$X3rd.GC.before),
                      mean(table.GC$X3rd.GC.after),
                      increase.third),
                    names=rep(c("%GC before","%GC after","GC increase"),3),
                    las=2,
                    ylim=c(-5,70),
                    col=rep(c("grey50","grey50","steelblue"),3),
                    axes=F,
                    main="Editing effect on the GC content \nper position")

axis(side=2,col="grey50",lwd=2)
axis(side=4,col="steelblue",lwd=2)
lines(x=c(barplot.pos[1],barplot.pos[3]),y=c(60,60))
text(x=barplot.pos[2], y=60, "1st Pos", pos=3, cex=1)
lines(x=c(barplot.pos[4],barplot.pos[6]),y=c(50,50))
text(x=barplot.pos[5], y=50, "2nd Pos", pos=3, cex=1)
lines(x=c(barplot.pos[7],barplot.pos[9]),y=c(40,40))
text(x=barplot.pos[8], y=40, "3rd Pos", pos=3, cex=1)

sd.first=sqrt(var(subset(table.types, position=="first")$A.to.C)+
                var(subset(table.types, position=="first")$A.to.G)+
                var(subset(table.types, position=="first")$T.to.C)+
                var(subset(table.types, position=="first")$T.to.G)+
                var(subset(table.types, position=="first")$C.to.A)+
                var(subset(table.types, position=="first")$C.to.T)+
                var(subset(table.types, position=="first")$G.to.A)+
                var(subset(table.types, position=="first")$G.to.T))/
  sqrt(nrow(subset(table.types, position=="first")))
sd.second=sqrt(var(subset(table.types, position=="second")$A.to.C)+
                 var(subset(table.types, position=="second")$A.to.G)+
                 var(subset(table.types, position=="second")$T.to.C)+
                 var(subset(table.types, position=="second")$T.to.G)+
                 var(subset(table.types, position=="second")$C.to.A)+
                 var(subset(table.types, position=="second")$C.to.T)+
                 var(subset(table.types, position=="second")$G.to.A)+
                 var(subset(table.types, position=="second")$G.to.T))/
  sqrt(nrow(subset(table.types, position=="second")))
sd.third=sqrt(var(subset(table.types, position=="third")$A.to.C)+
                var(subset(table.types, position=="third")$A.to.G)+
                var(subset(table.types, position=="third")$T.to.C)+
                var(subset(table.types, position=="third")$T.to.G)+
                var(subset(table.types, position=="third")$C.to.A)+
                var(subset(table.types, position=="third")$C.to.T)+
                var(subset(table.types, position=="third")$G.to.A)+
                var(subset(table.types, position=="third")$G.to.T))/
  sqrt(nrow(subset(table.types, position=="third")))


pos.int=c(mean(table.GC$X1st.GC.before),
          mean(table.GC$X1st.GC.after),
          increase.first,
          mean(table.GC$X2nd.GC.before),
          mean(table.GC$X2nd.GC.after),
          increase.second,
          mean(table.GC$X3rd.GC.before),
          mean(table.GC$X3rd.GC.after),
          increase.third)

conf.int=c(sd(table.GC$X1st.GC.before)/sqrt(nrow(table.GC)),
           sd(table.GC$X1st.GC.after)/sqrt(nrow(table.GC)),
           sd.first,
           sd(table.GC$X2nd.GC.before)/sqrt(nrow(table.GC)),
           sd(table.GC$X2nd.GC.after)/sqrt(nrow(table.GC)),
           sd.second,
           sd(table.GC$X3rd.GC.before)/sqrt(nrow(table.GC)),
           sd(table.GC$X3rd.GC.after)/sqrt(nrow(table.GC)),
           sd.third)*1.96

segments(barplot.pos,y0=pos.int-conf.int,y1=pos.int+conf.int,lwd=1)
arrows(barplot.pos, pos.int-conf.int, 
       barplot.pos, pos.int+conf.int, 
       lwd = 1, angle = 90,
       code = 3, length = 0.05)

mtext(side=4,text="Number of edits increasing GC content",col="steelblue",padj=4)
mtext(side=2,text="%GC content",col="grey50",padj=-4)


liste=c("ch","ht","km","kv","lp","pl","sm")

for (i in liste){

increase.first=mean(subset(table.types,organism==i & position=="first")$A.to.C)+
  mean(subset(table.types,organism==i & position=="first")$A.to.G)+
  mean(subset(table.types,organism==i & position=="first")$T.to.C)+
  mean(subset(table.types,organism==i & position=="first")$T.to.G)-
  mean(subset(table.types,organism==i & position=="first")$C.to.A)-
  mean(subset(table.types,organism==i & position=="first")$C.to.T)-
  mean(subset(table.types,organism==i & position=="first")$G.to.A)-
  mean(subset(table.types,organism==i & position=="first")$G.to.T)
increase.second=mean(subset(table.types,organism==i & position=="second")$A.to.C)+
  mean(subset(table.types,organism==i & position=="second")$A.to.G)+
  mean(subset(table.types,organism==i & position=="second")$T.to.C)+
  mean(subset(table.types,organism==i & position=="second")$T.to.G)-
  mean(subset(table.types,organism==i & position=="second")$C.to.A)-
  mean(subset(table.types,organism==i & position=="second")$C.to.T)-
  mean(subset(table.types,organism==i & position=="second")$G.to.A)-
  mean(subset(table.types,organism==i & position=="second")$G.to.T)
increase.third=mean(subset(table.types,organism==i & position=="third")$A.to.C)+
  mean(subset(table.types,organism==i & position=="third")$A.to.G)+
  mean(subset(table.types,organism==i & position=="third")$T.to.C)+
  mean(subset(table.types,organism==i & position=="third")$T.to.G)-
  mean(subset(table.types,organism==i & position=="third")$C.to.A)-
  mean(subset(table.types,organism==i & position=="third")$C.to.T)-
  mean(subset(table.types,organism==i & position=="third")$G.to.A)-
  mean(subset(table.types,organism==i & position=="third")$G.to.T)

barplot.pos=barplot(c(mean(subset(table.GC,organism==i)$X1st.GC.before),
                      mean(subset(table.GC,organism==i)$X1st.GC.after),
                      increase.first,
                      mean(subset(table.GC,organism==i)$X2nd.GC.before),
                      mean(subset(table.GC,organism==i)$X2nd.GC.after),
                      increase.second,
                      mean(subset(table.GC,organism==i)$X3rd.GC.before),
                      mean(subset(table.GC,organism==i)$X3rd.GC.after),
                      increase.third),
                    names=rep(c("%GC before","%GC after","GC increase"),3),
                    las=2,
                    ylim=c(-5,70),
                    col=rep(c("grey50","grey50","steelblue"),3),
                    axes=F,
                    main=paste0("Editing effect on the GC content \nper position within ",i))

axis(side=2,col="grey50",lwd=2)
axis(side=4,col="steelblue",lwd=2)
lines(x=c(barplot.pos[1],barplot.pos[3]),y=c(60,60))
text(x=barplot.pos[2], y=60, "1st Pos", pos=3, cex=1)
lines(x=c(barplot.pos[4],barplot.pos[6]),y=c(50,50))
text(x=barplot.pos[5], y=50, "2nd Pos", pos=3, cex=1)
lines(x=c(barplot.pos[7],barplot.pos[9]),y=c(40,40))
text(x=barplot.pos[8], y=40, "3rd Pos", pos=3, cex=1)

sd.first=sqrt(var(subset(table.types,organism==i & position=="first")$A.to.C)+
                var(subset(table.types,organism==i & position=="first")$A.to.G)+
                var(subset(table.types,organism==i & position=="first")$T.to.C)+
                var(subset(table.types,organism==i & position=="first")$T.to.G)+
                var(subset(table.types,organism==i & position=="first")$C.to.A)+
                var(subset(table.types,organism==i & position=="first")$C.to.T)+
                var(subset(table.types,organism==i & position=="first")$G.to.A)+
                var(subset(table.types,organism==i & position=="first")$G.to.T))/
  sqrt(nrow(subset(table.types,organism==i & position=="first")))
sd.second=sqrt(var(subset(table.types,organism==i & position=="second")$A.to.C)+
                 var(subset(table.types,organism==i & position=="second")$A.to.G)+
                 var(subset(table.types,organism==i & position=="second")$T.to.C)+
                 var(subset(table.types,organism==i & position=="second")$T.to.G)+
                 var(subset(table.types,organism==i & position=="second")$C.to.A)+
                 var(subset(table.types,organism==i & position=="second")$C.to.T)+
                 var(subset(table.types,organism==i & position=="second")$G.to.A)+
                 var(subset(table.types,organism==i & position=="second")$G.to.T))/
  sqrt(nrow(subset(table.types,organism==i & position=="second")))
sd.third=sqrt(var(subset(table.types,organism==i & position=="third")$A.to.C)+
                var(subset(table.types,organism==i & position=="third")$A.to.G)+
                var(subset(table.types,organism==i & position=="third")$T.to.C)+
                var(subset(table.types,organism==i & position=="third")$T.to.G)+
                var(subset(table.types,organism==i & position=="third")$C.to.A)+
                var(subset(table.types,organism==i & position=="third")$C.to.T)+
                var(subset(table.types,organism==i & position=="third")$G.to.A)+
                var(subset(table.types,organism==i & position=="third")$G.to.T))/
  sqrt(nrow(subset(table.types,organism==i & position=="third")))


pos.int=c(mean(subset(table.GC,organism==i)$X1st.GC.before),
          mean(subset(table.GC,organism==i)$X1st.GC.after),
          increase.first,
          mean(subset(table.GC,organism==i)$X2nd.GC.before),
          mean(subset(table.GC,organism==i)$X2nd.GC.after),
          increase.second,
          mean(subset(table.GC,organism==i)$X3rd.GC.before),
          mean(subset(table.GC,organism==i)$X3rd.GC.after),
          increase.third)

conf.int=c(sd(subset(table.GC,organism==i)$X1st.GC.before)/sqrt(nrow(subset(table.GC,organism==i))),
           sd(subset(table.GC,organism==i)$X1st.GC.after)/sqrt(nrow(subset(table.GC,organism==i))),
           sd.first,
           sd(subset(table.GC,organism==i)$X2nd.GC.before)/sqrt(nrow(subset(table.GC,organism==i))),
           sd(subset(table.GC,organism==i)$X2nd.GC.after)/sqrt(nrow(subset(table.GC,organism==i))),
           sd.second,
           sd(subset(table.GC,organism==i)$X3rd.GC.before)/sqrt(nrow(subset(table.GC,organism==i))),
           sd(subset(table.GC,organism==i)$X3rd.GC.after)/sqrt(nrow(subset(table.GC,organism==i))),
           sd.third)*1.96

segments(barplot.pos,y0=pos.int-conf.int,y1=pos.int+conf.int,lwd=1)
arrows(barplot.pos, pos.int-conf.int, 
       barplot.pos, pos.int+conf.int, 
       lwd = 1, angle = 90,
       code = 3, length = 0.05)

mtext(side=4,text="Number of edits increasing GC content",col="steelblue",padj=4)
mtext(side=2,text="%GC content",col="grey50",padj=-4)

}


dev.off()

