##############################
# BML BIOLOGICAL HYPPOTHESIS #
##############################

workingdirectory="/Users/Lucas/Documents/ENS/BIO_M1_stage/Report/Figures"
setwd(workingdirectory)

normal.lake=(matrix(0,5,18))
normal.lake[1,]=c(0.4,0.7,0.8,0.5,0.3,0.2,
                  0.2,0.2,0.2,0.2,0.2,0.2,
                  0.2,0.2,0.2,0.2,0.2,0.2)
normal.lake[2,]=c(0.05,0.2,0.5,0.6,0.7,0.3,
                  0.2,0.2,0.2,0.2,0.2,0.2,
                  0.2,0.2,0.2,0.2,0.2,0.2)
normal.lake[3,]=c(0,0.1,0.2,0.2,0.3,0.55,
                  0.6,0.3,0.2,0.2,0.2,0.2,
                  0.2,0.2,0.2,0.2,0.2,0.2)
normal.lake[4,]=c(0,0,0,0.05,0.05,0.1,
                  0.1,0.2,0.2,0.2,0.2,0.2,
                  0.2,0.2,0.2,0.2,0.2,0.2)
normal.lake[5,]=c(0,0,0,0,0,0.05,
                  0.07,0.1,0.1,0.2,0.2,0.2,
                  0.2,0.2,0.2,0.2,0.2,0.2)

row.names(normal.lake)=c("Photosynthesizer1","Photosynthesizer2","Photosynthesizer3","Grazer1","Grazer2")

pdf(file = "Figure12A.pdf", width=12, height=5)

plot(normal.lake["Photosynthesizer2",],type="n",ylab="Relative abundance",
     main="Expected relative abundance in Reservoir",
     xlab="Time",
     ylim=c(0,1))
lines(normal.lake["Photosynthesizer1",],col="lightgreen",lwd=2)
lines(normal.lake["Photosynthesizer2",],col="green",lwd=2)
lines(normal.lake["Photosynthesizer3",],col="green4",lwd=2)
lines(normal.lake["Grazer1",],col="blue",lwd=2)
lines(normal.lake["Grazer2",],col="blue4",lwd=2)
lines(x=c(6,6),y=c(0.8,0.9))
lines(x=c(12,12),y=c(0.8,0.9))
text(x=3, y=0.8, "Spring", pos=3, cex=0.8)
text(x=9, y=0.8, "Summer", pos=3, cex=0.8)
text(x=15, y=0.8, "Autumn", pos=3, cex=0.8)

legend("topleft", # places a legend at the appropriate place
       c("Photosynthesizer1","Photosynthesizer2","Photosynthesizer3","Grazer1","Grazer2"), # puts text in the legend
       lty=c(1,1,1,1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5,2.5,2.5),
       col=c("lightgreen","green","green4","blue","blue4"),
       ncol=5,bty="n") # gives the legend lines the correct color and width
dev.off()

pit.lake=(matrix(0,5,18))

pit.lake[1,]=c(0.05,0.4,0.6,0.2,0.1,0.05,
               0.05,0.05,0.05,0.05,0.05,0.05,
                  0,0,0,0,0,0)
pit.lake[2,]=c(0,0.1,0.2,0.5,0.6,0.4,
                  0.2,0.1,0.1,0.1,0.1,0.1,
                  0.1,0.1,0.1,0.1,0.1,0.1)
pit.lake[3,]=c(0.05,0.05,0.1,0.1,0.3,
                  0.55,0.7,0.6,0.4,0.1,0.1,
                  0.05,0,0,0,0,0,0)
pit.lake[4,]=c(0,0,0,0.05,0.05,0.1,
               0.1,0.3,0.55,0.7,0.6,0.4,
               0.1,0.1,0.05,0,0,0)
pit.lake[5,]=c(0.1,0.1,0.2,0.2,0.2,0.2,
                  0.2,0.2,0.2,0.3,0.3,0.4,
                  0.6,0.7,0.8,0.8,0.8,0.7)

row.names(pit.lake)=c("Grazer1","Grazer2","Grazer3","Grazer4","Fungi")

pdf(file = "Figure12B.pdf", width=12, height=5)

plot(pit.lake["Grazer1",],type="n",ylab="Relative abundance",
     main="Relative abundance in End Pit Lake",
     xlab="Time",
     ylim=c(0,1))
lines(pit.lake["Grazer1",],col="green",lwd=2)
lines(pit.lake["Grazer2",],col="lightblue",lwd=2)
lines(pit.lake["Grazer3",],col="blue",lwd=2)
lines(pit.lake["Grazer4",],col="blue4",lwd=2)
lines(pit.lake["Fungi",],col="red",lwd=2)
lines(x=c(6,6),y=c(0.8,0.9))
lines(x=c(12,12),y=c(0.8,0.9))
text(x=3, y=0.8, "Spring", pos=3, cex=0.8)
text(x=9, y=0.8, "Summer", pos=3, cex=0.8)
text(x=15, y=0.8, "Autumn", pos=3, cex=0.8)

legend("topleft", # places a legend at the appropriate place
       c("Grazer1","Grazer2","Grazer3","Grazer4","Fungi"), # puts text in the legend
       lty=c(1,1,1,1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5,2.5,2.5),
       col=c("green","lightblue","blue","blue4","red"),
       ncol=5,bty="n") # gives the legend lines the correct color and width
dev.off()

