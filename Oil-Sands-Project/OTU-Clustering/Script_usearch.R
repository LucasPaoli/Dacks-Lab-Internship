##------------------------------------------------------##
##    Rscript to finish the Script.sh, UPARSE script    ##
##------------------------------------------------------##

args <- commandArgs(T)
workingdirectory=as.character(args[1])

setwd(workingdirectory)


n=6 ## Corresponding to the colonne attributed to Class, might need to be checked.
blast_output = read.delim(file="otu_taxonomy.txt",header=F,sep=";")
blast_output_save = read.delim(file="otu_taxonomy_save.txt",header=F,sep=";")


check_doublons <- function(blast_output)
{
    doublons=c(0)
    for (i in 1:(nrow(blast_output)-1)){
        if (blast_output[i,1]==blast_output[i+1,1]){
            doublons=c(doublons,i+1)
        }
    }
    if (length(doublons)>1){
        doublons=doublons[2:length(doublons)]
        blast_output=blast_output[-doublons,]
    }
    return(blast_output)
}

blast_output=check_doublons(blast_output)
blast_output_save=check_doublons(blast_output_save)

OTU = read.delim(file="otutab.txt",row.names=1)

Otus_unknown = read.delim(file="otus_unknown.txt",header=F)
names(Otus_unknown)="V1"
Otus_unknown[,2:length(blast_output)]="Unknown"
blast_output=rbind(blast_output,Otus_unknown)

Otus_unknown = read.delim(file="otus_unknown.txt",header=F)
names(Otus_unknown)="V1"
Otus_unknown[,2:length(blast_output_save)]="Unknown"
blast_output_save=rbind(blast_output_save,Otus_unknown)

##Créationd des matrices de transitions
add_validated_OTU <- function(x,y,z,n)  ## Add Validated OTU to the abundance matrix
{
  temp=x
  names=unique(as.vector(y[,n]))
  for (i in 1:nrow(y)){
    for (j in 1:length(names)){
      if (y[i,n]==names[j]){
        temp[as.character(names[j]),][is.na(temp[as.character(names[j]),])] <- 0
        temp[as.character(names[j]),]=
          temp[as.character(names[j]),]+z[as.character(y[i,1]),]
      }
    }
  }
  return(temp)
}

Class=OTU[0,]
Class=add_validated_OTU(Class,blast_output,OTU,n)
Species=OTU[0,]
Species=add_validated_OTU(Species,blast_output_save,OTU,4)

write.table(Class, "Class.txt", sep="\t")
write.table(Species, "Species.txt", sep="\t")

