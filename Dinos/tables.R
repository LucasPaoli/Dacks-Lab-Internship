##################################################
# Statistical Analysis of Dinos ARN editing Data #
##################################################

library(xtable)

#########################
# EDITING GLOBAL TRENDS #
#########################

setwd("/Users/Lucas/Dropbox/Chris/Very Final Stuff/Pipeline_CSV_files_for_Lucas/Master_outfiles")
table=read.delim(file="global_master_edit.csv")

liste=c("ch","ht","km","kv","lp","pl","sm")


###########
# TABLE 1 #
###########

table_1=as.data.frame(matrix(0,7,10))
row.names(table_1)=liste

for (i in liste){
  table_1[i,1]=nrow(subset(table,organism==i))
  table_1[i,2]=sum(subset(table,organism==i)$nucleotide.length)
  table_1[i,3]=sum(subset(table,organism==i)$amino.acid.length)
  table_1[i,4]=sum(subset(table,organism==i)$number.edits)
  table_1[i,5]=mean((subset(table,organism==i)$number.edits)/(subset(table,organism==i)$nucleotide.length)*100)
  table_1[i,6]=sd((subset(table,organism==i)$number.edits)/(subset(table,organism==i)$nucleotide.length)*100)
  table_1[i,7]=mean(subset(table,organism==i)$percent.non.synonymous.edits)
  table_1[i,8]=sd(subset(table,organism==i)$percent.non.synonymous.edits)
  table_1[i,9]=mean((subset(table,organism==i)$number.amino.acid.edits)/(subset(table,organism==i)$amino.acid.length)*100)
  table_1[i,10]=sd((subset(table,organism==i)$number.amino.acid.edits)/(subset(table,organism==i)$amino.acid.length)*100)
  
}

row.names(table_1)=c("C. horridum",
                     "H. triquetra",
                     "K. mikimotoi",
                     "K. veneficum",
                     "L. polyedrum",
                     "P. lunula",
                     "S. minutum")
names(table_1)=c("# Genes","Nuc Length", "AA Length", "# Edits", "Avg % Edits","+-" ,"Avg % Non-Syn","+-%" , "Avg % AA Change","+-%")

print(xtable(table_1,digits=c(0,0,0,0,0,2,2,2,2,2,2)))

###########
# TABLE 3 #
###########

table_3=as.data.frame(matrix(0,7,5))
row.names(table_3)=liste

for (i in liste){
  table_3[i,1]=mean(subset(table,organism==i)$GC.before)
  table_3[i,2]=sd(subset(table,organism==i)$GC.before)
  table_3[i,3]=mean(subset(table,organism==i)$GC.after)
  table_3[i,4]=sd(subset(table,organism==i)$GC.after)
  if (i != "lp"){
    hist(subset(table,organism==i)$GC.before)
    hist(subset(table,organism==i)$GC.after)
    table_3[i,5]=t.test(subset(table,organism==i)$GC.before,subset(table,organism==i)$GC.after)$p.value}
}

row.names(table_3)=c("C. horridum",
                     "H. triquetra",
                     "K. mikimotoi",
                     "K. veneficum",
                     "L. polyedrum",
                     "P. lunula",
                     "S. minutum")
names(table_3)=c("Avg %GC Before","+-", "Avg %GC After", "+-%", "P value")

print(xtable(table_3,digits=c(0,2,2,2,2,5)))

###########
# TABLE 4 #
###########

table_4=as.data.frame(matrix(0,7,6))
row.names(table_4)=liste

for (i in liste){
  table_4[i,1]=sum(subset(table,organism==i)$first.position.edits)
  table_4[i,2]=sum(subset(table,organism==i)$second.position.edits)
  table_4[i,3]=sum(subset(table,organism==i)$third.position.edits)
  table_4[i,4]=mean(subset(table,organism==i)$percent.edits.in.first.two.positions)
  table_4[i,5]=sd(subset(table,organism==i)$percent.edits.in.first.two.positions)
  if (i != "lp"){
    table_4[i,6]=chisq.test(c(table_4[i,1],table_4[i,2],table_4[i,3]))$p.value}
}

row.names(table_4)=c("C. horridum",
                     "H. triquetra",
                     "K. mikimotoi",
                     "K. veneficum",
                     "L. polyedrum",
                     "P. lunula",
                     "S. minutum")
names(table_4)=c("# 1st Pos Edits","# 2nd Pos Edits", "# 3rd Pos Edits","Avg % 1st Two Pos", "+-", "P value")

print(xtable(table_4,digits=c(0,0,0,0,2,2,5)))

###########
# TABLE 5 #
###########

table_5=as.data.frame(matrix(0,7,10))
row.names(table_5)=liste

for (i in liste){
  table_5[i,1]=mean(subset(table,organism==i)$fraction.polyT.before)
  table_5[i,2]=sd(subset(table,organism==i)$fraction.polyT.before)
  table_5[i,3]=mean(subset(table,organism==i)$fraction.polyT.after)
  table_5[i,4]=sd(subset(table,organism==i)$fraction.polyT.after)
  if (i != "lp"){
    hist(subset(table,organism==i)$fraction.polyT.before)
    hist(subset(table,organism==i)$fraction.polyT.after)
    table_5[i,5]=t.test(subset(table,organism==i)$fraction.polyT.before,
                        subset(table,organism==i)$fraction.polyT.after)$p.value}
  table_5[i,6]=mean(subset(table,organism==i)$fraction.70.percent.polyT.before)
  table_5[i,7]=sd(subset(table,organism==i)$fraction.70.percent.polyT.before)
  table_5[i,8]=mean(subset(table,organism==i)$fraction.70.percent.polyT.after)
  table_5[i,9]=sd(subset(table,organism==i)$fraction.70.percent.polyT.after)
  if (i != "lp"){
    hist(subset(table,organism==i)$fraction.70.percent.polyT.before)
    hist(subset(table,organism==i)$fraction.70.percent.polyT.after)
    table_5[i,10]=t.test(subset(table,organism==i)$fraction.70.percent.polyT.before,
                         subset(table,organism==i)$fraction.70.percent.polyT.after)$p.value}
}

row.names(table_5)=c("C. horridum",
                     "H. triquetra",
                     "K. mikimotoi",
                     "K. veneficum",
                     "L. polyedrum",
                     "P. lunula",
                     "S. minutum")

names(table_5)=c("% polyT before","+-","% polyT after","+-%", "P value","% 70 polyT Before", "+-%%","% 70 polyT After", "+-%%%", "P value")

print(xtable(table_5,digits=c(0,2,2,2,2,5,2,2,2,2,5)))

###########
# TABLE 2 #
###########

setwd("/Users/Lucas/Dropbox/Chris/Very Final Stuff/Pipeline_CSV_files_for_Lucas/Master_edit_type_outfiles")
table.edit=read.delim(file="Global_master_edit_editing_types.csv")

table.edit_2=as.data.frame(matrix(0,7,12))
row.names(table.edit_2)=liste

for (i in liste){
  table.edit_2[i,1]=sum(subset(table.edit,organism==i)$A.to.T)
  table.edit_2[i,2]=sum(subset(table.edit,organism==i)$A.to.G)
  table.edit_2[i,3]=sum(subset(table.edit,organism==i)$A.to.C)
  table.edit_2[i,4]=sum(subset(table.edit,organism==i)$T.to.A)
  table.edit_2[i,5]=sum(subset(table.edit,organism==i)$T.to.G)
  table.edit_2[i,6]=sum(subset(table.edit,organism==i)$T.to.C)
  table.edit_2[i,7]=sum(subset(table.edit,organism==i)$G.to.A)
  table.edit_2[i,8]=sum(subset(table.edit,organism==i)$G.to.T)
  table.edit_2[i,9]=sum(subset(table.edit,organism==i)$G.to.C)
  table.edit_2[i,10]=sum(subset(table.edit,organism==i)$C.to.A)
  table.edit_2[i,11]=sum(subset(table.edit,organism==i)$C.to.T)
  table.edit_2[i,12]=sum(subset(table.edit,organism==i)$C.to.G)
}

table.edit_2=sweep(table.edit_2,1,rowSums(table.edit_2),"/")*100

row.names(table.edit_2)=c("C. horridum",
                          "H. triquetra",
                          "K. mikimotoi",
                          "K. veneficum",
                          "L. polyedrum",
                          "P. lunula",
                          "S. minutum")
names(table.edit_2)=c("A/T","A/G", "A/C", "T/A", "T/G","T/C" ,"G/A","G/T" ,"G/C","C/A","C/T","C/G")

print(xtable(table.edit_2,digits=rep(2,13)))

###########
# TABLE S1 #
###########

setwd("/Users/Lucas/Dropbox/Chris/Very Final Stuff/Pipeline_CSV_files_for_Lucas/Master_codon_pref_outfiles")
table.codon=read.delim(file="Global_master_edit_codons.csv")

table.codon_S1=as.data.frame(matrix(0,21,7))
liste.AA=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","STOP","T","V","W","Y")
row.names(table.codon_S1)=liste.AA
names(table.codon_S1)=liste

for (j in liste){
  for (i in liste.AA){
    table.codon_S1[i,j]=chisq.test(
      cbind(subset(table.codon,organism==j & amino.acid==i)$genome.usage,
            subset(table.codon,organism==j & amino.acid==i)$mRNA.usage),correct=F)$p.value
  }
}

names(table.codon_S1)=c("C. horridum",
                        "H. triquetra",
                        "K. mikimotoi",
                        "K. veneficum",
                        "L. polyedrum",
                        "P. lunula",
                        "S. minutum")

print(xtable(table.codon_S1,digits=c(0,rep(4,7))))

####################
# Table Simulation #
####################

setwd("/Users/Lucas/Dropbox/Chris/Very Final Stuff/Pipeline_CSV_files_for_Lucas/Simulation_stats")
table.simulation=read.delim(file="Global_simulation.csv")

setwd("/Users/Lucas/Dropbox/Chris/Very Final Stuff/Pipeline_CSV_files_for_Lucas/Master_edit_type_outfiles")
table.edit=read.delim(file="Global_master_edit_editing_types.csv")

table.edit[,3:length(table.edit)]=sweep(table.edit[,3:length(table.edit)],1,rowSums(table.edit[,3:length(table.edit)]),"/")*100
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
table.edit[is.nan(table.edit)] <- 0
table.edit[,3:length(table.edit)]=round(table.edit[,3:length(table.edit)], digits = 2)

setwd("/Users/Lucas/Dropbox/Chris/Very Final Stuff/Pipeline_CSV_files_for_Lucas/Master_outfiles")
table=read.delim(file="global_master_edit.csv")

table.position=table[,c("organism","gene","first.position.edits","second.position.edits","third.position.edits")]
table.position[,3:length(table.position)]=sweep(table.position[,3:length(table.position)],1,rowSums(table.position[,3:length(table.position)]),"/")*100
table.position[is.nan(table.position)] <- 0
table.position[,3:length(table.position)]=round(table.position[,3:length(table.position)], digits = 2)

table.stats=table.simulation[,1:4]
table.stats[,3]=rep(NA,nrow(table.stats))
table.stats[,4]=rep(NA,nrow(table.stats))
names(table.stats)[3:4]=c("edit position","edit type")

levels(table.simulation$organism)=unique(c(levels(table.simulation$organism),levels(table.position$organism)))
levels(table.stats$organism)=unique(c(levels(table.simulation$organism),levels(table.position$organism)))
levels(table.simulation$gene)=unique(c(levels(table.simulation$gene),levels(table.position$gene)))
levels(table.stats$gene)=unique(c(levels(table.simulation$gene),levels(table.position$gene)))


for (i in 1:nrow(table.stats)){
  
  transit1=cbind(c(subset(table.position,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$first.position.edits,
                   subset(table.position,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$second.position.edits,
                   subset(table.position,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$third.position.edits),
                 
                 c(subset(table.simulation,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$num.1st.pos,
                   subset(table.simulation,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$num.2nd.pos,
                   subset(table.simulation,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$num.3rd.pos))
  
  table.stats[i,3]=chisq.test(transit1[rowSums(transit1)>0,],correct=F)$p.value
  
  
  transit2=cbind(c(subset(table.edit,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$A.to.T,
                   subset(table.edit,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$A.to.G,
                   subset(table.edit,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$A.to.C,
                   subset(table.edit,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$T.to.A,
                   subset(table.edit,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$T.to.G,
                   subset(table.edit,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$T.to.C,
                   subset(table.edit,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$G.to.A,
                   subset(table.edit,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$G.to.T,
                   subset(table.edit,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$G.to.C,
                   subset(table.edit,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$C.to.A,
                   subset(table.edit,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$C.to.T,
                   subset(table.edit,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$C.to.G),
                 
                 c(subset(table.simulation,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$A.to.T,
                   subset(table.simulation,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$A.to.G,
                   subset(table.simulation,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$A.to.C,
                   subset(table.simulation,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$T.to.A,
                   subset(table.simulation,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$T.to.G,
                   subset(table.simulation,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$T.to.C,
                   subset(table.simulation,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$G.to.A,
                   subset(table.simulation,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$G.to.T,
                   subset(table.simulation,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$G.to.C,
                   subset(table.simulation,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$C.to.A,
                   subset(table.simulation,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$C.to.T,
                   subset(table.simulation,organism==table.stats[i,"organism"] & gene==table.stats[i,"gene"])$C.to.G))
  
  table.stats[i,4]=chisq.test(transit2[rowSums(transit2)>0,],correct=F)$p.value
}

setwd("/Users/Lucas/Dropbox/Chris/Very Final Stuff/Pipeline_CSV_files_for_Lucas/Simulation_stats")
write.table(table.stats,"stats_results.csv",sep=",",row.names=F)

