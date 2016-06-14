##################################################
# Statistical Analysis of Dinos ARN editing Data #
##################################################

## We focus on 4 organisms : kv, km, sm, pl
## And 4 genes family, atp, psa, psb, pet


##############################
# Loading, transforming Data #
##############################

# /!\
# Optional if already done.
# /!\

##################
# SLIDING WINDOW #
##################


names=c("name","percent_above_average_edits","number_obs",
        "DF_(N-2)","nucleotide_pearson_correlation","nucleotide_t_value",
        "amino_acid_pearson_correlation","amino_acid_t_value",
        "organism","plastid","reference")

ref_organism=c("Ptricornutum","Vbrassicaformis","Ctobin","Ehuxleyi","Amphidinium")
ref_org=c("Pt","Vb","Ct","Eh","Ac")

for (i in 1:5){
  
  setwd("/Users/Lucas/dropbox/Chris/Very\ Final\ Stuff/SW_New_Kv_and_Ht/Trimmed/Both")
  
  kv=read.delim(file=paste0("Kveneficum_New/",ref_organism[i],"/Kv_",ref_org[i],"_sliding_window_out.csv"),sep=",")
  kv=na.omit(kv) 
  kv=cbind(kv,rep("kv",length(kv$name)))
  kv=cbind(kv,rep("kv",length(kv$name)))
  kv=cbind(kv,rep(ref_org[i],length(kv$name)))
  names(kv)=names

  setwd("/Users/Lucas/dropbox/Chris/Very\ Final\ Stuff/Very_Final_Sliding_Window/Trimmed/Both")
  
  km=read.delim(file=paste0("Kmikimotoi/",ref_organism[i],"/Km_",ref_org[i],"_sliding_window_out_trim.csv"),sep=",")
  km=na.omit(km)
  km=cbind(km,rep("km",length(km$name)))
  km=cbind(km,rep("km",length(km$name)))
  km=cbind(km,rep(ref_org[i],length(km$name)))
  names(km)=names
  
  pl=read.delim(file=paste0("Plunula/",ref_organism[i],"/Pl_",ref_org[i],"_sliding_window_out_trim.csv"),sep=",")
  pl=na.omit(pl)
  pl=cbind(pl,rep("pl",length(pl$name)))
  pl=cbind(pl,rep("pl",length(pl$name)))
  pl=cbind(pl,rep(ref_org[i],length(pl$name)))
  names(pl)=names
  
  sm=read.delim(file=paste0("Sminutum/",ref_organism[i],"/Sm_",ref_org[i],"_sliding_window_out_trim.csv"),sep=",")
  sm=na.omit(sm)
  sm=cbind(sm,rep("sm",length(sm$name)))
  sm=cbind(sm,rep("sm",length(sm$name)))
  sm=cbind(sm,rep(ref_org[i],length(sm$name)))
  names(sm)=names
  
  table_trans=rbind(kv,km,pl,sm)
  
  table=cbind(table_trans[,1],table_trans[,1:length(table_trans)])
  names(table)=c("gene","gene_family","percent_above_average_edits","number_obs",
                  "DF_(N-2)","nucleotide_pearson_correlation","nucleotide_t_value",
                  "amino_acid_pearson_correlation","amino_acid_t_value",
                  "organism","plastid","reference")
 
  table[,1]=gsub(".*_atp","atp",table[,1])
  table[,1]=gsub(".*_psa","psa",table[,1])
  table[,1]=gsub(".*_psb","psb",table[,1])
  table[,1]=gsub(".*_pet","pet",table[,1])
  
  table[,1]=gsub("_trimmed","",table[,1])
  
  table[,2]=gsub(".*atp.*","atp",table[,2])
  table[,2]=gsub(".*psa.*","psa",table[,2])
  table[,2]=gsub(".*psb.*","psb",table[,2])
  table[,2]=gsub(".*pet.*","pet",table[,2])
  
  table[,length(table)-1]=gsub("kv","fuco",table[,length(table)-1])
  table[,length(table)-1]=gsub("km","fuco",table[,length(table)-1])
  table[,length(table)-1]=gsub("sm","pere",table[,length(table)-1])
  table[,length(table)-1]=gsub("pl","pere",table[,length(table)-1])
  
  
  table=subset(table,gene_family== "psa" | gene_family== "psb" | gene_family== "atp" | gene_family== "pet")
  table.house=subset(table,gene_family!= "psa" & gene_family!= "psb" & gene_family!= "atp" & gene_family!= "pet")
  
  setwd("/Users/Lucas/dropbox/Chris/Very\ Final\ Stuff/Very_Final_Sliding_Window/Lucas")

  write.table(table,paste0("SW_analysis_",ref_organism[i],"_lucas.csv"),sep=",",row.names=F)
  write.table(table.house,paste0("SW_analysis_housekeep_",ref_organism[i],"_lucas.csv"),sep=",",row.names=F)
}

setwd("/Users/Lucas/dropbox/Chris/Very\ Final\ Stuff/Very_Final_Sliding_Window/Lucas")

SW1=read.delim(file="SW_analysis_Ptricornutum_lucas.csv",sep=",")
SW2=read.delim(file="SW_analysis_Ctobin_lucas.csv",sep=",")
SW3=read.delim(file="SW_analysis_Ehuxleyi_lucas.csv",sep=",")
SW4=read.delim(file="SW_analysis_Vbrassicaformis_lucas.csv",sep=",")
SW5=read.delim(file="SW_analysis_Amphidinium_lucas.csv",sep=",")

SWH1=read.delim(file="SW_analysis_housekeep_Ptricornutum_lucas.csv",sep=",")
SWH2=read.delim(file="SW_analysis_housekeep_Ctobin_lucas.csv",sep=",")
SWH3=read.delim(file="SW_analysis_housekeep_Ehuxleyi_lucas.csv",sep=",")
SWH4=read.delim(file="SW_analysis_housekeep_Vbrassicaformis_lucas.csv",sep=",")
SWH5=read.delim(file="SW_analysis_housekeep_Amphidinium_lucas.csv",sep=",")

SW=rbind(SW1,SW2,SW3,SW4,SW5)
SWH=rbind(SWH1,SWH2,SWH3,SWH4,SWH5)

write.table(SW,"SW_analysis_global_lucas.csv",sep=",",row.names=F)
write.table(SWH,"SWH_analysis_global_lucas.csv",sep=",",row.names=F)

################################
# EDITING SCORE ON AMINO ACIDS #
################################

setwd("/Users/Lucas/dropbox/Chris/Very\ Final\ Stuff/Amino_Acid_Analysis_V3")

names=c("gene","num.AA.changes","num.identical.before",
        "num.identical.after","num.similar.before","num.similar.after",
        "average.edit.score.diff","organism","reference")

ref_organism=c("Ptricornutum","Vbrassicaformis","Ctobin","Ehuxleyi","Amphidinium")
ref_org=c("Pt","Vb","Ct","Eh","Ac")

for (i in 1:5){
  
  kv=read.delim(file=paste0("Kveneficum_New/",ref_organism[i],"/Kv_",ref_org[i],"_AA_outN.csv"),sep=",")
  kv=na.omit(kv)
  kv=cbind(kv,rep("kv",length(kv$gene)))
  kv=cbind(kv,rep(ref_org[i],length(kv$gene)))
  names(kv)=names
  
  km=read.delim(file=paste0("Kmikimotoi/",ref_organism[i],"/Km_",ref_org[i],"_AA_out.csv"),sep=",")
  km=na.omit(km)
  km=cbind(km,rep("km",length(km$gene)))
  km=cbind(km,rep(ref_org[i],length(km$gene)))
  names(km)=names
  
  pl=read.delim(file=paste0("Plunula/",ref_organism[i],"/Pl_",ref_org[i],"_AA_out.csv"),sep=",")
  pl=na.omit(pl)
  pl=cbind(pl,rep("pl",length(pl$gene)))
  pl=cbind(pl,rep(ref_org[i],length(pl$gene)))
  names(pl)=names
  
  sm=read.delim(file=paste0("Sminutum/",ref_organism[i],"/Sm_",ref_org[i],"_AA_out.csv"),sep=",")
  sm=na.omit(sm)
  sm=cbind(sm,rep("sm",length(sm$gene)))
  sm=cbind(sm,rep(ref_org[i],length(sm$gene)))
  names(sm)=names
  
  table_trans=rbind(kv,km,pl,sm)
  table=cbind(table_trans[,1],table_trans[,1:2],table_trans[,7:8],table_trans[,8],table_trans[,9])
  names(table)=c("gene","gene_family","num.AA.changes",
                 "average.edit.score.diff","organism","plastid","reference")
  
  table[,2]=gsub("atp.*","atp",table[,2])
  table[,2]=gsub("psa.*","psa",table[,2])
  table[,2]=gsub("psb.*","psb",table[,2])
  table[,2]=gsub("pet.*","pet",table[,2])
  
  table[,6]=gsub("kv","fuco",table[,6])
  table[,6]=gsub("km","fuco",table[,6])
  table[,6]=gsub("sm","pere",table[,6])
  table[,6]=gsub("pl","pere",table[,6])
  
  
  table=subset(table,gene_family== "psa" | gene_family== "psb" | gene_family== "atp" | gene_family== "pet")
  table.house=subset(table,gene_family!= "psa" & gene_family!= "psb" & gene_family!= "atp" & gene_family!= "pet")
  
  write.table(table,paste0("table_analysis_",ref_organism[i],"_lucas.csv"),sep=",",row.names=F)
  write.table(table.house,paste0("table_analysis_housekeep_",ref_organism[i],"_lucas.csv"),sep=",",row.names=F)
  
}


table1=read.delim(file="table_analysis_Ptricornutum_lucas.csv",sep=",")
table2=read.delim(file="table_analysis_Ctobin_lucas.csv",sep=",")
table3=read.delim(file="table_analysis_Ehuxleyi_lucas.csv",sep=",")
table4=read.delim(file="table_analysis_Vbrassicaformis_lucas.csv",sep=",")
table5=read.delim(file="table_analysis_Amphidinium_lucas.csv",sep=",")

tableH1=read.delim(file="table_analysis_housekeep_Ptricornutum_lucas.csv",sep=",")
tableH2=read.delim(file="table_analysis_housekeep_Ctobin_lucas.csv",sep=",")
tableH3=read.delim(file="table_analysis_housekeep_Ehuxleyi_lucas.csv",sep=",")
tableH4=read.delim(file="table_analysis_housekeep_Vbrassicaformis_lucas.csv",sep=",")
tableH5=read.delim(file="table_analysis_housekeep_Amphidinium_lucas.csv",sep=",")

table=rbind(table1,table2,table3,table4,table5)
tableH=rbind(tableH1,tableH2,tableH3,tableH4,tableH5)


write.table(table,"table_analysis_global_lucas.csv",sep=",",row.names=F)
write.table(tableH,"table_analysis_housekeep_global_lucas.csv",sep=",",row.names=F)


