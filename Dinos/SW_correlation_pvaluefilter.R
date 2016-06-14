################################
# SLIDING CORRELATION ANALYSIS #
################################

setwd("/Users/Lucas/dropbox/Chris/Very\ Final\ Stuff/Very_Final_Sliding_Window/Lucas")
SW_sub=read.delim(file="SW_analysis_global_lucas.csv",sep=",")
SW_sub=subset(SW_sub,dt(SW_sub$amino_acid_t_value,df=SW_sub$DF_.N.2.)<=0.05)

setwd("/Users/Lucas/dropbox/Chris/Very\ Final\ Stuff/Amino_Acid_Analysis_V3")
table_sub=read.delim(file="table_analysis_global_lucas.csv",sep=",")

SW_sub$gene=as.character(SW_sub$gene)
table_sub$gene=as.character(table_sub$gene)

# Adding the pearson correlation to the table_sub with edit scores

for (i in 1:nrow(table_sub)){
  if (length(subset(SW_sub, gene == table_sub[i,"gene"] &
                      organism == table_sub[i,"organism"] &
                      reference == table_sub[i,"reference"])$amino_acid_pearson_correlation)==0){
    table_sub[i,8]=NA
  }else{
    table_sub[i,8]=subset(SW_sub, gene == table_sub[i,"gene"] &
                            organism == table_sub[i,"organism"] &
                            reference == table_sub[i,"reference"])$amino_acid_pearson_correlation
  }
}

table_sub=na.omit(table_sub)

## CORRELATION

# ALL

cor.test(table_sub[,"average.edit.score.diff"],table_sub[,"V8"])
xyplot(table_sub[,"average.edit.score.diff"]~table_sub[,"V8"])
# gene
table_sub.gene=subset(table_sub, gene_family == "psb")
xyplot(table_sub.gene[,"average.edit.score.diff"]~table_sub.gene[,"V8"],groups=table_sub.gene$reference)
cor.test(table_sub.gene[,"average.edit.score.diff"],table_sub.gene[,"V8"])

# PEREDININ

table_sub.P=subset(table_sub, organism== "sm" | organism== "pl")
cor.test(table_sub.P[,"average.edit.score.diff"],table_sub.P[,"V8"])
xyplot(table_sub.P[,"average.edit.score.diff"]~table_sub.P[,"V8"],groups=table_sub.P$reference)

pdf(file = "correlation_SW_score_peredinin_pvalue.pdf", width=8, height=8)
xyplot(table_sub.P[,"average.edit.score.diff"]~table_sub.P[,"V8"],groups=table_sub.P$reference,
       main="Correlation for Peredinin plastids",
       ylab="Average editing score",
       xlab="Sliding window pearson coefficient",
       sub="corr. coeff. = 0.1792564, p-value = 0.07738",
       par.settings = list(superpose.symbol = list(pch = 16)),
       auto.key=list(title="Reference organism :",column=5,cex=0.8)
)
dev.off()

# gene family
table_sub.P.gene=subset(table_sub.P, gene_family == "psb")
xyplot(table_sub.P.gene[,"average.edit.score.diff"]~table_sub.P.gene[,"V8"],groups=table_sub.P.gene$reference)
cor.test(table_sub.P.gene[,"average.edit.score.diff"],table_sub.P.gene[,"V8"])

xyplot(table_sub.P[,"average.edit.score.diff"]~table_sub.P[,"V8"],groups=table_sub.P$gene_family,
       main="Correlation for Peredinin plastids",
       ylab="Average editing score",
       xlab="Sliding window pearson coefficient",
       sub="psb : corr. coeff. = 0.020885, p-value = 0.8808",
       par.settings = list(superpose.symbol = list(pch = 16)),
       auto.key=list(title="gene family :",column=4,cex=0.8)
)

# FUCOXANTHIN

table_sub.F=subset(table_sub, organism== "kv" | organism== "km")
cor.test(table_sub.F[,"average.edit.score.diff"],table_sub.F[,"V8"])
xyplot(table_sub.F[,"average.edit.score.diff"]~table_sub.F[,"V8"],groups=table_sub.F$reference)

pdf(file = "correlation_SW_score_fucoxanthin_pvalue.pdf", width=8, height=8)
xyplot(table_sub.F[,"average.edit.score.diff"]~table_sub.F[,"V8"],groups=table_sub.F$reference,
       main="Correlation for Fucoxanthin plastids",
       ylab="Average editing score",
       xlab="Sliding window pearson coefficient",
       sub="corr. coeff. = -0.0888516, p-value = 0.2879",
       par.settings = list(superpose.symbol = list(pch = 16)),
       auto.key=list(title="Reference organism :",column=5,cex=0.8)
)
dev.off()

# gene family
table_sub.gene=subset(table_sub, gene_family == "atp")
xyplot(table_sub.gene[,"average.edit.score.diff"]~table_sub.gene[,"V8"],groups=table_sub.gene$reference)
cor.test(table_sub.gene[,"average.edit.score.diff"],table_sub.gene[,"V8"])

xyplot(table_sub.F[,"average.edit.score.diff"]~table_sub.F[,"V8"],groups=table_sub.F$gene_family,
       main="Correlation for Fucoxanthin plastids",
       ylab="Average editing score",
       xlab="Sliding window pearson coefficient",
       sub="atp : corr. coeff. = -0.4698892 , p-value = 0.0003372",
       par.settings = list(superpose.symbol = list(pch = 16)),
       auto.key=list(title="gene family :",column=4,cex=0.8)
)