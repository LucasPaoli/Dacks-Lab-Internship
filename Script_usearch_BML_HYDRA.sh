#!/bin/sh

##These commands set up the Grid Environment for your job:
#PBS -N Usearch_BML2
#PBS -l nodes=1,walltime=1000:00:00
#PBS -q default
#PBS -M paoli@ualberta.ca
#PBS -m abe
##Print the Host name of Node
echo "Hostname is " $HOSTNAME
##print the time and date
##Print  working Direcitory
pwd
##Print Date
date
##send command to Cluster nodes


#  Script.sh
#
#  Created by Lucas Paoli on 15/03/2016. Contact : lucas.paoli@ens.fr | lucas.paoli@gmail.com
#  On OS X 10.9.5
#
# Purpose : Write the Uparse pipeline
#
# Requires usearch, blastn (BLAST is assumed to be in the PATH) and seqtk.
# Tested with BLAST: 2.2.31+, usearch: v8.1.1861_i86osx32, seqtk: 1.0-r77-dirty
#
# Usearch : http://www.drive5.com/usearch/ (binary to downoald, no install needed)
# Seqtk on github : https://github.com/lh3/seqtk (download and use "make" to compile) \
#   On MasOS another option is "brew install seqtk"

#############
# Variables #
#############

# TO BE CHANGED :

PATHDATA="/home/paoli/BML/Data"
OUTPUT="/home/paoli/BML/usearch/output_BML2"
N_THREADS=8
DB_PR2="/home/paoli/opt/OTU_clustering_LP/DB/gb203_pr2.fasta"
DB_SILVA="/home/paoli/opt/OTU_clustering_LP/DB/SILVA_123_SSURef_Nr99_tax_silva.fasta"
SAMPLES="7C	7D	7E	7F	7G	7H	8A	8B	8C	8D	8E	8F	8G	8H	9A	9B	9C	9D	3A	3B	3C	3D	3E	3F	3G	3H	4A	4B	4C	4D	4E	4F	4G	4H	5A	5B	5C	5D	5E	5F	5G	5H	6A	6B	6C	6D"
# String, space separated with sample identifiers. Each string characterizing \
# a pair of files (fwd and rev) corresponding to a sample.
LOOP="M02031"
# First characteristic string (berfore the first ":" usually), usually corresponding \
# to the name of the sequencing machine

# MIGHT NEED TO BE CHANGED

PATHSEQTK="/home/paoli/opt/seqtk"
USEARCH="/home/paoli/opt/lotus/usearch"

# DO NOT CHANGE

LOG="${OUTPUT}/log"

######################
# Custom bash script #
######################

# SAMPLE NAMING

sample_number=0

for i in $(echo "${SAMPLES}")
do
sed -i "s/@""${LOOP}""/@SMPL"${sample_number}":""${LOOP}""/g" "${PATHDATA}"/${i}*.fastq
# M02031 is to be replaced, depending on the files you have.
#rm "${PATHDATA}"/*.fastq.txt
let "sample_number=${sample_number}+1"
"${PATHSEQTK}"/seqtk trimfq "${PATHDATA}"/${i}*R1*.fastq > "${PATHDATA}"/${i}_R1.fastq
"${PATHSEQTK}"/seqtk trimfq "${PATHDATA}"/${i}*R2*.fastq > "${PATHDATA}"/${i}_R2.fastq
done

##################
# OTU CLUSTERING #
##################

mkdir "${OUTPUT}"
mkdir "${LOG}"

## Amelioration keypoints
##########################
## Improve merging results
# Remove low quality ends DONE
# Filtering step before merging 50 bp window truncating \
# end of sequences if quality to low \
# might improve the number of merged hits
## Less unknown sequences after blast
# Filter Bacteria with blast against SILVA using -perc_identity 0.97 DONE
## Add Swarm clustering

# Merging the paired end sequences
"${USEARCH}" \
-fastq_mergepairs "${PATHDATA}"/*_R1.fastq \
-fastqout "${OUTPUT}"/merged.fq \
-fastq_minovlen 30 \
-log "${LOG}"/merge.log

rm "${PATHDATA}"/*_R1.fastq
rm "${PATHDATA}"/*_R2.fastq

# Showing some stats on the merged sequences
"${USEARCH}" \
-fastq_stats "${OUTPUT}"/merged.fq \
-log "${LOG}"/stats_merged.log

# Filtering the sequences (could use maxee 1 and minlen 200/250 to be ess stringent)
"${USEARCH}" \
-fastq_filter "${OUTPUT}"/merged.fq \
-fastq_maxee 0.75 \
-fastq_minlen 300 \
-fastq_maxns 0 \
-fastaout "${OUTPUT}"/filtered.fa \
-log "${LOG}"/filt.log

# Doesn't work because there is no fasta_stats, idealy, to be corrected with another function
#"${USEARCH}" \
#-fastq_stats "${OUTPUT}"/filtered.fa \
#-log "${LOG}"/stats_filtered.log

# Filtering chimeras from the sequences (useful for the abundance matrix)
#"${USEARCH}" \
#-uchime_ref "${OUTPUT}"/filtered.fa \
#-db "${DB_PR2}" \
#-strand plus \
#-minh 1.0 \
#-nonchimeras "${OUTPUT}"/filtered_nonch.fa \
#-uchimeout "${OUTPUT}"/filtered.uchime \
#-uchimealns "${OUTPUT}"/filtered.aln \
#-log "${LOG}"/uchime_filtered.log

# Comment the uchime_ref and uncomment this line if you want to do fast tests on the script
cp "${OUTPUT}"/filtered.fa "${OUTPUT}"/filtered_nonch.fa

# Doesn't work because there is no fasta_stats, idealy, to be corrected with another function
#"${USEARCH}" \
#-fastq_stats "${OUTPUT}"/filtered.fa \
#-log "${LOG}"/stats_filter_nonch.log

# Dereplication of the filtered sequences
"${USEARCH}" \
-derep_fulllength "${OUTPUT}"/filtered_nonch.fa \
-sizeout \
-fastaout "${OUTPUT}"/uniques.fa \
-log "${LOG}"/derep.log

# OTU clustering, with the UPARSE algorithm
"${USEARCH}" \
-cluster_otus "${OUTPUT}"/uniques.fa \
-minsize 2 \
-sizeout \
-sizein \
-otus "${OUTPUT}"/otus.fa \
-relabel Otu \
-uparseout "${OUTPUT}"/uparseout.txt \
-log "${LOG}"/cluster_otus.log

# Add Swarm clustering here

# Doesn't work because there is no fasta_stats, idealy, to be corrected with another function
#"${USEARCH}" \
#-fastq_stats "${OUTPUT}"/otus.fa \
#-log "${LOG}"/stats_otus.log

# Producing the abundance matrix for the OTUs
"${USEARCH}" \
-usearch_global "${OUTPUT}"/filtered_nonch.fa \
-db "${OUTPUT}"/otus.fa \
-strand plus \
-id 0.97 \
-otutabout "${OUTPUT}"/otutab.txt \
-biomout "${OUTPUT}"/otutab.json \
-log "${LOG}"/usearch_global.log

#######################
# Taxonomy assignment #
#######################

# Blast against SILVA to identify Bacterium sequences
blastn \
-db "${DB_SILVA}" \
-query "${OUTPUT}"/otus.fa \
-out "${OUTPUT}"/otu_bacteria_output.txt \
-max_target_seqs 1 \
-num_threads "${N_THREADS}" \
-outfmt 6

sed -i $'s/\t/;/g' "${OUTPUT}"/otu_bacteria_output.txt
rm "${OUTPUT}"/otu_bacteria_output.txt.txt

# List of the OTUs with a Bacteria/Archae match
grep -e 'Bacteria;' -e 'Archaea;' "${OUTPUT}"/otu_bacteria_output.txt \
| cut -f1 -d$';' > "${OUTPUT}"/otu_bacteria.txt
sed -i '/^$/d' "${OUTPUT}"/otu_bacteria.txt
rm "${OUTPUT}"/otu_bacteria.txt.txt

# Delete the bacteria from all the OTUs as they are contamination
cp "${OUTPUT}"/otutab.txt "${OUTPUT}"/otutab_save.txt
for i in $(seq 1 $(grep -c "." "${OUTPUT}"/otu_bacteria.txt))
do
y=$(sed "${i}q;d" "${OUTPUT}"/otu_bacteria.txt)
sed -i "/"$y$'\t/d' "${OUTPUT}"/otutab.txt
done
rm "${OUTPUT}"/otutab.txt.txt

# Blast against PR2
# The argument "-perc_identity 95" can be added for more stringency
blastn \
-db "${DB_PR2}" \
-query "${OUTPUT}"/otus.fa \
-out "${OUTPUT}"/otu_taxonomy.txt \
-max_target_seqs 1 \
-num_threads "${N_THREADS}" \
-outfmt 6

# Slight modifications to the blast Output to be better read by R afterwards
sed -i $'s/\t/;/g' "${OUTPUT}"/otu_taxonomy.txt
cp "${OUTPUT}"/otu_taxonomy.txt "${OUTPUT}"/otu_taxonomy_save.txt
sed -i 's/|Eukaryota|/;Eukaryota|/g' "${OUTPUT}"/otu_taxonomy_save.txt
rm "${OUTPUT}"/otu_taxonomy_save.txt.txt
sed -i 's/|/;/g' "${OUTPUT}"/otu_taxonomy.txt
rm "${OUTPUT}"/otu_taxonomy.txt.txt

# Delete the bacteria from the blast output
for i in $(seq 1 $(grep -c "." "${OUTPUT}"/otu_bacteria.txt))
do
y=$(sed "${i}q;d" "${OUTPUT}"/otu_bacteria.txt)
sed -i "/"$y';/d' "${OUTPUT}"/otu_taxonomy.txt
sed -i "/"$y';/d' "${OUTPUT}"/otu_taxonomy_save.txt
done
rm "${OUTPUT}"/otu_taxonomy*.txt.txt

# List of all the OTUs beside bacteria contamination
cut -f1 -d$'\t' "${OUTPUT}"/otutab.txt > "${OUTPUT}"/otus_unknown.txt
sed -i '/^$/d' "${OUTPUT}"/otus_unknown.txt
rm "${OUTPUT}"/otus_unknown.txt.txt

# List of OTUs with PR2 match
cut -f1 -d$';' "${OUTPUT}"/otu_taxonomy.txt > "${OUTPUT}"/otus_known_list.txt
sed -i '/^$/d' "${OUTPUT}"/otus_known_list.txt
rm "${OUTPUT}"/otus_known_list.txt.txt

# Delete the known OTUs from all the OTUs to have the unknown OTUs.
for i in $(seq 1 $(grep -c "." "${OUTPUT}"/otus_known_list.txt))
do
y=$(sed "${i}q;d" "${OUTPUT}"/otus_known_list.txt)
sed -i "/"$y"$/d" "${OUTPUT}"/otus_unknown.txt
done
rm "${OUTPUT}"/otus_unknown.txt.txt

# The end is done with R
# Rscript /home/paoli/opt/OTU_clustering_LP/Script_usearch.R "${OUTPUT}"

# Let's do the analysis!
# Rscript Analysis.R "${sample_number}






