#!/bin/sh

#  Script.sh
#
#  Created by Lucas Paoli on 15/03/2016. Contact : lucas.paoli@ens.fr | lucas.paoli@gmail.com
#  On OS X 10.9.5
#
# Purpose : Filter out Streptophytes and Metazoa

cd $1

# CHECK THAT STREPTOPHYTA ONLY CONTAINS EMBRYOPHYTA

grep -v "Metazoa" Species.txt | grep -v "Streptophyta" > Species_filtered.txt

grep -v "Metazoa" Class.txt | grep -v "Streptophyta"  > Class_filtered.txt






