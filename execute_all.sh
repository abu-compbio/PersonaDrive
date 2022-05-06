#!/usr/bin/env bash

# This bash script can be used to run all python files
# sequentially for the PersonaDrive algorithm. The required
# python libraries are numpy, scipy and networkx. The
# input files are provided within the data folder.

cancer_type="LUAD"
dataset="CCLE" #TCGA or CCLE
network="ST11" # where ST11-> String 11.5 | ST10 -> String 10.5 | DW -> DawnRank network


#########################################################
#	Construct PBNs
# The 'constructing_PBNs.py' script will construct the personalized bipartite networks (PBNs).
#########################################################

printf "################################################\n"
printf "    1. Personalized Bipartite Networks (PBNs).... \n"
printf "################################################\n\n\n"
python constructing_PBNs.py -d $dataset -c $cancer_type -n $network

#########################################################
#
#	Rank Mutated Genes
#
#########################################################

printf "################################################\n"
printf "    2 - Rank Mutated Genes ...\n"
printf "################################################\n\n\n"

python PersonDrive.py -d $dataset -c $cancer_type -n $network

#########################################################
#
#	Evaluation
#
#########################################################
printf "################################################\n"
printf "    3 - Evaluation ...\n"
printf "################################################\n\n\n"
python evaluation.py -d $dataset -c $cancer_type -n $network
