import argparse
import numpy as np
import os
import pandas as pd
import statistics
import matplotlib.pyplot as plt
from tqdm import tqdm as tqdm
from numpy import trapz
from sklearn import metrics
from matplotlib.ticker import MaxNLocator
from Data import *
class Data:

    def __init__(self, cancer,dataset,netowrk,mut_ppi_filter=True):

        print("---- Data ----")
        print("A) LOAD Data")
        self.cancer=cancer
        self.dataset=dataset
        self.netowrk=netowrk
        self.mut_ppi_filter=mut_ppi_filter
        #ppi networks
        if "COAD" in self.cancer:
            self.cancer_="COAD"
        elif "LUAD" in self.cancer:
            self.cancer_="LUAD"
        elif "BRCA" in self.cancer:
            self.cancer_="BRCA"
        networks={
            "ST11":"data/string_v11.5_network.csv",
            "ST10":"data/string_v10.5_network.csv",
            "DW":"data/dawnrank_network.csv"
        }
        ppi_network = networks[netowrk]
        #load the ppi network
        self.ppi_network=self.load_PPI(ppi_network)
        self.mut_dic=self.load_mutation_data(self.cancer,self.dataset,self.ppi_network,self.mut_ppi_filter)

    ##construct mutation dictionary
    def load_mutation_data(self,cancer,dataset,ppi,filter=True):

        mutation_matrix = pd.read_csv("data/"+cancer.split("_")[0]+"_"+dataset+"/MUT.csv", index_col=0,sep=",")

        # filter mutation matrix to contain only mutated genes exist in the PPI netowrk
        if filter:
            mutation_matrix=mutation_matrix.drop(index=[g for g in mutation_matrix.index if g not in ppi])

        #get sample IDs
        samples = mutation_matrix.columns.tolist()


        #get the set of mutated genes for each sample
        mutations = mutation_matrix.isin([1])

        # save the mutated genes of each sample in a dictionary
        mut_dic ={}
        for s in samples:
            s=s.replace(".",'-')
            mut_dic[s]= mutations.index[mutations[s] == True].tolist()
        print("   1) Mutation data loaded for ",len(samples)," samples. ",mutation_matrix.shape)
        return mut_dic

    def get_mutation_data(self):
        return self.mut_dic

    ## Load PPI
    def load_PPI(self,file):
        ppi_matrix = pd.read_csv(file,sep="\t")
        ppi = {}
        for i in range(len(ppi_matrix.index)):
            ge1 = ppi_matrix.loc[i,'genes1']
            ge2 = ppi_matrix.loc[i,'genes2']
            if ge1 not in ppi:
                ppi[ge1]=set()
                ppi[ge1].add(ge2)
            else:
                ppi[ge1].add(ge2)
            if ge2 not in ppi:
                ppi[ge2]=set()
                ppi[ge2].add(ge1)
            else:
                ppi[ge2].add(ge1)
        return ppi

    #load drug tqrgets reference set data
    def load_gdsc_depmap_reference_sets(self,cancer):


        data="data/reference_sets/{}_DepMap_GDSC2.txt".format(cancer.split("_")[0])
        gdsc_depmap={}
        with open(data) as ifile:
            for line in ifile.readlines():
                gdsc_depmap[line.split("\t")[0]]= line.strip().split("\t")[1:]

        print("   3)" ,"Reference sets from  [GDSC U DepMap] data loaded")
        return gdsc_depmap

    #load CGC reference set data
    def load_cgc_reference_sets(self,cancer):

        cgc_data = pd.read_csv("data/reference_sets/CGC.csv", index_col=0)
        cgc = cgc_data.index.tolist()
        print("   1) Loaded ",len(cgc)," drivers from CGC data")

        #identify cgc cancer type specific drivers
        cgc_cancer_type=set()
        cgc_keywords={
            "COAD":["colon","colorectal"],
            "LUAD":["lung"],
            "BRCA":["breast"]
        }

        for driver in cgc_data.index:
            if pd.isnull(cgc_data["Tumour Types(Somatic)"][driver]):
                continue
            for keyword in cgc_keywords[self.cancer_]:
                if any(value in cgc_data["Tumour Types(Somatic)"][driver].lower() for value in cgc_keywords[self.cancer_]):
                    #if keyword in cgc_data["Tumour Types(Somatic)"][driver].lower():
                    cgc_cancer_type.add(driver)

        print("   2) Identified ",len(cgc_cancer_type)," ",cancer," drivers from CGC data")
        return cgc,cgc_cancer_type

    #load NCG reference set data
    def load_ncg_reference_sets(self,cancer,cgc):
        #NCG
        with open ("data/reference_sets/NCG6_cancergenes.tsv","r") as ifile:
            ncg=set()
            ncg_cgc=set()

            ncg_keywords={
                "COAD":"colorectal",
                "LUAD":"lung",
                "BRCA":"breast"
            }

            for line in ifile.readlines():
                line=line.split("\t")
                if ncg_keywords[self.cancer_] in line[4]:
                    if line[1] in cgc:
                        ncg_cgc.add(line[1])
                    ncg.add(line[1])

        print("   3) Loaded ",len(ncg)," drivers from NCG data")
        print("   4) Identified ",len(ncg_cgc)," drivers from [NCG ∩ CGC] data")
        return ncg,ncg_cgc

    #load CancerMine reference set data
    def load_cancermine_reference_sets(self,cancer,cgc):
        #CancerMine
        cancermine = set()
        cancermine_cgc=set()
        cancermine_keywords={
            "COAD":"coad",
            "LUAD":"luad",
            "BRCA":"breast"
        }
        with  open('data/reference_sets/cancermine_{}.txt'.format(cancermine_keywords[self.cancer_])) as ifile:
            for line in ifile:
                gene = line.split()
                if gene[0] in cgc:
                    cancermine_cgc.add(gene[0])
                cancermine.add(gene[0])


        print("   5) Loaded ",len(cancermine)," drivers from CancerMine data")
        print("   6) Identified ",len(cancermine_cgc)," drivers from [CancerMIne ∩ CGC] data")
        return cancermine,cancermine_cgc
