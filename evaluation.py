#**************************************************************************#
# * Run the command 'python evaluation.py' in the Terminal to implement IMC. #
# Or directly run this file in the IDE.
# * This function calculates mean precision, recall and F1 of the results
# provided by PersonaDrive against reference sets
# * Make sure you select the correct cancer name and data type               #
#**************************************************************************#

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


class evaluation:
    def __init__(self, mut_dic, reference_sets,cancer,dataset,results,cgc=None):
        print("C) Evaluation Started")
        self.cancer=cancer
        self.mut_dic=mut_dic    #mutation dictionary
        self.reference_sets=reference_sets #dictionary of reference sets
        self.results=results    #path to methods results
        self.dataset=dataset #ccle or tcga
        self.cgc=cgc            #cgc reference set used only if the data type is CCLE

        #list the result files
        methods=[]
        for m in os.listdir(self.results):
            if ".txt" in m :
                methods.append(self.results+'/'+m.split('.txt')[0])


        #calculate the top N value
        for ref_version in self.reference_sets:
            reference = self.reference_sets[ref_version] #select the current reference set
            reference_sizes=[] # list to save the personalized ref. set sizes
            selected_samples=[] # a list to save the set of samples with ref. set size >=3
            for i in mut_dic: # for each sample in the mutation dictionary
                if self.dataset=="TCGA":
                    reference_i = set(mut_dic[i]).intersection(set(reference)) # if tcga data, the personalized ref. set is the overlap of mutations of sample i with ref. genes
                elif self.dataset=="CCLE": # if ccle data, the personalized ref. set is the overlap of mutations of sample i with targets of sample i filtered by CGC
                    if i not in reference:
                        continue

                    reference_i = set(mut_dic[i]).intersection(set(reference[i])).intersection(set(cgc))

                if len(reference_i) >= 3: # add the current sample to the selected samples for evaluation, if the size of personalized ref set is greater than 2
                    reference_sizes.append(len(reference_i))
                    selected_samples.append(i)

            N=2*int((statistics.median(reference_sizes))) # top N is twice the median of peronalized reference set sizes
            print( "   -)Selected top N for: ",ref_version," is: ",N)

            title = self.dataset+"_"+self.cancer+"_["+ref_version+"]_"+'[max_N={}]'.format(N)
            #run the evaluation
            precision_scores,recall_scores,f1_scores = self.exclusive_evaluation( reference=reference, methods = methods, top_N=N, evaluation_ttile=title, mut_dic=mut_dic, selected_samples=selected_samples,dataset=self.dataset)

            #plot the results
            self.plot(f1_scores,precision_scores,recall_scores,title)

    def calculate_precision(self,drivers,reference,mut_dic,sample_i,k,dataset):
        #calculate the precision score for sample_i
        if dataset == "TCGA":
            ref_i = set(mut_dic[sample_i]).intersection(set(reference))
        elif dataset == "CCLE":

            ref_i = set(mut_dic[sample_i]).intersection(set(reference[sample_i])).intersection(set(self.cgc))
        drivers_ref = set(drivers).intersection(set(ref_i))
        precision = len(drivers_ref)/float(k)
        return precision

    def calculate_recall(self,drivers,reference,mut_dic,sample_i,total_samples,dataset):
        if dataset == "TCGA":
            ref_i = set(mut_dic[sample_i]).intersection(set(reference))
        elif dataset == "CCLE":
            ref_i = set(mut_dic[sample_i]).intersection(set(reference[sample_i])).intersection(set(self.cgc))
        drivers_ref = set(drivers).intersection(set(ref_i))

        if len(ref_i) == 0:
            recall = 0
            total_samples += 1
        else:
            recall = len(drivers_ref)/float(len(ref_i))
            total_samples += 1
        return recall, total_samples

    def load_method_results(self,model):
        drivers={}
        with open(str(model)+'.txt') as ifile:

            for line in ifile.readlines():
                drivers[line.strip().split("\t")[0].replace(".","-")]=line.strip().split("\t")[1:]

        return drivers

    def exclusive_evaluation(self,methods,top_N,evaluation_ttile,reference,mut_dic,selected_samples,dataset):

        #methods results dictionaries
        methods_precision_results = {}
        methods_recall_results = {}
        methods_f1_results = {}

        #progress bar
        pbar = tqdm(total=len(methods))

        for method in methods:
            pbar.update(1)
            #get method name
            method_name = method.split('/')[-1]

            #intitialize the final result array lists
            methods_precision_results[method_name]= []
            methods_recall_results[method_name]= []
            methods_f1_results[method_name]= []


            #intitialize the kth results arrays
            method_k_f1 = []
            method_k_precision = []
            method_k_recall = []

            prec_dic = {}
            rec_dic = {}


            #laod method results
            method_drivers = self.load_method_results(method)
            for k in range(1,top_N+1):

                #initialize the precision,f1, and recall scores for the current k
                recall = 0
                precision = 0
                f1_score = 0
                total_samples = 0

                for i in method_drivers: # for each sample i in S
                    if i not in selected_samples: # check that the current sample is selected (which satifies the critiria of selecing a sample to be evaluaated base on its reference set size)
                        continue
                    if len(method_drivers[i])>=k:
                        #get the the first k drivers for sample i
                        drivers = method_drivers[i][:k]
                        #calculate the samples's precision at top k
                        prec = self.calculate_precision(drivers,reference,mut_dic,i,k   ,dataset)
                        precision +=prec
                        #calculate the samples's recall at top k
                        rec, total_samples = self.calculate_recall(drivers,reference,mut_dic,i,total_samples,dataset)
                        recall +=rec

                if total_samples == 0:
                    continue
                # calculate the average precision, recall, and F1 scores
                prec_dic[k] = precision/total_samples
                rec_dic[k] = recall/total_samples
                method_k_recall.append(recall/total_samples)
                method_k_precision.append(precision/total_samples)

                if prec_dic[k] ==0 and rec_dic[k] == 0:
                    method_k_f1.append(f1_score)
                    continue
                else:
                    f1_score = 2 * ((prec_dic[k] * rec_dic[k]) / (prec_dic[k] + rec_dic[k]))
                    method_k_f1.append(f1_score)


            methods_f1_results[method_name]=method_k_f1
            methods_recall_results[method_name]=method_k_recall
            methods_precision_results[method_name]=method_k_precision
        return methods_precision_results ,methods_recall_results,methods_f1_results

    def plot(self,f1_results,precision_results,recall_results,title):
        fig, (ax1, ax2,ax3) = plt.subplots(1, 3,figsize=(20,6))
        fontsize=11
        lw=2.0 #markersize
        ms=6 #linewidth
        methods=[m for m in f1_results]
        for m in methods:
            x_axis=[x for x in range(1,len(f1_results[m])+1)]
            ax1.plot(x_axis, precision_results[m],'-o',markersize=ms,linewidth=lw) #precision
            ax2.plot(x_axis, recall_results[m],'-o',markersize=ms,linewidth=lw) #recall
            ax3.plot(x_axis, f1_results[m],'-o',markersize=ms,linewidth=lw) #f1

        ##################################################################
        #                           legends                              #
        ##################################################################
        #precision
        legend_precision = []
        for m in precision_results:
            auc=metrics.auc(x_axis,precision_results[m])
            legend_precision.append('{} (AUC: {})'.format(m, str(round(auc, 2))))
        legend = ax1.legend(legend_precision, loc=8,fancybox=True, fontsize= fontsize, framealpha=0,
                          edgecolor = 'b', ncol= 1, bbox_to_anchor=(0.69,0.90),borderaxespad=0.)

        #recall
        legend_recall = []
        for m in precision_results:
            auc=metrics.auc(x_axis,recall_results[m])
            legend_recall.append('{} (AUC: {})'.format(m, str(round(auc, 2))))
        legend = ax2.legend(legend_recall, loc=8,fancybox=True, fontsize= fontsize, framealpha=0,
                          edgecolor = 'b', ncol= 1, bbox_to_anchor=(0.69,0.0),borderaxespad=0.)

        #f1
        legend_f1 = []
        for m in f1_results:
            auc=metrics.auc(x_axis,f1_results[m])
            legend_f1.append('{} (AUC: {})'.format(m, str(round(auc, 2))))
        legend = ax3.legend(legend_f1, loc=8,fancybox=True, fontsize= fontsize, framealpha=0,
                          edgecolor = 'b', ncol= 1, bbox_to_anchor=(0.69,0.0),borderaxespad=0.)

        ##################################################################
        #                         end [legends]                          #
        ##################################################################

        fig.suptitle(title, fontsize=16)
        #x label
        ax1.set_xlabel('Top N genes',fontsize=fontsize)#,
        ax2.set_xlabel('Top N genes',fontsize=fontsize)#,
        ax3.set_xlabel('Top N genes',fontsize=fontsize)#,

        #y label
        ax1.set_ylabel('Precision',fontsize=fontsize)#,
        ax2.set_ylabel('Recall',fontsize=fontsize)#,
        ax3.set_ylabel('F1 Score',fontsize=fontsize)#,

        ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax3.xaxis.set_major_locator(MaxNLocator(integer=True))


        plt.show()
        plt.close()




class ccle_evaluation:

    def __init__(self, cancer,dataset,netowrk):
        self.cancer = cancer
        self.netowrk = netowrk
        self.dataset=dataset
        DataC=Data(cancer,dataset,netowrk)
        print("A) Load Refernce Data")
        self.cgc,self.cgc_specific=DataC.load_cgc_reference_sets(cancer)
        drug_targets=DataC.load_gdsc_depmap_reference_sets(cancer)


        #laod mutation dataset
        mut_dic=DataC.get_mutation_data()
        self.reference_sets={"GDSC U DepMap":drug_targets}
        self.results="results/"+self.dataset +"/"+self.cancer+"_"+self.netowrk

        evaluation(mut_dic,self.reference_sets,self.cancer,self.dataset,self.results,cgc=self.cgc)



class tcga_evaluation:

    def __init__(self, cancer,dataset,netowrk):

        self.cancer = cancer
        self.netowrk=netowrk
        if "COAD" in self.cancer:
            self.cancer_="COAD"
        elif "LUAD" in self.cancer:
            self.cancer_="LUAD"
        elif "BRCA" in self.cancer:
            self.cancer_="BRCA"
        DataC=Data(cancer,dataset,netowrk)
        self.dataset=dataset
        print("B) Load Refernce Data")
        cgc,cgc_cancer_type=DataC.load_cgc_reference_sets(cancer)
        ncg,ncg_cgc=DataC.load_ncg_reference_sets(cancer,cgc)
        cancermine,cancermine_cgc=DataC.load_cancermine_reference_sets(cancer,cgc)

        #laod mutation dataset
        mut_dic=DataC.get_mutation_data()
        self.reference_sets={ "CGC":cgc_cancer_type,
                            "NCG":ncg_cgc,
                            "CancerMine":cancermine_cgc,
                            "All_CGC":cgc,
                            "All_NCG":ncg,
                            "All_cancermine":cancermine}
        self.results="results/"+self.dataset +"/"+self.cancer+"_"+self.netowrk

        evaluation(mut_dic,self.reference_sets,self.cancer,self.dataset,self.results)



if __name__ == "__main__":

    description="Evaluation Framewok"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-d', "--dataset", type=str, required=True, default='TCGA', help="Dataset: TCGA or CCLE")
    parser.add_argument('-c', '--cancer', type=str, required=True, default="COAD", help="Cancer Type")
    parser.add_argument('-n', '--network', type=str, required=False, default="ST11", help="PPI Network")

    args = parser.parse_args()
    dataset = args.dataset
    cancer = args.cancer
    network = args.network
    if dataset=="TCGA":
        tcga_evaluation(cancer,dataset,network)
    elif dataset=="CCLE":
        ccle_evaluation(cancer,dataset,network)
