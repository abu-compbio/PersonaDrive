#***********************************
# Run the command 'python constructing_personalized_networks.py' in the Terminal to start the
# the construction of the personalized networkxs
# Or directly run this file in your chosen IDE.
# This may take several hours.
# When is done, the files of prepared data are saved in the subfolders of graphs/[dataset]/[cancer].
#***********************************
import argparse
from tqdm import tqdm as tqdm
import os
import sys
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite

## Load DEGs matrix
def load_DEGs(file):
    DEGs_matrix = pd.read_csv(file, index_col=0,sep=",")
    DEGs_matrix=DEGs_matrix.transpose()
    DEGs_matrix=DEGs_matrix.rename(columns = {i:i.replace('.','-') for i in DEGs_matrix})
    return DEGs_matrix

## Load MUTATION matrix
def load_mutations(file):
    M_matrix = pd.read_csv(file, index_col=0)
    M_matrix=M_matrix.rename(columns = {i:i.replace('.','-') for i in M_matrix})
    return M_matrix

## Load PPI
def load_ppi(file):
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

## construct personalized networks
def construct_PBNs(M,O,PPI,cancer,dataset):

    #create the output folder
    if not os.path.exists("graphs/"+dataset+"/"+cancer):
        os.mkdir("graphs/"+dataset+"/"+cancer)

    #####################################
    #   Starting with mutation data     #
    #####################################
    #list samples in M matrix
    samples = M.columns.tolist()
    M = M.isin([1])

    #dictionary that contains samples ids and their mutations
    M_dictionary ={}
    for i in samples:
        m_i = M.index[M[i] == True].tolist()
        M_dictionary[i]= list(set(m_i))

    #####################################
    #   Pre-processing outliers data    #
    #####################################
    pbar = tqdm(total=len(samples))
    #genes is a set of genes in outliers dataset
    genes = O.index.tolist()
    #list of outliers
    D_set= []
    for i in samples:
        pbar.update(1)
        for gene in genes:
            if i in O.columns:
                if O[i][gene] == True:
                    # outliers are in form of: CCLEXXXX_TP53
                    D_set.append(str(i+'_'+gene))
    print("There are: ",len(D_set)," outliers.")
    pbar.close()


    #####################################
    #  Constructing Bipartite networks  #
    #####################################
    pbar = tqdm(total=len(samples))

    for i in samples:

        if 'BP_'+i+'.gml' not in os.listdir("graphs/"+dataset+"/"+cancer):
            pbar.update(1)

            # initialize the bipartite networks B_i
            Bi = nx.Graph()
            Bi.add_nodes_from(M_dictionary[i], bipartite=0)
            Bi.add_nodes_from(D_set, bipartite=1)


            pbar2=tqdm(total=len(M_dictionary[i]),leave=False)
            Mi=M_dictionary[i]
            for u in Mi:
                pbar2.update(1)
                if u in PPI:
                    for v in D_set:
                        v_d = v.split('_')[1]
                        # check if the selected mutated gene is mutated in the sample of current outlier
                        if u in M_dictionary[v.split('_')[0]]:
                            if v_d in PPI:
                                #check if outlier is a neighbor of the mutation
                                if v_d in PPI[u]:
                                    if v_d !=u:
                                        Bi.add_edge(u,v)


            bipartite0_nodes = {n for n, d in Bi.nodes(data=True) if d['bipartite']==0}
            bipartite1_nodes = {n for n, d in Bi.nodes(data=True) if d['bipartite']==1}

            #deleting nodes with 0 degree
            for v in bipartite1_nodes:
                if int(Bi.degree(v))==0:
                    Bi.remove_node(v)

            for u in bipartite0_nodes:
                if int(Bi.degree(u))==0:
                    Bi.remove_node(u)


            bipartite0_nodes = {n for n, d in Bi.nodes(data=True) if d['bipartite']==0}
            bipartite1_nodes = {n for n, d in Bi.nodes(data=True) if d['bipartite']==1}

            print("Mutations:  ",len(bipartite0_nodes),"   Outliers:",len(bipartite1_nodes))
            nx.write_gml(Bi, "graphs/"+dataset+"/"+cancer + '/BP_'+i+'.gml')
        pbar.close()


if __name__ == "__main__":
    data_folder="data/"

    description="Personalized Bipartite Networks (PBNs)"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-d', "--dataset", type=str, required=False, default='TCGA', help="Dataset: TCGA or CCLE")
    parser.add_argument('-c', '--cancer', type=str, required=False, default="COAD", help="Cancer Type")

    args = parser.parse_args()
    dataset = args.dataset
    cancer = args.cancer
    #~~~~~~~~~~~~~Step 1：~~~~~~~~~~~~~~~~~~
    # load PPI data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ppi_file = data_folder+"ppi_dawnrank.csv"
    ppi = load_ppi(ppi_file)
    print('PPI network loaded...')

    #~~~~~~~~~~~~~Step 2：~~~~~~~~~~~~~~~~~~
    # load outliers data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEGs_matrix = load_DEGs(data_folder+cancer+"_"+dataset+"/DEGs_data.csv")
    print('DEGs loaded...')

    #~~~~~~~~~~~~~Step 3：~~~~~~~~~~~~~~~~~~
    # load mutation data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    M_matrix = load_mutations(data_folder+cancer+"_"+dataset+"/mutation_data.csv")
    print('Mutations loaded...')

    #~~~~~~~~~~~~~Step 4：~~~~~~~~~~~~~~~~~~
    # pre-process the matrices to contain only genes in the PPI
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEGs_matrix=DEGs_matrix.drop(index=[n for n in DEGs_matrix.index if n not in ppi])
    M_matrix=M_matrix.drop(index=[n for n in M_matrix.index if n not in ppi])

    #~~~~~~~~~~~~~Step 5：~~~~~~~~~~~~~~~~~~
    # construct the personalized networks
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    construct_PBNs(M_matrix,DEGs_matrix,ppi,cancer,dataset)
