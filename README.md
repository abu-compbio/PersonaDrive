# PersonaDrive: A Method for the Identification and Prioritization of Personalized Cancer Drivers
### _This is the original repository for the PersonaDrive paper_

**Installing the dependencies**
```
pip install -r requirements.txt
```

## **Input**

### 1. Personalized Bipartite Networks (PBNs)

There are three input files for the Personalized Bipartite Networks' (PBNs) construction step: the Protein Protein Interaction (PPI) network edges file (STRING netowkrs (v11.5 or v10.5) orDawnRank netowrk), a binary matrix of dysregulated genes (DEGs), and a binary matrix of mutated genes (MUT). Files for each cacner type are located in data folder.

#### 1.1. PPI Networks:
We employ three different interaction networks in our evaluations; STRING network v11.5, the STRING network v10.5 employed in
_(Dinstag and Shamir, 2020)_  and the DawnRank gene interaction network of _(Hou  and  Ma,  2014)_.
The files are located at data/[version]_network.csv

```
gene1 gene2 score
g1 g2 0.9
g3 g8 0.6
```
#### 1.2. Biological pathways:
We employ two different KEGG versions (Kanehisa et al., 2020) for the input set
of biological pathways, the KEGG Release 101 (denoted as v1) and the KEGG pathways used in _Dinstag and Shamir (2020)_ (denoted as v2).
The files are located at data/kegg_pathways_[version].csv

#### 1.3. Mutation Data:

The file is located at data/[cancer]/MUT.csv
```
        p1  p2  ... pn
g1      0   1   ... 1
g2      1   1   ... 0
gx      0   0   ... 1
```
#### 1.4 DEGs Data:

The file is located at data/[cancer]/DEGs.csv
```
        g1      g2     ...  gy
p1      False   True   ...  False
p2      True    True   ...  False
...
```
**Note**: we use the R [code](https://github.com/shahcompbio/drivernet/blob/master/R/getPatientOutlierMatrix.R) from Bashashati et al., (2012) to generate the set of DEGs.

### 2. PersonaDrive Framework (Prioritizing Mutated Genes)
There are two input data for the PersonaDrive framework to prioritize mutated genes in _Bi_ network: the geneerated _.gml_ PBNs' files, and KEGG pathways data retrieved from the supplementary material of _Dinstag and Shamir, (2020)_. The constructed PBNs will be located at graphs/[dataset]/[cancer]_[network]/.

### 3. Evaluation Framewok
#### 3.1 Evaluations with Reference Sets Relevant for Cohort Studies
The personalized reference sets are constructed with respect to several relevant reference sets of known cancer genes: Cancer Gene Census (CGC), Network of Cancer Genes (NCG), and CancerMine. Files are located at data/reference_sets/.

#### 3.2 Evaluations Based on Cell Line Data
For this type of evaluations, for each available cell line we define a novel reference gene set by compiling the target genes of drugs that are found to be sensitive based on data from GDSC _(Yanget al., 2013)_ and DepMap databases for that cell line. Files are located at data/reference_sets/.

#### 3.2 Evaluations Based on Enrichment Analysis
For this type of evaluations, we evaluate the methods based on KEGG and Reactome (Fabregat et al., 2018) enrichment analysis by checking the amounts of overlaps between the pathways enriched signiﬁcantly in the genes output by some personalized prioritization method and those that are enriched in cell line reference sets constructed from drug sensitivity data.
## **Run**

For more details on the execution parameters please refer to the python files.

1. Constructing PBNs:

```
python constructing_PBNs.py -d TCGA -c COAD -n ST11
```

1. Rank Mutated Genes:

```
python PersonDrive.py -d TCGA -c COAD -n ST11
```

2. Evaluation

```
python evaluation.py -d TCGA -c COAD -n ST11
```


## **Outputs**
- The 'constructing_PBNs.py' script will construct the personalized bipartite networks (PBNs).

- The 'PersonDrive.py' script will output the personalized ranking for each sample in the chosen cancer type and dataset.  

- The 'evaluation.py' script will compute the mean precision, recall and F1 scores and plot them.

## **Data Availability**
The used dataset is available from the authors upon request (aissa.houdjedj@antalya.edu.tr).
