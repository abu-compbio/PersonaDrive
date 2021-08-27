# PersonaDrive: A Method for the Identification and Prioritization of Personalized Cancer Drivers
### _This is the original repository for the PersonaDrive paper_

**Installing the dependencies**
```
pip install -r requirements.txt
```

## **Input**

### 1. Personalized Bipartite Networks (PBNs)

There are three input files for the Personalized Bipartite Networks' (PBNs) construction step: the Protein Protein Interaction (PPI) network edges file (DawnRank or STRING netowkrs), a binary matrix of dysregulated genes, and a binary matrix of mutated genes. Files for each cacner type are located in data folder.

#### 1.1. PPI Networks:
We employ two different interaction networks in our evaluations; the STRING network employed in
_(Dinstag and Shamir, 2020)_ and the DawnRank gene interaction network of _(Hou  and  Ma,  2014)_.
The files are located at data/ppi_[version].csv

```
gene1 gene2 score
g1 g2 0.9
g3 g8 0.6
```

#### 1.2. Mutation Data:

The file is located at data/[cancer]/mutation_data.csv
```
        p1  p2  ... pn
g1      0   1   ... 1
g2      1   1   ... 0
gx      0   0   ... 1
```
#### 1.3 DEGs Data:

The file is located at data/[cancer]/DEGs_data.csv
```
        g1      g2     ...  gy
p1      False   True   ...  False
p2      True    True   ...  False
...
```

### 2. PersonaDrive Framework (Prioritizing Mutated Genes)
There are two input data for the PersonaDrive framework to prioritize mutated genes in _Bi_ network: the geneerated _.gml_ PBNs' files, and KEGG pathways data retrieved from the supplementary material of (Dinstag and Shamir, 2020). The constructed PBNs will be located at graphs/[dataset]/[cancer]/.

### 3. Evaluation Framewok
#### 3.1 Evaluations with Reference Sets Relevant for Cohort Studies
The personalized reference sets are constructed with respect to several relevant reference sets of known cancer genes: Cancer Gene Census (CGC), Network of Cancer Genes (NCG), and CancerMine. Files are located at data/reference_sets/.

#### 3.2 Evaluations Based on Cell Line Data
For this type of evaluations, for each available cell line we define a novel reference gene set by compiling the target genes of drugs that are found to be sensitive based on data from GDSC _(Yanget al., 2013)_ and DepMap databases for that cell line. Files are located at data/reference_sets/.

#### 3.3 Evaluations Based on Enriched Pathway Overlap _EPO_
For this type of evaluations, we find the set of enriched KEGG or Reactome pathwaysusing the g:GOSt tool (the core of g:Profiler tool) (Raudvereet al., 2019).

## **Run**

For more details on the execution parameters please refer to the python files.

1. Constructing PBNs:

```
python constructing_PBNs.py -d TCGA -c COAD
```

1. Rank Mutated Genes:

```
python PersonDrive.py -d TCGA -c COAD
```

2. Evaluation

```
python evaluation.py -d TCGA -c COAD
```


## **Outputs**
- The 'constructing_PBNs.py' script will construct the personalized bipartite networks (PBNs).

- The 'PersonDrive.py' script will output the personalized ranking for each sample in the chosen cancer type and dataset.  

- The 'evaluation.py' script will compute the mean precision, recall and F1 scores and plot them.

## **Data Availability**
The used dataset is available from the authors upon request (aissa.houdjedj@antalya.edu.tr).
