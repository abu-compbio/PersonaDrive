# PersonaDrive: A Method for the Identification and Prioritization of Personalized Cancer Drivers
### _This is the original repository for the PersonaDrive paper_

The Data Folder contains all the necessary input data to reproduce the results discussed in our main paper.

## **Somatic Mutations and Gene Expressions**

Files are located at data/[cancer]_[dataset]/. Where dataset is defined as TCGA or CCLE.

### Data Sources
   - TCGA: We employed the TCGAbiolinks R package to compile the relevant data pertaining to the samples from COAD, LUAD, and BRCA datasets (we retrieved both FPKM and raw read counts datasets).

   - DepMap: We employed the DepMap database version 20Q1 (Ghandi et al., 2019) to gather the mutations and gene expression data of the cancer cell line samples.

   Regarding the expression data of normal cell line samples, we employ the normal human bronchial epithelial (NHBE) cell line from the Expression Atlas database of Petryszak et al. (2016) for LUAD and the CCD-18Co human normal colon myofibroblasts data from (Ferrer-Mayorga et al. 2019) for COAD.

### Data Preprocessing
   1 - The gene expression values in the original dataset are in terms of FPKM values, as an input to our method, we convert the FPKM values to TPM values as follows:  TPMi= (FPKMi / ∑FPKMj) * 10^6

   2 - We use TPM values as an input to PersonaDrive and DawnRank. However, raw read counts values for both normal and tumor samples are used as input for PRODIGY to perform differential expression analysis using DeSEQ2. Among the evaluated methods, only SCS requires expression data from paired normal and tumor samples. To evaluate all the employed methods on a larger dataset, we extract differentially expressed genes for SCS from unpaired gene expression data by retrieving the log-fold change values using DeSEQ2 R package.

   Note that PersonaDrive and DawnRank assume that TPM values are in log domain  to perform differential expression analysis.
   3 - Since RNA-seq TPM expression values for CCLE datasets retrieved from DepMap are already in log domain, we only apply the log2 transformation for TCGA expression values and expression data of normal cell line samples using a pseudo-count of 1.

   4 - For CCLE dataset, we pre-process both normal and tumor data matrices to contain only genes G_{T,N} that exist in both the tumor and normal data as follows: G_{T,N}= G_T ∩ G_N.

   5 - For both the TCGA and CCLE datasets, we filter out the silent mutations from the somatic mutation data.



## **PPI Networks**
We employ three different interaction networks in our evaluations: STRING network v11.5, the STRING network v10.5 employed in (Dinstag and Shamir, 2020), and the DawnRank gene interaction network of (Hou and Ma, 2014). The files are located at data/[version]_network.csv

STRING network v11.5: string_v11.5_network.csv
STRING network v10.5: string_v10.5_network.csv
DawnRank network: dawnrank_network.csv


## **Biological pathways**
We employ two different KEGG versions (Kanehisa et al., 2020) for the input set of biological pathways, the KEGG Release 101 (denoted as v1) and the KEGG pathways used in Dinstag and Shamir (2020) (denoted as v2). The files are located at data/kegg_pathways_[version].csv.


## **Reference sets**
The personalized reference sets are constructed with respect to several relevant reference sets of known cancer genes.

### Reference sets relevant for cohort studies
  1 - CGC: (the file is located at data/reference_sets/CGC.csv).
  2 - NCG: The Network of Cancer Genes (the file is located at data/reference_sets/NCG6_cancergenes.tsv).
  3 - CancerMine, uses text-mining to catalogue cancer-associated genes where it also extracts information about the type of the cancer (the files are located at data/reference_sets/cancermine_[cancer].txt)

  *** The number of genes in each reference is provided in Supplementary Table 5.

  *** (Please read the main paper for more details)

### Evaluations based on cell line data
  We construct potential personalized reference sets based on data specific to each sample, we propose the use of cell line data coupled with drug sensitivity data. For this type of evaluation, for each available cell line we define a novel reference gene set by compiling the target genes of drugs that are found to be sensitive based on data from GDSC   (Yang et al., 2013) and DepMap databases (the files are located at data/reference_sets/[cancer]_DepMap_GDSC2.txt).

  *** The statistical information of the targets in GDSC and the DepMap reference sets for each cancer type under study are available in Supplementary Tables 2-4.

  *** (Please read the main paper for more details)

## **Enrichment analysis**

  We transformed the KEGG Release 101 (denoted as v1) to a ".gmt" matrix to be used as an input to g:GOSt tool (the core of g:Profiler tool).

  file is located at data/enrichment_analysis/kegg_pathways_v1.gmt.

  *** (Please read the main paper for more details)


## **References**

  - Ghandi,   M.  et  al.  (2019).    Next-generation  characterization  of  the  cancer  cell  line encyclopedia.  Nature,  569(7757),  503–508.

  - Petryszak, R. et al. (2016). Expression Atlas update—an integrated database of gene and protein expression in humans, animals and plants. Nucleic Acids Research, 44(Database issue), D746–D752.

  - Ferrer-Mayorga, G. et al. (2019). Vitamin D and Wnt3A have additive and partially overlapping modulatory effects on gene expression and phenotype in human colon ﬁbroblasts. Scientiﬁc Reports, 9(1), 8085.

  - Dinstag, G. and Shamir, R. (2020). PRODIGY: personalized prioritization of driver genes. Bioinformatics, 36(6), 1831–1839.

  - Hou, J. P. and Ma, J. (2014). DawnRank: discovering personalized driver genes in cancer. Genome Medicine, 6(7), 56.

  - Kanehisa, M. et al. (2020). KEGG: integrating viruses and cellular organisms. Nucleic Acids Research, 49(D1), D545–D551.

  - Yang, W. et al. (2013). Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells. Nucleic Acids Research, 41(Database issue), D955–961.
