# Host transcriptional response to SARS-Cov2 by sex and age.

## Motivation and Background

**Observations**: SARS-CoV-2, a novel coronavirus has rapidly spread around the world and has shown a wide range of clinical manifestations among different population groups. It has been previously shown that the entry of SARS-CoV-2 into the host cells depends on its binding to angiotensin-converting enzyme 2 (ACE2) receptors [1]. It is further induced upon exposure to interferon which suggests that SARS-CoV-2 exploits the host's antiviral responses [2]. Based on these findings it has been hypothesized that SARS-CoV-2 employs mechanisms similar to the infection of bronchial epithelial cells at low multiplicity of infection (MOI) which does not result in extensive transcription of interferon-stimulated genes at 24 hours post infection [3]. An important consequence of this is that the viral load and transmissibility peaks at the time of symptom onset [4]. It has also been shown that age and sex play a role is shaping the clinical outcomes. Older individuals and males tend to be more susceptible. The current study that forms the basis of our project delineates the host gene expression profiles based on viral load in a time dependent manner while accounting for sex and age but does not show how the interactive effects of sex and age shape the gene expression profiles. Additionally there is a lack of understanding of how these expression profiles can be used to predict disease outcomes over time. 

**Research Question**: Is there a difference in gene expression profiles among SARS-CoV-2 infected patients based on age and sex? Can these differences be leveraged in machine learning models to predict strength of host immune response and disease outcome over time? 

# Division of Labour 

| Name | Background | Degree | Affiliations | Job Assignment | Projected Contributions |
| :-------------: | :-------------: | :-------------: | :-------------: | :-------------: | :-------------: |
| Ekpereka Amutaigwe | Medical Microbiology | M.Sc. | GSAT, Hoodless Lab, Terry Fox Lab, BC Cancer Research Center | Gene Ontology and Machine learning | 25%  |
| Credo Casmil | Biochemistry | M.Sc. | GSAT, Blakney Lab, Michael Smith Labs | QC of data, DE analysis | 25% |
| Dollina Dodani | Computer Science; Molecular biology and biochemistry |  | Bioinformatics,  OVCARE, Vancouver General Hospital | Machine Learning and Data Cleaning | 25% |
| Aditi Nagaraj Nallan | B.Sc in Biotechnology, Zoology, Chemistry  |  Masters  | Bioinformatics, LSI, Steven Hallam lab | Data download, DE analysis in R |  25% |

# Datasets
We will work with the raw counts that were generated as follows by the researchers:
(i) RNA sequencing using the Illumina NextSeq 500 platform.
(ii) Single end reads were adapter and quality trimmed by Trimmomatic v0.39.
(iii) Alignement to Ensembl v96 human transcriptome by Kallisto v0.46.
(iv) Gene level annotations by Biomart.

The count matrices will be supplemented with phenotypic and meta data from the GSE152075 expression set. Though the GSE file can be obtained, the assayData contains 0 features (exprs) hence the use of the raw counts.

A subdirectory named `Datasets` has been added. This contains the count matrices and phenoData. It contains the following:
1. The `GSE152075_raw_counts_GEO.txt` which was obtained as a supplemental file contains 37784 rows and 484 columns. The rows represent individual gene counts while the columns represent each sample.
2. The `PhenoData` was obtained from the GSE152075. It has 484 rows and 43 columns that describe the metadata of the study and the following can be obtained from it:
  - It shows that total RNA was isolated from nasopharyngeal swabs using the Roche MagNAPure or Qiagen BioRobot.
  - Researchers determined if host-specific gene expression differences   were correlated to SARS-CoV-2 infection status, host age, sex, and viral load in NP swabs from 430 SARS-CoV-2-infected individuals and 54 negative controls.

## Specific questions and Methodology
In order to achieve the final goals of our project, we aim to answer:

**1) What genes are significantly differentially expressed between the two groups (healthy and infected) while accounting for factors such as age, sex?**

  Computational and statistical methods used:
  
  * Reading data into R using the GEOquery package
  * Exploratory data analysis and sanity checks using Tidyverse tools
  * Leveraging DESeq2 model to infer differential gene expression
 
**2) Are the significantly differentially expressed genes directly connected to immune response?**

We would first identify the genes that are associated with immune functions in humans by carrying out a literature review and using databases such as DAVID [5] and this will be followed by comparing these genes to our list of DE genes obtained from step 1.

**3) Could we leverage the information from the genes that are differentially expressed to predict the level of immune response in infected individuals, given their viral load, age, and sex?**

To answer this question, we would develop a machine learning model to classify patients based on the differentially expressed genes and associated metadata. Briefly, we hope to first use dimensionality reduction techniques such as Principal Component Analysis (PCA), we would then apply a supervised learning model (Support Vector Machines) to classify patients. 

## Reason for switch from inital proposal

Our original proposal was based on analyzing Type II diabetes (T2D) in three groups of patients (healthy, insulin-sensitive, and insulin-resistant). While we had access to raw reads (in .fastq format), we were unable to download the prepossessed data or raw counts (in an expression matrix format). We therefore moved to the currently proposed data set that was easily accessible in the desired format and would allow us to ask similar statistical questions that we initially proposed. 

## References

[1] Hoffmann M, Kleine-Weber H, Schroeder S, Krüger N, Herrler T, Erichsen S, et al. SARS-CoV-2 Cell Entry Depends on ACE2 and TMPRSS2 and Is Blocked by a Clinically Proven Protease Inhibitor. Cell. 2020. April 16;181(2):271–280.e8. 10.1016/j.cell.2020.02.052

[2] Hou YJ, Okuda K, Edwards CE, Martinez DR, Asakura T, Dinnon KH, et al. SARS-CoV-2 Reverse Genetics Reveals a Variable Infection Gradient in the Respiratory Tract. Cell. 2020. May;S0092867420306759.

[3] Blanco-Melo D, Nilsson-Payant BE, Liu W-C, Uhl S, Hoagland D, Møller R, et al. Imbalanced Host Response to SARS-CoV-2 Drives Development of COVID-19. Cell. 2020. 28;181(5):1036–1045.e9. 10.1016/j.cell.2020.04.026

[4] Zou L, Ruan F, Huang M, Liang L, Huang H, Hong Z, et al. SARS-CoV-2 Viral Load in Upper Respiratory Specimens of Infected Patients. N Engl J Med. 2020. March 19;382(12):1177–9. 10.1056/NEJMc2001737 

[5] Huang DW, Sherman BT, Lempicki RA. Systematic and integrative analysis of large gene lists using DAVID Bioinformatics Resources. Nature Protoc. 2009;4(1):44-57.


