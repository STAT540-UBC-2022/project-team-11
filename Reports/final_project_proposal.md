# Host transcriptional response to SARS-Cov2 by sex and age.

## Motivation and background

SARS-CoV-2, a novel coronavirus, has rapidly spread worldwide and has a wide range of clinical manifestations among different population groups. After induction in the host, the pathogen exploits the host's antiviral responses with viral load, and transmissibility peaked at symptom onset [1, 2]. Age and sex play a role in shaping clinical outcomes. In particular, older individuals and males tend to be more susceptible. Lieberman et al. [3] first described the dataset used in this project. Their study aimed to use RNA sequencing profiles from nasopharyngeal (NB) swabs to examine host response gene expression across infected status, viral load, age and sex in a longitudinal data set (time-dependent manner).  

We aim to extend the original research by exploring the interactive effects of age and sex on gene expression profiles. Specifically, we want to examine any statistical difference in gene expression profiles among SARS-CoV-2 infected patients based on age and sex. Additionally, we would like to leverage these differences in machine learning models to predict the strength of host immune response and disease outcome over time. While several studies have used machine learning techniques to predict the onset [4] and mortality risk [5] of COVID-19, few have extended such models to estimate a natural immune response in the host using RNA-Seq longitudinal data.

## Division of labour 

| Name | Background | Degree | Affiliations | Job Assignment | Projected Contributions |
| :-------------: | :-------------: | :-------------: | :-------------: | :-------------: | :-------------: |
| Ekpereka Amutaigwe | Medical Microbiology | M.Sc. | GSAT, Hoodless Lab, Terry Fox Lab, BC Cancer Research Center | Gene Ontology and Machine learning | 25%  |
| Credo Casmil | Biochemistry | M.Sc. | GSAT, Blakney Lab, Michael Smith Labs | QC of data, DE analysis | 25% |
| Dollina Dodani | Computer Science; Molecular biology and biochemistry | M.Sc. | Bioinformatics,  OVCARE, Vancouver General Hospital, Talhouk Lab | Machine Learning and Data Cleaning | 25% |
| Aditi Nagaraj Nallan | B.Sc in Biotechnology, Zoology, Chemistry  |  M.Sc.  | Bioinformatics, LSI, Steven Hallam lab | Data download, DE analysis in R |  25% |

## Datasets
We will work with the raw counts that were generated as follows by the researchers:<br/>
(i) RNA sequencing using the Illumina NextSeq 500 platform.<br/>
(ii) Single end reads were adapter and quality trimmed by Trimmomatic v0.39. <br/>
(iii) Alignement to Ensembl v96 human transcriptome by Kallisto v0.46. <br/>
(iv) Gene level annotations by Biomart.<br/>

The count matrices were supplemented with phenotypic and metadata from the GSE152075 expression set. Though the GSE file can be obtained, the assayData contains 0 features (exprs), and hence we will use the raw counts.

We added a subdirectory named `Datasets` containing the count matrices and phenoData. Specifically, it comprises of:<br/>
1.`GSE152075_raw_counts_GEO.txt` obtained as a supplemental file contains 37784 rows and 484 columns. Rows represent individual gene counts, and the columns represent each sample.<br/>
2.`PhenoData` was obtained from the GSE152075. It has 484 rows and 43 columns that describe study metadata. It shows that total RNA was isolated from nasopharyngeal swabs using the Roche MagNAPure or Qiagen BioRobot. Researchers determined if host-specific gene expression differences correlated to SARS-CoV-2 infection status, host age, sex, and viral load in NP swabs from 430 SARS-CoV-2-infected individuals and 54 negative controls.<br/>

## Specific questions and methodology
To achieve the final goals of our project, we aim to answer:

**1) What genes are significantly differentially expressed between the two groups (healthy and infected) while accounting for factors such as age and sex?**

  Computational and statistical methods used:
  
  * Reading data into R using the GEOquery package
  * Exploratory data analysis and sanity checks using Tidyverse tools
  * Leveraging DESeq2 model to infer differential gene expression
 
**2) Are the significantly differentially expressed genes directly connected to immune response?**

We would first identify the genes associated with immune functions in humans by carrying out a literature review and using databases such as DAVID [5]. We would then compare these genes to our list of DE genes obtained from step 1.

**3) Could we leverage the information from the differentially expressed genes to predict the level of immune response in infected individuals, given their viral load, age, and sex?**

We would develop a machine learning model to classify patients based on the differentially expressed genes and associated metadata. Briefly, we would first use dimensionality reduction techniques such as Principal Component Analysis (PCA) and then apply a supervised learning model (Support Vector Machines) to classify patients. 

## Reason for a switch from initial proposal

We initially proposed studying Type II diabetes (T2D) in three groups of patients (healthy, insulin-sensitive, and insulin-resistant). While we had access to raw reads (in .fastq format), we could not download the prepossessed data or raw counts (in an expression matrix format). We, therefore, moved to the currently proposed data set that was easily accessible in the desired format and would allow us to ask similar statistical questions that we initially proposed. 

## References

[1] Hoffmann M, Kleine-Weber H, Schroeder S, Krüger N, Herrler T, Erichsen S, et al. SARS-CoV-2 Cell Entry Depends on ACE2 and TMPRSS2 and Is Blocked by a Clinically Proven Protease Inhibitor. Cell. 2020. April 16;181(2):271–280.e8. 10.1016/j.cell.2020.02.052

[2] Hou YJ, Okuda K, Edwards CE, Martinez DR, Asakura T, Dinnon KH, et al. SARS-CoV-2 Reverse Genetics Reveals a Variable Infection Gradient in the Respiratory Tract. Cell. 2020. May; S0092867420306759.

[3] Blanco-Melo D, Nilsson-Payant BE, Liu W-C, Uhl S, Hoagland D, Møller R, et al. Imbalanced Host Response to SARS-CoV-2 Drives Development of COVID-19. Cell. 2020. 28;181(5):1036–1045.e9. 10.1016/j.cell.2020.04.026

[4] Zoabi, Yazeed, Shira Deri-Rozov, and Noam Shomron. "Machine learning-based prediction of COVID-19 diagnosis based on symptoms." npj digital medicine 4.1 (2021): 1-5.

[5] Mahdavi, Mahdi, et al. "A machine learning based exploration of COVID-19 mortality risk." Plos one 16.7 (2021): e0252384.

[6] Huang DW, Sherman BT, Lempicki RA. Systematic and integrative analysis of large gene lists using DAVID Bioinformatics Resources. Nature Protoc. 2009;4(1):44-57.

