# Technical Report - Analysis and Results 

## Project aims 

## Data exploration and manipulation 

## RNA-Seq analysis

Before creating, the DGEList, `Combat_seq()` is utilized to adjust the count matrix for batch effects. This adjust count matrix as well the modified metadata (including the imputed sex annotations) are used to build the DGEList. 

After TMM normalization and filtering of lowly expressed genes (cpm > 10 for at least 2 samples), data visualization is perfomed. 

The boxplot show the distribution of gene expression levels of 20 random samples. Only a small number of samples have different distributions i.e. sample 167 and 169.
![Box plot](Results/Plots/box_plot.jpg)  

Additional pairwise correlation between samples was ploted on a heat map. Based on the clustering, higher correlation was present in both negative and positive SARS tested samples as long as they were within the senior and adult age category. This was similar to what liberman et al. (2020) observed.
![Heat map](Results/Plots/heatmap.png)

To test our 1st hypothesis, we examined if interaction between infection status, age category and sex would result in differentially expressed genes. To acomplish this, 3 design matrices were made:
Matrix 1: `Age * infection status`
Matrix 2: `Infection status * age`
Matrix 3: `Sex * Age`

Applying a cutoff of |1| for the log fold change and adjust p-value of 0.05: 

- Interaction between infection status and age which resuted in 126 and 15 down and up DEGs, respectively.

![](Results/Plots/MD_1.png)

- Interaction between infection status and sex which resulted in 8 and 7 down and up DEGs, respectively.

![](Results/Plots/MD_2.png)

- Interaction between age and sex which resulted in 0 and 2 down and up DEGs, respectively. 

![](Results/Plots/MD_3.png)


## GO ontology analysis

## Machine learning analysis

Our objective for the machine learning analysis was to be able to predict infection status of a patient, given their transcriptomic data. In order to accomplish this, we used design matrix 3 (described in RNA-Seq section above) that examines the interactive effects of age and sex. We filtered the DEGs to only include the top 100 genes (sorted by logFC) irrespective of their significance. This is because design matrix 3 provided only 2 significantly DEGs which is not enough to train a model. We picked a value of 100 to adhere to our computational limitations. 

We trained two models for our classification problem:

1) K Nearest Neighbors (KNN): 

2) Logistic regression 



## Conclusions 
