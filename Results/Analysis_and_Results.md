# Technical Report - Analysis and Results 

## Project aims 

## Data exploration and manipulation 

## RNA-Seq analysis

## GO ontology analysis

## Machine learning analysis

Our objective for the machine learning analysis was to be able to predict infection status of a patient, given their transcriptomic data. In order to accomplish this, we used design matrix 3 (described in RNA-Seq section above) that examines the interactive effects of age and sex. We filtered the DEGs to only include the top 100 genes (sorted by logFC) irrespective of their significance. This is because design matrix 3 provided only 2 significantly DEGs which is not enough to train a model. We picked a value of 100 to adhere to our computational limitations. 

We trained two models for our classification problem:

1) K Nearest Neighbors (KNN): 

2) Logistic regression 



## Conclusions 
