# Technical Report - Analysis and Results 

## Project aims 

## Data exploration and manipulation 

## RNA-Seq analysis

## GO ontology analysis

## Machine learning analysis

Our objective for the machine learning analysis was to be able to predict infection status of a patient, given their transcriptomic data. In order to accomplish this, we used design matrix 3 (described in RNA-Seq section above) that examines the interactive effects of age and sex. We filtered the DEGs to only include the top 100 genes (sorted by logFC) irrespective of their significance. This is because design matrix 3 provided only 2 significantly DEGs which is not enough to train a model. We picked a value of 100 to adhere to our computational limitations. 

Our data set consisting of 467 samples, 100 genes, 2 metadata variables (age and sex), and the predictor variabe (infection status). We first split our data set into a training set (80%) and a validation set (20%).

We trained two models for our classification problem:

1) K Nearest Neighbors (KNN): This model examines the k nearest points (based on Euclidean distance) to make a classification. We first scaled all explanatory variables of our data set K is user-defined and must be fine-tuned before fitting the model. We tested 20 values of k (1 to 20) and noticed that the classification error on the test set was minimum when k is 5, 7 or between 8 to 15. We decided to fit our model with k=15 using the `knn` function from `class` package to avoid over fitting of the model. 

![](../Results/Plots/ML.png)

2) Logistic regression: This model assigns a probability (0 to 1) to an unknown sample of belonging to a class (infection status positive or negative). Since we had 102 explanatory variables, we were not able to cope with the computational burden. We then decided to perform feature selection analysis using step wise logistic regression. We build a null model and a full model consisting of a constant term and all genes, respectively. Using this, we identified 3 genes (EFS, MTCO2P22 and EXOC5P1) that significantly contributed to the classification (explained most of the variance in the data set). We used the `glm` function to fit our model. 

**Model selection**
The two models acquired the follwing accuracy metrics: 

| Model         | Sensitivity   | Specificity  | Accuracy |
| ------------- |:-------------:| :-----:| :--------------: |
| KNN           | 0.943 | 1.00 | 0.947 |
| Logistic regression      | 0.965 | 1.00 | 0.968 |

We can see that while both models have a specificity of 0 (that is, no COVID-19 positive patient was misclassified as negative), they differed in accuracy scores. Logistic regression is more accurate and sensitive while making a prediction based on only three genes (compared to KNN that uses 100 genes). Thus, logistic regression would be a slightly better model.

## Conclusions 
