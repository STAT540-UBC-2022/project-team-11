# Technical Report - Analysis and Results 

## Project aims 
Using data published by Liberman et al, that consists of RNA-Seq as well as associated metadata, we wanted to answer:

**1) What genes are significantly differentially expressed between the two groups (healthy and infected) while accounting for factors such as age and sex?**

**2) Are the differentially expressed genes (identified in aim 1) directly connected to immune response?**

**3) Could we leverage RNA-Seq data to predict infection status of an individual?**

## Data exploration and manipulation 

## RNA-Seq analysis

Before creating, the DGEList, `Combat_seq()` is utilized to adjust the count matrix for batch effects. This adjust count matrix as well the modified metadata (including the imputed sex annotations) are used to build the DGEList. 

After TMM normalization and filtering of lowly expressed genes (cpm > 10 for at least 2 samples), data visualization is perfomed. 

The boxplot show the distribution of gene expression levels of 20 random samples. Only a small number of samples have different distributions i.e. sample 167 and 169.

![](https://github.com/STAT540-UBC-2022/project-team-11/blob/580decab757e78a930b6954cbb1988a954c6f572/Results/Plots/box_plot.jpg)

Additional pairwise correlation between samples was ploted on a heat map. Based on the clustering, higher correlation was present in both negative and positive SARS tested samples as long as they were within the senior and adult age category. This was similar to what liberman et al. (2020) observed.

![](https://github.com/STAT540-UBC-2022/project-team-11/blob/dc53ae687dd40eb342cfcfafd5c49cee6f15f213/Results/Plots/heatmap.png)

To test our 1st hypothesis, we examined if interaction between infection status, age category and sex would result in differentially expressed genes. To acomplish this, 3 design matrices were made:
Matrix 1: `Age * infection status`.
Matrix 2: `Infection status * age`.
Matrix 3: `Sex * Age`.

Applying a cutoff of |1| for the log fold change and adjust p-value of 0.05: 

- Interaction between infection status and age which resuted in 126 and 15 down and up DEGs, respectively.

![](https://github.com/STAT540-UBC-2022/project-team-11/blob/263abeff2f45447b7568219d6215f820a1164090/Results/Plots/MD_1.png)

- Interaction between infection status and sex which resulted in 8 and 7 down and up DEGs, respectively.

![](https://github.com/STAT540-UBC-2022/project-team-11/blob/263abeff2f45447b7568219d6215f820a1164090/Results/Plots/MD_2.png)

- Interaction between age and sex which resulted in 0 and 2 down and up DEGs, respectively. 

![](https://github.com/STAT540-UBC-2022/project-team-11/blob/263abeff2f45447b7568219d6215f820a1164090/Results/Plots/MD_3.png)


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

Using the data first described by Liberman et al., we were able to accomplish our three goals. Using EdgeR we identified significantly DEGs while examining interactive effects of [infection status and age], [infection status and sex], and [age and sex]. We further leveraged GO ontology analysis to identify genes that were directly related to the immune response system in humans. We concluded that there were no significant associations to immune function using the interaction of age and sex. Lastly, we successfully build a classification model (with an accuracy of 97%) that utilizes three genes (EFS, MTCO2P22 and EXOC5P1) to predict if an individual has COVID-19.
As future work, we would like to explore the biological mechanisms of why our model (specially the logistic regression generated such high accuracy scores). This could be done by a more detailed pathway analysis (using GO) or wet lab techniques such a knock out models. 

