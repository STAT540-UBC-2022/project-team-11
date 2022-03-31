Host transcriptional response to SARS-Cov2 by Sex and Age
========================================================
author: Ekpereka Amutaigwe, Credo Casmil, Dollina Dodani, Aditi Nallan
date: 30 March 2022
autosize: true

Background
========================================================

- SARS-CoV-2 has shown a wide range of clinical manifestations among different population groups.
- Clinically, COVID-19 cases tend to be more severe for older adults and males
- [Lieberman et al. (Plos Biology 2020)](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000849) examined host response gene expression across infection status and viral load from 430 individuals with SARS-CoV-2 and 54 negative controls.

![Figure from Lieberman et al. (Plos Biology 2020) paper](Heatmap_paper.png)

Project Aim
========================================================

- Examine the statistical difference in gene expression profiles among SARS-CoV-2 infected patients by extending the original research to explore the interactive effects of age and sex.

- Leverage these differences in a Machine Learning model to predict infection state.


Dataset
========================================================



- Obtained from Lieberman et al. (Plos Biology 2020) paper as:
  - **Main data file**: expression values - one per gene, per sample
  - **Metadata file**: several covariates/experimental design variables for each sample. 


```r
# Expression matrix
dim(count_mat)
```

```
[1] 35784   484
```

```r
# Metadata
dim(m)
```

```
[1] 484   6
```


Metadata
========================================================


```r
options(height = 100)
summary(m)
```

```
    Sample               Age        Gender      Viral_load    Sars_test
 Length:484         Min.   : 2.0   F   :231   Min.   : 0.00   neg: 54  
 Class :character   1st Qu.:39.5   M   :200   1st Qu.:17.73   pos:430  
 Mode  :character   Median :54.0   NA's: 53   Median :20.44            
                    Mean   :54.6              Mean   :18.76            
                    3rd Qu.:70.5              3rd Qu.:23.54            
                    Max.   :90.0              Max.   :30.54            
                    NA's   :17                NA's   :17               
     Batch    
 Q      : 57  
 P      : 54  
 L      : 43  
 O      : 28  
 J      : 26  
 F      : 24  
 (Other):252  
```


Visualization of missing data
========================================================

.pull-left[
![plot of chunk unnamed-chunk-4](test-figure/unnamed-chunk-4-1.png)

]

.pull-right[

```r
#Check for percentage of missing data
pct_miss(m)
```

```
[1] 2.995868
```

```r
#Check for percent rows that have missing data
pct_miss_case(m)
```

```
[1] 10.95041
```
]
