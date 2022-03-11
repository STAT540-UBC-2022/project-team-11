Data Cleaning and initial DE analysis
================
Credo Casmil, Aditi Nagaraj Nallan, Ekpereka Amutaigwe and Dollina
Dodani
2022-03-10

# Load packages

``` r
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(statmod))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggbiplot))
suppressPackageStartupMessages(library(ggplot2))
```

``` r
eset <- getGEO("GSE152075", getGPL = FALSE)[[1]]
getGEOSuppFiles("GSE152075")

#Read in the count matrix
raw_counts <- read.csv("GSE152075/GSE152075_raw_counts_GEO.txt", sep = "")
count_mat <- as.matrix(raw_counts,row.names="gene_id")
```

``` r
#Read in the metadata
pdata <- eset@phenoData@data %>% as_tibble()

#Clean the metadata to include only columns of interest for the initial analysis
pdata_clean = pdata %>%
  select(title, `age:ch1`, characteristics_ch1.3, characteristics_ch1.4, `sars-cov-2 positivity:ch1`)
colnames(pdata_clean) = c("Title", "Age", "Gender", "Batch", "Sars_test")

#Modify the cleaned data set as per further DE analysis
pdata_mod <- pdata_clean %>%
  subset(Age != "Unknown") %>%
  transform(Age = as.numeric(Age))
```

    ## Warning in eval(substitute(list(...)), `_data`, parent.frame()): NAs introduced
    ## by coercion

``` r
colnames(pdata_mod) = c("Title", "Age", "Gender", "Batch", "Sars_test")
pdata_mod[is.na(pdata_mod)] <- 90

pdata_mod = pdata_mod %>%
  mutate(Age_category =  case_when(
    Age < 18 ~ "Child",
    Age >= 18 & Age < 35 ~ "Young Adult",
    Age >= 35 & Age < 65  ~ "Adult",
    TRUE ~ "Senior"
  ))

data = pdata %>%
       select(title, `age:ch1`) %>% 
       subset(`age:ch1` == "Unknown")

drop <- data$title
count_mat_final  = count_mat[,!colnames(count_mat) %in% drop]
```

# Using EdgeR

``` r
#Creating DGE List 
dge <- DGEList(counts = count_mat, samples = pdata_clean, group = pdata_clean$Sars_test)
dim(dge)
```

    ## [1] 35784   484

``` r
head(apply(dge$counts, 2, sum)) # total gene counts per sample
```

    ## POS_001 POS_002 POS_003 POS_004 POS_005 POS_006 
    ## 1230730 2519860 2849945 2130923 4445875 7467886

``` r
#removing lowly expressed genes
keep_edge <- rowSums(cpm(dge)>100) >= 50
dge_mod <- dge[keep_edge,]
dim(dge_mod)
```

    ## [1] 699 484

``` r
dge$samples$lib.size <- colSums(dge$counts) # Reset library sizes
```

``` r
#Calculating TMM normalization factors and directly adding them to the DGEList
dge_norm = calcNormFactors(dge_mod, method = "TMM")

cpm = cpm(dge_norm, log = FALSE, normalized.lib.sizes = TRUE)
log2cpm = log2(cpm + 1)

#Transforming object from wide to long format for plotting after randomly subsetting the data
ran_samp <- subset(log2cpm[,150:170])

longExpr = ran_samp %>% 
           as.data.frame() %>% 
           rownames_to_column("gene") %>%
           pivot_longer(cols = !gene,
                        values_to = "Expression",
                        names_to = "sample_ID")

#Density plot for a random subset of data for ease of visualization 
longExpr  %>% 
  ggplot(aes(x = Expression, color = sample_ID)) +
  geom_density() +
  labs( x = "Expression", y = "Density", title = "Density plot showing distribution of gene expression across 20 random samples")
```

![](DataCleanup_DE_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
#Box plot for a random subset of data for ease of visualization 
longExpr %>% 
  ggplot(aes(x = sample_ID, y = Expression)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  labs( x = "Sample ID", y = "Gene Expression", title = "Box plot showing distribution of gene expression across 20 random samples")
```

![](DataCleanup_DE_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
#Calculate correlation using log2 transformed random subset CPM values from earlier
cormat = round(cor(ran_samp), 2)

#Plot heatmap
pheatmap(cormat, border_color = NA, cluster_rows = TRUE, cellheight=9, cellwidth = 9)
```

![](DataCleanup_DE_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
#Setting up model matrix: Batch corrected and Batch not corrected
designMatrix1 = model.matrix(~Sars_test , data = dge_norm$samples)
designMatrix2 = model.matrix(~Sars_test + Batch , data = dge_norm$samples)
```

``` r
#Calculation of variance weights and generation of mean-variance trend plot
de_model1 = voom(dge_norm, design = designMatrix1, plot = TRUE)
```

![](DataCleanup_DE_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
de_model2 = voom(dge_norm, design = designMatrix2, plot = TRUE)
```

![](DataCleanup_DE_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

``` r
dge_disp1 <- estimateDisp(dge_norm, designMatrix1, robust = TRUE)
range(dge_disp1$prior.df)
```

    ## [1] 3.648507 7.814706

``` r
plotBCV(dge_disp1,  cex=0.5)
```

![](DataCleanup_DE_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
dge_disp2 <- estimateDisp(dge_norm, designMatrix2, robust = TRUE)
range(dge_disp2$prior.df)
```

    ## [1] 5.108595 9.173370

``` r
plotBCV(dge_disp2,  cex=0.5)
```

![](DataCleanup_DE_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

``` r
lfit1 <- glmFit(dge_disp1, designMatrix1)
lrt1 <- glmLRT(lfit1, coef = "Sars_testpos")
toptags1 <- topTags(lrt1)$table %>% signif(3)

lfit2 <- glmFit(dge_disp2, designMatrix2)
lrt2 <- glmLRT(lfit2, coef = "Sars_testpos")
toptags2 <- topTags(lrt2)$table %>% signif(3)

#Filtering for up and downregulated genes
lrt2_filter_upreg <- filter(lrt2$table, PValue < 0.05, logFC >1)
dim(lrt2_filter_upreg)
```

    ## [1] 45  4

``` r
lrt2_filter_downreg <- filter(lrt2$table, PValue < 0.05, logFC < -1)
dim(lrt2_filter_downreg)
```

    ## [1] 21  4

``` r
de1 <- decideTestsDGE(lrt1, adjust.method="BH", p.value = 0.05)
de1tags1 <- rownames(dge_disp1)[as.logical(de1)]
plotSmear(lrt1, de.tags=de1tags1)
```

![](DataCleanup_DE_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
de1_summary <- summary(de1)
de1_summary
```

    ##        Sars_testpos
    ## Down            261
    ## NotSig          226
    ## Up              212

``` r
de2 <- decideTestsDGE(lrt2, adjust.method="BH", p.value = 0.05)
de1tags2 <- rownames(dge_disp2)[as.logical(de2)]
plotSmear(lrt2, de.tags=de1tags2)
```

![](DataCleanup_DE_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

``` r
de2_summary <- summary(de2)
de2_summary
```

    ##        Sars_testpos
    ## Down             21
    ## NotSig          635
    ## Up               43

# Principal Component Analysis

``` r
# Convert rownames to column
lrt2_upreg_mod <- lrt2_filter_upreg %>% 
                  rownames_to_column("gene")

# Filter data to retain only the upregulated genes
log2cpm_dat <-  log2cpm %>% 
                as.data.frame() %>% 
                rownames_to_column("gene") %>% 
                filter(gene %in% lrt2_upreg_mod$gene) %>% 
  # Change data to a long form
                pivot_longer(cols = !gene,
                             values_to = "Expression",
                             names_to = "sample_ID")

# Rename column to join by
pdata_clean <- dplyr::rename(pdata_clean, sample_ID = Title)

# Join expression and metadata datasets
DEG_new <- log2cpm_dat %>% 
           left_join(pdata_clean, 
                     by = "sample_ID") 

# Transform data back to a wide format
DEG_new_trans <- pivot_wider(DEG_new, 
                             id_cols = c(sample_ID, Age, Gender, Batch, Sars_test), 
                             names_from = gene, 
                             values_from = Expression)

# Use to drop missing values incase there's any before PCA
DEG_new_trans2 <- DEG_new_trans %>% 
                  drop_na() 

# Perform PCA
pca_DEG_new_trans2 <- prcomp(DEG_new_trans2[,-c(1, 2, 3, 4, 5, 6), 
                                            center = TRUE, 
                                            scale = TRUE])
# See what PCA result looks like
summary(pca_DEG_new_trans2)
```

    ## Importance of components:
    ##                           PC1    PC2     PC3     PC4     PC5     PC6     PC7
    ## Standard deviation     8.4766 3.8123 2.72281 2.20963 1.71975 1.64250 1.58904
    ## Proportion of Variance 0.5314 0.1075 0.05483 0.03611 0.02187 0.01995 0.01867
    ## Cumulative Proportion  0.5314 0.6389 0.69368 0.72979 0.75166 0.77161 0.79029
    ##                            PC8     PC9   PC10    PC11    PC12   PC13    PC14
    ## Standard deviation     1.38595 1.32163 1.2846 1.26467 1.20597 1.1917 1.10633
    ## Proportion of Variance 0.01421 0.01292 0.0122 0.01183 0.01076 0.0105 0.00905
    ## Cumulative Proportion  0.80449 0.81741 0.8296 0.84144 0.85220 0.8627 0.87175
    ##                           PC15    PC16    PC17   PC18    PC19    PC20    PC21
    ## Standard deviation     1.07343 1.03128 1.01882 1.0000 0.97823 0.95343 0.93600
    ## Proportion of Variance 0.00852 0.00787 0.00768 0.0074 0.00708 0.00672 0.00648
    ## Cumulative Proportion  0.88027 0.88814 0.89581 0.9032 0.91028 0.91701 0.92349
    ##                           PC22    PC23    PC24   PC25    PC26    PC27    PC28
    ## Standard deviation     0.89360 0.87567 0.84877 0.8385 0.80795 0.79137 0.77922
    ## Proportion of Variance 0.00591 0.00567 0.00533 0.0052 0.00483 0.00463 0.00449
    ## Cumulative Proportion  0.92939 0.93506 0.94039 0.9456 0.95042 0.95505 0.95954
    ##                           PC29    PC30    PC31    PC32    PC33    PC34    PC35
    ## Standard deviation     0.76670 0.73707 0.71068 0.66445 0.63508 0.59826 0.59407
    ## Proportion of Variance 0.00435 0.00402 0.00374 0.00327 0.00298 0.00265 0.00261
    ## Cumulative Proportion  0.96389 0.96790 0.97164 0.97490 0.97789 0.98053 0.98314
    ##                           PC36    PC37    PC38    PC39    PC40    PC41   PC42
    ## Standard deviation     0.58516 0.55085 0.54100 0.53512 0.50332 0.48730 0.4511
    ## Proportion of Variance 0.00253 0.00224 0.00216 0.00212 0.00187 0.00176 0.0015
    ## Cumulative Proportion  0.98567 0.98792 0.99008 0.99220 0.99407 0.99583 0.9973
    ##                          PC43    PC44
    ## Standard deviation     0.4351 0.41358
    ## Proportion of Variance 0.0014 0.00126
    ## Cumulative Proportion  0.9987 1.00000

``` r
# What are the components of the PCA result
str(pca_DEG_new_trans2)
```

    ## List of 5
    ##  $ sdev    : num [1:44] 8.48 3.81 2.72 2.21 1.72 ...
    ##  $ rotation: num [1:44, 1:44] -0.244 -0.162 -0.186 -0.138 -0.17 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:44] "CXCL10" "CXCL9" "DDX58" "DDX60" ...
    ##   .. ..$ : chr [1:44] "PC1" "PC2" "PC3" "PC4" ...
    ##  $ center  : Named num [1:44] 4.25 3.93 5.17 4.63 5.36 ...
    ##   ..- attr(*, "names")= chr [1:44] "CXCL10" "CXCL9" "DDX58" "DDX60" ...
    ##  $ scale   : logi FALSE
    ##  $ x       : num [1:484, 1:44] -6.31 -9.52 5.08 -6.87 -3.02 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : NULL
    ##   .. ..$ : chr [1:44] "PC1" "PC2" "PC3" "PC4" ...
    ##  - attr(*, "class")= chr "prcomp"

``` r
# Visualize the PCA result by gender
ggbiplot(pca_DEG_new_trans2, 
         ellipse = TRUE, 
         var.axes = FALSE, 
         groups = DEG_new_trans2$Gender) +
  ggtitle("PCA of SARS_Cov2 gene expression by gender") +
  theme_bw()
```

![](DataCleanup_DE_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
# Visualize PCA result by infection status
ggbiplot(pca_DEG_new_trans2, 
         ellipse = TRUE, 
         groups = DEG_new_trans2$Sars_test) + 
  ggtitle("PCA of SARS_Cov2 gene expression by infection status") +
  theme_bw()
```

![](DataCleanup_DE_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

``` r
# Determine total variance explained by each principal component
variance_expl <- pca_DEG_new_trans2$sdev^2 / sum(pca_DEG_new_trans2$sdev^2)

# create an elbow plot
qplot(c(1:44), variance_expl) +
  geom_line() +
  xlab("Principal Component") +
  ylab("Variance explained") +
  ggtitle("Elbow Plot of Principal Components") +
  ylim(0, 1) +
  theme_bw()
```

![](DataCleanup_DE_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->
