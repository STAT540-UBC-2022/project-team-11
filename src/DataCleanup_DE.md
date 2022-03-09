Data Cleaning and initial DE analysis
================
Credo Casmil and Aditi Nagaraj Nallan
2022-03-09

# Load packages

``` r
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("GEOquery"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("statmod"))
suppressPackageStartupMessages(library("ggrepel"))
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
