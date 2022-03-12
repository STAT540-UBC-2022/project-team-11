
# Progress Report

## Deviation from the final proposal.

As described in our final proposal, we are still using the supplemental
file `GSE152075_raw_counts_GEO.txt` and `PhenoData` from the
`GSE152075`as a source of the raw counts and metadata, respectively.
Analysis of the aforementioned data sets for differential gene
expression, we initially proposed to leverage DESeq2. However, we have
opted for EdgeR because the authors of the paper used DESeq2 and based
on the premise that no statistical modelling can fully capture all
biological phenomena, we wanted to see if a different pipeline of
analysis would be comparable to what was originally found. All group
members contributed to the current progress of our project as described
in the final proposal.

## Progress of Analyses

The following steps describe the methodology used so far in out
analyses: 1. Preliminary data exploration was done using `dplyr`
functions. Using `ggplot2`and `pheatmap` packages, boxplots and heatmaps
were generated. 2. Using functions from the `EdgeR` package, data a
`DGEList` was created amd the count data was normalized using the “TMM”
method. 3. The `decideTestsDGE()` function was used to calculate
differentially expressed genes. 4. From the list of up-regulated and
down-regulated genes obtained from step 4, principal component analyses
(PCA) was performed using functions from the `stat` package. 5. A
logistic regression model was generated to predict infection status
based on RNAseq data.

Based on the above steps, aims 1 and 3 of our final proposal were
partially answered. For aim 1, we have characterized the differentially
expressed genes based on the infection status, i.e. positive vs
negative. As for aim 3, the infection status was predicted after
sub-setting our PCA data into training and validation sets. Other than
the switch from DESeq2 to EdgeR for expression analyses, our goals as
stipulated in our final proposal remain the same.

## Results

Differential gene expression was observed between individuals with
negative and positive COVID infection status. 45 genes were up-regulated
while 21 were down-regulated when sequencing batch was accounted for.
When compared to what the authors of the paper found, we can determine
that we have positive results though we have \~20 less genes than they
got. This discrepancy can be attributed to the genes that we dropped
based on our high CPM and sample thresholds. The following is a table
showing the upregulated genes:

    ##  [1] "|Genes    |    logFC|   logCPM|        LR|    PValue|"
    ##  [2] "|:--------|--------:|--------:|---------:|---------:|"
    ##  [3] "|CMPK2    | 1.899225| 5.432461| 18.249487| 0.0000194|"
    ##  [4] "|CXCL10   | 4.187512| 6.435500| 19.410655| 0.0000105|"
    ##  [5] "|CXCL9    | 4.110525| 5.641219| 24.617031| 0.0000007|"
    ##  [6] "|DDX58    | 2.044954| 5.981685| 20.962829| 0.0000047|"
    ##  [7] "|DDX60    | 1.702926| 5.201112| 19.307998| 0.0000111|"
    ##  [8] "|DDX60L   | 2.049598| 6.437360| 16.327786| 0.0000533|"
    ##  [9] "|FYB1     | 1.526096| 5.593276|  8.180510| 0.0042343|"
    ## [10] "|GBP1     | 2.322803| 6.965568| 26.941275| 0.0000002|"

The following is a table showing the downregulated genes:

    ##  [1] "|Genes   |     logFC|   logCPM|        LR|    Pvalue|"
    ##  [2] "|:-------|---------:|--------:|---------:|---------:|"
    ##  [3] "|BPIFB1  | -1.138585| 6.659925|  6.195652| 0.0128065|"
    ##  [4] "|CBX5    | -1.187887| 5.366671| 10.225757| 0.0013849|"
    ##  [5] "|CCNI    | -1.051094| 5.817436| 11.889177| 0.0005646|"
    ##  [6] "|EEF2    | -1.089281| 6.478402|  6.680745| 0.0097460|"
    ##  [7] "|GAPDH   | -1.097878| 7.635841| 10.205172| 0.0014005|"
    ##  [8] "|MUC5B   | -1.427564| 5.343531|  5.055832| 0.0245432|"
    ##  [9] "|PCSK1N  | -3.299301| 7.753679|  6.870801| 0.0087615|"
    ## [10] "|PIGR    | -1.045433| 7.024179| 10.235317| 0.0013778|"
