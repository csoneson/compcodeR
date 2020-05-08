
# Comparison report, differential expression of RNAseq data
Created by the compcodeR package, version 1.25.1

Date: Fri May  8 11:50:21 2020

Data set:

```
mydata
```
Number of samples per condition:

```
5
```
Included replicates (for repeated simulated data sets):

```
1
```
Differential expression methods included in the comparison:

```
voom.3.45.0.limma.TMM
edgeR.3.31.0.exact.TMM.movingave.tagwise
```
Parameter values:

```
                   value
fdr.threshold       0.05
tpr.threshold       0.05
typeI.threshold     0.05
mcc.threshold       0.05
ma.threshold        0.05
fdc.maxvar          1500
overlap.threshold   0.05
fracsign.threshold  0.05
nbrtpfp.threshold   0.05
signal.measure      mean
```
---

<a name='contents'></a>
## Contents
- [Area under the ROC curve](#auc)

- [False discovery rate](#fdr)

- [True positive rate](#tpr)

- [Spearman correlation between scores, single replicate](#correlation)

---

---
<a name='auc'></a>
## AUC [(Contents)](#contents)
A receiver operating characteristic (ROC) curve is a way to summarize the ability of a test or ranking procedure to rank truly positive (i.e., truly differentially expressed) genes ahead of truly negative (i.e., truly non-differentially expressed). To create the ROC curve for a given differential expression method, the genes are ranked in decreasing order by the score, which is assigned to them during the differential expression analysis and quantifies the degree of statistical significance or association with the predictor (the condition). For a given threshold, all genes with scores above the threshold are classified as 'positive' and all genes with scores below the threshold are classified as 'negative'. Comparing these assignments to the true differential expression status, a true positive rate and a false positive rate can be computed and marked in a plot. As the threshold is changed, these pairs of values trace out the ROC curve. A good test procedure gives a ROC curve which passes close to the upper left corner of the plot, while a poor test corresponds to a ROC curve close to the diagonal. The area under the ROC curve (AUC) summarizes the performance of the ranking procedure. A good method gives an AUC close to 1, while a poor method gives an AUC closer to 0.5. Each boxplot below summarizes the AUCs across all data set replicates included in the comparison. 

<img src="/private/var/folders/24/8k48jl6d249_n_qfxwsl6xvm0000gn/T/RtmpURItjk/compcodeR_figure/auc-1.png" title="plot of chunk auc" alt="plot of chunk auc" style="display: block; margin: auto auto auto 0;" />
---
<a name='fdr'></a>
## FDR [(Contents)](#contents)
The false discovery rate (FDR) indicates the fraction of truly non-differentially expressed genes that we expect to find among the genes that we consider to be differentially expressed. For high-dimensional problems, where many statistical tests are performed simultaneously (such as gene expression studies) it is more relevant to attempt to control the FDR than to control the gene-wise type I error rate, since it is almost certain that at least one gene will show a low nominal p-value even if the null hypothesis is true. To control the FDR, typically, the nominal p-values are adjusted for the large number of tests that are performed. The figures below indicate the observed rate of false discoveries (i.e., the fraction of truly non-differentially expressed genes among the genes that are considered significant) at an adjusted p-value threshold of 0.05. Only methods returning corrected p-values or FDR estimates are included. Each boxplot summarizes the values obtained across all data set replicates that are included in the comparison. For a good method, the observed FDR should not be too high above the imposed adjusted p-value threshold (indicated by a dashed vertical line). If the observed FDR is much larger than the imposed adjusted p-value threshold, the fraction of false discoveries is not controlled at the claimed level. 

<img src="/private/var/folders/24/8k48jl6d249_n_qfxwsl6xvm0000gn/T/RtmpURItjk/compcodeR_figure/fdr-1.png" title="plot of chunk fdr" alt="plot of chunk fdr" style="display: block; margin: auto auto auto 0;" />
---
<a name='tpr'></a>
## TPR [(Contents)](#contents)
The true positive rate (TPR) indicates the fraction of truly non-differentially expressed genes that are indeed considered significant by a method at a given significance threshold. A good method gives a high true positive rate, while at the same time keeping the false discovery rate under control. The figures below show the observed rate of true positives at an adjusted p-value threshold of 0.05. Only methods returning corrected p-values or FDR estimates are included. Each boxplot summarizes the values obtained across all data set replicates included in the comparison.

<img src="/private/var/folders/24/8k48jl6d249_n_qfxwsl6xvm0000gn/T/RtmpURItjk/compcodeR_figure/tpr-1.png" title="plot of chunk tpr" alt="plot of chunk tpr" style="display: block; margin: auto auto auto 0;" />
---
<a name='correlation'></a>
## Spearman correlation between scores [(Contents)](#contents)
The table below shows, for each pair of compared differential expression methods, the Spearman correlation between the scores that they assign to the genes. The value of the correlation is always between -1 and 1, and a high positive value of the Spearman correlation indicates that the compared methods rank the genes in a similar fashion. The results are also shown in a 'heatmap', where the color indicates the Spearman correlation. Finally, the methods are clustered using hierarchical clustering, with a dissimilarity measure defined as 1 - Spearman correlation. This visualizes the relationships among the compared differential expression methods, and groups together methods that rank the genes similarly.

#### 5  samples/condition [(Contents)](#contents)

```
                                         voom.3.45.0.limma.TMM
voom.3.45.0.limma.TMM                                1.0000000
edgeR.3.31.0.exact.TMM.movingave.tagwise             0.8163897
                                         edgeR.3.31.0.exact.TMM.movingave.tagwise
voom.3.45.0.limma.TMM                                                   0.8163897
edgeR.3.31.0.exact.TMM.movingave.tagwise                                1.0000000
```

<img src="/private/var/folders/24/8k48jl6d249_n_qfxwsl6xvm0000gn/T/RtmpURItjk/compcodeR_figure/correlation-1-1.png" title="plot of chunk correlation-1" alt="plot of chunk correlation-1" style="display: block; margin: auto auto auto 0;" />
---
