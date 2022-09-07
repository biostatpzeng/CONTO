# CONTO: Composite null hypothesis test for trans-ethnic genetic overlap
# Introduction
**CONTO** is a novel gene-centric genetic overlap detection method. The identification of population-common genes across the whole genome can be effectively handled under the high-dimensional framework of composite null hypothesis testing by borrowing the idea of joint significance test (JST). Unlike previous studies which analyzed individual SNPs, it instead focuses on a set of SNPs located within a gene simultaneously and assess their joint significance with the trait of interest.

Like JST, we take Pmax=max(P1, P2) as our test statistic for the detection of for trans-ethnic genetic overlap. However, in contrast to JST which uses the zero-one uniform distribution as its null distribution, we directly build the null distribution of Pmax to correct the conservativeness of JST by borrowing the idea given in, which was proposed under the context of high-dimensional epigenetic mediation analysis. 

It needs to highlight that the proposed method above implicitly assumes that the two P values are uncorrelated with each other. When implementing our method, we first decorrelate test statistics for each gene across populations by multiplying Z scores by the inverse of a correlation matrix. The cross-population correlation coefficient is calculated with Z scores of null genes (e.g., those with P1>0.05 and P2>0.05). The uncorrelated Z scores can be in turn transformed into two-sided P values based on the normal approximation. Theoretically, this decorrelation strategy maximizes the transformed test statistics and the original ones; therefore, it has the minimal influence on identifying shared associations.

CONTO is implemented in the R statistical environment.
# Required input data
CONTO requires two types of input data:

GENE-level summary ststistics of two populations in terms of effect size and their standard error are required as inputs.

Summary statistics, e.g.,
```ruby
               P.x         P.y
A1BG          0.824555    0.042210
A1CF          0.661079    0.403511
A3GALT2       0.712489    0.203850
...

```
The summary statistics file should be at least four columns (i.e. P.x and P.y). P.x is the p-value of GENEs in one population; P.y is the p-value of GENEs in the other population.

# Example
```ruby
library(fdrtool)
source("CONTO.R")
input_pvalues = read.table("scz.txt",sep="\t",header=T)
ax=conto_estimation(input_pvalues,lambda=0.5)

X1         X2         efdr11     efdr10     efdr01
0.8713924  0.5784122  0.8126806  0.6138184  0.7875082

```
# Cite
Jiahao Qiao<sup>$</sup>, Zhonghe Shao<sup>$</sup>, Yuxuan Wu, Ting Wang<sup>#</sup> and Ping Zeng<sup>#</sup> (2022). Detecting associated genes for complex traits shared across East Asian and European populations under the framework of composite null hypothesis testing.

# Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Ping Zeng](https://github.com/biostatpzeng) via zpstat@xzhmu.edu.cn.
