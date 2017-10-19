# pathwayPCA
A Bioconductor package for extracting principal components from expressed pathways



## Introduction
We aim to write a package to collect, organize, and document a suite of existing `R` scripts and files. The purpose of this is to ensure that biologists and bioinformaticians will be easily able to apply our work to their existing data. Our core values for this project are as follows:

  - Rely on as few external packages as possbile. This will require more development work at the beginning, but it will make future development, bug fixes, patches, and enhancements much easier.
- Document *everything*. Once again, this will require more up-front work, but it will yield more informed end-users and make transitioning between development teams seamless.
- Simplify. We want object names to have structure and files to be organized. This is more for future developers, but we expect this will help us as well.

To see the current work on the project, please visit [our package site](https://github.com/gabrielodom/pathwayPCA) on GitHub.



## Overview
The `pathwayPCA` package exists to extract principal components (PCs) from pre-defined gene or protein sets (called pathways) expressed in the users' data. As shown in [Chen et al (2008)](https://academic.oup.com/bioinformatics/article/24/21/2474/191290), [Chen et al (2010)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3480088/), and [Chen (2011)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3215429/), modelling outcomes based on PCs extracted from pathways yields superior results to modelling outcomes on all genes. This package will enable users to extract pathway PCs for outcome association or prediction via three principal component analysis (PCA) modifications

1. Vanilla PCA ([Pearson, 1901](http://stat.smmu.edu.cn/history/pearson1901.pdf))
2. Adaptive, Elastic-Net, Sparse PCA ([Zou et al, 2006](http://www.it.dtu.dk/projects/manifold/Papers/sparsepc.pdf); [Zou and Zhang, 2009](https://projecteuclid.org/euclid.aos/1245332831); and [Chen, 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3215429/))
3. Supervised PCA ([Bair et al, 2006](http://amstat.tandfonline.com/doi/abs/10.1198/016214505000000628#.Wej7y1tSxhE); [Chen et al, 2008](https://academic.oup.com/bioinformatics/article/24/21/2474/191290); and [Chen et al, 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3480088/))



## User Supplied Inputs
We expect the end-user to input the following information:

### Predictor Matrix
We require a predictor matrix, $\textbf{X}$. This data frame or matrix is expected to be in [tidy form](https://www.jstatsoft.org/article/view/v059i10): $\textbf{X}$ has $n$ rows, one for each observation, and $p$ columns, one for each gene or protein measured. Moreover, because the point of this package is to extract PCs from pathways, we need to know the names of the genes or proteins. Therefore, $\textbf{x}$ must include the names of the genes or proteins as column names.

For data input into an `R` environment, we recommend the end-user consider the [`readr::`](https://cran.r-project.org/web/packages/readr/readr.pdf), [`readxl::`](https://cran.r-project.org/web/packages/readxl/readxl.pdf) (both part of the [Tidyverse](https://www.tidyverse.org/)) and [`data.table::`](https://cran.r-project.org/web/packages/data.table/data.table.pdf) packages. See this [white paper](https://csgillespie.github.io/efficientR/5-3-importing-data.html#fast-data-reading) for a comparison of speed and robustness of these packages. See this [Stack Overflow discussion](https://stackoverflow.com/questions/21435339/data-table-vs-dplyr-can-one-do-something-well-the-other-cant-or-does-poorly) for comments on the strengths and weaknesses of the `data.table::` and Tidyverse packages. An anecdotal comparison could be that SQL users prefer `data.table::`, while most other users find the Tidyverse syntax easier to learn.

### Response Matrix
We request a response matrix, $\textbf{Y}$. This data frame or matrix should have one column if the analysis expected is regression or classification and two columns (event time and censoring indicator) for survival. If $\textbf{Y}$ is `NULL`, then we assume that the end user requests PCs extracted from the pathways without any analysis. We further require that the row order of the response matrix matches the row order of the predictor matrix.

### Gene Set List
We require a list of pathway information. This list must be composed of the following three elements

1. `pathways`: a list of the hundreds or thousands of pathways in the gene set. Each element of the list will be a character vector containing the names of the genes or proteins. For example, a list entry could be
```{r}
"ALDH1B1" "ALDH2" "UGDH" "ALDH9A1" "ALDH3A2" "ALDH7A1" "UGT2B10" "UGT2A3" "UGT2B15"
```
which is the name of nine genes in a particular pathway as a character vector.
2. `TERMS`: a character vector of the names of the list elements in `pathways`. For example, the name of the pathway considered in the previous example is `KEGG_ASCORBATE_AND_ALDARATE_METABOLISM`





