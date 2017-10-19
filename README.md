# pathwayPCA
A Bioconductor package for extracting principal components from expressed pathways

Authors: Gabriel J. Odom, Yuguang Ban, and Xi Chen.
Date: 2017-10-19



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

2. `TERMS`: a character vector of the names of the list elements in `pathways`. For example, the name of the pathway considered in the previous example is `KEGG_ASCORBATE_AND_ALDARATE_METABOLISM`.

3. `setsize`: a named integer vector of the number of genes or proteins in each pathway. For example, the entry for the `KEGG_ASCORBATE_AND_ALDARATE_METABOLISM` pathway is 9.

While we encourage the end-users to supply a pathway, this argument will default to either the "GO" or "MSigDB" pathways. This default gene set is a topic of further discussion. In future versions of this package, we will attempt to have these pathways automatically update from a web source as the package is loaded or attached via the `.onLoad()` or `.onAttach()` namespace hooks, respectively. 

### Additional Developer Arguments
We allow additional arguments to be passed to internal function calls via `...` ("dots") functionality. These additional arguments are available to developers and data scientists who have the requisite knowledge of this package's internal function design to make modifications to internal function defaults. We do not expect the common end-user to have any need for this.

One exception would be the choice of PCA variant. While I strongly believe this should be hidden or defaulted (at minimum, it should appear after the `...` in the arguments list, preculding any partial match possibilities), it should *at least* be addressed in the user-facing help documentation.



## Main Function(s) Outputs
Our aim in this package is for the end-user to have little knowledge of this code or `R`'s inner workings and still be able to analyze his/her data.

### A New Class for Returned Objects
With user ease in mind, we plan to have our user-facing functions return objects of a specific class, say `pathwayPCA_object` for instance. [To keep with our core value of simplicity, this class name choice obviously needs some work.] We will then develop methods for common functions which apply to objects of our class via the `UseMethod()` function. For example, we may want to return a model summary from a function output, so we will create a method for the common `summary()` function in `base::` `R`. This ensures that end-users can call functions to which they are already comfortable, such as `plot()`, `print()`, and others, while we already prescribe the behavior of these functions when called on objects of our class.

### Summary Graphics
We will create plotting functions that produce annotated horizontal bar charts showing the significant pathways related to the response vector. The vertical axis will be specific pathways ranked by their negative log adjusted $p$-values, and the horizontal axis will be these $p$-values. We will annotate the bars themselves with summary statistics concerning the number of genes in that pathway.

The above graphics can be created irrespective of the type of response. However, we will also have a response-specific suite of graphics. This is a topic of further discussion

### Summary Tables
Because of the flexibility of the `UseMethod()` function, we can return a single object that is handled differently by the `plot()` and `print()` functions. One main benefit to this is that we can create extractor functions to print summary tables from the pathway analysis. This table would include information such as the pathway name, raw $p$-value, adjusted $p$-value, significane ranking, number of genes in the pathway, number of genes expressed in the data, and potentially many more measurements.

### List of PC Matrices
We also plan to return pathway-specific extracted PCs, either as a named list or as a column-concatenated matrix. I think that the class object would contain all this information in it anyway, but we would have an extractor function to access this list or matrix.

### Developer Specifics
We can also return model calls, model return specifics, and a report on algorithm convergence. This information would not ever be seen by the common end-user, but it would be available to data scientists or developers interested in debugging or troubleshooting an analysis pipeline.



## The Internal Function Outline
This is where I'm picking up next. I plan to write an analysis pipeline based on the requirements listed above. My next steps are: draft outlines for each of the necessary external and internal functions, piece out the existing code into these outlines, modify the outline documentation and existing code to match, and workflow test.
