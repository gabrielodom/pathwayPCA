# `pathwayPCA` Development Notes
A Bioconductor package for extracting principal components from expressed pathways

Authors: Gabriel J. Odom, Yuguang Ban, and Xi Chen.
Date: 2017-10-19


<br>

## Introduction
*******************************************************************************
We aim to write a package to collect, organize, and document a suite of existing `R` scripts and files. The purpose of this is to ensure that biologists and bioinformaticians will be easily able to apply our work to their existing data. After a discussion with Prof. Chen, we will not include any support for prediction; this package will address pathway to response attribution only. Our core values for this project are as follows:

  - Rely on as few external packages as possbile. This will require more development work at the beginning, but it will make future development, bug fixes, patches, and enhancements much easier.
- Document *everything*. Once again, this will require more up-front work, but it will yield more informed end-users and make transitioning between development teams seamless.
- Simplify. We want object names to have structure and files to be organized. This is more for future developers, but we expect this will help us as well.

To see the current work on the project, please visit [our package site](https://github.com/gabrielodom/pathwayPCA) on GitHub.



### A Vocabulary Primer
For a tutorial on genomics, proteomics, lipidomics, and metabolomics, collectively known as "-omics" or "MS-omics" (mass spectrometry), from a computer scientist's persepctive, see [this BMC Bioinformatics paper](https://doi.org/10.1186/1471-2105-15-S7-S9).

  - Proteomics: the study of biological processes via the analysis of protein expression or state in cells or tissue. Proteins are built on peptides, which are amino acid chains constructed by mRNA translation.
  - Lipidomics: the systems-level study of lipids and their interactions. Lipids are often clustered into [eight categories](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3995129/) sharing certain chemical or physical properties: fatty acyls, glycerolipids, glycerophospholipids, sphingolipids, saccharolipids, polyketides, sterol lipids, and prenol lipids. According to [Feng and Prestwich (2005)](https://www.crcpress.com/Functional-Lipidomics/Feng-Prestwich/p/book/9781574444674), lipids that occur rarely or in small quantities are often the most effectual lipids in biological processes, meaning they are particularly important in disease diagnostics and in understanding pathology.
  - Metabolomics: the study of small molecular end products (called metabolomes) from cellular regulatory pathways. 
  - Gemomics: the study of the complete set of DNA within a single cell of a living organism. As of 2013, this field of study was far and away the largest and best developed field of study of all the -omics branches.

*******************************************************************************
  


<br>

## Overview
The `pathwayPCA` package exists to extract principal components (PCs) from pre-defined gene or protein sets (called pathways) expressed in the users' data. As shown in [Chen et al (2008)](https://academic.oup.com/bioinformatics/article/24/21/2474/191290), [Chen et al (2010)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3480088/), and [Chen (2011)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3215429/), modelling outcomes based on PCs extracted from pathways yields superior results to modelling outcomes on all genes. This package will enable users to extract pathway PCs for outcome association or prediction via three principal component analysis (PCA) modifications

1. Vanilla PCA ([Pearson, 1901](http://stat.smmu.edu.cn/history/pearson1901.pdf))
2. Adaptive, Elastic-Net, Sparse PCA ([Zou et al, 2006](http://www.it.dtu.dk/projects/manifold/Papers/sparsepc.pdf); [Zou and Zhang, 2009](https://projecteuclid.org/euclid.aos/1245332831); and [Chen, 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3215429/))
3. Supervised PCA ([Bair et al, 2006](http://amstat.tandfonline.com/doi/abs/10.1198/016214505000000628#.Wej7y1tSxhE); [Chen et al, 2008](https://academic.oup.com/bioinformatics/article/24/21/2474/191290); and [Chen et al, 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3480088/))

*******************************************************************************


<br>

## User Supplied Inputs
We expect the end-user to input the following information:

### Gene / Protein Input Matrix
We require an input matrix, $\textbf{X}$. This data frame or matrix is expected to be in [tidy form](https://www.jstatsoft.org/article/view/v059i10): $\textbf{X}$ has $n$ rows, one for each observation, and $p$ columns, one for each gene or protein measured. Moreover, because the point of this package is to extract PCs from pathways, we need to know the names of the genes or proteins. Therefore, $\textbf{x}$ must include the names of the genes or proteins as column names.

For data input into an `R` environment, we recommend the end-user consider the [`readr::`](https://cran.r-project.org/web/packages/readr/readr.pdf), [`readxl::`](https://cran.r-project.org/web/packages/readxl/readxl.pdf) (both part of the [Tidyverse](https://www.tidyverse.org/)) and [`data.table::`](https://cran.r-project.org/web/packages/data.table/data.table.pdf) packages. See this [white paper](https://csgillespie.github.io/efficientR/5-3-importing-data.html#fast-data-reading) for a comparison of speed and robustness of these packages. See this [Stack Overflow discussion](https://stackoverflow.com/questions/21435339/data-table-vs-dplyr-can-one-do-something-well-the-other-cant-or-does-poorly) for comments on the strengths and weaknesses of the `data.table::` and Tidyverse packages. An anecdotal comparison could be that SQL users prefer `data.table::`, while most other users find the Tidyverse syntax easier to learn.

### Response Matrix
We request a response matrix, $\textbf{Y}$. This data frame or matrix should have one column if the analysis expected is regression or classification and two columns (event time and censoring indicator) for survival. If $\textbf{Y}$ is `NULL`, then we assume that the end user requests PCs extracted from the pathways without any analysis or attribution ranking. We further require that the row order of the response matrix matches the row order of the gene / protein matrix.

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

*******************************************************************************


<br>

## Main Function(s) Outputs
Our aim in this package is for the end-user to have little knowledge of this code or `R`'s inner workings and still be able to analyze his/her data.

### A New Class for Returned Objects
With user ease in mind, we plan to have our user-facing functions return objects of a specific class, say `pathwayPCA_object` for instance. [To keep with our core value of simplicity, this class name choice obviously needs some work.] We will then develop methods for common functions which apply to objects of our class via the `UseMethod()` function. For example, we may want to return a model summary from a function output, so we will create a method for the common `summary()` function in `base::` `R`. This ensures that end-users can call functions to which they are already comfortable, such as `plot()`, `print()`, `table()`, and others, while we already prescribe the behavior of these functions when called on objects of our class.

### Summary Graphics
We will create plotting functions that produce annotated horizontal bar charts showing the significant pathways related to the response vector. The vertical axis will be specific pathways ranked by their negative log adjusted $p$-values, and the horizontal axis will be these $p$-values. We will annotate the bars themselves with summary statistics concerning the number of genes in that pathway.

The above graphics can be created irrespective of the type of response. However, we will also have a response-specific suite of graphics. This is a topic of further discussion.

### Summary Tables
Because of the flexibility of the `UseMethod()` function, we can return a single object that is handled differently by the `plot()` and `print()` functions. One main benefit to this is that we can create extractor functions to print summary tables from the pathway analysis. This table would include information such as the pathway name, raw $p$-value, adjusted $p$-value, significance ranking, number of genes in the pathway, number of genes expressed in the data, and potentially many more measurements (including measurements of model-specific importance: $p$-values of independence tests on the Schoenfeld residuals of a Cox Proportional Hazards model when $\textbf{Y}$ is a survival-type response).

### List of PC Matrices
We also plan to return pathway-specific extracted PCs, either as a named list or as a column-concatenated matrix. I think that the class object would contain all this information in it anyway, but we would have an extractor function to access this as a list or matrix.

### Developer Specifics
We can also return model calls, model return specifics, and a report on algorithm convergence. This information would not ever be seen by the common end-user, but it would be available to data scientists or developers interested in debugging or troubleshooting an analysis pipeline.

*******************************************************************************


<br>

## Current Work
This is where I'm picking up next. I plan to write an analysis pipeline based on the requirements listed above. My next steps are: draft outlines for each of the necessary external and internal functions, piece out the existing code into these outlines, modify the outline documentation and existing code to match, and workflow test.

### The Function Outline
This is where we will document the inputs and outputs of the internal and external functions. We will decide which functions will be internal later.

[SPOKE WITH JAMES] An aside: after some testing (see the `Testing_Passing_Function_as_Argument.R` file in the `inst/` folder), we can pass functions as an argument to another function. Therefore, our PC extraction and pathway testing functions can take in functions appropriate to their methods: `pca()` / `svd()`, `aes.pca()`, or `suprPCA()` for PC extraction and `lm()`, `glm()`, `survival::coxph()` and `nnet::multinom()` for pathway testing. Because of this, we can integrate development modularity and add other PCA variants or interesting models later, and the existing code can remain unchanged.

[SPOKE WITH JAMES 2] On applying a model function to each PC matrix in a pathway list: we will generate a model summary within the list apply, then create model-specific extraction functions. For example, assume we have 1000 pathways, and assume we ran a linear model testing the association between $\textbf{Y}$ and each pathway. Then, we can create a linear model summary accessor function to extract and tabulate all the relavent information from the output of `summary(lm(y ~ x))`. Creating statistic-specific extraction functions can wait for the next version of the package.

Functions will be as follows:

  - A function to create an input object: take in $\textbf{X}$, $\textbf{Y}$, and the geneset list. Perform necessary compatability checks. Return an S4 object with class called (for the sake of argument) `pathwayExpression`. This is the function that will issue errors and warnings of the data or geneset has the wrong form, including a warning for a `NULL` response matrix.
  - Split the `pathwayExpression` object into a list of matrices by pathway. Remove all the pathway list matrices with fewer than $k = 3$ columns. Return (silently?) a character vector of all the annihilated gene pathways (this should probably be accessible to developers). Honestly, this function should probably be called internally within the function that creates the `pathwayExpression` objects.
  - Extract PCs from each pathway matrix in the `pathwayExpression` object: take in the list of matrices, the number of PCs to extract (default to $k = 3$), and the extraction method (possibly as a function? How did Ben do that with the `slidR` package? \textbf{See paragraph above.}). Return a named list of pathway PC matrices, each with $k$ columns and $N$ rows. This function will require the methods (PCA, AES-PCA, and Supervised PCA) functions to be written seperately and either called internally or passed as an argument. Passed as an argument allows us to write new PCA methods in the future, without rewriting this function. These functions should all support parallel / cluster computing arguments, as this will probably bottleneck the code.
    + Vanilla PCA extraction: take in a list of matrices and the number of PCs to extract (defaults to $k = 3$). Return a list of matrices with $k$ columns as the first $k$ eigenvectors of the pathway matrix or pathway covariance matrix. NOTE: inspect the computational requirement for PC extraction when $N > path_p$ and $N < path_p$. Some pathways may have more features than observations, which makes calculating and decomposing the covariance matrix inefficient. Also, for $N > path_p$, we still have to multiply the $N \times path_p$ pathway matrix by the $p \times k$ eigenvectors. Taking the SVD directly may be more efficient still, if $N$ isn't too large.
    + AES-PCA extraction: same as above, but we also can take in tuning parameters? Check into this.
    + Supervised PCA: same as above, but check for additional necessary arguments to build in or default
  - Test each pathway in the `pathwayExpression` object for significance with the response, $\textbf{Y}$. This function should also probably parallelized as well (but I don't think running 1500 $N \times 3$ linear models is nearly as expensive as running 1500 SVDs or AES-PCAs).
    + Inputs: the list of pathway-extracted PCs from the "extract PCs" function (above), the response vector or matrix $\textbf{Y}$, and the model function to apply across the list of extracted PCs. Once again, this should support a functional, not character, input, so that we can add attribution models in the future. \textbf{See [SPOKE WITH JAMES] paragraph above.}
    + Output: oh man, this is where things get crazy. Do we return a data frame with one row for each pathway? Do we return a list of model summaries and try to do `purrr::` list extraction after the fact? That would be more difficult, but it would allow us to return a lot more stuff than the hard-coded data frame return option. Maybe we return a list of summaries, but then write seperate extraction functions based on the model choice? That does make more sense I suppose... [TALK TO JAMES ABOUT THIS.] Ok, here's what we're going to do: we will return a list of model summaries as a `genesetSummary` object.
  - Tabulate `genesetSummary` results: this will be a suite of functions (probably called via `UseMethod()`). Each function will take in a `genesetSummary` object and return a data frame of the pathway names, unadjusted model $p$-values, the number of expressed measures in the pathway, the total number of measures in the pathway, and other model-specific statistics / measurements. We will have an internal function for each of the four model types (`lm()`, `glm()`, `survival::coxph()` and `nnet::multinom()`).
  - Adjust $p$-values and sort the results data frame: this will take in a data frame returned by the function above, adjust the model $p$-values, sort the columns by the adjusted $p$-values, and return the sorted data frame.
  - Graphics: take in a sorted table from above, and create the following graphics
    + Horizontal Bar Chart: $x$-axis will be $-\log(p_{adj})$, $y$-axis will be the pathways with the smallest 20-30 adjusted $p$-values with labels as the pathway names, the box colour will be the percentage of genes expressed out of genes included in the pathway; mark each horizontal bar with an $a / b$ fraction, where $a$ is the number of expressed genes / proteins / lipids in the pathway, and $b$ is the total number of genes / proteins / lipids in the pathway.
    + Some other graphs?
    
#### Notes on Documentation

  - If you document classes that extend another class or depend on another file, make sure to list those dependencies in the `@include` field of the file documentation. Multiple files are listed with space delimiter. This field then updated the `Collates:` field of the `DESCRIPTION` file.
  - When you add a dataset to your package via `devtools::use_data()`, the name of the object should match the name of the file in the data/ directory.
  - If you document multiple functions in one file, follow the template in `validClass_Omics.R`. Use the `@rdname` to link the multiple functions together (all with the same function name after the tag). Also use colons to seperate section headers from section bodies -- a line break isn't enough.
  
#### Completed and Documented Functions

  - `createClass_Omics.*`: class declarations
  - `validClass_Omics.*()`: validity functions
  - `create_Omics.*()`: generation functions for each class
  - `expressedOmes()`: given a pathway list, extract the -ome names from the data matrix which are in this pathway list. Then extract the matching columns from the MS design matrix.  For supervised PCA, the minimum number of genes per pathway (given by `trim`) should match the `min.features` argument in the `pathway_tScores()` and `pathway_tControl()` functions. Otherwise, you will get this error:
  ```
  Error in quantile.default(abs(cur.tt), 1 - (min.features/p)) : 
  'probs' outside [0,1]
  ```
  - `extract_aesPCs()`: given the list of pathway gene expression matrices, extract the AES-PCA principal components from each pathway. For supervised PCA, the minimum number of genes per pathway (given by `trim`) should match the `min.features` argument in the `pathway_tScores()` and `pathway_tControl()` functions. Otherwise, you will get this error:
  ```
  Error in quantile.default(abs(cur.tt), 1 - (min.features/p)) : 
  'probs' outside [0,1]
  ```
  - `permTest_OmicsSurv()`: given the list of pathway PCs, fit each pathway PC matrix to a given survival response via the Cox Proportional Hazards model, and record the AIC of the true model. Then, permute the response some thousands of times, fit Cox PH models to each, and record the AICs of these permuted-response models. Compare the true model AIC to the permuted-response model AIC, and record the proportion of models for which the the AIC of the permuted-response model is less than the AIC of the true model. This proportion is the $p$-value for the specific pathway. The `permTest_OmicsSurv()` function returns a named vector of all pathway permuted $p$-values.
  - `permTest_OmicsReg()`: ~~IN PROGRESS~~ DONE.
  - `permTest_OmicsClassif()`: ~~IN PROGRESS~~ DONE.
  - The `aespca()` functions and utilities:
    + `matrixRoot()`: in the `calculate_matrixRoot.R` file. Take the $r^{th}$ root of a matrix, for $r > 0$.
    + `normalize()`: in the `unknown_matrixNorm.R` file. I am still at a loss for what this function does, but it modifies the output of the `lars.lsa()` function.
    + `lars.lsa()`: in the `calculate_LARS.R` file. Least Angle and LASSO Regression (modified from the `lars::` package.)
    + `aespca()`: in the `calculate_AESPCA.R` file. Perform adaptive, elastic-net PCA on a pathway design matrix.
  - The Cox PH Supervised PCA functions (in the `superPC_CoxPH_model.R` file):
    + `coxTrain_fun()`: Takes in data, a vector of response times, and a censor vector. Return a scaled score vector.
    + `.coxscor()`: internal - returns score vectors
    + `.coxvar()`: internal - returns variances
    + `.coxstuff()`: internal
    + `mysvd()`: Takes in "tall" -omics data matrix, and returns its centered SVD and centering vector.
  - The Regression Supervised PCA function: `olsTrain_fun()`
  - The Classification Supervised PCA function: `glmTrain_fun()`
  - Supervised PCA Wrapper Functions:
    + `superpc.train()`: Train the supervised PCA model.
    + `superpc.st()`: Given a fit object returned by `superpc.train()`, extract and test the PCs for significance.
    + `pathway_tScores()`: a wrapper function that allows the user to apply both the training and extract/test functions across a list of pathways. For survival data, no `NA`s are allowed.
    + `pathway_tScores()`: Find the pathway-specific Student's-$t$scores.
    + The `randomControlSample` functions: `sample_Survivalresp()`, `sample_Regresp()`, and `sample_Classifresp()`. These functions either take a parametric bootstrap sample or permuted sample of the response input. For survival data, the response input is a set of event time and censoring vectors, while the response vector is "as is" for regression and classification.
    + `pathway_tControl()`: Find the pathway-specific control for the Student's-$t$ scores. This function calls the `randomControlSample` functions to find a comparative null-distribution Student's-$t$ score for each pathway. For survival data, no `NA`s are allowed.
    + `weibullMix_optimParams()`: Given control the $t$ scores and the number of parameters in each pathway, find the estimated parameter values that minimize the Weibull mixture log-likelihood function.
    + `pathway_pValues()`: Given a set of optimal Weibull mixture parameters and a vector of control $t$ scores, calculated the associated $p$-values for each pathway. This function will also ajdust these $p$-values for multiple comparisons if requested (FDR correction by the Benjamini & Hochberg (1995) step-up FDR-controlling procedure is the default).
    + `adjustRaw_pVals()`: The $p$-value adjustment function. This is a direct copy of `multtest::mt.rawp2adjp()` from Bioconductor. I ported it to our package because I was having some `devel` and `bioc-devel` build issues on Travis. Unfortunately, those problems persisted even after this port, but I think I'll leave it in anyway (one less package dependency).
  - Adjustment functions:
    + `adjust_and_sort()`: this function is a nice wrapper around the `adjustRaw_pVals()` function that allows us to simply take in a named vector as returned by the `permTest_OmicsSurv()`, `permTest_OmicsReg()`, or `permTest_OmicsCateg()` functions (for Supervised PCA) or as returned by the `permTest_OmicsCateg()` function (for AES-PCA).
    + `superPCA_pathway_pvals()`: ~~BUILT BUT NOT YET DOCUMENTED. Currently, this function is still in `Test_supervisedPCA_wrapper.R`.~~ This is completed, and in the file `superPC_wrapper.R`.
    
#### Functions Still to Build
These are functions called in functions I have built so far, but I have yet to build them. I still need to find an example for regression data and logistic data to test these functions to find out what I'm missing.

  - The additions to `glmTrain_fun()`, `superpc.train()`, and `superpc.st()` for n-ary and ordered responses. Currently I have a `stop()` command for these two.
  - Export the second `calculate_pathway_pvalues()` function (the one with the `if()` statements and error checks) from the `Test_supervisedPCA_wrapper.R` file to a proper function file in `R/`. NOTE: This function should call the `adjust_and_sort()` function, rather than have that code defined internally (or maybe have the sorting as an external last step?). Ok, so actually export the third version: `superPCA_pathway_pvals()`. This function calls the `adjust_and_sort()` function. DONE.
  - Graphs: check out the link between `Cytoscape` and `R`. Apparently the `Cytoscape` grapics are nice.


### Piecing in Existing Code
On 22 Feb, I met with James to discuss this progress. I still need to:

1. Create a wrapper function for the AES-PCA workflow (similar to the wrapper for Supervised PCA). DONE: this function is in `aesPC_wrapper.R`.
2. Finish the Supervised PCA vignette. DONE.
3. Find out why `browseVignettes()` can't find the vignette I've written. Because when you install a package from GitHub, it doesn't automatically build the vignette. Building vignettes can take a while. If you want this to happen, use `devtools::install_github("gabrielodom/pathwayPCA", build_vignettes = TRUE)`.
4. Write the AES-PCA vignette. DONE.
5. Curate a small, example data set and gene pathway set for our examples and vignettes. DONE

Also after that meeting, I renamed many files to fit my `superPC`, `aesPC`, and `createClass` file groups.

### Still to Do

1. Check the "Functions Still to Build" section: specifically, we have graphics functions to build and we need to add functionality for more GLM variants.
2. Examples for all the main functions.

#### AES-PCA
We have adapted and documented code in the `aes.pca.R` file to our package. We still need to modify this code heavily to make it more efficient.

*******************************************************************************


<br>
<br>

# Testing
James sent me comments as a Word file. I am adding these comments to the file `James_Tests_20180228.Rmd` in the `testing_files/` directory. I will repeat each of his errors / comments here <span style="color:blue"> in  blue </span>, and my fixes / responses below them.

## Installing Package Error
When James installs the package using `devtools::install_github("gabrielodom/pathwayPCA", build_vignettes = TRUE)`, <span style="color:blue"> he sees the following error: </span>
```
Error: processing vignette 'AES-PCA_Walkthrough.Rmd' failed with diagnostics:
there is no package called 'tidyverse'
```
This error is due to the package `tidyverse` being unavailable on his machine. There are three workarounds, as described in [this GitHub thread](https://github.com/lme4/lme4/issues/233):

1. Do not build the vignettes: `devtools::install_github("gabrielodom/pathwayPCA", build_vignettes = FALSE)`
2. Install the dependencies of the vignettes first, then install the package
```
install.packages(c("tidyverse", "reshape2"))
devtools::install_github("gabrielodom/pathwayPCA", build_vignettes = TRUE)
```
3. Install the dependencies in the function call:
```
devtools::install_github("gabrielodom/pathwayPCA",
                         build_vignettes = TRUE,
                         dependencies = TRUE)
```

While I do not believe with absolute certainty that all three methods will work, I do believe that *at least one* of the three methods will work.

*******************************************************************************

<br>

## Data Format and Structure:

### Input

<span style="color:blue">
I think the user input data should be a more general format, e.g. vector, list, data frame, matrix. </span>

1. <span style="color:blue"> less people are familiar with tidy form: can it be data frame or list for input if it would not give us any data cleaning, converting, or computing efficiency issue later in the analysis? Internally use “tibble” to pass data </span>
2. <span style="color:blue"> can we convert the following input formats internally and provide the message such as “numeric values are coerced to factor...” </span>

I strongly disagree with automatic data wrangling. There is really no way to know that the data we would use to create an `Omics*` object is what the user intended without the user confirming the data creation. At the most crucial juncture of the data analysis, we should force the user to check that their inputs are *exactly* what they intend. Once again, the "garbage in, garbage out" rule applies heavily.

However, I do believe that the creation and inclusion of a *data wrangling* vignette would most likely be beneficial to the average user. Furthermore, we can edit the error messages in the `Omics*` object creation to include a link to the help files, and even some pointers on common problems in the data creation step. Because the data structure of the `Omics*` objects are so crucial to the proper execution of the remainder of the analysis, I should beef up the errors, warnings, and messages at this step.

### Gene Set List
<span style="color:blue"> Do we need “setsize” at all for input? I think we can calculate this one, and return this plus the ratio of number of available assayed features/genes to input gene set. </span>

This is correct. We do not need the user to supply the `setsize` object. We can request a list of the gene pathways and a vector of the pathway names (if more detail is necessary than the list names -- GOXXXXXX vs "Kegg glycolosis etc pathway"). If the `TERMS` element of the pathway list is `NULL`, then we will copy the list names to the `TERMS` entry. Further, we will calculate the number of genes in each pathway at the outset of analysis, and store this number as a (named?) vector element of the pathway list called `setsize`.

TO DO: the above. Also, finish the documentation in the `Omics*_class` and `create_Omics*` function files. They currently have no information on the restrictions of the `pathwaySet` object.
