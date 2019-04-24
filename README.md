### `pathwayPCA`:  an R package for integrative pathway analysis with modern PCA methodology and gene selection

Initial Date: 2017-10-19


<br>

## Introduction
*******************************************************************************
With the advance in high-throughput technology for molecular assays, multi-omics datasets have become increasingly available. However, most currently available pathway analysis software do not provide estimates on sample-specific pathway activities, and provide little or no functionalities for analyzing multiple types of omics data simultaneously. To address these challenges, we present pathwayPCA, a unique integrative pathway analysis software that utilizes modern statistical methodology on principal component analysis (PCA) and gene selection. 

The main features of pathwayPCA include: 

1.	Performing pathway analysis for datasets with binary, continuous, or survival outcomes with computational efficiency. 
2.	Extracting relevant genes from pathways using the SuperPCA and AESPCA approaches.
3.	Computing PCs based on the selected genes. These estimated latent variables represent pathway activity for individual subjects, which can be used to perform integrative pathway analysis, such as multi-omics analysis, or predicting survival time.
4.	Can be used to analyze studies with complex experimental designs that include multiple covariates and/or interaction effects. For example, testing whether pathway associations with clinical phenotype are different between male and female subjects.
5.	Performing analyses with enhanced computational efficiency with parallel computing and enhanced data safety with S4-class data objects.

<br>

## Installing the Package
*******************************************************************************
`pathwayPCA` is a package for R, so you need [R](https://cloud.r-project.org/) first. We also strongly recommend the [RStudio](https://www.rstudio.com/products/rstudio/download/) integrated development environment as a user-friendly graphical wrapper for R.

### Stable Build
The stable build of our package will be available on Bioconductor in May of 2019. To access Bioconductor packages, first install BiocManager, then use BiocManager to install this package:
```
install.packages("BiocManager")
BiocManager::install("pathwayPCA")
```

### Development Build
Because we are currently in the development phase for version 2 of this package, you can install the package from GitHub. In order to install a package from GitHub, you will need the `devtools::` package (https://github.com/r-lib/devtools) and either [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (for Windows) or [Xcode](https://developer.apple.com/xcode/) (for Mac). Then you can install the development version of the [`pathwayPCA` package](https://github.com/gabrielodom/pathwayPCA) from [GitHub](https://github.com/):
```
devtools::install_github("gabrielodom/pathwayPCA")
```

If you are using R version 3.5 or later, and want access to the frozen build for this version, please use
```
devtools::install_github("gabrielodom/pathwayPCA", ref = "stable_3_5")
```

<br>

## Help Tickets
*******************************************************************************
To see the current work on the project, please visit [our package development site](https://github.com/gabrielodom/pathwayPCA) on GitHub, or [our package website](https://gabrielodom.github.io/pathwayPCA/).

If you find bugs in our code, or you feel that some functionality is poorly explained, please submit an issue ticket here: https://github.com/gabrielodom/pathwayPCA/issues. Helpful issue tickets give a [minimum working](https://www.jaredknowles.com/journal/2013/5/27/writing-a-minimal-working-example-mwe-in-r) and [reproducible](http://adv-r.had.co.nz/Reproducibility.html) example whenever possible. Please peruse the included links for advice on writing good help ticket requests.

<br>

## Development Principles
*******************************************************************************

We aim to write a package to collect, organize, and document a suite of existing `R` scripts and files. The purpose of this is to ensure that biologists and bioinformaticians will be easily able to apply our work to their existing data. This package will address pathway to response attribution only. Our core values for this project are as follows:

  - Rely on as few external packages as possbile. This will require more development work at the beginning, but it will make future development, bug fixes, patches, and enhancements much easier.
  - Document *everything*. Once again, this will require more up-front work, but it will yield more informed end-users and make transitioning between development teams seamless.
  - Simplify. We want object names to have structure and files to be organized. This is more for future developers, but we expect this will help us as well.
