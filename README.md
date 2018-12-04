---
title: "README"
author: "Gabriel J. Odom, Yuguang James Ban, Lily Wang, & Xi Steven Chen"
date: "March 20, 2018"
output: html_document
---

### `pathwayPCA`: A Bioconductor package for extracting principal components from expressed pathways

Initial Date: 2017-10-19


<br>

## Introduction
*******************************************************************************

This package enables users to test for significant relationships between phenotype data and principal components extracted from pathway-specific expression data subsets. The principal components can be extracted with adaptive, elastic-net, sparse (AES) and Supervised principal component analysis (PCA) methods. Please see our User Guides for instruction on using the functions in this package. We plan to release this package on Bioconductor with the Fall 2018 release.

<br>

## Help Tickets
*******************************************************************************
To see the current work on the project, please visit [our package development site](https://github.com/gabrielodom/pathwayPCA) on GitHub, or [our package website](https://gabrielodom.github.io/pathwayPCA/).

If you find bugs in our code, or you feel that some functionality is poorly explained, please submit an issue ticket here: https://github.com/gabrielodom/pathwayPCA/issues. Helpful issue tickets give a [minimum working](https://www.jaredknowles.com/journal/2013/5/27/writing-a-minimal-working-example-mwe-in-r) and [reproducible](http://adv-r.had.co.nz/Reproducibility.html) example whenever possible. Please peruse the included links for advice on writing good help ticket requests.

<br>

## Development Principles
*******************************************************************************

We aim to write a package to collect, organize, and document a suite of existing `R` scripts and files. The purpose of this is to ensure that biologists and bioinformaticians will be easily able to apply our work to their existing data. After a discussion with Prof. Chen, we will not include any support for prediction; this package will address pathway to response attribution only. Our core values for this project are as follows:

  - Rely on as few external packages as possbile. This will require more development work at the beginning, but it will make future development, bug fixes, patches, and enhancements much easier.
  - Document *everything*. Once again, this will require more up-front work, but it will yield more informed end-users and make transitioning between development teams seamless.
  - Simplify. We want object names to have structure and files to be organized. This is more for future developers, but we expect this will help us as well.
