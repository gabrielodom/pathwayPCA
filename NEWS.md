# pathwayPCA 1.1.4
### 2020-02-26

* Exported the `aespca()` function for use in the `rnaEditr::` package (currently in development)
* Reduced R's version requirement to 3.1 for legacy user compatibility
* Cleaned up formatting of examples
* Updated website
* Cleaned up examples in vignette 5
* Reformatted NEWS file to put newest news first (dumb mistake on Gabriel's part)



# pathwayPCA 1.1.1
### 2019-06-06

Our build on Bioconductor 3.9 devel fails for the second vignette. This patch
resolves this issue.



# pathwayPCA 0.99.5
### 2019-04-12

We have been accepted to Bioconductor!
See <https://github.com/Bioconductor/Contributions/issues/1000>



# pathwayPCA 0.99.1
### 2019-02-01

We are submitting to Bioconductor soon, so we are resolving as many of the
`BiocCheck()` ERRORs, WARNINGs, and NOTEs. See
<https://github.com/gabrielodom/pathwayPCA/issues/64>



# pathwayPCA 0.98.0
### 2018-12-13

See the issues within the "Bioconductor Submission" and "Vignette Work"
milestones for detailed descriptions of these changes, their discussion, and
motivation: <https://github.com/gabrielodom/pathwayPCA/issues>

* Standardized naming conventions to UpperCamelCase for consistency
* Changed class name pathwaySet -> pathwayCollection
* Added functions `SubsetPathwayData`, `getPathPCLs`, `LoadOntoPCs`,
`write_gmt`, `getSampleIDs`, `CheckSampleIDs`, `getTrimPathwayCollection`,
`getPathwayCollection`, `CheckPwyColl`, `CheckAssay`, and print and subset
methods for `pathwayCollection` objects
* Added sample ID requirements for input assays and phenotype data. Changed the
phenotype required class from anything to a data frame.
* Made the lars implementation a touch more robust. If the algorithm fails to
converge, then we default back to regular SVD
* The AESPCA function can also return pathway-specific vanilla PCA. Also, we
inspected the parametric p-values compared to the permutation-based p-values.
They congruence is nearly perfect for all of the data sets we tested. Thus, we
changed the default p-value estimation method for AESPCA to be parametric.
* Method functions (`*_pVals()`) now return a list with the p-values data frame,
list of PC vectors (in a data frame), and list of loading vectors (in a data
frame). The `getPathPCLs` function will subset this list to return the PCs,
loadings, and administrivia of the method PCA output for a single pathway
* added the June *homo sapiens* pathway collection from WikiPathways:
`wikipwsHS_Entrez_pathwayCollection`
* Lily wrote a vignette geared to show off the diverse functionality of the
package, so this is the new main vignette. I broke off the plots from chapter 4
into their own chapter. The five vignettes I wrote are now supplemental chapters
* updated vignettes and website


# pathwayPCA 0.0.0.9000
### 2018-03-20


* Added a `NEWS.md` file to track changes to the package.
* Built website.