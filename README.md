# Optodroplets
Small script to analyze fluorescence microscopy images of light-dependent formation of optodroplets in living cells.<br/>

## What does this script do?
* simple segmentation of the nucleus
* exclusion of larger aggregates that are already present before blue light illumination
* quantification of:
  * concentration of the protein of interest (i. e. mean nuclear fluorescence intensity)
  * optodroplet abundance (i. e. normalized coefficient of variation of nuclear fluorescence intensity)
* fit of the data to compare different proteins of interest
* outputs:
  * text file listing the above-mentioned parameters for each cell
  * pdf with a basic droplet abundance vs. concentration plot
  * fit parameters (only in R command line)
  * optional: an image with the created nuclear and aggregate masks for visual inspection (for each cell)

## How is this script used?

### Requirements
* R version 4.0 (tested with 4.0.1 & 4.0.5)
* R packages:
  * [EBImage](https://bioconductor.org/packages/release/bioc/html/EBImage.html) version 4.29.2
  * [this.path](https://cran.r-project.org/package=this.path) version 0.4.4

### Running the analysis on the test data set
* Install the listed requirements via [Bioconductor](https://bioconductor.org/) (for [EBImage](https://bioconductor.org/packages/release/bioc/html/EBImage.html)) or directly in R (for [this.path](https://cran.r-project.org/package=this.path) via `install.packages("this.path")`)
* Download the [functions](https://github.com/AnneRademacher/Optodroplets/tree/main/functions) and [example_data](https://github.com/AnneRademacher/Optodroplets/tree/main/example_data) folders as well as the R script [optodroplets_induction.R](https://github.com/AnneRademacher/Optodroplets/blob/main/optodroplets_induction.R) into one folder on your local computer (note that [example_data](https://github.com/AnneRademacher/Optodroplets/tree/main/example_data) is about 250 MB and contains the full data set from Fig. 4 of Rademacher *et al.* (2021) listed below)
* Run the R script [optodroplets_induction.R](https://github.com/AnneRademacher/Optodroplets/blob/main/optodroplets_induction.R) either in R Studio or from the R console using `source("optodroplets_induction.R")`

## Associated scientific publications
This is associated with the following publications:
* Rademacher A, Erdel F, Weinmann R and Rippe K (2021) Assessing the phase separation propensity of proteins in living cells through optodroplet formation. Methods in Molecular Biology. *In Preparation*.
* Erdel F, Rademacher A, Vlijm R, Tünnermann J, Frank L, Weinmann R, Schweigert E, Yserentant K, Hummert J, Bauer C, Schumacher S, Al Alwash A, Normand C, Herten DP, Engelhardt J, Rippe K. Mouse Heterochromatin Adopts Digital Compaction States without Showing Hallmarks of HP1-Driven Liquid-Liquid Phase Separation. *Mol Cell*. 2020 Apr 16;78(2):236-249.e7. PMID: 3210700. PMCID: [PMC7163299](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc7163299/). DOI: [10.1016/j.molcel.2020.02.005](https://doi.org/10.1016/j.molcel.2020.02.005).
