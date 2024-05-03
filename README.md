# Identifying Differentially Expressed Genes in Response to Exercise Training in Mice

> **_NOTE:_**  For a Markdown version of the analysis, refer to the [index.md](index.md) file.

## Introduction

In this study, we downloaded a publicly available gene expression dataset [GSE227516](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227516) from Gene Expression Omnibus (GEO) which investigated the effects of exercise on pancreatic islet function in mice. Our aim is to analyze this dataset using DESeq2 to identify genes that show altered expression in response to exercise training.


## Dataset

The gene expression dataset used in this analysis is  [GSE227516](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227516) from Gene Expression Omnibus (GEO).

## Plots

All plots generated during the analysis are available in the [/results/plots](/results/plots) directory.

## Usage

1. **Install R**: If you haven't already, install R from [https://www.r-project.org/](https://www.r-project.org/).
2. **Install Packages**: Open R and install required packages:
    ```R
    install.packages(c("knitr", "rmarkdown"), repos = "https://cloud.r-project.org/")
    library("knitr")
    library("rmarkdown")
    ```
3. **Run Analysis**:
    - Download the provided Rmd file and open it in RStudio or any text editor.
    - Render the Rmd file to generate the output HTML file:
        ```R
        rmarkdown::render("DGE_Analysis.Rmd", output_file="index.html")
        ```
4. **Retrieve R Code**:
    - To retrieve the R code from the Rmd file, run:
        ```R
        knitr::purl("DGE_Analysis.Rmd")
        ```