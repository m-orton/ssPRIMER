# ssPRIMER
ssPRIMER (or species-specific PRIMER) is a shiny based software tool that can be used to design species-specific primer sets for qPCR assays. A multiple sequence alignment can be imported in by a user, and the tool will then guide the user through the process of designing and evaluating species-specific primer sets (and in futue iterations of the tool: Taqman probes). The tool is designed to create primer sets that maximize amplification efficiency for the target species (sensitivity) but minimize amplification efficiency for non-target species (specificity). This tool is designed to benefit the users of eDNA technology, including field biologists, ecologists, conservation researchers, and environmental consultants and could contribute to environmental biomonitoring using molecular methods.   



## Installation

### Currently you may experience issues running the tool locally but you can run the tool online here:

### [ssPRIMER](https://www.mattortonapps.com/shiny/ssPRIMER/)


Please ensure the following packages are installed in R/RStudio if you wish to run locally:

```
# Please ensure R version: 3.6.0 is installed prior to running this tool

# install.packages("plyr")
library(plyr)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("shiny")
library(shiny)
# install.packages("shinyjs")
library(shinyjs)
# install.packages("shinyalert")
library(shinyalert)
# install.packages("shinyWidgets")
library(shinyWidgets)
# install.packages("jsonlite")
library(jsonlite)
# install.packages("rhandsontable")
library(rhandsontable)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("Biostrings")
# BiocManager::install("DECIPHER")
library(Biostrings)
library(DECIPHER)
# install.packages("RSQLite")
library(RSQLite)
# install.packages("foreach")
library(foreach)
# install.packages("TmCalculator")
library(TmCalculator)
# install.packages("stringi")
library(stringi)
# install.packages("d3heatmap")
library(d3heatmap)
```
To run locally, please ensure OligoArrayLux is installed first:

[OligoArrayLux](http://unafold.rna.albany.edu/?q=DINAMelt/OligoArrayAux)

If using the test alignment ensure this file is in your working directory:

[Sample Alignment](sampleAlignment.fas)

Then simply use this command in RStudio (ensuring that you have R version 3.6.0 installed):
```
shiny::runGitHub('ssPRIMER', 'm-orton')
```

## Authors of Shiny App
Matthew Orton

Contributions by Dr. Sally Adamowicz, Kamil Chatila-Amos, Samantha Majoros and Alexandra Albin.

This tool relies greatly on the DECIPHER, Biostrings and TmCalculator R packages for the design of primer sets:

[DECIPHER Package](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html)

[Biostrings Package](http://bioconductor.org/packages/release/bioc/html/Biostrings.html)

[TmCalculator](https://cran.r-project.org/web/packages/TmCalculator/index.html)

It also relies on the OligoArrayLux software:

[OligoArrayLux](http://unafold.rna.albany.edu/?q=DINAMelt/OligoArrayAux)
