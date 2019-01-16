# ssPRIMER
ssPRIMER (or species-specific PRIMER) is a shiny based software tool that can be used to design species-specific primer sets for qPCR assays. A multiple sequence alignment can be imported in by a user, and the tool will then guide the user through the process of designing and evaluating species-specific primer sets (and in futue iterations of the tool: Taqman probes). The tool is designed to create primer sets that maximize amplification efficiency for the target species (sensitivity) but minimize amplification efficiency for non-target species (specificity). This tool is designed to benefit the users of eDNA technology, including field biologists, ecologists, conservation researchers, and environmental consultants and could contribute to environmental biomonitoring using molecular methods.                   

## Installation

Please ensure the following packages are installed in R/RStudio:

```
source("https://bioconductor.org/biocLite.R")
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("DECIPHER")
install.packages("ggplot2")
install.packages("shinyjs")
install.packages("rhandsontable")
install.packages("jsonlite")
install.packages("plyr")
```

Also ensure that the sample multiple sequence alignment is present in your working directory:


## Authors of Pipeline
Matthew Orton

Contributions by Dr. Sally Adamowicz and Kamil Chatila-Amos
