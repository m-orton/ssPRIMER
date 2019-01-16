# ssPRIMER
ssPRIMER (or species-specific PRIMER) is a shiny based software tool that can be used to design species-specific primer sets for qPCR assays. A multiple sequence alignment can be imported in by a user, and the tool will then guide the user through the process of designing and evaluating species-specific primer sets (and in futue iterations of the tool: Taqman probes). The tool is designed to create primer sets that maximize amplification efficiency for the target species (sensitivity) but minimize amplification efficiency for non-target species (specificity). This tool is designed to benefit the users of eDNA technology, including field biologists, ecologists, conservation researchers, and environmental consultants and could contribute to environmental biomonitoring using molecular methods.   Eventually this tool will be made available as a free online Shiny app.

## Installation

Please ensure the following packages are installed in R/RStudio:

```
source("https://bioconductor.org/biocLite.R")
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("DECIPHER")
install.packages("RSQLite")
install.packages("shiny")
install.packages("ggplot2")
install.packages("shinyjs")
install.packages("rhandsontable")
install.packages("jsonlite")
install.packages("plyr")
```

To run simply use the command in RStudio:
```
shiny::runGitHub('ssPRIMER', 'm-orton')
```

Alternatively you can run by ensuring the app.R script and following files are in your current working
R directory in RStudio:

[Sample Alignment](sampleAlignment.fas)

[Hybrid-min](hybrid-min.exe)

[Hybrid-ss-min](hybrid-ss-min.exe)


## Authors of Shiny App
Matthew Orton

Contributions by Dr. Sally Adamowicz and Kamil Chatila-Amos

This tool relies greatly on the DECIPHER and Biostrings R packages for design of primer sets:

[DECIPHER Package](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html)

[Biostrings Package](http://bioconductor.org/packages/release/bioc/html/Biostrings.html)
