# Please ensure R version: 3.6.0 is installed prior to running this tool

# Change max file upload to 20 Mb
options(shiny.maxRequestSize = 100*1024^2)
options(warn=-1)

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

# Read in of an example alignment for testing
exampleAlignment <- readDNAStringSet("sampleAlignment.fas")

ui <- tagList(useShinyalert(),
              useShinyjs(),
              tags$head(
                tags$style(
                  HTML(".shiny-notification {
                       height: 100px;
                       width: 600px;
                       position:fixed;
                       top: calc(50% - 50px);;
                       left: calc(50% - 400px);;
                       }
                       "
                  )
                )
              ),
              tags$head(
                HTML("<script>
                    var socket_timeout_interval
                    var n = 0
                    $(document).on('shiny:connected', function(event) {
                    socket_timeout_interval = setInterval(function(){
                    Shiny.onInputChange('count', n++)
                    }, 15000)
                    });
                    $(document).on('shiny:disconnected', function(event) {
                    clearInterval(socket_timeout_interval)
                    });
                    </script>"
                )
                
              ),
              navbarPage("ssPRIMER - A Web-Based Tool for Species-Specific Primer Design", id = "mainPage",
                         tabPanel("1. Upload Alignment", id='panel1', value="panel1",
                                  textOutput("keepAlive"),
                                  sidebarLayout(
                                    sidebarPanel(
                                      p("Welcome to the ssPRIMER homepage! ssPRIMER is a GUI based tool that provides a straightforward 
                                        process to designing species-specific primers and Taqman probes for qPCR Assays.
                                        To start please upload a mulitple sequence alignment in fasta format. 
                                        More options will then be presented to begin the design process.
                                        You can also run a test alignment to test out the functionality of the tool. 
                                        Once uploaded, an interactive visual representation of the alignment will be shown to the 
                                        right. This designs primers and probes primarily with integration from the DECIPHER, Biostrings and Tm Calculator packages."),
                                        
                                      p(strong("***Please note this tool is still in beta testing and we are working hard to make it as 
                                                   polished and functional as possible but it is currently running on a small server and 
                                                   can't handle many users simultaneously. If you are running into issues with the tool 
                                                   please contact us here: morton01@uoguelph.ca***")),
                                      
                                      tags$hr(),
                                      
                                      div(id="alignmentLoad",
                                          
                                          actionButton('exampleAlign', "Load Test Alignment", icon("info-circle"), 
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                          
                                          tags$hr(),
                                          
                                          fileInput('alignmentFile', 'Choose an alignment to upload (.fas or .fa format), 
                                                    max file size = 20 Mb. Ensure all sequences are the same length before 
                                                    uploading, this may require an alignment trimming step before upload. Also ensure that 
                                                    the sequences for the target species are named in exactly the same format.
                                                    Ex: Limnephilus_hyalinus and Limnephilus hyalinus would not be considered the same target group.',
                                                    multiple = F,
                                                    accept = c(
                                                      '.fa',
                                                      '.fas'
                                                    )
                                          )),
                                      
                                      tags$hr(),
                                      
                                      conditionalPanel(
                                        
                                        condition = "output.fileUploaded",
                                        
                                        div(id="target",
                                            p("First select a target species/OTU from your alignment that you would like to be 
                                              selectively amplified."),
                                            
                                            selectInput("inSelectTarget", "Target Species/OTU",
                                                        choices = "Pending Upload"),
                                            
                                            tags$hr(),
                                            
                                            p("You can also restrict to a specific region of the alignment you would like targeted 
                                              by your primer set. It is recommended you do not select a region smaller than 75 bp 
                                              as this will greatly limit the number of potential primer sets and will greatly lower 
                                              the probability that a suitable probe can be found. The alignment visualization will 
                                              regenerate with the region specified automatically."),
                                            
                                            sliderInput("inSlider1", "Restrict Amplification Region",
                                                        min = 1, max = 100, value = c(1,100))
                                        ),
                                        
                                        a(id = "toggleAlignSettings", "View more alignment settings", href = "#"),
                                        shinyjs::hidden(   
                                          div(id="alignSettings",    
                                              sliderInput("inSliderWidth", "Alignment Display Width",
                                                          min = 20, max = 200, value = 60, step = 20, post = "bp")
                                          )
                                        ),
                                        tags$hr(),
                                        
                                        div(id="stepRun",
                                            actionButton('run1', "Start Design Process!", icon("angle-double-right"), 
                                                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                        ),
                                        width = 3
                                      )),
                                    
                                    mainPanel(
                                      
                                      tags$style(type="text/css",
                                                 ".shiny-output-error { visibility: hidden; }",
                                                 ".shiny-output-error:before { visibility: hidden; }"
                                      ),
                                      
                                      conditionalPanel(
                                        
                                        condition = "output.fileUploaded",
                                        
                                        htmlOutput("inc"),
                                        
                                        width = 7)))
                         ),
                         
                         tabPanel("2. Set qPCR Reaction Conditions", id='panel2', value="panel2",
                                  sidebarLayout(
                                    sidebarPanel(
                                      
                                      # Chose Reaction Conditions for qPCR
                                      div(id="qPCRConditions",
                                          p("You can either choose a pre-defined option from the dropdown menu or adjust manually. 
                                            IDT qPCR parameters will load by default for easy ordering of primers from IDT. 
                                            More reaction presets will be added in the future."),
                                          
                                          selectInput("inSelect2", "Choose a Master Mix for qPCR",
                                                      choices = c("Default IDT qPCR Conditions","Chai Master Mix (Hot Start 2x)")),
                                          
                                          p("Set conditions manually:"),
                                          sliderInput("inSlider2", "[Mg]",
                                                      min = 0, max = 10, value = 3, step = 1, post = "mM"),
                                          
                                          sliderInput("inSlider3", "[K]",
                                                      min = 0, max = 250, value = 0, step = 10, post = "mM"),
                                          
                                          sliderInput("inSlider4", "[Na]",
                                                      min = 0, max = 250, value = 50, step = 10, post = "mM"),
                                          
                                          sliderInput("inSlider5", "[dNTPs]",
                                                      min = 0, max = 2, value = 0.8, step = 0.1, post = "mM"),
                                          
                                          sliderInput("inSlider6", "[Primers]",
                                                      min = 25, max = 500, value = 200, step = 25, post = "nM"),
                                          
                                          selectInput("inSelectPolymerase", "DNA Polymerase",
                                                      choices = c("Taq polymerase", "High-fidelity polymerase")),
                                          
                                          tags$hr(),
                                          
                                          actionButton('back1', "Back", icon("angle-double-left"), 
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                          
                                          actionButton('reset1', "Reset to Default", icon("undo"), 
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                          
                                          actionButton('run2', "Proceed to Primer and Probe Constraints", icon("angle-double-right"), 
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                      ), width = 5),
                                    mainPanel(div(),
                                              width=7))),
                         
                         tabPanel("3. Set Primer and Probe Constraints", id='panel3', value="panel3",
                                  sidebarLayout(
                                    sidebarPanel(
                                      
                                      div(id="Primers",   
                                          
                                          tabsetPanel(id="panels",     
                                                      
                                                tabPanel("Primary Constraints", id='panel3a', value="panel3a",
                                                               
                                                         tags$hr(),      
                                                               
                                                         p("Next, set the contraints on your primers and probes before designing them. 
                                                            It is highly recommended you start with the default if you are new to 
                                                            primer design. Default ranges and values are based on the IDT primer 
                                                            design guidelines for Taqman based qPCR as well as suggested default 
                                                            values from the Decipher package. More experienced users can navigate 
                                                            to the secondary contraints tab, the contraints shown below are the 
                                                            most critical to the design process."),
                                                               
                                                               
                                                          # Amplicon Length
                                                          uiOutput("slider9"),
                                                               
                                                          # Minimum Coverage 
                                                          # h5("Minimum fraction of the target species sequences that must be covered with the 
                                                          #     primer+probe set."),
                                                          sliderInput("inSlider12", "Min Target Coverage (minimum fraction of target 
                                                                      sequences from the alignment that must be covered with the 
                                                                      primers and probe)",
                                                                      min = 20, max = 100, value = 90, post = "%"),
                                                               
                                                          sliderInput("inSlider31", "Min Group Coverage (minimum fraction of target sequences
                                                                      that must have sequence data in the region specified by the user, no gaps)",
                                                                      min = 20, max = 100, value = 90, post = "%"),
                                                               
                                                          # Minimum Amplification Efficiency
                                                          #h5("Minimum efficiency of hybridization desired for the primer set."),
                                                          sliderInput("inSlider13", "Min Target Hybridization Efficiency (minimum fraction 
                                                                       of target amplicons that will be amplified with the specified primer 
                                                                       set each PCR cycle, ***Very high hybridization efficiency of 95-100 will lower 
                                                                       specificity to the target species but low efficiency of 50 or below will greatly 
                                                                       lower the sensitivity of the primer set***)",
                                                                       min = 20, max = 100, value = 80, post = "%"),
                                                               
                                                          sliderInput("inSlider14", "Max Non-Target Hybridization Efficiency 
                                                                      (maximum fraction of non-target amplicons that will be amplified 
                                                                       with the specified primer set each PCR cycle. Ideally this should 
                                                                       be set as low as possible however setting it lower will limit the 
                                                                       number of potential primer sets that can be generated)",
                                                                       min = 0, max = 70, value = 10, post = "%"),
                                                               
                                                          # Primer set Tm
                                                          sliderInput("inSlider15", "Optimal Primer Set Annealing Temperature (C) (the tool will try and find a primer set
                                                                      that matches as closely as possible to this temperature)",
                                                                      min = 40, max = 70, value = 60, step = 0.1),
                                                         
                                                          sliderInput("inSlider17", "Primer Set Annealing Temperature Range (C) (generated primer sets
                                                                       must fall within this range of annealing temperatures)",
                                                                       min = 40, max = 70, value = c(55, 65), step = 0.1),
                                                               
                                                          sliderInput("inSlider29", "Max Annealing Temperature Difference (C) between FP and RP
                                                                       (ideally should be within 0-3 degrees difference)",
                                                                       min = 0, max = 10, value = 3, step = 0.1),
                                                               
                                                          sliderInput("inSlider26", "Min Increase in Probe Annealing Temperature (C)
                                                                       (ideally should be 7-10 degrees higher than primer set)",
                                                                       min = 1, max = 20, value = 7, step = 0.1),
                                                         
                                                          tags$hr(),

                                                          actionButton('back2', "Back", icon("angle-double-left"), 
                                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                                               
                                                          actionButton('reset2a', "Reset to Default", icon("undo"), 
                                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                                               
                                                          actionButton('advSettings', "Go to Secondary Constraints", icon("cogs"), 
                                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                                               
                                                          actionButton('run3a', "Design Primer and Probe Sets!", icon("angle-double-right"), 
                                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                                      ),
                                                      
                                                      tabPanel("Secondary Constraints", id='panel3b', value="panel3b",
                                                               
                                                          tags$hr(),
                                                               
                                                          sliderInput("inSlider16", "Primer Length Range",
                                                                      min = 18, max = 30, value = c(18, 23), step = 1, post = "bp"),
                                                               
                                                          sliderInput("inSlider30", "Probe Length Range",
                                                                      min = 20, max = 40, value = c(24, 28), step = 1, post = "bp"),
                                                               
                                                          sliderInput("inSlider27", "Primer GC Percentage Range",
                                                                      min = 30, max = 80, value = c(30, 70), step = 0.1, post = "%"), 
                                                               
                                                          sliderInput("inSlider28", "Probe GC Percentage Range",
                                                                      min = 30, max = 80, value = c(30, 70), step = 0.1, post = "%"), 
                                                               
                                                          sliderTextInput("inSlider23", "Primer-Dimer Hybridization Efficiency 
                                                                          (less likely to form dimers < 1e-07 > more likely to form dimers. 
                                                                           Ideally you want this set as low as possible however setting it 
                                                                           lower will also limit the number of potential primer sets.)",
                                                                           choices = c("1e-10", "1e-09", "1e-08", "1e-07", "1e-06", "1e-05", "1e-04"), 
                                                                           selected = c("1e-07")),
                                                               
                                                          p("Long stretches of runs or repeats are generally unfavorable for your primers due to 
                                                            mispriming. Run and repeat settings are set higher by default to allow for a greater 
                                                            number of potential primer and probe sets but can be set lower if needed."),
                                                               
                                                          sliderInput("inSlider20", "Max Run Length for Primers 
                                                                      (ex: AAA would be a run length of 3)",
                                                                       min = 0, max = 10, value = 4, step = 1, post = "bp"),
                                                               
                                                          sliderInput("inSlider21", "Max Repeat Length for Primers
                                                                      (ex: CGCG would be a repeat length of 2)",
                                                                       min = 0, max = 6, value = 3, step = 1, post = "bp"),
                                                               
                                                          sliderInput("inSlider22", "Max Run Length for Probes 
                                                                      (ex: AAA would be a run length of 3)",
                                                                       min = 0, max = 10, value = 5, step = 1, post = "bp"),
                                                               
                                                          sliderInput("inSlider24", "Max Repeat Length for Probes
                                                                      (ex: CGCG would be a repeat length of 2)",
                                                                      min = 0, max = 6, value = 4, step = 1, post = "bp"),
                                                               
                                                          tags$hr(),
                                                               
                                                          actionButton('back3', "Back", icon("angle-double-left"), 
                                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                                               
                                                          actionButton('reset2b', "Reset to Default", icon("undo"), 
                                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                                               
                                                      ))), width = 5),
                                    mainPanel(div(),
                                              width=7))),
                         
                         tabPanel("4. Evaluate and Download Primer Sets", id='panel5', value="panel5",
                                  sidebarLayout(
                                    sidebarPanel(
                                      
                                      p("Download primer and probe sets
                                        in .csv"), 
                                      
                                      div(id="PrimerTableDownload",
                                          downloadButton('downloadPrimers', 'Download')
                                      ), width = 2),
                                    
                                    mainPanel(
                                      tabsetPanel(type = "tabs",
                                                  
                                                  tabPanel("Potential Primer and Probe Sets", id='panelPrimerTable', value="panelPrimerTable",
                                                           
                                                           tags$hr(), 
                                                           
                                                           p("Here the top 5 primer and probe sets will be shown and ranked according to 
                                                              their target coverage (fraction of target sequences from the alignment that 
                                                              are covered with the primers and probe). "),
                                                           
                                                           tags$hr(), 
                                                           
                                                           div(rHandsontableOutput('primerTable')),
                                                  
                                                           tags$hr(), 
                                                          
                                                           p("Not sure what a parameter means? Go to the help menu:"),
                                                          
                                                           actionButton('help', "Help Menu", icon("angle-double-right"), 
                                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                                           ),
                                                          
                                                  tabPanel("Primer Binding Visualization",  id='panelPrimerBind', value="panelPrimerBind",
                                                  
                                                           tags$hr(),      
                                                          
                                                           p("Select one of the primer sets from the dropdown menu and a primer binding visualization 
                                                             will be generated based on your target species. Identical base positions are marked by dots 
                                                             and mismatched bases will be highlighted to show differences between target and non-target 
                                                             at regions where the primers and probe bind. All target species sequences will be sorted to 
                                                             the top of the alignment."),
                                                          
                                                           tags$hr(), 
                                                          
                                                           selectInput("inSelectPrimer", "Choose a Primer and Probe Set to Visualize",
                                                                       choices = c("Waiting for Primer Sets")),
                                                            
                                                           actionButton('run4a', "Generate Visualization", icon("angle-double-right"), 
                                                                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                                            
                                                           tags$hr(),      
          
                                                           htmlOutput("inc2")),
                                                  
                                                  tabPanel("Primer Set Comparison (Bar Plot)",  id='panelPrimerGraphs', value="panelPrimerGraphs",
                                                           
                                                           tags$hr(), 
                                                           
                                                           p("Here you can visualize in a bar plot the differences between each primer and probe
                                                              set based on many different parameters selectable from the dropdown menu shown below."),
                                                           
                                                           tags$hr(), 
                                                  
                                                           selectInput("inSelectPrimer2", "Choose a Parameter to Compare across Sets",
                                                                       choices = c("Waiting for Primer Sets")),
                                                           
                                                           tags$hr(),
                                                           
                                                           plotOutput('plot')),
                                                  
                                                  tabPanel("Primer Set Comparison (Heat Map)",  id='panelPrimerGraphs', value="panelPrimerGraphs",
                                                           
                                                           tags$hr(), 
                                                           
                                                           p("Here you can visualize in a heat map the differences between each primer and probe
                                                              set based on their individual parameters. Max and min parameters set during the 
                                                              primer design process are also presented for comparison"),
                                                           
                                                           tags$hr(),
                                                           
                                                           actionButton('heatmap', "Generate Heat Map", icon("angle-double-right"), 
                                                                        style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                                           
                                                           tags$hr(),
                                                           
                                                           d3heatmapOutput("heatmap")),
                                                  
                                                  tabPanel("Dimerization Check",  id='panelPrimerDimer', value="panelPrimerDimer",
                                                           
                                                           tags$hr(), 
                                                           
                                                           p("To check for possible dimerization, it is recommended that you head over to this website:"),
                                                           a(href="http://www.primer-dimer.com/", "PrimerDimer", target="_blank"),
                                                           
                                                           tags$hr(), 
                                                           
                                                           p("You can simply copy the text output shown below for the generated primer sets, paste them into this tool using 
                                                              the link provided above and the tool will find the most likely candidates for dimerization amongst the 
                                                              generated primer sets. Select the multiplex analysis option. Look for dimer products less than -9 kcal/mol as these 
                                                              are the most problematic. Cross-dimerization products between primer and probe sets can be ignored if multiple primer 
                                                              sets are inputted at once."),
            
                                                           tags$hr(),
                                                           
                                                           htmlOutput("primertext"),
                                                  
                                                           tags$hr(),  
                                                            
                                                           p("PrimerDimer is a web-based, condition-independent primer dimerization tool and has been shown to 
                                                              outperform other software tools in the accuracy of its dimerization predictions:"),
                                                           p("Lu, J et al. (2017) PrimerSuite: A High-Throughput Web-Based Primer Design Program for Multiplex Bisulfite PCR. 
                                                              Sci. Rep.7, 41328; doi: 10.1038/srep41328."))
                                                           
                                                           
                                    )))),
                         tabPanel("About", id='About', value="panel6",
                                  sidebarLayout(
                                    sidebarPanel(
                                      # Chose Reaction Conditions for qPCR
                                      div(id="about1",
                                          p("ssPRIMER (or species-specific PRIMER) is a web-based software tool that can be used to 
                                             design species-specific primer sets and Taqman probes for qPCR assays. A multiple sequence alignment can 
                                             be imported in by a user, and the tool will then guide the user through the process of 
                                             designing and evaluating species-specific primer sets and Taqman probes. 
                                             The tool is designed to create primer sets that maximize amplification 
                                             efficiency for the target species (sensitivity) but minimize amplification efficiency 
                                             for non-target species (specificity). This tool is designed to benefit the users of eDNA 
                                             technology, including field biologists, ecologists, conservation researchers, and 
                                             environmental consultants and could contribute to environmental biomonitoring using 
                                             molecular methods."),
                                          
                                          tags$hr(),
                                          
                                          p("Author of tool: Matthew G. Orton, morton01@uoguelph.ca"), 
                                          p("With contributions from Dr. Sally J. Adamowicz, Kamil Chatila-Amos, Alexandra Albin and Samantha Majoros."),
                                          p("If you have any feedback to give on the tool, please contact me at my email above."),
                                          p("You can also view the source code here:"),
                                          a(href="https://github.com/m-orton/ssPRIMER", "ssPRIMER source code", target="_blank"), 
                                          
                                          tags$hr(), 
                                          
                                          p("Issues with the tool? You can also report bugs here:"),
                                          a(href="https://github.com/m-orton/ssPRIMER/issues", "Report a bug", target="_blank"), 
                                          
                                          tags$hr(),
                                          
                                          p("ssPRIMER is licensed under the GNU General Public License v3.0:"),
                                          a(href="https://github.com/m-orton/ssPRIMER/blob/master/LICENSE", "GPLv3.0 License for ssPRIMER", target="_blank"),
                                          
                                          tags$hr(),
                          
                                          p("This tool relies greatly on the DECIPHER, Biostrings and TmCalculator R packages and OligoArrayAux 
                                             software for design of primer sets:"),
                                          a(href="https://bioconductor.org/packages/release/bioc/html/DECIPHER.html", "DECIPHER, ", target="_blank"),
                                          a(href="http://bioconductor.org/packages/release/bioc/html/Biostrings.html", "Biostrings, ", target="_blank"),
                                          a(href="https://cran.r-project.org/web/packages/TmCalculator/index.html", "TmCalculator, ", target="_blank"),
                                          a(href="http://unafold.rna.albany.edu/?q=DINAMelt/OligoArrayAux", "OligoArrayAux", target="_blank"),
                                          
                                          tags$hr(),
                                          
                                          p("Citation for DECIPHER:"),
                                          p("Wright, ES (2016) Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in R. The R Journal, 8(1), 352-359."),
                                          
                                          tags$hr(),
                                          
                                          p("Citation for Biostrings:"),
                                          p("Pagès H, Aboyoun P, Gentleman R, DebRoy S (2019) Biostrings: Efficient manipulation of biological strings. 
                                             R package version 2.52.0."),
                                          
                                          tags$hr(),
                                          
                                          p("Citations for OligoArrayAux:"), 
                                          p("Markham, NR & Zuker, M (2005) DINAMelt web server for nucleic acid melting prediction. 
                                             Nucleic Acids Res., 33, W577-W581."),
                                          p("Markham, NR & Zuker, M (2008) UNAFold: software for nucleic acid folding and hybridization. 
                                             In Keith, J. M., editor, Bioinformatics, Volume II. Structure, Function and Applications, 
                                             number 453 in Methods in Molecular Biology, chapter 1, pages 3–31. Humana Press, 
                                             Totowa, NJ. ISBN 978-1-60327-428-9."), 
                                          
                                          tags$hr(),
                                          
                                          p("Citation for PrimerDimer:"), 
                                          p("Lu, J et al. (2017) PrimerSuite: A High-Throughput Web-Based Primer Design Program for Multiplex Bisulfite PCR. 
                                             Sci. Rep.7, 41328; doi: 10.1038/srep41328."),
                                          
                                          tags$hr(),
                                          
                                          p("Links for other packages used:"), 
                                          a(href="https://cran.r-project.org/web/packages/shiny/index.html", "shiny, ", target="_blank"),
                                          a(href="https://github.com/rstudio/d3heatmap", "d3HeatMap, ", target="_blank"),
                                          a(href="https://github.com/jrowen/rhandsontable", "rhandsontable, ", target="_blank"),
                                          a(href="https://github.com/daattali/shinyjs", "shinyjs, ", target="_blank"),
                                          a(href="https://github.com/daattali/shinyalert", "shinyalert, ", target="_blank"),
                                          a(href="https://github.com/dreamRs/shinyWidgets", "shinyWidgets, ", target="_blank"),
                                          a(href="https://cran.r-project.org/web/packages/foreach/index.html", "foreach, ", target="_blank")

                                      ), width = 5),
                                    mainPanel(div(),
                                              width=6))),
                         tabPanel("Help", id='Help', value="panel7",
                                  sidebarLayout(
                                    sidebarPanel(
                                      # Gives definitions of each parameter to help users
                                      div(id="help1",
                                          p(strong("setID:")),
                                          p("A numerical identifier assigned for each generated primer and probe set."), 
                                          tags$hr(),
                                          p(strong("identifier:")),
                                          p("The target species/OTU that was selected on the initial webpage."), 
                                          tags$hr(),
                                          p(strong("targetCoverage:")),
                                          p("Represented as a % value, this is an average of the target coverages of the forward primer, reverse primer and probe. 
                                             Target coverage is the fraction of target sequences (either belonging to the same target species or target OTU) from 
                                             the alignment that the primers and probe will recognize. Generated primer and probe sets by default are sorted according to 
                                             this value."), 
                                          tags$hr(),
                                          p(strong("ampliconLength:")),
                                          p("The length in bp of the amplicon generated from the primer set."), 
                                          tags$hr(),
                                          p(strong("forward_primer, reverse_primer, probe_seq:")),
                                          p("The sequence in 5' to 3' for either the forward primer, reverse primer or probe."), 
                                          tags$hr(),
                                          p(strong("lengthFP, lengthRP, length_probe:")),
                                          p("Length in bp for either the forward primer, reverse primer or probe."), 
                                          tags$hr(),
                                          p(strong("gcFP, gcRP, gcProbe:")),
                                          p("% gc content for either the forward primer, reverse primer or probe."), 
                                          tags$hr(),
                                          p(strong("Hybridization Efficiency (forward, reverse and average):")),
                                          p("Represented as % values, these values represent the predicted fraction of target amplicons that will be amplified every cycle of pcr.
                                             This value is presented individually for both forward and reverse primers and as an average. Non-target hybridization efficiency as
                                             seen in the primary constraints tab represents the same concept but for non-target amplicons. Non-target efficiency represents the level of 
                                             specificity (as higher non-target efficiency would indicate lower specificity) of the assay whearas target efficiency represents the level 
                                             of sensitivity of the assay. On the secondary constaints tab there is also a setting for primer-dimer hybridization efficiency which would 
                                             represent the predicted fraction of primers (either forward or reverse) that will anneal to each other each cycle of pcr 
                                             (which should be minimized as much as possible)."), 
                                          tags$hr(),
                                          p(strong("annealAverage:")),
                                          p("Represented in degrees C and averaged between forward and reverse primers, this value represents the predicted annealing temperature 
                                             at which both the forward and reverse primers should anneal optimally to the target amplicon."), 
                                          tags$hr(),
                                          p(strong("annealDiffprimer:")),
                                          p("Represented in degrees C, the absolute difference in temperature between the forward and reverse primers."), 
                                          tags$hr(),
                                          p(strong("annealProbe:")),
                                          p("Represented in degrees C, the predicted annealing temperature at which the probe should anneal optimally to the target amplicon."), 
                                          tags$hr(),
                                          p(strong("annealDiffProbe:")),
                                          p("Represented in degrees C, the increase in annealing temperature of the probe subtracted from the annealing average of the primers."), 
                                          tags$hr(),
                                          p(strong("startPosFP, endPosFP, etc. :")),
                                          p("Start and end positions in the alignment for the primers and probe."), 
                                          tags$hr(),
                                          p(strong("ampliconSeq:")),
                                          p("The amplicon sequence generated from the first target sequence of the alignment for each generated primer and probe set. 
                                             This can be useful if for instance a gBlock needs to be synthesized for standard curve experiments"),
                                          
                                          tags$hr(),
                                          
                                          actionButton('back4', "Back to Homepage", icon("angle-double-left"), 
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                          
                                          actionButton('back5', "Back to Primer Results", icon("angle-double-left"), 
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                          
                                      ), width = 5),
                                    mainPanel(div(),
                                              width=6)))
                         

                         
              )
)

server <- function(input, output, session){
  
  # Setting conditional tabs so will only appear with click of an action button
  observe({
    toggle(condition = input$run1, selector = "#mainPage li a[data-value=panel2]")
  })
  
  observe({
    toggle(condition = input$run2, selector = "#mainPage li a[data-value=panel3]")
  })
  
  observe({
    toggle(condition = input$run3a, selector = "#mainPage li a[data-value=panel5]")
  })
  
  # Hide alignment settings
  shinyjs::onclick("toggleAlignSettings",
                   shinyjs::toggle(id = "alignSettings", anim = TRUE))    
  
  reactUpload <- reactiveValues(upload=NULL, upload2=NULL, alignment=NULL)
  
  observeEvent(input$exampleAlign, {
    reactUpload$upload <- exampleAlignment
  })
  
  # Reactive variable for alignment upload
  contentsrea <- reactive({

    reactUpload$upload2 <- input$alignmentFile
    if (is.null(reactUpload$upload2) && is.null(reactUpload$upload)){
      return(NULL)
    }
    if (is.null(reactUpload$upload)){
      reactUpload$alignment <- readDNAStringSet(reactUpload$upload2$datapath)
    }
    if (is.null(reactUpload$upload2)){
      reactUpload$alignment <- reactUpload$upload
    }
    splitNames <- as.character(unlist(lapply(strsplit((names(reactUpload$alignment)), "[*]", fixed=TRUE), function(x) return(x[1]))))
  })
  
  # Reactive variable for alignment length
  contentsrea2 <- reactive({
    dnaLength <- as.numeric(nchar(reactUpload$alignment[[1]]))
  })
  
  # qPCR Reaction Pre-sets selectable in the UI                     
  observeEvent(input$inSelect2, {
    if(input$inSelect2 == 'Chai Master Mix (Hot Start 2x)'){
      updateSliderInput(session, "inSlider2", value = 6)
      updateSliderInput(session, "inSlider3", value = 100)
      updateSliderInput(session, "inSlider4", value = 0)
      updateSliderInput(session, "inSlider5", value = 0.6)
      updateSliderInput(session, "inSlider6", value = 200)
    } 
    if(input$inSelect2 == 'Default IDT qPCR Conditions'){
      updateSliderInput(session, "inSlider2", value = 3)
      updateSliderInput(session, "inSlider3", value = 0)
      updateSliderInput(session, "inSlider4", value = 50)
      updateSliderInput(session, "inSlider5", value = 0.8)
      updateSliderInput(session, "inSlider6", value = 200)
    }
  })
  
  reactSlider1 <- reactiveValues(sliderInput1=NULL, sliderInput2=NULL, sliderInput3=NULL, sliderInput4=NULL)
  
  # Observe slider for target range of alignment
  observeEvent(input$inSlider1, {
    reactSlider1$sliderInput1 <- as.numeric(input$inSlider1[2]-input$inSlider1[1])
    reactSlider1$sliderInput2 <- input$inSlider1[1]
    reactSlider1$sliderInput3 <- input$inSlider1[2]
  })
  
  # Variable for display of alignment - width adjustment
  observeEvent(input$inSliderWidth, {
    reactSlider1$sliderInput4 <- input$inSliderWidth
  })
  
  # Generation of the alignment visualization from the uploaded alignment
  contentsrea3 <- reactive({
    if (is.null(reactUpload$alignment)){
      return(NULL)
    } else {
      alignment <- subseq(reactUpload$alignment, reactSlider1$sliderInput2, reactSlider1$sliderInput3)
      BrowseSeqs(alignment, htmlFile = paste(tempdir(), "/myAlignment.html", sep = ""), openURL = FALSE,
                 patterns=c("-", "A", "C", "G", "T"), 
                 colors=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2"), colWidth=reactSlider1$sliderInput4)
    }
  })
  
  # Display html file for alignment visualization
  getPage<-function() {
    return(includeHTML( file (contentsrea3() )))
  }
  output$inc<-renderUI({
    if (is.null(reactUpload$alignment)){
      return(NULL)
    } else {
      getPage()
    }
  })
  
  # Only show a portion of the UI before upload and then reveal after alignment uplaod
  output$fileUploaded <- reactive({
    return(!is.null(contentsrea()))
  })
  
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  
  # Update sliders based on alignment properties
  observe({
    updateSelectInput(session, "inSelectTarget", choices = contentsrea())
    updateSliderInput(session, "inSlider1", value = c(1, contentsrea2()), max = contentsrea2())
  })
  
  reactTarget <- reactiveValues(target=NULL)
  
  # Calculation of amplicon length range slider based on target range amplification region
  output$slider9 <- renderUI({
    if(reactSlider1$sliderInput1 < 75){
      sliderInput("inSlider9", "Amplicon Length Range", min = (reactSlider1$sliderInput1-10), 
                  max = reactSlider1$sliderInput1, value = c((reactSlider1$sliderInput1-5), reactSlider1$sliderInput1), post="bp")
    } else {
      sliderInput("inSlider9", "Amplicon Length Range", min = 75, max = reactSlider1$sliderInput1, value = c(75, 150), post="bp")
    }
  })
  
  observeEvent(input$inSlider9, {
    if(as.numeric(input$inSlider9[2]-input$inSlider9[1]) < 75){
      shinyalert("Oh no!", "It seems your maximum amplicon length is shorter than 75 bp. You may have difficulty finding primer sets with such a short amplicon length.", type = "warning")
    }
  })
  
  # Observe which target is selected
  observeEvent(input$inSelectTarget, {
    reactTarget$target <- as.character(input$inSelectTarget)
  })
  
  # Primer design Variables
  primerVar <- reactiveValues(startPos=NULL, endPos=NULL, minLength=NULL, maxLength=NULL, primer=NULL,
                              probe=NULL, mono=NULL, div=NULL, dNTP=NULL, annealTempProbe=NULL, minProductSize=NULL, 
                              maxProductSIze=NULL, annealTempPrimer=NULL, minEfficiency=NULL, numPrimerSets=NULL, 
                              minCoverage=NULL, numPrimerSets=NULL, primerPosStart=NULL, primerPosEnd=NULL, gcPrimerMin=NULL,
                              gcPrimerMax=NULL, gcProbeMin=NULL, gcProbeMax=NULL, annealDiffPrimerMin=NULL, annealDiffPrimerMax=NULL,
                              probeMinLength=NULL, probeMaxLength=NULL, annealTempProbe=NULL, monovalentNa=NULL, monovalentK=NULL,
                              divalentMg=NULL, primerConc=NULL, dNTP=NULL, dimers=NULL, maxNonTarget=NULL, setID=NULL, startPosFP=NULL,
                              endPosFP=NULL, startPosProbe=NULL, endPosProbe=NULL, startPosRP=NULL, endPosRP=NULL, chosenSetId=NULL, 
                              repeatL=NULL, runL=NULL, repeatLP=NULL, runLP=NULL, polymerase=NULL, annealRangeMin=NULL, annealRangeMax=NULL)
  
  # Observe primer design choices selected by user
  observeEvent(input$inSelectPolymerase, {
    if(input$inSelectPolymerase == 'Taq polymerase'){
        primerVar$polymerase <- TRUE
    }
    if(input$inSelectPolymerase == 'High-fidelity polymerase'){
        primerVar$polymerase <- FALSE
    }
  })
  
  observeEvent(input$inSlider1[1], {
    primerVar$startPos <- as.numeric(input$inSlider1[1])
  })
  
  observeEvent(input$inSlider1[2], {
    primerVar$endPos <- as.numeric(input$inSlider1[2])
  })
  
  observeEvent(input$inSlider2, {
    primerVar$divalentMg <- as.numeric(input$inSlider2)
  })
  
  observeEvent(input$inSlider3, {
    primerVar$monovalentK <- as.numeric(input$inSlider3)
  })
  
  observeEvent(input$inSlider4, {
    primerVar$monovalentNa <- as.numeric(input$inSlider4)
  })
  
  observeEvent(input$inSlider5, {
    primerVar$dNTP <- as.numeric(input$inSlider5)
  })
  
  observeEvent(input$inSlider6, {
    primerVar$primerConc <- as.numeric(input$inSlider6)
  })
  
  observeEvent(input$inSlider9[1], {
    primerVar$minProductSize <- as.numeric(input$inSlider9[1])
  })
  
  observeEvent(input$inSlider9[2], {
    primerVar$maxProductSize <- as.numeric(input$inSlider9[2])
  })
  
  observeEvent(input$inSlider12, {
    primerVar$minCoverage <- as.numeric(input$inSlider12) / 100
  })
  
  observeEvent(input$inSlider13, {
    primerVar$minEfficiency <- as.numeric(input$inSlider13) / 100
  })
  
  observeEvent(input$inSlider14, {
    primerVar$maxNonTarget <- as.numeric(input$inSlider14)
  })
  
  observeEvent(input$inSlider15, {
    if(as.numeric(input$inSlider15) < as.numeric(input$inSlider17[1]) | as.numeric(input$inSlider15) > as.numeric(input$inSlider17[2])){
      shinyalert("Oh no!", "Please ensure your optimal annealing temperature falls within the annealing temperature range you have set.", type = "warning")
    }
    primerVar$annealTempPrimer <- as.numeric(input$inSlider15)
  })
  
  observeEvent(input$inSlider16[1], {
    primerVar$minLength <- as.numeric(input$inSlider16[1])
  })
  
  observeEvent(input$inSlider16[2], {
    primerVar$maxLength <- as.numeric(input$inSlider16[2])
  })
  
  observeEvent(input$inSlider17[1], {
    primerVar$annealRangeMin <- as.numeric(input$inSlider17[1])
  })
  
  observeEvent(input$inSlider17[2], {
    primerVar$annealRangeMax <- as.numeric(input$inSlider17[2])
  })
  
  observeEvent(input$inSlider20, {
    primerVar$runL <- as.numeric(input$inSlider20) + 1
  })
  
  observeEvent(input$inSlider21, {
    primerVar$repeatL <- as.numeric(input$inSlider21) + 1
  })
  
  observeEvent(input$inSlider22, {
    primerVar$runLP <- as.numeric(input$inSlider22) + 1
  })
  
  observeEvent(input$inSlider23, {
    primerVar$dimers <- as.numeric(input$inSlider23)
  })
  
  observeEvent(input$inSlider24, {
    primerVar$repeatLP <- as.numeric(input$inSlider24) + 1
  })
  
  observeEvent(input$inSlider26, {
    primerVar$annealTempProbe <- as.numeric(input$inSlider26)
  })

  observeEvent(input$inSlider27[1], {
    primerVar$gcPrimerMin <- as.numeric(input$inSlider27[1]) / 100
  })
  
  observeEvent(input$inSlider27[2], {
    primerVar$gcPrimerMax <- as.numeric(input$inSlider27[2]) / 100
  })

  observeEvent(input$inSlider28[1], {
    primerVar$gcProbeMin <- as.numeric(input$inSlider28[1]) / 100
  })
  
  observeEvent(input$inSlider28[2], {
    primerVar$gcProbeMax <- as.numeric(input$inSlider28[2]) / 100
  })

  observeEvent(input$inSlider29, {
    primerVar$annealDiffPrimerMax <- as.numeric(input$inSlider29)
  })
  
  observeEvent(input$inSlider30[1], {
    primerVar$probeMinLength <- as.numeric(input$inSlider30[1]) - 1
  })
  
  observeEvent(input$inSlider30[2], {
    primerVar$probeMaxLength <- as.numeric(input$inSlider30[2]) + 1
  })
  
  observeEvent(input$inSlider31, {
    primerVar$minGCoverage <- as.numeric(input$inSlider31) / 100
  })
  
  # Vars for table, visualizations and graphing outputs
  dfPrimers <- reactiveValues(primerTab=NULL, primerBindTab=NULL, primerBindTab2=NULL, primerGraphTab=NULL)
  
  # Tab navigation
  observeEvent(input$run1, {
    updateNavbarPage(session, "mainPage", selected = 'panel2')
  })
  
  observeEvent(input$run2, {
    updateNavbarPage(session, "mainPage", selected = 'panel3')
  })
  
  observeEvent(input$run3a, {
    updateNavbarPage(session, "mainPage", selected = 'panel5')
  })
  
  observeEvent(input$advSettings, {
    updateTabsetPanel(session, "panels", selected = 'panel3b')
  })
  
  observeEvent(input$back1, {
    updateTabsetPanel(session, "mainPage", selected = 'panel1')
  })
  
  observeEvent(input$back2, {
    updateTabsetPanel(session, "mainPage", selected = 'panel2')
  })
  
  observeEvent(input$back3, {
    updateTabsetPanel(session, "panels", selected = 'panel3a')
  })
  
  observeEvent(input$help, {
    updateTabsetPanel(session, "mainPage", selected = 'panel7')
  })
  
  observeEvent(input$back4, {
    updateTabsetPanel(session, "mainPage", selected = 'panel1')
  })
  
  observeEvent(input$back5, {
    updateTabsetPanel(session, "mainPage", selected = 'panel5')
  })
  
  # Generating primer and probe sets 
  observeEvent(input$run3a, {
    
    withProgress(message = "Starting Design Process...", value = 0, {
      
      dbConn <- dbConnect(SQLite(), ":memory:")
      Seqs2DB(reactUpload$alignment, "XStringSet", dbConn, "userAlignment")
      desc <- dbGetQuery(dbConn, "select description from Seqs")
      desc <- unlist(lapply(strsplit(desc$description, "userAlignment", fixed=TRUE),
                            function(x) return(x[length(x)])))
      Add2DB(data.frame(identifier=desc), dbConn)
      
      incProgress(0.1, detail = "Organizing sequences into tiles...")
      
      tiles <- TileSeqs(dbConn)
      
      incProgress(0.1, detail = "Starting design of primers...")
      
      # Generate Initial Primer and Probe Sets that will then be filtered by GC content further down
      primers <- try(DesignPrimers(tiles, identifier=reactTarget$target, start=primerVar$startPos, end=primerVar$endPos,
                                   minLength=primerVar$minLength, maxLength=primerVar$maxLength, minCoverage=primerVar$minCoverage,
                                   minGroupCoverage=primerVar$minGCoverage, annealingTemp=primerVar$annealTempPrimer, 
                                   P=(primerVar$primerConc/1000000000), monovalent=((primerVar$monovalentNa+primerVar$monovalentK)/1000), 
                                   divalent=(primerVar$divalentMg/1000), dNTPs=(primerVar$dNTP/1000), minEfficiency=primerVar$minEfficiency, 
                                   numPrimerSets=50, minProductSize=primerVar$minProductSize, maxProductSize=primerVar$maxProductSize, 
                                   primerDimer=primerVar$dimers, taqEfficiency=primerVar$polymerase), silent = TRUE)
      
      if(dim(primers)[1] == 0){ 
        # Error for initial design of the primer sets
        shinyalert("Oh no!", "Sorry no primer sets met your specified constraints. Please check the following constraints: amplicon length, 
                    target amplification region, qPCR reaction settings, min and max primer length, primer-dimer hybridization efficiency,
                    min target and min group coverage, min target hybridization efficiency and primer annealing temperature.
                    The most probable causes are either your target hybridization efficiency value is set too high or your min target 
                    coverage value is set too high.", type = "error")
        updateNavbarPage(session, "mainPage", selected = 'panel3')
        
      } else {
        
        incProgress(0.2, detail = "Filtering primer sets based on settings...")
        
        primers$setID <- row.names(primers)
        primers$length_FP <- nchar(primers$forward_primer)
        primers$length_RP <- nchar(primers$reverse_primer)
        
        # Subsetting alignment by target amplification region and using that alignment for calculation of primer positions
        dna <- SearchDB(dbConn)
        alignment <- subseq(dna, reactSlider1$sliderInput2, reactSlider1$sliderInput3)
        names(alignment) <- desc
        u_alignment <- unique(alignment)
        names(u_alignment) <- names(alignment)[match(u_alignment, alignment)]
        alignment <- u_alignment
        # move the target group to the top
        w <- which(names(alignment)==reactTarget$target)
        alignment <- c(alignment[w], alignment[-w])
        targetSequence <- alignment[1]
        
        # Forwar primer positioning
        forwardPrimers <- foreach(i=1:nrow(primers)) %do% DNAString(primers$forward_primer[i])
        matchPatternFP <- foreach(i=1:length(forwardPrimers)) %do% vmatchPattern(forwardPrimers[[i]], targetSequence, max.mismatch=2)
        matchPatternFP_End <- foreach(i=1:length(forwardPrimers)) %do% matchPatternFP[[i]]@ends
        matchPatternFP_End <- as.numeric(unlist(matchPatternFP_End))
        matchPatternFP_Start <- foreach(i=1:length(forwardPrimers)) %do% matchPatternFP[[i]]@width0
        matchPatternFP_Start <- matchPatternFP_End - as.numeric(unlist(matchPatternFP_Start)) + 1
        primers$startPosFP <- matchPatternFP_Start
        primers$endPosFP <- matchPatternFP_End
        
        # Reverse primer positioning
        reversePrimers <- foreach(i=1:nrow(primers)) %do% DNAString(primers$reverse_primer[i])
        reversePrimers <- foreach(i=1:length(reversePrimers)) %do% reverseComplement(reversePrimers[[i]])
        matchPatternRP <- foreach(i=1:length(reversePrimers)) %do% vmatchPattern(reversePrimers[[i]], targetSequence, max.mismatch=2)
        matchPatternRP_End <- foreach(i=1:length(reversePrimers)) %do% matchPatternRP[[i]]@ends
        matchPatternRP_End <- as.numeric(unlist(matchPatternRP_End))
        matchPatternRP_Start <- foreach(i=1:length(reversePrimers)) %do% matchPatternRP[[i]]@width0
        matchPatternRP_Start <- matchPatternRP_End - as.numeric(unlist(matchPatternRP_Start)) + 1
        primers$startPosRP <- matchPatternRP_Start
        primers$endPosRP <- matchPatternRP_End
        
        # Generating an amplicon sequence from positions
        ampliconSeqs <- foreach(i=1:nrow(primers)) %do% subseq(targetSequence, primers$startPosFP[i], primers$endPosRP[i])
        ampliconSeqs <- foreach(i=1:nrow(primers)) %do% as.character(ampliconSeqs[[i]])
        primers$ampliconSeq <- unlist(ampliconSeqs)
        
        # Run detection in primer sets (ex: AAA is a run length of 3 in a primer)
        regexRun <- paste("(A{",primerVar$runL,"}|G{", primerVar$runL, "}|C{", primerVar$runL, "}|T{", primerVar$runL, "})", sep="")
        runCheckFP <- foreach(i=1:nrow(primers)) %do% stri_count_regex(primers$forward_primer[i], regexRun, opts_regex = list())
        runCheckFP_2 <- which(runCheckFP == 0)
        runCheckRP <- foreach(i=1:nrow(primers)) %do% stri_count_regex(primers$reverse_primer[i], regexRun, opts_regex = list())
        runCheckRP_2 <- which(runCheckRP == 0)
        runCheck <- intersect(runCheckFP_2, runCheckRP_2)
        primers <- primers[runCheck,]
        
        if(dim(primers)[1] == 0){ 
          # Error for initial design of the primer sets
          shinyalert("Oh no!", "It seems you have set your max run length setting too low for your primers 
                      and thus no primers could be found. You will have to increase this value.", type = "error")
          updateNavbarPage(session, "panels", selected = 'panel3b')
          
        } else {
        
          # Repeat detection in primer sets (ex: CGCG would be a repeat length of 2 in a primer)
          regexRepeat <- paste("(AA{",primerVar$repeatL, "}|AT{", primerVar$repeatL, "}|AG{", primerVar$repeatL, "}|AC{", primerVar$repeatL,
                              "}|CA{",primerVar$repeatL, "}|CT{", primerVar$repeatL, "}|CG{", primerVar$repeatL, "}|CC{", primerVar$repeatL, 
                              "}|GA{",primerVar$repeatL, "}|GT{", primerVar$repeatL, "}|GG{", primerVar$repeatL, "}|GC{", primerVar$repeatL, 
                              "}|TA{",primerVar$repeatL, "}|TT{", primerVar$repeatL, "}|TG{", primerVar$repeatL, "}|TC{", primerVar$repeatL, "})", sep="")
          repeatCheckFP <- foreach(i=1:nrow(primers)) %do% stri_count_regex(primers$forward_primer[i], regexRepeat, opts_regex = list())
          repeatCheckFP_2 <- which(repeatCheckFP == 0)
          repeatCheckRP <- foreach(i=1:nrow(primers)) %do% stri_count_regex(primers$reverse_primer[i], regexRepeat, opts_regex = list())
          repeatCheckRP_2 <- which(repeatCheckRP == 0)
          repeatCheck <- intersect(repeatCheckFP_2, repeatCheckRP_2)
          primers <- primers[repeatCheck,]
          
          if(dim(primers)[1] == 0){ 
            # Error for initial design of the primer sets
            shinyalert("Oh no!", "It seems you have set your max repeat length setting set too low for your primers 
                        and thus no primers could be found. You will have to increase this value.", type = "error")
            updateNavbarPage(session, "panels", selected = 'panel3b')
            
          } else {
          
            # Amplicon length calculation
            primers$ampliconLength <- primers$endPosRP - primers$startPosFP + 1
            
            # Hybr efficiency rounding
            primers$FP_hybr_efficiency <- round(primers$forward_efficiency, 2)
            primers$RP_hybr_efficiency <- round(primers$reverse_efficiency, 2)
            
            # Calculating GC content
            gcFP <- foreach(i=1:nrow(primers)) %do% letterFrequency(DNAString(primers$forward_primer[i]), "GC", as.prob=TRUE) 
            gcRP <- foreach(i=1:nrow(primers)) %do% letterFrequency(DNAString(primers$reverse_primer[i]), "GC", as.prob=TRUE)
                
            # Filtering by GC Content Criteria set by the user
            gcCheckFP <- which(gcFP >= primerVar$gcPrimerMin & gcFP <= primerVar$gcPrimerMax)
            gcCheckRP <- which(gcRP >= primerVar$gcPrimerMin & gcRP <= primerVar$gcPrimerMax)
            gcCheck <- intersect(gcCheckFP, gcCheckRP)
            
            if(length(gcCheck) == 0){ 
              shinyalert("Oh no!","Sorry no primer sets met your specified constraints. Please check the following constraint: 
                          % gc content for primers. You will have to increase the acceptable % gc range of your primers.", type = "error")
              updateNavbarPage(session, "panels", selected = 'panel3b')
               
            } else {
              
              # Assigning GC content to primers dataframe
              primers <- primers[gcCheck, ]
              gcFP <- foreach(i=1:nrow(primers)) %do% letterFrequency(DNAString(primers$forward_primer[i]), "GC", as.prob=TRUE)
              gcFP <- unlist(gcFP)
              gcRP <- foreach(i=1:nrow(primers)) %do% letterFrequency(DNAString(primers$reverse_primer[i]), "GC", as.prob=TRUE)
              gcRP <- unlist(gcRP)
              primers$gcFP <- round(gcFP, 3)
              primers$gcRP <- round(gcRP, 3)
    
              # Finding non-target hybridization values and filtering according to non-target hybridization criteria set by the user
              nonTFP <- regmatches(primers$mismatches_forward, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", primers$mismatches_forward))
              nonTFP <- foreach(i=1:nrow(primers)) %do% as.numeric(nonTFP[[i]])
              nonTRP <- regmatches(primers$mismatches_reverse, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", primers$mismatches_reverse))
              nonTRP <- foreach(i=1:nrow(primers)) %do% as.numeric(nonTRP[[i]])
              nonTargetCheck <- foreach(i=1:nrow(primers)) %do% append(nonTFP[[i]], nonTRP[[i]])
              nonTargetCheck2 <- foreach(i=1:nrow(primers)) %do% which(nonTargetCheck[[i]] > 50)
              nonTargetCheck <- foreach(i=1:nrow(primers)) %do% nonTargetCheck[[i]][-nonTargetCheck2[[i]]]
              nonTargetCheck <- foreach(i=1:nrow(primers)) %do% max(nonTargetCheck[[i]])
              primers$nonTarget <- as.numeric(nonTargetCheck)
              nonTargetCheck <- which(primers$nonTarget > primerVar$maxNonTarget)
              
              if(length(nonTargetCheck) > 0){
                 primers <- primers[-nonTargetCheck,]
              }
              
              if(dim(primers)[1] == 0){ 
                 shinyalert("Oh no!", "Sorry no primer sets met your specified constraints. Please check the following constraint: 
                             Max non-target hybridization efficiency. You will have to increase this value to allow for more inclusivity 
                             of non-target species/OTUs.", type = "error")
                 updateNavbarPage(session, "mainPage", selected = 'panel3')
                 
              } else {
                
                primers <- (primers[,c("setID", "identifier", "ampliconLength", "startPosFP", "endPosFP", "length_FP","start_reverse", 
                                       "startPosRP", "endPosRP", "length_RP", "forward_primer", "reverse_primer", "FP_hybr_efficiency", 
                                       "RP_hybr_efficiency", "gcFP", "gcRP", "forward_coverage", "reverse_coverage", "ampliconSeq")])
                
                primers <- data.frame(primers, stringsAsFactors = FALSE)
                primers <- as.data.frame(lapply(primers, function(y) gsub(",", "", y)))
                primers <- primers[complete.cases(primers), ]
                primers[, ] <- lapply(primers[, ], as.character)
                primers$setID <- 1:nrow(primers)
                
                # Averaging of target coverage across forward and reverse primers
                # this value will then get averaged with the probe target coverage to give a single metric of coverage
                primers$forward_coverage <- as.numeric(primers$forward_coverage)
                primers$reverse_coverage <- as.numeric(primers$reverse_coverage)
                primers$coverage <- (primers$forward_coverage + primers$reverse_coverage) / 2
                
                # Tm using R package TmCalculator and annealing temperature calculation
                predictedTm_FP <- foreach(i=1:nrow(primers)) %do% Tm_GC(primers$forward_primer[i], ambiguous = FALSE, userset = NULL, variant = "Primer3Plus",
                                                                        Na = primerVar$monovalentNa, K = primerVar$monovalentNa, Tris = 0, 
                                                                        Mg = primerVar$divalentMg, dNTPs = primerVar$dNTP, saltcorr = 0, mismatch = TRUE)
                
                primers$Annealing_FP <- round(as.numeric(unlist(predictedTm_FP)) - 3, 2)
                
                predictedTm_RP <- foreach(i=1:nrow(primers)) %do% Tm_GC(primers$reverse_primer[i], ambiguous = FALSE, userset = NULL, variant = "Primer3Plus",
                                                                        Na = primerVar$monovalentNa, K = primerVar$monovalentNa, Tris = 0, 
                                                                        Mg = primerVar$divalentMg, dNTPs = primerVar$dNTP, saltcorr = 0, mismatch = TRUE)
                
                primers$Annealing_RP <- round(as.numeric(unlist(predictedTm_RP)) - 3, 2)
                
                primers$annealDiffprimer <- round(abs(primers$Annealing_FP - primers$Annealing_RP), 2)
                
                annealDiffCheck <- which(primers$annealDiffprimer >= primerVar$annealDiffPrimerMax)
                primers <- primers[-annealDiffCheck,]
                
                if(dim(primers)[1] == 0){ 
                  shinyalert("Oh no!", "Sorry no primer sets met your specified constraints. Please check the following constraint: 
                              Max anneal temperature difference between FP and RP. You will have to increase this value.", type = "error")
                  updateNavbarPage(session, "mainPage", selected = 'panel3')
                  
                } else {
                  
                  # Calculating an average of the annealing temperature and hybr efficiency of each primer
                  annealAverage <- (as.numeric(primers$Annealing_FP) + as.numeric(primers$Annealing_RP)) / 2
                  primers$annealAverage <- annealAverage
                  
                  # Determining which primers fall within the annealing range
                  annealAverageCheck <- which(primers$annealAverage >= primerVar$annealRangeMin & primers$annealAverage <= primerVar$annealRangeMax)
                  primers <- primers[annealAverageCheck,]
                  
                  if(dim(primers)[1] == 0){ 
                    shinyalert("Oh no!", "Sorry no primer sets met your specified constraints. Please check the following constraint: 
                              Annealing temperature range of primer set. You will have to increase the annealing temperature range of your primers.", type = "error")
                    updateNavbarPage(session, "mainPage", selected = 'panel3')
                    
                    } else {
                    
                    # Calculating an average of hybr efficiency of each primer
                    primerHybrAverage <- (as.numeric(primers$FP_hybr_efficiency) + as.numeric(primers$RP_hybr_efficiency)) / 2
                    primers$hybrEfficiencyAverage <- primerHybrAverage
                    
                    # Calculating probe start and stop positions that must be met 
                    probeStart <- as.numeric(primers$endPosFP) 
                    probeEnd <- as.numeric(primers$startPosRP) 
                    
                    incProgress(0.2, detail = "Starting design of probes...")
                    
                    # Generates tileseqs between a range of probe lengths set by user
                    probeTiles <- foreach(i=1:(primerVar$probeMaxLength-primerVar$probeMinLength)) %do% TileSeqs(dbConn, identifier = reactTarget$target, minLength = primerVar$probeMinLength,  maxLength = (primerVar$probeMinLength + i))
                    dfTarget <- do.call("rbind", probeTiles)
                    
                    # Find probes within target amplification region
                    probeFind1 <- foreach(i=1:length(probeStart)) %do%  which(dfTarget$start_aligned > probeStart[i])
                    probeFind2 <- foreach(i=1:length(probeStart)) %do%  which(dfTarget$end_aligned < probeEnd[i])
                    probeFind <- foreach(i=1:length(probeStart)) %do% intersect(probeFind1[[i]], probeFind2[[i]])
                    listProbeFind <- foreach(i=1:length(probeFind)) %do% dfTarget[probeFind[[i]],]
                    
                    if(length(listProbeFind[[1]]) == 0){ 
                      shinyalert("Oh no!", "Sorry no probes were found from within your target amplification region. You will have to increase either your target 
                                  amplification region or your amplicon length range. Also look at reducing the lengths of your probe if your amplicon
                                  length range is small.", type = "error")
                      updateNavbarPage(session, "mainPage", selected = 'panel1')
                      
                    } else {
                      
                      incProgress(0.2, detail = "Filtering probes based on settings...")
                      
                      names(listProbeFind) <- 1:nrow(primers)
                      probes <- do.call("rbind", listProbeFind)
                      probes$setID <- round(as.numeric(rownames(probes)), 1)
                      probes$setID <- round(probes$setID, 0)
                      probes$length_probe <- nchar(probes$target_site)
                      
                      # Match to target species using group coverage - must be highly similar
                      groupCoverage <- which(probes$groupCoverage >= primerVar$minCoverage)
                      probes <- probes[groupCoverage,]
                      
                      if(dim(probes)[1] == 0){ 
                        shinyalert("Oh no!", "Sorry no probes were found with such a high target coverage. You will have to lower the value for this setting.", type = "error")
                        updateNavbarPage(session, "mainPage", selected = 'panel3')
                        
                      } else {
                        
                        # Run detection in probe sets (ex: AAA is a run length of 3 in a primer)
                        regexRun2 <- paste("(A{",primerVar$runLP,"}|G{", primerVar$runLP, "}|C{", primerVar$runLP, "}|T{", primerVar$runLP, "})", sep="")
                        runCheckProbe <- foreach(i=1:nrow(probes)) %do% stri_count_regex(probes$target_site[i], regexRun2, opts_regex = list())
                        runCheckProbe2 <- which(runCheckProbe == 0)
                        probes <- probes[runCheckProbe2,]
                        
                        if(dim(probes)[1] == 0){ 
                          shinyalert("Oh no!", "Sorry no probes were found with the run length that was set. You will have to increase the run length of your probes.", type = "error")
                          updateNavbarPage(session, "mainPage", selected = 'panel3b')
                          
                        } else {
                        
                          # Repeat detection in probe sets (ex: CGCG would be a repeat length of 2 in a primer)
                          regexRepeat2 <- paste("(AA{",primerVar$repeatLP, "}|AT{", primerVar$repeatLP, "}|AG{", primerVar$repeatLP, "}|AC{", primerVar$repeatLP,
                                               "}|CA{",primerVar$repeatLP, "}|CT{", primerVar$repeatLP, "}|CG{", primerVar$repeatLP, "}|CC{", primerVar$repeatLP, 
                                               "}|GA{",primerVar$repeatLP, "}|GT{", primerVar$repeatLP, "}|GG{", primerVar$repeatLP, "}|GC{", primerVar$repeatLP, 
                                               "}|TA{",primerVar$repeatLP, "}|TT{", primerVar$repeatLP, "}|TG{", primerVar$repeatLP, "}|TC{", primerVar$repeatLP, "})", sep="")
                          repeatCheckProbe <- foreach(i=1:nrow(probes)) %do% stri_count_regex(probes$target_site[i], regexRepeat2, opts_regex = list())
                          repeatCheckProbe2 <- which(repeatCheckProbe == 0)
                          probes <- probes[repeatCheckProbe2,]
                          
                          if(dim(probes)[1] == 0){ 
                            shinyalert("Oh no!", "Sorry no probes were found with the repeat length that was set. You will have to increase the repeat length of your probes.", type = "error")
                            updateNavbarPage(session, "panels", selected = 'panel3b')
                            
                          } else {
                          
                            # GC Check for the probe
                            gcProbe <- foreach(i=1:nrow(probes)) %do% letterFrequency(DNAString(probes$target_site[i]), "GC", as.prob=TRUE) 
                            gcProbe <- unlist(gcProbe)
                            probes$gcProbe <- round(gcProbe, 3) 
                            probes$probe_seq <- probes$target_site  
                            
                            gcCheck <- which(probes$gcProbe >= primerVar$gcProbeMin & probes$gcProbe <= primerVar$gcProbeMax)
                            probes <- probes[gcCheck,]
                            
                            if(dim(probes)[1] == 0){ 
                              shinyalert("Oh no!", "Sorry no probes were found in the % gc range that was set. You will have to increase the range of acceptable % gc content for the probes.", type = "error")
                              updateNavbarPage(session, "panels", selected = 'panel3b')
                              
                            } else {
                              
                              # Annealing Check for the probe
                              tmProbe <- foreach(i=1:nrow(probes)) %do% Tm_GC(probes$probe_seq[i], ambiguous = FALSE, userset = NULL, variant = "Primer3Plus",
                                                                              Na = primerVar$monovalentNa, K = primerVar$monovalentNa, Tris = 0, 
                                                                              Mg = primerVar$divalentMg, dNTPs = primerVar$dNTP, saltcorr = 0, mismatch = TRUE)
  
                              probes$annealProbe <- round(as.numeric(unlist(tmProbe)) - 3, 2)
                              
                              probes <- (probes[,c("setID", "length_probe", "start_aligned", "end_aligned", "groupCoverage", "annealProbe", "gcProbe", "probe_seq")])
                              
                              probeList <- lapply(unique(probes$setID), function(x) 
                                probes[probes$setID == x,])
                              
                              annealCheck <- foreach(i=1:length(probeList)) %do% which(probeList[[i]]$annealProbe < (primers$annealAverage[i] + primerVar$annealTempProbe))
                              
                              if(length(annealCheck[[1]]) == 0){ 
                                shinyalert("Oh no!", "Sorry no probes were found with such a high annealing temperature. You will have to lower the min annealing temperature of the probe.", type = "error")
                                updateNavbarPage(session, "mainPage", selected = 'panel3')
                                
                              } else {
                                
                                # Subsetting to select probes by ordering by target coverage
                                probeList <- foreach(i=1:length(probeList)) %do% probeList[[i]][-annealCheck[[i]],]
                                probeList <- foreach(i=1:length(probeList)) %do% probeList[[i]][order(probeList[[i]]$groupCoverage, decreasing = TRUE),]
                                probeList <- foreach(i=1:length(probeList)) %do% probeList[[i]][1,]
                                
                                probes <- do.call("rbind", probeList)
                                probes <- probes[1:nrow(primers),]
                                
                                # Merging primer sets to probes and modifying/organizing columns
                                primerprobes <- merge(primers, probes, by.x="setID", by.y="setID")
                                primerprobes$setID <- 1:nrow(primerprobes)
                                primerprobes$annealDiffProbe <- round(primerprobes$annealProbe - primerprobes$annealAverage, 2)
                                primerprobes$startPosFP <- as.numeric(primerprobes$startPosFP)
                                primerprobes$endPosFP <- as.numeric(primerprobes$endPosFP)
                                primerprobes$startPosRP <- as.numeric(primerprobes$startPosRP)
                                primerprobes$endPosRP <- as.numeric(primerprobes$endPosRP)
  
                                # Probe positioning
                                probeDNA <- foreach(i=1:nrow(primerprobes)) %do% DNAString(primerprobes$probe_seq[i])
                                matchPatternProbe <- foreach(i=1:length(probeDNA)) %do% vmatchPattern(probeDNA[[i]], targetSequence, max.mismatch=2)
                                matchPatternProbe_End <- foreach(i=1:length(probeDNA)) %do% matchPatternProbe[[i]]@ends
                                matchPatternProbe_End <- as.numeric(unlist(matchPatternProbe_End))
                                matchPatternProbe_Start <- foreach(i=1:length(probeDNA)) %do% matchPatternProbe[[i]]@width0
                                matchPatternProbe_Start <- matchPatternProbe_End - as.numeric(unlist(matchPatternProbe_Start)) + 1
                                primerprobes$startPosProbe <- matchPatternProbe_Start
                                primerprobes$endPosProbe <- matchPatternProbe_End
                                
                                # Double check to make sure probe annealing is above the min and subset these out of the results
                                annealProbeCheck2 <- which(primerprobes$annealDiffProbe > primerVar$annealTempProbe)
                                primerprobes <- primerprobes[annealProbeCheck2,]
                                
                                if(dim(primerprobes)[1] == 0){ 
                                  shinyalert("Oh no!", "Sorry no probes were found with such a high annealing temperature. You will have to lower the min annealing temperature of the probe.", type = "error")
                                  updateNavbarPage(session, "mainPage", selected = 'panel3')
                                  
                                } else {
                                
                                  # Check to make sure probes dont overlap with primers
                                  probePosCheck1 <- which(primerprobes$endPosProbe < primerprobes$startPosRP)
                                  primerprobes <- primerprobes[probePosCheck1,]
                                  probePosCheck2 <- which(primerprobes$startPosProbe > primerprobes$endPosFP)
                                  primerprobes <- primerprobes[probePosCheck2,]
                                  
                                  if(dim(primerprobes)[1] == 0){ 
                                    
                                    shinyalert("Oh no!", "Sorry no probes were found from within your target amplification region. You will have to increase either your target amplification region or your amplicon length range. 
                                    Also look at reducing the lengths of your probe if your amplicon length range is small.", type = "error")
                                    updateNavbarPage(session, "mainPage", selected = 'panel1')
                                    
                                  } else {
                                    
                                    incProgress(0.2, detail = "Preparing primer and probe sets...")
                                  
                                    # Target coverage calculation
                                    primerprobes$coverage <- as.numeric(primerprobes$coverage)
                                    primerprobes$groupCoverage <- as.numeric(primerprobes$groupCoverage)
                                    primerprobes$targetCoverage <- round((primerprobes$coverage + primerprobes$groupCoverage) / 2, 4)
                                    primerprobes <- primerprobes[order(primerprobes$targetCoverage, decreasing = TRUE),]
      
                                    # Selecting the top five matches (ordered by target coverage)
                                    primerprobes <- primerprobes[1:5,]
                                    primerprobes$setID <- 1:nrow(primerprobes)
                                    primerVar$setID <- primerprobes$setID
                                    updateSelectInput(session, "inSelectPrimer", choices = primerVar$setID)
                                    
                                    dfPrimers$primerBindTab <- primerprobes
                                    
                                    primerprobes <- (primerprobes[,c("setID","identifier","targetCoverage", "ampliconLength", "forward_primer","reverse_primer", "probe_seq", 
                                                                     "length_FP", "length_RP", "length_probe", "gcFP", "gcRP", "gcProbe", "FP_hybr_efficiency", 
                                                                     "RP_hybr_efficiency", "hybrEfficiencyAverage", "annealAverage", "annealDiffprimer", 
                                                                     "annealProbe", "annealDiffProbe", "startPosFP", "endPosFP", "startPosProbe", "endPosProbe",
                                                                     "startPosRP", "endPosRP", "ampliconSeq")])
                                    
                                    # Converting from decimal to percentage so its more clear to the user
                                    primerprobes[,3] <- primerprobes[,3] * 100
                                    primerprobes[,11:16] <- lapply(primerprobes[,11:16], as.numeric)
                                    primerprobes[,11:16] <- primerprobes[,11:16] * 100
                                    
                                    # Then converting to character again so can be easily imported to rhanstontable
                                    primerprobes[, ] <- lapply(primerprobes[, ], as.character)
                                    dfPrimers$primerTab <- primerprobes
                                    
                                    dfPrimers$primerGraphTab <- dfPrimers$primerTab
                                        
                                    dfPrimers$primerGraphTab$annealAveragePrimer <- dfPrimers$primerGraphTab$annealAverage
                                        
                                    dfPrimers$primerGraphTab <- (dfPrimers$primerGraphTab[,c("setID","targetCoverage", "ampliconLength", "gcFP", "gcRP", "gcProbe", 
                                                                                             "FP_hybr_efficiency", "RP_hybr_efficiency", "hybrEfficiencyAverage", 
                                                                                             "annealAveragePrimer", "annealDiffprimer", "annealProbe", "annealDiffProbe")])
                                    
                                    dfPrimers$primerGraphTab[, ] <- lapply(dfPrimers$primerGraphTab[, ], as.numeric)
                                    
                                    updateSelectInput(session, "inSelectPrimer2", choices = colnames(dfPrimers$primerGraphTab[,2:13]), selected = "targetCoverage")
                                    
                                    incProgress(0.2, detail = "Finished design of primer and probe sets!")
                                    
                                    Sys.sleep(0.5)
                                    
                                    # Primer Table using RHandsontable functionality 
                                    output$primerTable <- renderRHandsontable({
                                      rhandsontable(dfPrimers$primerTab) %>% hot_cols(columnSorting = TRUE, readOnly = FALSE, manualColumnResize = TRUE)
                                    })
                                    
                                    # Render text for dimerization
                                    p1 <- foreach(i=1:nrow(dfPrimers$primerTab)) %do% paste(">", dfPrimers$primerTab$setID[i], "_FP", sep="")
                                    p2 <- foreach(i=1:nrow(dfPrimers$primerTab)) %do% paste(dfPrimers$primerTab$forward_primer[i]) 
                                    p3 <- foreach(i=1:nrow(dfPrimers$primerTab)) %do% paste(">", dfPrimers$primerTab$setID[i], "_RP", sep="")
                                    p4 <- foreach(i=1:nrow(dfPrimers$primerTab)) %do% paste(dfPrimers$primerTab$reverse_primer[i])
                                    p5 <- foreach(i=1:nrow(dfPrimers$primerTab)) %do% paste(">", dfPrimers$primerTab$setID[i], "_Probe", sep="")
                                    p6 <- foreach(i=1:nrow(dfPrimers$primerTab)) %do% paste(dfPrimers$primerTab$probe_seq[i])
                                    p7 <- rep(paste(""), nrow(dfPrimers$primerTab))
                                    p8 <- rep(paste(""), nrow(dfPrimers$primerTab))
                                                                  
                                    output$primertext = renderUI({                              
                                      HTML(paste(p1, p2, p3, p4, p5, p6, p7, p8, sep='<br/>'))                          
                                    })
                                    
                                    # Switch to last panel with primer and probe table
                                    updateNavbarPage(session, "mainPage", selected = 'panel5')
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      
    })
    
  })
  
  
  # Prevent timeouts 
  output$keepAlive <- renderText({
    req(input$count)
    paste("keep alive ", input$count)
  })
  
  # Primer binding visualization code
  
  # Subset by selected primer and probe set
  observeEvent(input$inSelectPrimer, {
    dfPrimers$primerBindTab2 <- dfPrimers$primerBindTab[dfPrimers$primerBindTab$setID %in% input$inSelectPrimer,]
  })
  
  observeEvent(input$run4a, {
    
    dbConn <- dbConnect(SQLite(), ":memory:")
    Seqs2DB(reactUpload$alignment, "XStringSet", dbConn, "userAlignment")
    desc <- dbGetQuery(dbConn, "select description from Seqs")
    desc <- unlist(lapply(strsplit(desc$description, "userAlignment", fixed=TRUE),
                             function(x) return(x[length(x)])))
    Add2DB(data.frame(identifier=desc), dbConn)
    
    # Subset alignment by target amplification region and reorganize alignment with target species as first sequences shown  
    dna <- SearchDB(dbConn)
    alignment <- subseq(dna, reactSlider1$sliderInput2, reactSlider1$sliderInput3)
    names(alignment) <- desc
    u_alignment <- unique(alignment)
    names(u_alignment) <- names(alignment)[match(u_alignment, alignment)]
    alignment <- u_alignment
    # move the target group to the top
    w <- which(names(alignment)==reactTarget$target)
    alignment <- c(alignment[w], alignment[-w])
    
    # Code for generating the alignment visualization
    BrowseSeqs(alignment, htmlFile = paste(tempdir(), "/myAlignment2.html", sep = ""), openURL = FALSE, colorPatterns=c(dfPrimers$primerBindTab2$startPosFP, 
               dfPrimers$primerBindTab2$endPosFP, dfPrimers$primerBindTab2$startPosProbe, dfPrimers$primerBindTab2$endPosProbe, 
               dfPrimers$primerBindTab2$startPosRP, dfPrimers$primerBindTab2$endPosRP), highlight=1, colors=c("0072B2", "#E69F00", "#56B4E9", "#009E73", "#999999"))
    
    getPage2 <-function() {
      return(includeHTML( file ( paste(tempdir(), "/myAlignment2.html", sep = "") )))
    }
    
    output$inc2 <-renderUI({getPage2()})
      
  })
  
  # Primer and probe bar graph visualization code
  observeEvent(input$inSelectPrimer2, {
    
    dfPrimers$primerGraphTab$setID <- as.character(dfPrimers$primerGraphTab$setID)
    
    # Code for plotting bar graph according to primer and probe variable
    output$plot <- renderPlot({
      ggplot(dfPrimers$primerGraphTab, aes_string(x="setID", y=input$inSelectPrimer2, fill="setID")) + geom_bar(stat="identity", width=0.30) +
        theme(axis.title=element_text(size=15), axis.text.x=element_text(face="bold", size=12), axis.text.y=element_text(face="bold", size=12), 
              legend.position = "none", legend.title = element_blank()) + 
        scale_fill_manual(values = c("#0072B2", "#E69F00", "#56B4E9", "#009E73", "#999999"))
    })
  })
  
  # Primer and probe heat map visualization code
  observeEvent(input$heatmap, {
    
    dfHeatMap <- dfPrimers$primerGraphTab

    # Adding max and min parameters set during the design process for comparison on the heatmap
    dfHeatMap[nrow(dfHeatMap) + 1,] = list(setID="Max Parameter", targetCoverage="100", ampliconLength=primerVar$maxProductSize, gcFP=primerVar$gcPrimerMax * 100, 
                                           gcRP=primerVar$gcPrimerMax * 100, gcProbe=primerVar$gcProbeMax * 100, FP_hybr_efficiency="100", RP_hybr_efficiency="100", 
                                           hybrEfficiencyAverage="100", annealAveragePrimer=primerVar$annealRangeMax, 
                                           annealDiffprimer=primerVar$annealDiffPrimerMax, annealProbe="80", annealDiffProbe="20")
    
    dfHeatMap[nrow(dfHeatMap) + 1,] = list(setID="Min Parameter", targetCoverage="20", ampliconLength=primerVar$minProductSize, gcFP=primerVar$gcPrimerMin * 100, 
                                           gcRP=primerVar$gcPrimerMin * 100, gcProbe=primerVar$gcProbeMin * 100, FP_hybr_efficiency="20", RP_hybr_efficiency="20", 
                                           hybrEfficiencyAverage="20", annealAveragePrimer=primerVar$annealRangeMin, annealDiffprimer="0",
                                           annealProbe="40", annealDiffProbe=primerVar$annealTempProbe)
    dfHeatMap <- rbind(dfHeatMap[6,], dfHeatMap[1:5,], dfHeatMap[7,])
    
    dfHeatMap2 <- dfHeatMap[,2:13]

    m <- as.matrix(sapply(dfHeatMap2, as.numeric))  
    colnames(m) <- colnames(dfHeatMap2)
    rownames(m) <- dfHeatMap$setID
    
    # Code for plotting heatmap according to primer and probe variables using d3heatmap
    output$heatmap <- renderD3heatmap({
      d3heatmap(m, colors = c("#d0d1e6", "#a6bddb", "#74a9cf", "#3690c0", "#0570b0", "#045a8d", "#023858"), scale = "column", dendrogram = "none", 
                xaxis_font_size = 15, yaxis_font_size = 15, xaxis_height = 150)
    })
  })

  # Shiny Download handler for download of generated primer sets
  output$downloadPrimers <- downloadHandler(
    filename = function() { 
      paste(gsub("\\s", "", reactTarget$target), "_", "PrimerProbeSet", ".csv", sep="")
    },
    content = function(file) {
      write.csv(dfPrimers$primerTab, file)
    })
  
  
  # Reset to defaults for each Tab
  observeEvent(input$reset1, {
    shinyjs::reset("qPCRConditions")
  })
  
  observeEvent(input$reset2a, {
    shinyjs::reset("Primers")
  })
  
  observeEvent(input$reset2b, {
    shinyjs::reset("Primers")
  })
  
  observeEvent(input$reset3, {
    shinyjs::reset("globalConstraints")
  })
  
}

shinyApp(ui, server)
