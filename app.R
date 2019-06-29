# Please ensure R version: 3.6.0 is installed prior to running this tool

# Change max file upload to 20 Mb
options(shiny.maxRequestSize = 20*1024^2)
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

exampleAlignment <- readDNAStringSet("sampleAlignment.fas")

ui <- tagList(useShinyalert(),
              useShinyjs(),
              tags$head(
                tags$style(
                  HTML(".shiny-notification {
                       height: 100px;
                       width: 800px;
                       position:fixed;
                       top: calc(50% - 50px);;
                       left: calc(50% - 400px);;
                       }
                       "
                  )
                )
              ),
              navbarPage("ssPRIMER - A Web-Based Tool for Species-Specific Primer Design", id = "mainPage",
                         tabPanel("1. Upload Alignment", id='panel1', value="panel1",
                                  sidebarLayout(
                                    sidebarPanel(
                                      p("Welcome to the ssPRIMER homepage! ssPRIMER is a GUI based tool that provides a straightforward 
                                        process to designing species-specific primers and Taqman probe for qPCR Assays.
                                        To start please upload a mulitple sequence alignment in fasta format. 
                                        More options will then be presented to begin the design process.
                                        You can also run a test alignment to test out the functionality of the tool. 
                                        Once uploaded, an interactive visual representation of the alignment will be shown to the 
                                        right using integration primarily with the DECIPHER, Biostrings and Tm Calculator packages."),
                                        
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
                                              by your primer set. It is recommeded you do not select a region smaller than 75 bp 
                                              as this will greatly limit the number of potential primer sets. The alignment 
                                              visualization will regenerate with the region specified automatically."),
                                            
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
                                                                       with the specified primer set each PCR cycle)",
                                                                       min = 0, max = 30, value = 10, post = "%"),
                                                               
                                                          # Primer set Tm
                                                          #h5("Target annealing temperature for the primer set."),
                                                          sliderInput("inSlider15", "Optimal Primer Set Annealing Temperature",
                                                                      min = 40, max = 70, value = 60, step = 0.1, post = "C"),
                                                               
                                                          sliderInput("inSlider29", "Max Primer Set Tm Difference (ideally should be within 0-3 degrees difference)",
                                                                      min = 0, max = 6, value = 3, step = 0.1, post = "C"),
                                                               
                                                          sliderInput("inSlider26", "Min Increase in Probe Annealing Temperature (ideally should be 7-10 degrees higher than primers)",
                                                                      min = 1, max = 20, value = 7, step = 0.1, post = "C"),

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
                                                                           min = 20, max = 40, value = c(25, 28), step = 1, post = "bp"),
                                                               
                                                               sliderInput("inSlider27", "Primer GC Percentage Range",
                                                                           min = 30, max = 80, value = c(30, 70), step = 0.1, post = "%"), 
                                                               
                                                               sliderInput("inSlider28", "Probe GC Percentage Range",
                                                                           min = 30, max = 80, value = c(30, 70), step = 0.1, post = "%"), 
                                                               
                                                               sliderTextInput("inSlider23", "Primer-Dimer Hybridization Efficiency (less likely to form dimers < 1e-07 > more likely to form dimers)",
                                                                               choices = c("1e-09", "1e-08", "1e-07", "1e-06", "1e-05"), 
                                                                               selected = c("1e-07")),
                                                               
                                                               p("***These settings are still being tested and are not currently functional***"),
                                                               
                                                               sliderInput("inSlider20", "Max Run Length for Primers and Probes (ex: AAA would be a run length of 3)",
                                                                           min = 0, max = 10, value = 4, step = 1, post = "bp"),
                                                               
                                                               sliderInput("inSlider21", "Max Repeat Length for Primers and Probes (ex: CGCG would be a repeat length of 2)",
                                                                           min = 0, max = 5, value = 2, step = 1, post = "bp"),
                                                               
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
                                                           
                                                           div(rHandsontableOutput('primerTable'))),
                                                  
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
                                                  
                                                  tabPanel("Primer Set Comparison",  id='panelPrimerGraphs', value="panelPrimerGraphs",
                                                           
                                                           tags$hr(), 
                                                           
                                                           p("Here you can visualize in bar graph form the differences between each primer and probe
                                                             set based on many different parameters selectable from the dropdown menu shown below."),
                                                           
                                                           tags$hr(), 
                                                  
                                                           selectInput("inSelectPrimer2", "Choose a Parameter to Compare across Sets",
                                                                       choices = c("Waiting for Primer Sets")),
                                                           
                                                           tags$hr(),
                                                           
                                                           plotOutput('plot')
                                                  )
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
                                          p("Wright ES (2016). “Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in R.” The R Journal, 8(1), 352-359."),
                                          
                                          tags$hr(),
                                          
                                          p("Citation for Biostrings:"),
                                          p("Pagès H, Aboyoun P, Gentleman R, DebRoy S (2019). Biostrings: Efficient manipulation of biological strings. 
                                             R package version 2.52.0."),
                                          
                                          tags$hr(),
                                          
                                          p("Citations for OligoArrayAux:"), 
                                          p("Markham, N. R. & Zuker, M. (2005) DINAMelt web server for nucleic acid melting prediction. 
                                             Nucleic Acids Res., 33, W577-W581."),
                                          p("Markham, N. R. & Zuker, M. (2008) UNAFold: software for nucleic acid folding and hybridization. 
                                             In Keith, J. M., editor, Bioinformatics, Volume II. Structure, Function and Applications, 
                                             number 453 in Methods in Molecular Biology, chapter 1, pages 3–31. Humana Press, 
                                             Totowa, NJ. ISBN 978-1-60327-428-9.") 

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
                 patterns=c("A", "C", "G", "T", "-"), 
                 colors=c("#999999", "#E69F00", "#56B4E9", "#009E73","#0072B2"), colWidth=reactSlider1$sliderInput4)
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
                              endPosFP=NULL, startPosProbe=NULL, endPosProbe=NULL, startPosRP=NULL, endPosRP=NULL, chosenSetId=NULL)
  
  # Observe primer design choices selected by user
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
    primerVar$annealTempPrimer <- as.numeric(input$inSlider15)
  })
  
  observeEvent(input$inSlider16[1], {
    primerVar$minLength <- as.numeric(input$inSlider16[1])
  })
  
  observeEvent(input$inSlider16[2], {
    primerVar$maxLength <- as.numeric(input$inSlider16[2])
  })
  
  observeEvent(input$inSlider23, {
    primerVar$dimers <- as.numeric(input$inSlider23)
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
    primerVar$probeMinLength <- as.numeric(input$inSlider30[1])
  })
  
  observeEvent(input$inSlider30[2], {
    primerVar$probeMaxLength <- as.numeric(input$inSlider30[2])
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
  
  # Generating primer and probe sets 
  observeEvent(input$run3a, {
    
    withProgress(message = "Starting Design Process...", value = 0, {
      
      Sys.sleep(0.5)
      
      updateNavbarPage(session, "mainPage", selected = 'panel5')
      
      dbConn <- dbConnect(SQLite(), ":memory:")
      Seqs2DB(reactUpload$alignment, "XStringSet", dbConn, "userAlignment")
      desc <- dbGetQuery(dbConn, "select description from Seqs")
      desc <- unlist(lapply(strsplit(desc$description, "userAlignment", fixed=TRUE),
                            function(x) return(x[length(x)])))
      Add2DB(data.frame(identifier=desc), dbConn)
      
      tiles <- TileSeqs(dbConn)
      
      incProgress(0.3, detail = "Starting design of primers...")
      
      Sys.sleep(0.5)
      
      # Generate Initial Primer and Probe Sets that will then be filtered by GC content further down
      primers <- try(DesignPrimers(tiles, identifier=reactTarget$target, start=primerVar$startPos, end=primerVar$endPos,
                                   minLength=primerVar$minLength, maxLength=primerVar$maxLength, minCoverage=primerVar$minCoverage,
                                   minGroupCoverage=primerVar$minGCoverage, annealingTemp=primerVar$annealTempPrimer, 
                                   P=(primerVar$primerConc/1000000000), monovalent=((primerVar$monovalentNa+primerVar$monovalentK)/1000), 
                                   divalent=(primerVar$divalentMg/1000), dNTPs=(primerVar$dNTP/1000), minEfficiency=primerVar$minEfficiency, 
                                   numPrimerSets=30, minProductSize=primerVar$minProductSize, maxProductSize=primerVar$maxProductSize, 
                                   primerDimer=primerVar$dimers), silent = TRUE)
      
      if(dim(primers)[1] == 0){ 
        # Error for initial design of the primer sets
         shinyalert("Oh no!", "Sorry no primer sets met your specified constraints. Please check the following constraints: amplicon length, 
                    target amplification region, qPCR reaction settings, min and max primer length, primer-dimer hybridization efficiency,
                    min target and min group coverage, min target hybridization efficiency and primer annealing temperature.
                    The most probable causes are either your target hybridization efficiency value is set too high or your min target 
                    coverage value is set too high.", type = "error")
        updateNavbarPage(session, "mainPage", selected = 'panel3')
      } else {
        primers$setID <- row.names(primers)
        primers$length_FP <- nchar(primers$forward_primer)
        primers$length_RP <- nchar(primers$reverse_primer)
        
        # Subsetting alignment by target amplification region and using that alignment for calculation of primers
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
        
        # Amplicon length calculation
        primers$ampliconLength <- primers$endPosRP - primers$startPosFP + 1
        
        primers$FP_hybr_efficiency <- round(primers$forward_efficiency, 2)
        primers$RP_hybr_efficiency <- round(primers$reverse_efficiency, 2)
        
        # Calculating GC content
        gcFP <- foreach(i=1:nrow(primers)) %do% letterFrequency(DNAString(primers$forward_primer[i]), "GC", as.prob=TRUE) 
        gcRP <- foreach(i=1:nrow(primers)) %do% letterFrequency(DNAString(primers$reverse_primer[i]), "GC", as.prob=TRUE)
            
        # Filtering by GC Content Criteria set by the user
        gcCheckFP <- which(gcFP > primerVar$gcPrimerMin & gcFP < primerVar$gcPrimerMax)
        gcCheckRP <- which(gcRP > primerVar$gcPrimerMin & gcRP < primerVar$gcPrimerMax)
        gcCheck <- intersect(gcCheckFP, gcCheckRP)
        
        if(length(gcCheck) == 0){ 
           shinyalert("Oh no!","Sorry no primer sets met your specified constraints. Please check the following constraint: 
                      % gc content for primers. You will have to increase the % gc content range of your primers.", type = "error")
           updateNavbarPage(session, "mainPage", selected = 'panel3b')
        } else {
          # Assigning GC content to primers dataframe
          primers <- primers[gcCheck, ]
          gcFP <- foreach(i=1:nrow(primers)) %do% letterFrequency(DNAString(primers$forward_primer[i]), "GC", as.prob=TRUE)
          gcFP <- unlist(gcFP)
          gcRP <- foreach(i=1:nrow(primers)) %do% letterFrequency(DNAString(primers$reverse_primer[i]), "GC", as.prob=TRUE)
          gcRP <- unlist(gcRP)
          primers$gcFP <- round(gcFP, 2)
          primers$gcRP <- round(gcRP, 2)

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
                         Max non-target hybridization efficiency. You will have to increase this value to allow more inclusivity 
                         of non-target species/OTUs.", type = "error")
             updateNavbarPage(session, "mainPage", selected = 'panel3')
          } else {
            primers <- (primers[,c("setID", "identifier", "ampliconLength", "startPosFP", "endPosFP", "length_FP","start_reverse", 
                                   "startPosRP", "endPosRP", "length_RP", "forward_primer", "reverse_primer", "FP_hybr_efficiency", 
                                   "RP_hybr_efficiency", "gcFP", "gcRP", "forward_coverage", "reverse_coverage")])
            
            primers <- data.frame(primers, stringsAsFactors = FALSE)
            primers <- as.data.frame(lapply(primers, function(y) gsub(",", "", y)))
            primers <- primers[complete.cases(primers), ]
            primers[, ] <- lapply(primers[, ], as.character)
            primers$setID <- 1:nrow(primers)
            
            # Averaging of target coverage across forward and reverse primers
            primers$forward_coverage <- as.numeric(primers$forward_coverage)
            primers$reverse_coverage <- as.numeric(primers$reverse_coverage)
            primers$coverage <- (primers$forward_coverage + primers$reverse_coverage) / 2
            
            # Annealing temperature calculation
            predictedAnnealing_FP <- foreach(i=1:nrow(primers)) %do% Tm_GC(primers$forward_primer[i], ambiguous = FALSE, userset = NULL, variant = "Primer3Plus",
                                                                           Na = primerVar$monovalentNa, K = primerVar$monovalentNa, Tris = 0, 
                                                                           Mg = primerVar$divalentMg, dNTPs = primerVar$dNTP, saltcorr = 0, mismatch = TRUE)
            primers$Annealing_FP <- round(as.numeric(unlist(predictedAnnealing_FP)) - 3, 2)
            
            predictedAnnealing_RP <- foreach(i=1:nrow(primers)) %do% Tm_GC(primers$reverse_primer[i], ambiguous = FALSE, userset = NULL, variant = "Primer3Plus",
                                                                           Na = primerVar$monovalentNa, K = primerVar$monovalentNa, Tris = 0, 
                                                                           Mg = primerVar$divalentMg, dNTPs = primerVar$dNTP, saltcorr = 0, mismatch = TRUE)
            primers$Annealing_RP <- round(as.numeric(unlist(predictedAnnealing_RP)) - 3, 2)
            
            primers$annealDiffprimer <- round(abs(primers$Annealing_FP - primers$Annealing_RP), 2)
            
            annealDiffCheck <- which(primers$annealDiffprimer > primerVar$annealDiffPrimerMax)
            primers <- primers[-annealDiffCheck,]
            
            if(dim(primers)[1] == 0){ 
              shinyalert("Oh no!", "Sorry no primer sets met your specified constraints. Please check the following constraint: 
                          Max anneal temperature difference between FP and RP. You will have to increase this value.", type = "error")
              updateNavbarPage(session, "mainPage", selected = 'panel3')
            } else {
              annealAverage <- (as.numeric(primers$Annealing_FP) + as.numeric(primers$Annealing_RP)) / 2
              primers$annealAverage <- annealAverage
              
              primerHybrAverage <- (as.numeric(primers$FP_hybr_efficiency) + as.numeric(primers$RP_hybr_efficiency)) / 2
              primers$hybrEfficiencyAverage <- primerHybrAverage
              
              # Calculating probe start and stop positions that must be met
              probeStart <- as.numeric(primers$startPosFP) + as.numeric(primers$length_FP) + 1
              probeEnd <- as.numeric(primers$endPosRP) - as.numeric(primers$length_RP) - 1
              
              incProgress(0.4, detail = "Starting design of probes...")
              
              # Four different tileseqs at different probe lengths within probe length range set by user
              probesT1 <- TileSeqs(dbConn, identifier = reactTarget$target, minLength = (primerVar$probeMinLength)-1, maxLength = primerVar$probeMinLength)
              
              probesT2 <- TileSeqs(dbConn, identifier = reactTarget$target, minLength = (primerVar$probeMinLength), maxLength = (primerVar$probeMinLength)+1)
              
              probesT3 <- TileSeqs(dbConn, identifier = reactTarget$target, minLength = (primerVar$probeMaxLength)-2, maxLength = (primerVar$probeMaxLength)-1)
              
              probesT4 <- TileSeqs(dbConn, identifier = reactTarget$target, minLength = (primerVar$probeMaxLength)-1, maxLength = primerVar$probeMaxLength)
              
              dfTarget <- rbind(probesT1, probesT2, probesT3, probesT4)
              
              # Find probes within target amplification region
              probeFind <- foreach(i=1:length(probeStart)) %do%  which(dfTarget$start_aligned > probeStart[i] & dfTarget$end_aligned < probeEnd[i])
              listProbeFind <- foreach(i=1:length(probeFind)) %do% dfTarget[probeFind[[i]],]
              
              if(length(listProbeFind[[1]]) == 0){ 
                shinyalert("Oh no!", "Sorry no probes were found from within your target amplification region. You will have to increase either your target 
                            amplification region or your amplicon length range. Also look at reducing the length range of your probe if your amplicon
                            length range is small.", type = "error")
                updateNavbarPage(session, "mainPage", selected = 'panel1')
              } else {
                names(listProbeFind) <- 1:nrow(primers)
                probes <- do.call("rbind", listProbeFind)
                probes$setID <- round(as.numeric(rownames(probes)), 1)
                probes$setID <- round(probes$setID, 0)
                probes$length_probe <- nchar(probes$target_site)
                
                # Match to target species using group coverage - must be highly similar
                groupCoverage <- which(probes$groupCoverage >= primerVar$minCoverage)
                probes <- probes[groupCoverage,]
                
                if(dim(probes)[1] == 0){ 
                  shinyalert("Oh no!", "Sorry no probes were found with such a high target coverage. You will have to lower the value for this setting.",
                             type = "error")
                  updateNavbarPage(session, "mainPage", selected = 'panel3')
                } else {
                  # GC Check for the probe
                  gcProbe <- foreach(i=1:nrow(probes)) %do% letterFrequency(DNAString(probes$target_site[i]), "GC", as.prob=TRUE) 
                  gcProbe <- unlist(gcProbe)
                  probes$gcProbe <- round(gcProbe, 2) 
                  probes$probe_seq <- probes$target_site  
                  
                  gcCheck <- which(probes$gcProbe >= primerVar$gcProbeMin & probes$gcProbe <= primerVar$gcProbeMax)
                  probes <- probes[gcCheck,]
                  
                  if(dim(probes)[1] == 0){ 
                    shinyalert("Oh no!", "Sorry no probes were found in the gc % range that was set. You will have to increase the range of 
                                acceptable % gc content for the probes.", type = "error")
                    updateNavbarPage(session, "mainPage", selected = 'panel3b')
                  } else {
                    # Annealing Check for the probe
                    annealProbe <- foreach(i=1:nrow(probes)) %do% Tm_GC(probes$probe_seq[i], ambiguous = FALSE, userset = NULL, variant = "Primer3Plus",
                                                                        Na = primerVar$monovalentNa, K = primerVar$monovalentNa, Tris = 0, 
                                                                        Mg = primerVar$divalentMg, dNTPs = primerVar$dNTP, saltcorr = 0, mismatch = TRUE)
                    annealProbe <- unlist(annealProbe)
                    probes$annealProbe <- round(annealProbe, 2)
                    
                    probes <- (probes[,c("setID", "length_probe", "start_aligned", "end_aligned", "groupCoverage", "annealProbe", "gcProbe", "probe_seq")])
                    
                    probeList <- lapply(unique(probes$setID), function(x) 
                      probes[probes$setID == x,])
                    
                    annealCheck <- foreach(i=1:length(probeList)) %do% which(probeList[[i]]$annealProbe >= (annealAverage[i] + primerVar$annealTempProbe))
                    
                    if(length(annealCheck[[1]]) == 0){ 
                      shinyalert("Oh no!", "Sorry no probes were found with such a high annealing temperature. You will have to lower the 
                                  min annealing temperature of the probe.", type = "error")
                      updateNavbarPage(session, "mainPage", selected = 'panel3')
                    } else {
                      probeList <- foreach(i=1:length(probeList)) %do% probeList[[i]][annealCheck[[i]],]
                      probeList <- foreach(i=1:length(probeList)) %do% probeList[[i]][order(probeList[[i]]$groupCoverage),]
                      probeList <- foreach(i=1:length(probeList)) %do% probeList[[i]][1,]
                      
                      probes <- do.call("rbind", probeList)
                      probes <- probes[1:nrow(primers),]
                      probes$setID <- 1:nrow(primers)
                      
                      # Merging primer sets to probes and modifying/organizing probe columns
                      primerprobes <- merge(primers, probes, by.x="setID", by.y="setID")
                      primerprobes$setID <- 1:nrow(primerprobes)
                      primerprobes$annealDiffProbe <- round(primerprobes$annealProbe - primerprobes$annealAverage, 2)
                      primerprobes$startPosFP <- as.numeric(primerprobes$startPosFP)
                      primerprobes$endPosFP <- as.numeric(primerprobes$endPosFP)
                      primerprobes$startPosProbe <- as.numeric(primerprobes$start_aligned)
                      primerprobes$endPosProbe <- as.numeric(primerprobes$end_aligned)
                      primerprobes$startPosRP <- as.numeric(primerprobes$startPosRP)
                      primerprobes$endPosRP <- as.numeric(primerprobes$endPosRP)
                      primerprobes$coverage <- as.numeric(primerprobes$coverage)
                      primerprobes$groupCoverage <- as.numeric(primerprobes$groupCoverage)
                      primerprobes$targetCoverage <- round((primerprobes$coverage + primerprobes$groupCoverage) / 2, 4)
                      primerprobes <- primerprobes[order(primerprobes$targetCoverage, decreasing = TRUE),]
                      primerprobes <- primerprobes[1:5,]
                      primerprobes$setID <- 1:nrow(primerprobes)
                      
                      primerVar$setID <- primerprobes$setID
                      updateSelectInput(session, "inSelectPrimer", choices = primerVar$setID)
                      
                      dfPrimers$primerBindTab <- primerprobes
                      
                      primerprobes <- (primerprobes[,c("setID","identifier","targetCoverage", "ampliconLength", "forward_primer","reverse_primer", "probe_seq", 
                                                       "length_FP", "length_RP", "length_probe", "gcFP", "gcRP", "gcProbe", "FP_hybr_efficiency", 
                                                       "RP_hybr_efficiency", "hybrEfficiencyAverage", "annealAverage", "annealDiffprimer", 
                                                       "annealProbe", "annealDiffProbe", "startPosFP", "endPosFP", "startPosProbe", "endPosProbe",
                                                       "startPosRP", "endPosRP")])
                      primerprobes[, ] <- lapply(primerprobes[, ], as.character)
                      dfPrimers$primerTab <- primerprobes
                      
                      dfPrimers$primerGraphTab <- dfPrimers$primerTab
                          
                      dfPrimers$primerGraphTab$annealAveragePrimer <- dfPrimers$primerGraphTab$annealAverage
                          
                      dfPrimers$primerGraphTab <- (dfPrimers$primerGraphTab[,c("setID","targetCoverage", "ampliconLength", "gcFP", "gcRP", "gcProbe", 
                                                                               "FP_hybr_efficiency", "RP_hybr_efficiency", "hybrEfficiencyAverage", 
                                                                               "annealAveragePrimer", "annealDiffprimer", "annealProbe", "annealDiffProbe")])
                      
                      dfPrimers$primerGraphTab[, ] <- lapply(dfPrimers$primerGraphTab[, ], as.numeric)
                      
                      updateSelectInput(session, "inSelectPrimer2", choices = colnames(dfPrimers$primerGraphTab[,2:13]), selected = "targetCoverage")
                      
                      incProgress(0.3, detail = "Finished design of primer and probe sets!")
                      
                      Sys.sleep(0.5)
                      
                      # Primer Table using RHandsontable functionality 
                      output$primerTable <- renderRHandsontable({
                        rhandsontable(dfPrimers$primerTab) %>% hot_cols(columnSorting = TRUE, readOnly = FALSE, manualColumnResize = TRUE)
                      })
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
               dfPrimers$primerBindTab2$startPosRP, dfPrimers$primerBindTab2$endPosRP), highlight=1, colors=c("#999999", "#E69F00", "#56B4E9", "#009E73","0072B2"))
    
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
              legend.title = element_blank()) + 
        scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73","#0072B2")) 
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
