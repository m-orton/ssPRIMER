# Change max file upload is 30 Mb
options(shiny.maxRequestSize = 30*1024^2)

# source("https://bioconductor.org/biocLite.R")
# biocLite("DECIPHER")
# biocLite("Biostrings")
library(Biostrings)
library(DECIPHER)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("shiny")
library(shiny)
# install.packages("shinyjs")
library(shinyjs)
# install.packages("rhandsontable")
library(rhandsontable)
# install.packages("jsonlite")
library(jsonlite)
# install.packages("plyr")
library(plyr)

exampleAlignment <- readDNAStringSet("sampleAlignment.fas")

ui <- tagList(useShinyjs(),
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
                               process to designing species-specific primers for Taqman probe based qPCR Assays.
                               To start please upload a mulitple sequence alignment in fasta format. 
                               More options will then be presented to begin the design process.
                               You can also run a test alignment to test out the functionality of the tool. 
                               Once uploaded, an interactive visual representation of the alignment will be shown to the 
                               right using integration primarily with the Decipher R and Bioconductor packages. 
                               Design of probes along with primer sets is currently being tested and will be added soon."),
                            
                               tags$hr(),
                            
                                  div(id="alignmentLoad",
                                          
                                     actionButton('exampleAlign', "Load Test Alignment", icon("info-circle"), 
                                                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                          
                                     tags$hr(),
                                          
                                     p("Please ensure that the sequences for the target species are named in exactly the same format.
                                       Ex: Limnephilus_hyalinus and Limnephilus hyalinus would not be considered the same target group."),
                                          
                                        fileInput('alignmentFile', 'Choose an alignment to upload (.fas or .fa format), 
                                                   max file size = 30 Mb. Enusre all sequences are the same length before 
                                                   uploading, this may require an alignment trimming step before upload.',
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
                                            p("First select a target group from your alignment that you would like to be 
                                               selectively amplified."),
                                        
                                            selectInput("inSelectTarget", "Target Species",
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
                                        div(id="stepRun",
                                            actionButton('run1', "Start Design Process!", icon("angle-double-right"), 
                                                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                        ),
                                        width = 3
                                      )),
                                    
                                    mainPanel(
                                      
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
                         
                         tabPanel("3. Set Primer Constraints", id='panel3', value="panel3",
                             sidebarLayout(
                                   sidebarPanel(
                                      
                                    div(id="Primers",   
                                          
                                       tabsetPanel(id="panels",     
                                                      
                                          tabPanel("Primary Constraints", id='panel3a', value="panel3a",
                                                               
                                             tags$hr(),      
                                                               
                                             p("Next, set the contraints on your primers before designing them. 
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
                                                          min = 20, max = 100, value = 80, post = "%"),
                                                               
                                              # Minimum Amplification Efficiency
                                              #h5("Minimum efficiency of hybridization desired for the primer set."),
                                              sliderInput("inSlider13", "Min Target Hybridization Efficiency (minimum fraction 
                                                          of target amplicons that will be amplified with the specified primer 
                                                          set each PCR cycle)",
                                                          min = 50, max = 100, value = 80, post = "%"),
                                                               
                                              sliderInput("inSlider14", "Max Non-Target Hybridization Efficiency 
                                                         (maximum fraction of non-target amplicons that will be amplified 
                                                          with the specified primer set each PCR cycle)",
                                                          min = 0, max = 20, value = 10, post = "%"),
                                                               
                                              # Primer set Tm
                                              #h5("Target annealing temperature for the primer set."),
                                              sliderInput("inSlider15", "Optimal Primer Set Tm",
                                                          min = 50, max = 70, value = 60, step = 0.1, post = "C"),
                                                               
                                              sliderInput("inSlider16", "Primer Length Range",
                                                          min = 18, max = 40, value = c(18, 23), post = "bp"),
                                                               
                                              sliderInput("inSlider18", "Number of Potential Primer Sets to Generate",
                                                          min = 1, max = 15, value = 5, step = 1),
                                                               
                                              actionButton('back2', "Back", icon("angle-double-left"), 
                                                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                                               
                                              actionButton('reset2a', "Reset to Default", icon("undo"), 
                                                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                                               
                                              actionButton('advSettings', "Go to Secondary Constraints", icon("cogs"), 
                                                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                                               
                                              actionButton('run3a', "Design Primer Sets!", icon("angle-double-right"), 
                                                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                              ),
                                                      
                                              tabPanel("Secondary Constraints", id='panel3b', value="panel3b",

                                              tags$hr(),
                                                               
                                              p("***These settings are still being tested and are not currently functional***"),
                                                               
                                              sliderInput("inSlider17", "Primer GC Ratio Range",
                                                          min = 35, max = 75, value = c(40, 60), step = 0.1, post = "%"), 
                                                               
                                              sliderInput("inSlider19", "GC Clamp Length",
                                                          min = 0, max = 5, value = c(1, 3), step = 1, post = "bp"),
                                                               
                                              sliderInput("inSlider20", "Max Run Length",
                                                          min = 0, max = 10, value = c(0, 4), step = 1, post = "bp"),
                                                               
                                              sliderInput("inSlider21", "Max Repeat Length",
                                                          min = 0, max = 10, value = c(0, 3), step = 1, post = "bp"),

                                              sliderInput("inSlider22", "Secondary Structure (Delta G)",
                                                          min = -10, max = 0, value = c(-9, 0), step = 1, post = " kcal/mol"),
                                                               
                                              sliderInput("inSlider23", "Self-Dimerization (Delta G)",
                                                          min = -10, max = 0, value = c(-9, 0), step = 1, post = " kcal/mol"),
                                                               
                                              sliderInput("inSlider24", "Cross-Dimerization (Delta G)",
                                                          min = -10, max = 0, value = c(-9, 0), step = 1, post = " kcal/mol"),
                                                               
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
                                      
                                      p("Download primer sets in .csv"), 
                                      
                                      div(id="PrimerTableDownload",
                                          downloadButton('downloadPrimers', 'Download')
                                      ), width = 2),
                                    
                                    mainPanel(
                                      tabsetPanel(type = "tabs",
                                         tabPanel("Primer Sets", div(rHandsontableOutput('primerTable')))
                                         # ***Tab section in progress***         
                                         # tabPanel("Primer Binding Visualizations",  id='panelPrimerBind', value="panelPrimerBind",
                                                           
                                         # tags$hr(),      
                                                           
                                         # p("Select one of the primer sets from the dropdown menu and a primer binding visualization will be generated based on your target species. The visualization will indicate 
                                         #   the presence of diagnostic bases present in the primer set that can be used to distinguish it from non-target species."),
                                                           
                                         # div(id="primerSelect", 
                                         # selectInput("inSelectPrimer", "Primer Set:",
                                         # choices = "Pending Upload"))
                                         # )
                                                           
                                      )))),
                         tabPanel("About", id='About', value="panel6",
                                  sidebarLayout(
                                    sidebarPanel(
                                      # Chose Reaction Conditions for qPCR
                                      div(id="about1",
                                          p("ssPRIMER (or species-specific PRIMER) is a web-based software tool that can be used to 
                                             design species-specific primer sets for qPCR assays. A multiple sequence alignment can 
                                             be imported in by a user, and the tool will then guide the user through the process of 
                                             designing and evaluating species-specific primer sets (and in futue iterations of the tool: 
                                             Taqman probes). The tool is designed to create primer sets that maximize amplification 
                                             efficiency for the target species (sensitivity) but minimize amplification efficiency 
                                             for non-target species (specificity). This tool is designed to benefit the users of eDNA 
                                             technology, including field biologists, ecologists, conservation researchers, and 
                                             environmental consultants and could contribute to environmental biomonitoring using 
                                             molecular methods."),
                                          
                                          tags$hr(),
                                          
                                          p("This tool relies greatly on the DECIPHER and Biostrings R packages for design of 
                                             primer sets:
                                             https://bioconductor.org/packages/release/bioc/html/DECIPHER.html
                                             http://bioconductor.org/packages/release/bioc/html/Biostrings.html"),
                                          
                                          tags$hr(),
                                          
                                          p("Author of tool: Matthew Orton, morton01@uoguelph.ca, 
                                             With contributions from Dr. Sally Adamowicz and Kamil Chatila-Amos") 
                                          
                                          ), width = 4),
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
    splitNames <- unlist(lapply(strsplit((names(reactUpload$alignment)), "[*]", fixed=TRUE), function(x) return(x[1])))
  })
  
  # Reactive varibale for alignment length
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
                 colors=c("#1E90FF", "#32CD32", "#9400D3", "#000000", "#EE3300"), colWidth=reactSlider1$sliderInput4)
    }
  })
  
  # Display html file for alignment visualization
  getPage<-function() {
    return(includeHTML( file (contentsrea3() )))
  }
  output$inc<-renderUI({
    getPage()
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
  
  reactTarget <- reactiveValues(
    target=NULL)
  
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
                              minCoverage=NULL, numPrimerSets=NULL, primerPosStart=NULL, primerPosEnd=NULL, gc=NULL)
  
  # Observe primer design choices
  observeEvent(input$inSlider1[1], {
    primerVar$startPos <- as.numeric(input$inSlider1[1])
  })
  
  observeEvent(input$inSlider1[2], {
    primerVar$endPos <- as.numeric(input$inSlider1[2])
  })
  
  observeEvent(input$inSlider16[1], {
    primerVar$minLength <- as.numeric(input$inSlider16[1])
  })
  
  observeEvent(input$inSlider16[2], {
    primerVar$maxLength <- as.numeric(input$inSlider16[2])
  })
  
  observeEvent(input$inSlider9[1], {
    primerVar$minProductSize <- as.numeric(input$inSlider9[1])
  })
  
  observeEvent(input$inSlider9[2], {
    primerVar$maxProductSize <- as.numeric(input$inSlider9[2])
  })
  
  observeEvent(input$inSlider15, {
    primerVar$annealTempPrimer <- as.numeric(input$inSlider15 - 5)
  })
  
  observeEvent(input$inSlider12, {
    primerVar$minCoverage <- as.numeric(input$inSlider12/100)
  })
  
  observeEvent(input$inSlider13, {
    primerVar$minEfficiency <- as.numeric(input$inSlider13/100)
  })
  
  observeEvent(input$inSlider18, {
    primerVar$numPrimerSets <- as.numeric(input$inSlider18)
  })
  
  observeEvent(input$inSlider17[1], {
    primerVar$gc <- as.numeric(input$inSlider17[1])
  })
  
  observeEvent(input$inSlider17[2], {
    primerVar$gc <- as.numeric(input$inSlider17[2])
  })
  
  # Var for table output
  dfPrimers <- reactiveValues(primerTab=NULL, primerBindTab=NULL)
  
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
  
  # Generating primer sets using the Decipher package
  observeEvent(input$run3a, {
    
    withProgress(message = "Starting Design Process...this may take a few mins...", value = 0, {
      
      Sys.sleep(1)
      
      updateNavbarPage(session, "mainPage", selected = 'panel5')
      
      dbConn <- dbConnect(SQLite(), ":memory:")
      Seqs2DB(reactUpload$alignment, "XStringSet", dbConn, "userAlignment")
      desc <- dbGetQuery(dbConn, "select description from Seqs")
      desc <- unlist(lapply(strsplit(desc$description, "userAlignment", fixed=TRUE),
                            function(x) return(x[length(x)])))
      #desc <- unlist(lapply(strsplit(desc, " ", fixed=TRUE), function(x) return(x[1])))
      Add2DB(data.frame(identifier=desc), dbConn)
      
      tiles <- TileSeqs(dbConn)
      
      incProgress(0.05, detail = "Finished loading design parameters...")
      incProgress(0.35, detail = "Starting design of primers...")
      
      Sys.sleep(1)
      
      primers <- DesignPrimers(tiles, identifier=reactTarget$target, start=primerVar$startPos, end=primerVar$endPos,
                               minLength=primerVar$minLength, maxLength=primerVar$maxLength, minCoverage=primerVar$minCoverage,
                               annealingTemp=primerVar$annealTempPrimer, minEfficiency=primerVar$minEfficiency, 
                               numPrimerSets=primerVar$numPrimerSets, minProductSize=primerVar$minProductSize,
                               maxProductSize=primerVar$maxProductSize)
      primers$setID <- row.names(primers)
      primers$length_FP <- nchar(primers$forward_primer)
      primers$length_RP <- nchar(primers$reverse_primer)
      primers$ampLength <- primers$start_reverse - primers$start_forward
      primers$annealTemp_FP <- primerVar$annealTempPrimer
      primers$annealTemp_RP <- primerVar$annealTempPrimer
      primers$forward_efficiency <- round(primers$forward_efficiency, 2)
      primers$reverse_efficiency <- round(primers$reverse_efficiency, 2)
      
      primers <- (primers[,c("setID", "identifier", "ampLength", "start_forward", "length_FP","start_reverse","length_RP",
                             "forward_primer", "reverse_primer", "forward_efficiency", "reverse_efficiency", "annealTemp_FP", 
                             "annealTemp_RP")])
      primers <- data.frame(primers, stringsAsFactors = FALSE)
      primers <- as.data.frame(lapply(primers, function(y) gsub(",", "", y)))
      primers <- primers[complete.cases(primers), ]
      
      dfPrimers$primerTab <- primers
      
      primerVar$primerPosStart <- dfPrimers$primerTab$start_forward
      primerVar$primerPosEnd <- dfPrimers$primerTab$start_reverse

      incProgress(0.6, detail = "Finished design of primer sets!")
      
      Sys.sleep(1)
    })
    
    # Primer Table using RHandsontable functionality 
    output$primerTable <- renderRHandsontable({
      rhandsontable(dfPrimers$primerTab) %>% hot_cols(columnSorting = TRUE, readOnly = FALSE, manualColumnResize = TRUE)
    })
    
  })
  
  # ***ANother section still in progress for primer binding visualization                         
  # observe({
  #   updateSelectInput(session, "inSelectPrimer", choices = as.character(primerVar$setID))
  # })
  
  # observeEvent(input$inSelectPrimer, {
  #  alignment <- subseq(reactUpload$alignment, primerVar$primerPosStart, primerVar$primerPosEnd)
  # })
  
  # Shiny Download handler for download of generated primer sets
  output$downloadPrimers <- downloadHandler(
    filename = function() { 
      paste("Primers", Sys.Date(), ".csv", sep="")
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
