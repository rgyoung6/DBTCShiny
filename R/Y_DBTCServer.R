# Written by Rob Young at the University of Guelph in Ontario Canada, Sept, 2023
#********************************************Main program section***********************************************
################# Server Function #############################################

#server <- function(input, output, session) {
shinyAppServer <- function(input, output, session) {

  ############## Initialize variables ############################################

  #dada implement reactive values
  primerFile <- reactiveValues(data = NA)
  dadaLocation <- reactiveValues(data = NA)

  #dada combine reactive values
  dadaCombineLoc <- reactiveValues(data = NA)

  #Make a BLASTable data base reactive values
  makeBlastDBFileLoc <- reactiveValues(data = NA)
  makeblastdbPath <- reactiveValues(data = NA)
  makeBlastTaxaDBLoc <- reactiveValues(data = NA)

  #Make a seqBLAST reactive values
  BLASTDatabasePath <- reactiveValues(data = NA)
  blastnPath <- reactiveValues(data = NA)
  querySeqPath <- reactiveValues(data = NA)

  # Taxon Assign
  taxaAssignFileLoc <- reactiveValues(data = NA)
  taxaAssignDBLoc <- reactiveValues(data = NA)

  #Combine Taxa Files
  combineTaxaFileLoc <- reactiveValues(data = NA)

  #Reduce Taxa Files
  reduceTaxaFileLoc <- reactiveValues(data = NA)

  #Combine Reduce Taxa Files
  combineReducedTaxaFileLoc <- reactiveValues(data = NA)

  #Other Variables
  ASVFile <- reactiveValues(data = NA)
  metaDataFile <- reactiveValues(data = NA)
  output$ASVFile_out <- shiny::renderText({as.character("No data file selected")})
  output$metaDataFile_out <- shiny::renderText({as.character("No data file selected")})
  ASVFileTable<-reactiveValues(data = NA)
  metaDataFileTable<-reactiveValues(data = NA)

  #Filtering variables
  markerColumns <- NULL #this is the number of columns for the combine markers output file
  maxReads <- reactiveValues(data = NA)
  mergedTable <- reactiveValues(data = NA)
  workingMergedTable <- reactiveValues(data = NA)
  firstLoadFlag <- reactiveValues(data = 0)

  #Dada files
  dadaDirectoryDisplayString<-reactiveValues(data = NA)
  primerFileDisplayString<-reactiveValues(data = NA)
  #Dada Combine
  dadaCombineFileDisplayString<-reactiveValues(data = NA)
  #Make DB
  makeBlastDBFileLocDisplayString<-reactiveValues(data = NA)
  makeblastdbPathDisplayString<-reactiveValues(data = NA)
  makeBlastTaxaDBLocDisplayString<-reactiveValues(data = NA)
  #BLAST
  BLASTDatabasePathDisplayString<-reactiveValues(data = NA)
  blastnPathDisplayString<-reactiveValues(data = NA)
  querySeqPathDisplayString<-reactiveValues(data = NA)
  #Taxon Assign
  taxaAssignFileLocDisplayString<-reactiveValues(data = NA)
  taxaAssignDBLocDisplayString<-reactiveValues(data = NA)
  #Combine Taxa Assign
  combineTaxaFileLocDisplayString<-reactiveValues(data = NA)
  #Reduce taxa assign
  reduceTaxaFileLocDisplayString<-reactiveValues(data = NA)
  #Combine Reduced Taxa assign
  combineReducedTaxaFileLocDisplayString<-reactiveValues(data = NA)

  ###################### Initialize the Map ######################################
  #Building the initial map and using the
  #Default static leaflet map before filtering parameters are applied.


  #Define the colour Palette for leaflet map
  palette_map <- leaflet::colorFactor(palette = c("#fbd300",
                                                  "#ff8600",
                                                  "#ea5f94",
                                                  "#9d02d7",
                                                  "#0000ff",
                                                  "#00008b"),
                                         domain = c("10", "100", "1000", "10000", "100000", "1000000"))

  output$mymap <- leaflet::renderLeaflet({
    leaflet::leaflet() %>%
      leaflet::addTiles() %>%
      leaflet::addLegend(position = "topright",
                pal = palette_map,
                values = c("10", "100", "1000", "10000", "100000", "1000000"),
                title = 'Reads',
                opacity = 0.6) %>%

      #View map full screen (note: only works in web browser)
      leaflet.extras::addFullscreenControl() %>%

      #Default map view --> Change to Guelph
      leaflet::setView(lng = -80.2262, lat = 43.5327, zoom = 3) %>%

      #Change leaflet map tiles
      leaflet::addProviderTiles(providers$Esri.WorldStreetMap)

  })

  ################## mappingDashboard Click observe event #######################
  shiny::observe({

    if (input$tab_being_displayed == "mappingDashboard"){
      if(is.na(ASVFileTable$data) && is.na(metaDataFileTable$data)){

        shiny::showModal(shiny::modalDialog(
          title = "No Data Loaded",
          "There are no data loaded in this instance of DBTCShine. Please go to the GPS/Grouping Import page and upload the data first."
        ))
      }else{

# Place an update here to get the data points on the map based on the information contained in the Mapped Data Table tab







        workMergedTable <- mergedTable$data
        workMergedTable <- workMergedTable[workMergedTable$Abundance >= input$abundance[1] & workMergedTable$Abundance <= input$abundance[2],]
        workMergedTable <- workMergedTable[workMergedTable$Final_Rank %in% input$finalRank,]
        workMergedTable <- workMergedTable[workMergedTable$Final_Taxa %in% input$finalTaxa,]
        workMergedTable <- workMergedTable[workMergedTable$superkingdom %in% input$kingdomFilterInput,]
        workMergedTable <- workMergedTable[workMergedTable$phylum %in% input$phylumFilterInput,]
        workMergedTable <- workMergedTable[workMergedTable$class %in% input$classFilterInput,]
        workMergedTable <- workMergedTable[workMergedTable$order %in% input$orderFilterInput,]
        workMergedTable <- workMergedTable[workMergedTable$family %in% input$familyFilterInput,]
        workMergedTable <- workMergedTable[workMergedTable$genus %in% input$genusFilterInput,]
        workMergedTable <- workMergedTable[workMergedTable$species %in% input$speciesFilterInput,]
        workMergedTable <- workMergedTable[workMergedTable$Sample %in% input$sampleFilterInput,]
        workMergedTable <- workMergedTable[workMergedTable$Run %in% input$runFilterInput,]
        workMergedTable <- workMergedTable[workMergedTable$Lab %in% input$labFilterInput,]
        workMergedTable <- workMergedTable[workMergedTable$Type %in% input$typeFilterInput,]
        workMergedTable <- workMergedTable[workMergedTable$Marker %in% input$markerFilterInput,]
        workMergedTable <- workMergedTable[workMergedTable$Date >= input$dateInput[1] & workMergedTable$Date <= input$dateInput[2],]

        leaflet::leafletProxy("mymap", data = as.data.frame(workMergedTable)) %>%
          clearMarkers() %>%
          clearMarkerClusters() %>%
          clearPopups() %>%
          #Adding labels to markers on map
          addCircleMarkers(lng = ~North,
                           lat = ~West,
                           color = ~palette_map(workMergedTable$AbundanceCategory),
                           clusterOptions = markerClusterOptions(spiderfyDistanceMultiplier=1.5),
                           popup= paste("<strong><h5>Species:", workMergedTable$Final_Taxa, "</strong>",
                                        "<br><h6>Sample:", workMergedTable$Sample,
                                        "<br><h6>Marker:", workMergedTable$Marker,
                                        "<br><h6>Date:", workMergedTable$Date,
                                        "<br><h6>Run:", workMergedTable$Run,
                                        "<br><h6>Type:", workMergedTable$Type,
                                        "<br><h6>Lab:", workMergedTable$Lab,
                                        "<h6>Coord(Lat, Lon):", workMergedTable$North,",", workMergedTable$West))

      }
    }
  })

  ################## Dada Submit Function #####################################
  # Get the path where all of the folders containing the fastq files are located
    volumes = shinyFiles::getVolumes()
    #Connect this to the shinyChooseButton
    shinyFiles::shinyFileChoose(input, "dadaDirectory", roots = volumes, session = session)

    # Get the file with the primer data for this analysis
    shiny::observeEvent(input$dadaDirectory, {

      if(!is.null(input$dadaDirectory)){
        tryCatch(
          expr = {

            dadaDirectory <- shinyFiles::parseFilePaths(volumes, input$dadaDirectory)
            dadaDirectoryDisplayString$data <- as.character(substr(dadaDirectory$datapath, 2, nchar(dadaDirectory$datapath)))
            output$dadaDirectoryDisplay <- shiny::renderText({as.character(dadaDirectoryDisplayString$data)})

         },
         error = function(e){
           print("Error - Dada Location Button choose file cancelled")
           dadaLocation$data <- NA
         },
         warning = function(w){
           print("Warning - Dada Location Button choose file cancelled")
           dadaLocation$data <- NA
         }
       )
      }
    },ignoreInit = TRUE)

    #Connect this to the shinyChooseButton
    shinyFiles::shinyFileChoose(input, "primerFile", roots = volumes, session = session)

    # Get the file with the primer data for this analysis
    shiny::observeEvent(input$primerFile, {
      if(!is.null(input$primerFile)){
        tryCatch(
          expr = {

            primerFile <- shinyFiles::parseFilePaths(volumes, input$primerFile)
            primerFileDisplayString$data <- as.character(substr(primerFile$datapath, 2, nchar(primerFile$datapath)))
            output$primerFileDisplay <- shiny::renderText({as.character(primerFileDisplayString$data)})

          },
          error = function(e){
            print("Error - Dada Location Button choose file cancelled")
            dadaLocation$data <- NA
          },
          warning = function(w){
            print("Warning - Dada Location Button choose file cancelled")
            dadaLocation$data <- NA
          }
        )
      }
    },ignoreInit = TRUE)

  shiny::observeEvent(input$dadaSubmit, {

#   if(!is.na(dadaDirectoryDisplayString$data) && !is.na(primerFileDisplayString$data)){
   if (!is.na(dadaDirectoryDisplayString$data) && is.character(dadaDirectoryDisplayString$data) && length(dadaDirectoryDisplayString$data) != 0 &&
       !is.na(primerFileDisplayString$data) && is.character(primerFileDisplayString$data) && length(primerFileDisplayString$data) != 0){

      # Create variables to call the dada_implement so that there are no conflicts
      # with the multithreading and the shiny app

      runFolderLoc <- force(dadaDirectoryDisplayString$data)
      primerFile <- force(primerFileDisplayString$data)

       if (force(input$uniOrbidirectional) == "Unidirectional"){
         unidirectional = TRUE
         bidirectional = FALSE
         fwdIdent <- ""
         revIdent <- ""
       }else if(force(input$uniOrbidirectional) == "Bidirectional"){
         unidirectional = FALSE
         bidirectional = TRUE
         fwdIdent <- force(input$fwdIdent)
         revIdent <- force(input$revIdent)
       }else{
         unidirectional = TRUE
         bidirectional = TRUE
         fwdIdent <- force(input$fwdIdent)
         revIdent <- force(input$revIdent)
       }
       printQualityPdf <- force(input$printQualityPdf)
       maxPrimeMis <- force(input$maxPrimeMis)
       fwdTrimLen <- force(input$fwdTrimLen)
       revTrimLen <- force(input$revTrimLen)
       maxEEVal <- force(input$maxEEVal)
       truncQValue <- force(input$truncQValue)
       truncLenValueF <- force(input$truncLenValueF)
       truncLenValueR <- force(input$truncLenValueR)
       error <- force(input$error)
       nbases <- force(input$nbases)
       maxMismatchValue <- force(input$maxMismatchValue)
       minOverlapValue <- force(input$minOverlapValue)
       trimOverhang <- force(input$trimOverhang)
       minFinalSeqLen <- force(input$minFinalSeqLen)

       tryCatch(
       expr = {
        #Run the Dada function here.

        shiny::showModal(shiny::modalDialog(
         title = "Dada analysis is underway.",
         "Processing, please stand by...", footer=""

        ))

        dada_implement(runFolderLoc = runFolderLoc,
                     primerFile = primerFile,
                     fwdIdent = fwdIdent,
                     revIdent = revIdent,
                     unidirectional = unidirectional,
                     bidirectional = bidirectional,
                     printQualityPdf = printQualityPdf,
                     maxPrimeMis = maxPrimeMis,
                     fwdTrimLen = fwdTrimLen,
                     revTrimLen = revTrimLen,
                     maxEEVal = maxEEVal,
                     truncQValue = truncQValue,
                     truncLenValueF = truncLenValueF,
                     truncLenValueR = truncLenValueR,
                     error = error,
                     nbases = nbases,
                     maxMismatchValue = maxMismatchValue,
                     minOverlapValue = minOverlapValue,
                     trimOverhang = trimOverhang,
                     minFinalSeqLen = minFinalSeqLen)
         removeModal()
         shiny::showModal(shiny::modalDialog(
           title = "Dada analysis is complete",
           "Please see output files in the target
           directory."
         ))
       },
       error = function(e){
         removeModal()
         shiny::showModal(shiny::modalDialog(
           title = "ERROR",
           "Dada Location Button choose file cancelled. Please refer to the R
           consol for more information."
         ))
         print("Error - Dada Location Button choose file cancelled")
         dadaLocation$data <- NA
       },
       warning = function(w){
         removeModal()
         shiny::showModal(shiny::modalDialog(
           title = "ERROR",
           "Dada Location Button choose file cancelled. Please refer to the R
           consol for more information."
         ))
         print("Warning - Dada Location Button choose file cancelled")
         dadaLocation$data <- NA
       })
    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please select a primer file and try submitting again!"
      ))
    }
   },ignoreInit = TRUE)

  ################## Dada Combine Function ####################################
  #Get the location of the dada output files you would like to combine

    # Get the path where all of the folders containing the fastq files are located
    volumes = shinyFiles::getVolumes()
    #Connect this to the shinyChooseButton
    shinyFiles::shinyFileChoose(input, "dadaCombineFile", roots = volumes, session = session)

    # Get the file with the primer data for this analysis
    shiny::observeEvent(input$dadaCombineFile, {

      if(!is.null(input$dadaCombineFile)){

        tryCatch(
          expr = {
            dadaCombineFile <- shinyFiles::parseFilePaths(volumes, input$dadaCombineFile)
            dadaCombineFileDisplayString$data <- as.character(substr(dadaCombineFile$datapath, 2, nchar(dadaCombineFile$datapath)))
            output$dadaCombineDisplay <- shiny::renderText({as.character(dadaCombineFileDisplayString$data)})
         },
         error = function(e){
           print("Error - Dada Location Button choose file cancelled")
           dadaLocation$data <- NA
         },
         warning = function(w){
           print("Warning - Dada Location Button choose file cancelled")
           dadaLocation$data <- NA
         }
       )
      }
    })

  #Running the data combine
  shiny::observeEvent(input$dadaCombine, {

#    if(!is.na(dadaCombineFileDisplayString$data)){
    if (!is.na(dadaCombineFileDisplayString$data) && is.character(dadaCombineFileDisplayString$data) && length(dadaCombineFileDisplayString$data) != 0){
      # Create variables for the arguments to avoid conflicts between the multithreading
      # and the shiny
      fileLoc= force(dadaCombineFileDisplayString$data)
      minLen = force(input$dadaCombineMinLen)

      tryCatch(
        expr = {
          #Run the Dada function here.
          shiny::showModal(shiny::modalDialog(
            title = "Dada combine analysis results is underway.",
            "Processing, please stand by...", footer=""
          ))

        #Run the Dada combine function here.
          combine_dada_output(fileLoc = fileLoc, minLen = minLen)

        removeModal()
        shiny::showModal(shiny::modalDialog(
          title = "Dada combine analysis results is complete",
          "Please see output files in the target directory."
        ))
        },
        error = function(e){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        },
        warning = function(w){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        }
      )

    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please select the file location where all of the DBTC dada output files are located that you wish to combine!"
      ))
    }
  },ignoreInit = TRUE)

  ################## Make BLAST DB Function ###################################

  # Get the path where all of the folders containing the fastq files are located
  volumes = shinyFiles::getVolumes()
  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "makeBlastDBFileLoc", roots = volumes, session = session)

  # Get the fasta file you want to use to build your db
  shiny::observeEvent(input$makeBlastDBFileLoc, {
    tryCatch(
      expr = {

        makeBlastDBFileLoc <- shinyFiles::parseFilePaths(volumes, input$makeBlastDBFileLoc)
        makeBlastDBFileLocDisplayString$data <- as.character(substr(makeBlastDBFileLoc$datapath, 2, nchar(makeBlastDBFileLoc$datapath)))
        output$makeBlastDBFileLocDisplay <- shiny::renderText({as.character(makeBlastDBFileLocDisplayString$data)})

      },
      error = function(e){
        print("Error")
        makeBlastDBFileLoc$data <- NA
      },
      warning = function(w){
        print("Warning")
        makeBlastDBFileLoc$data <- NA
      }
    )
  },ignoreInit = TRUE)

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "makeblastdbPath", roots = volumes, session = session)

  # Select where the makeblastdb program is on your computer
  shiny::observeEvent(input$makeblastdbPath, {
    tryCatch(
      expr = {

        makeblastdbPath <- shinyFiles::parseFilePaths(volumes, input$makeblastdbPath)
        makeblastdbPathDisplayString$data <- as.character(substr(makeblastdbPath$datapath, 2, nchar(makeblastdbPath$datapath)))
        output$makeblastdbPathDisplay <- shiny::renderText({as.character(makeblastdbPathDisplayString$data)})

      },
      error = function(e){
        print("Error")
        makeblastdbPath$data <- NA
      },
      warning = function(w){
        print("Warning")
        makeblastdbPath$data <- NA
      }
    )
  },ignoreInit = TRUE)

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "makeBlastTaxaDBLoc", roots = volumes, session = session)

  # Point to the NCBI taxonomic data base
  shiny::observeEvent(input$makeBlastTaxaDBLoc, {
    tryCatch(
      expr = {

        makeBlastTaxaDBLoc <- shinyFiles::parseFilePaths(volumes, input$makeBlastTaxaDBLoc)
        makeBlastTaxaDBLocDisplayString$data <- as.character(substr(makeBlastTaxaDBLoc$datapath, 2, nchar(makeBlastTaxaDBLoc$datapath)))
        output$makeBlastTaxaDBLocDisplay <- shiny::renderText({as.character(makeBlastTaxaDBLocDisplayString$data)})

      },
      error = function(e){
        print("Error")
        makeBlastTaxaDBLoc$data <- NA
      },
      warning = function(w){
        print("Warning")
        makeBlastTaxaDBLoc$data <- NA
      }
    )
  },ignoreInit = TRUE)

  # Run the making BLAST db code.
  shiny::observeEvent(input$makeBlastDB, {
    if(is.na(makeblastdbPathDisplayString$data)){
      makeblastdbPathDisplayString$data <- "makeblastdb"
    }
    # if(!is.na(makeBlastDBFileLocDisplayString$data) &&
    #    !is.na(makeblastdbPathDisplayString$data) &&
    #    !is.na(makeBlastTaxaDBLocDisplayString$data) &&
    #    !is.na(input$dbName) &&
    #    !is.na(input$makeBLASTDBMinLen)){

print("makeBlastDBFileLocDisplayString$data")
print(makeBlastDBFileLocDisplayString$data)
print("makeblastdbPathDisplayString$data")
print(makeblastdbPathDisplayString$data)
print("makeBlastTaxaDBLocDisplayString$data")
print(makeBlastTaxaDBLocDisplayString$data)


    if (!is.na(makeBlastDBFileLocDisplayString$data) && is.character(makeBlastDBFileLocDisplayString$data) && length(makeBlastDBFileLocDisplayString$data) != 0 &&
      !is.na(makeBlastTaxaDBLocDisplayString$data) && is.character(makeBlastTaxaDBLocDisplayString$data) && length(makeBlastTaxaDBLocDisplayString$data)!= 0 ){

print("In the if checking the file locations")

      # Create local variables to avoid conflicts with shiny and multithread
      fileLoc = force(makeBlastDBFileLocDisplayString$data)

      if(is.character(makeblastdbPathDisplayString$data) && length(makeblastdbPathDisplayString$data) != 0){

        makeblastdbPath = force(makeblastdbPathDisplayString$data)

      } else {

        makeblastdbPath = "makeblastdb"

      }

      taxaDBLoc = force(makeBlastTaxaDBLocDisplayString$data)
      dbName = force(input$dbName)
      minLen = force(input$makeBLASTDBMinLen)

      tryCatch(
        expr = {
          #Run the Dada function here.
          shiny::showModal(shiny::modalDialog(
            title = "Make BLAST database is underway.",
            "Processing, please stand by...", footer=""

          ))
          #Run the function
          make_BLAST_DB(fileLoc = fileLoc,
                        makeblastdbPath = makeblastdbPath,
                        taxaDBLoc = taxaDBLoc,
                        dbName = dbName,
                        minLen = minLen)
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "Make BLAST database is complete",
            "Please see output files in the target
           directory."
          ))
        },
        error = function(e){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "There was an error in running the makeBLASTDB function. Make sure you have permissions for the target folder. Please see the R output for further details."
          ))
        },
        warning = function(w){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "There was an error in running the makeBLASTDB function. Make sure you have permissions for the target folder. Please see the R output for further details."
          ))
        })

    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please fill in all of the necessary fields and submit again!"
      ))
    }
  },ignoreInit = TRUE)


  ################## BLAST sequences Function #################################

  # Get the path where all of the folders containing the fastq files are located
  volumes = shinyFiles::getVolumes()
  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "BLASTDatabasePath", roots = volumes, session = session)

  # Point to the NCBI taxonomic data base
  shiny::observeEvent(input$BLASTDatabasePath, {
    tryCatch(
      expr = {
        BLASTDatabasePath <- shinyFiles::parseFilePaths(volumes, input$BLASTDatabasePath)
        BLASTDatabasePathDisplayString$data <- as.character(substr(BLASTDatabasePath$datapath, 2, nchar(BLASTDatabasePath$datapath)))
        output$BLASTDatabasePathDisplay <- shiny::renderText({as.character(BLASTDatabasePathDisplayString$data)})
      },
      error = function(e){
        print("Error")
        BLASTDatabasePath$data <- NA
      },
      warning = function(w){
        print("Warning")
        BLASTDatabasePath$data <- NA
      }
    )
  },ignoreInit = TRUE)

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "blastnPath", roots = volumes, session = session)

  # Get the path where all of the folders containing the fastq files are located
  shiny::observeEvent(input$blastnPath, {
    tryCatch(
      expr = {
        blastnPath <- shinyFiles::parseFilePaths(volumes, input$blastnPath)
        blastnPathDisplayString$data <- as.character(substr(blastnPath$datapath, 2, nchar(blastnPath$datapath)))
        output$blastnPathDisplay <- shiny::renderText({as.character(blastnPathDisplayString$data)})
      },
      error = function(e){
        print("Error")
        blastnPath$data <- NA
      },
      warning = function(w){
        print("Warning")
        blastnPath$data <- NA
      }
    )
  },ignoreInit = TRUE)

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "querySeqPath", roots = volumes, session = session)

  # Point to the location of the fasta files you want to BLAST
  shiny::observeEvent(input$querySeqPath, {
    tryCatch(
      expr = {
        querySeqPath <- shinyFiles::parseFilePaths(volumes, input$querySeqPath)
        querySeqPathDisplayString$data <- as.character(substr(querySeqPath$datapath, 2, nchar(querySeqPath$datapath)))
        output$querySeqPathDisplay <- shiny::renderText({as.character(querySeqPathDisplayString$data)})
      },
      error = function(e){
        print("Error")
        querySeqPath$data <- NA
      },
      warning = function(w){
        print("Warning")
        querySeqPath$data <- NA
      }
    )
  },ignoreInit = TRUE)

  shiny::observeEvent(input$blastSequences, {

    if(is.na(blastnPathDisplayString$data) | is.null(blastnPathDisplayString$data)){
      blastnPathDisplayString$data <- "blastn"
    }

    print(paste0("BLASTDatabasePathDisplayString$data - ",BLASTDatabasePathDisplayString$data))
    print(paste0("blastnPathDisplayString$data - ", blastnPathDisplayString$data))
    print(paste0("querySeqPathDisplayString$data - ", querySeqPathDisplayString$data))
    print(paste0("input$BLASTminLen - ", input$BLASTminLen))
    print(paste0("input$BLASTResults - ", input$BLASTResults))
    print(paste0("input$blastSeqNumCores - ", input$blastSeqNumCores))

    if (!is.na(BLASTDatabasePathDisplayString$data) && is.character(BLASTDatabasePathDisplayString$data) && length(BLASTDatabasePathDisplayString$data) != 0 &&
        !is.na(blastnPathDisplayString$data) && is.character(blastnPathDisplayString$data) && length(blastnPathDisplayString$data) != 0 &&
        !is.na(querySeqPathDisplayString$data) && is.character(querySeqPathDisplayString$data) && length(querySeqPathDisplayString$data) != 0) {

      # Create local variables to avoid conflicts with shiny and multithread
      databasePath = force(BLASTDatabasePathDisplayString$data)
      blastnPath = force(blastnPathDisplayString$data)
      querySeqPath = force(querySeqPathDisplayString$data)
      minLen = force(input$BLASTminLen)
      BLASTResults = force(input$BLASTResults)
      numCores = force(input$blastSeqNumCores)

      tryCatch(
        expr = {
          #Run the Dada function here.

          shiny::showModal(shiny::modalDialog(
            title = "Sequence BLAST is underway.",
            "Processing, please stand by...", footer=""

          ))

          #Run the function
          seq_BLAST(databasePath = databasePath,
                    blastnPath = blastnPath,
                    querySeqPath = querySeqPath,
                    minLen = minLen,
                    BLASTResults = BLASTResults,
                    numCores = numCores)

          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "Sequence BLAST is complete",
            "Please see output files in the target directory."
          ))
        },
        error = function(e){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        },
        warning = function(w){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        }
      )
    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please fill in all of the necessary fields and submit again!"
      ))
    }
  },ignoreInit = TRUE)

  ################## Taxon Assign Function ####################################

  # Get the path where all of the folders containing the fastq files are located
  volumes = shinyFiles::getVolumes()

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "taxaAssignFileLoc", roots = volumes, session = session)

  # Point to the BLAST output files to get taxonomic assignment
  shiny::observeEvent(input$taxaAssignFileLoc, {
    tryCatch(
      expr = {
        taxaAssignFileLoc <- shinyFiles::parseFilePaths(volumes, input$taxaAssignFileLoc)
        taxaAssignFileLocDisplayString$data <- as.character(substr(taxaAssignFileLoc$datapath, 2, nchar(taxaAssignFileLoc$datapath)))
        output$taxaAssignFileLocDisplay <- shiny::renderText({as.character(taxaAssignFileLocDisplayString$data)})
      },
      error = function(e){
        print("Taxa Assign Error")
        taxaAssignFileLoc$data <- NA
      },
      warning = function(w){
        print("Taxa Assign Warning")
        taxaAssignFileLoc$data <- NA
      }
    )
  },ignoreInit = TRUE)

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "taxaAssignDBLoc", roots = volumes, session = session)

  # Point to the NCBI taxonomic data base
  shiny::observeEvent(input$taxaAssignDBLoc, {
    tryCatch(
      expr = {

        taxaAssignDBLoc <- shinyFiles::parseFilePaths(volumes, input$taxaAssignDBLoc)
        taxaAssignDBLocDisplayString$data <- as.character(substr(taxaAssignDBLoc$datapath, 2, nchar(taxaAssignDBLoc$datapath)))
        output$taxaAssignDBLocDisplay <- shiny::renderText({as.character(taxaAssignDBLocDisplayString$data)})

      },
      error = function(e){
        print("Error")
        taxaAssignDBLoc$data <- NA
      },
      warning = function(w){
        print("Warning")
        taxaAssignDBLoc$data <- NA
      }
    )
  },ignoreInit = TRUE)

  shiny::observeEvent(input$taxonAssign, {
    if (!is.na(taxaAssignFileLocDisplayString$data) && is.character(taxaAssignFileLocDisplayString$data) && length(taxaAssignFileLocDisplayString$data) != 0 && !is.na(taxaAssignDBLocDisplayString$data) && is.character(taxaAssignDBLocDisplayString$data) && length(taxaAssignDBLocDisplayString$data) != 0) {

      # Create local variables to avoid conflicts with shiny and multithread
       fileLoc = force(taxaAssignFileLocDisplayString$data)

       taxaDBLoc = force(taxaAssignDBLocDisplayString$data)
       numCores = force(input$taxaAssignNumCores)
       coverage = force(input$coverage)
       ident = force(input$ident)
       propThres = force(input$propThres)
       coverReportThresh = force(input$coverReportThresh)
       identReportThresh = force(input$identReportThresh)
       includeAllDada = force(input$includeAllDada)

       tryCatch(
         expr = {
           #Run the function here.
           shiny::showModal(shiny::modalDialog(
             title = "Taxonomic assingment is underway.",
             "See the R terminal for estimated time to completion.
             Processing, please stand by...", footer=""

           ))
           #Run the function
           taxon_assign(fileLoc = fileLoc,
                        taxaDBLoc = taxaDBLoc,
                        numCores = numCores,
                        coverage = coverage,
                        ident = ident,
                        propThres = propThres,
                        coverReportThresh = coverReportThresh,
                        identReportThresh = identReportThresh,
                        includeAllDada = includeAllDada)
            removeModal()

           shiny::showModal(shiny::modalDialog(
             title = "Taxonomic assingment is complete",
             "Please see output files in the target directory."
           ))
         },
         error = function(e){
           removeModal()
           shiny::showModal(shiny::modalDialog(
             title = "ERROR",
             "Please refer to the R consol for more information."
           ))
         },
         warning = function(w){
           removeModal()
           shiny::showModal(shiny::modalDialog(
             title = "ERROR",
             "Please refer to the R consol for more information."
           ))
         }
       )

    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please fill in all of the necessary fields and submit again!"
      ))
    }
  })

  ################## Combine Taxa Assign Function #############################

  # Get the path where all of the folders containing the fastq files are located
  volumes = shinyFiles::getVolumes()

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "combineTaxaFileLoc", roots = volumes, session = session)

  # Point to the folder with the Taxon Assign records you would like to combine
  shiny::observeEvent(input$combineTaxaFileLoc, {
    tryCatch(
      expr = {

        combineTaxaFileLoc <- shinyFiles::parseFilePaths(volumes, input$combineTaxaFileLoc)
        combineTaxaFileLocDisplayString$data <- as.character(substr(combineTaxaFileLoc$datapath, 2, nchar(combineTaxaFileLoc$datapath)))
        output$combineTaxaFileLocDisplay <- shiny::renderText({as.character(combineTaxaFileLocDisplayString$data)})

      },
      error = function(e){
        print("Error")
        combineTaxaFileLoc$data <- NA
      },
      warning = function(w){
        print("Warning")
        combineTaxaFileLoc$data <- NA
      }
    )
  },ignoreInit = TRUE)

  shiny::observeEvent(input$combineTaxa, {
    if (!is.na(combineTaxaFileLocDisplayString$data) && is.character(combineTaxaFileLocDisplayString$data) && length(combineTaxaFileLocDisplayString$data) != 0) {
      # Create local variables to avoid conflicts with shiny and multithread
      fileLoc = force(combineTaxaFileLocDisplayString$data)
      numCores = force(input$combineTaxaNumCores)

print(paste0("Here is the value of the combine_assign_output fileLoc...", fileLoc))
print(paste0("Here is the value of the combine_assign_output numCores...", numCores))

      tryCatch(
        expr = {
          #Run the function here.
          shiny::showModal(shiny::modalDialog(
            title = "Combining taxa assingment files is underway.",
            "Processing, please stand by...", footer=""

          ))
          #Run the function
          combine_assign_output(fileLoc = fileLoc,
                                numCores = numCores)
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "Combining taxa assingment files is complete",
            "Please see output files in the target directory."
          ))
        },
        error = function(e){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        },
        warning = function(w){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        }
      )
    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please fill in all of the necessary fields and submit again!"
      ))
    }
  })

  ################## Reduce Taxa Assign Function ##############################

  # Get the path where all of the folders containing the fastq files are located
  volumes = shinyFiles::getVolumes()

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "reduceTaxaFileLoc", roots = volumes, session = session)

  # Point to the folder with the Taxon Assign records you would like to combine
  shiny::observeEvent(input$reduceTaxaFileLoc, {
    tryCatch(
      expr = {

        reduceTaxaFileLoc <- shinyFiles::parseFilePaths(volumes, input$reduceTaxaFileLoc)
        reduceTaxaFileLocDisplayString$data <- as.character(substr(reduceTaxaFileLoc$datapath, 2, nchar(reduceTaxaFileLoc$datapath)))
        output$reduceTaxaFileLocDisplay <- shiny::renderText({as.character(reduceTaxaFileLocDisplayString$data)})

      },
      error = function(e){
        print("Error")
        reduceTaxaFileLocDisplayString$data <- NA
      },
      warning = function(w){
        print("Warning")
        reduceTaxaFileLocDisplayString$data <- NA
      }
    )
  },ignoreInit = TRUE)

  shiny::observeEvent(input$reduceTaxa, {
    if (!is.na(reduceTaxaFileLocDisplayString$data) && is.character(reduceTaxaFileLocDisplayString$data) && length(reduceTaxaFileLocDisplayString$data) != 0) {
      # Create local variables to avoid conflicts with shiny and multithread
      fileLoc = force(reduceTaxaFileLocDisplayString$data)
      numCores = force(input$reduceTaxaNumCores)

print(paste0("Here is the value of the reduceTaxa fileLoc...", fileLoc))
print(paste0("Here is the value of the reduceTaxa numCores...", numCores))

      tryCatch(
        expr = {
          #Run the function here.

          shiny::showModal(shiny::modalDialog(
            title = "Reduce ASV results with taxonomic assignment to unique taxa is underway.",
            "Processing, please stand by...", footer=""

          ))

          #Run the function
          reduce_taxa(fileLoc = fileLoc,
                      numCores = numCores)

          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "Reduce ASV results with taxonomic assignment to unique taxa is complete",
            "Please see output files in the target directory."
          ))
        },
        error = function(e){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        },
        warning = function(w){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        }
      )
    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please fill in all of the necessary fields and submit again!"
      ))
    }
  })

  ################## Combine Reduce Taxa Assign Function ######################

  # Get the path where all of the folders containing the fastq files are located
  volumes = shinyFiles::getVolumes()

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "combineReducedTaxaFileLoc", roots = volumes, session = session)

  # Point to the folder with the Taxon Assign records you would like to combine
  shiny::observeEvent(input$combineReducedTaxaFileLoc, {
    tryCatch(
      expr = {

        combineReducedTaxaFileLoc <- shinyFiles::parseFilePaths(volumes, input$combineReducedTaxaFileLoc)
        combineReducedTaxaFileLocDisplayString$data <- as.character(substr(combineReducedTaxaFileLoc$datapath, 2, nchar(combineReducedTaxaFileLoc$datapath)))
        output$combineReducedTaxaFileLocDisplay <- shiny::renderText({as.character(combineReducedTaxaFileLocDisplayString$data)})

      },
      error = function(e){
        print("Error")
        combineReducedTaxaFileLocDisplayString$data <- NA
      },
      warning = function(w){
        print("Warning")
        combineReducedTaxaFileLocDisplayString$data <- NA
      }
    )
  },ignoreInit = TRUE)

  shiny::observeEvent(input$combineReduceTaxa, {

     if (!is.na(combineReducedTaxaFileLocDisplayString$data) && is.character(combineReducedTaxaFileLocDisplayString$data) && length(combineReducedTaxaFileLocDisplayString$data) != 0) {
      # Create local variables to avoid conflicts with shiny and multithread
      fileLoc = force(combineReducedTaxaFileLocDisplayString$data)
      presenceAbsence = force(input$presenceAbsence)

      tryCatch(
        expr = {
          #Run the function here.

          shiny::showModal(shiny::modalDialog(
            title = "Combine reduced taxonomic results for multiple markers for the same samples is underway.",
            "Processing, please stand by...", footer=""

          ))

          #Run the function
          combine_reduced_output(fileLoc = fileLoc, presenceAbsence = presenceAbsence)

          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "Combine reduced taxonomic results for multiple markers for the same samples is complete",
            "Please see output files in the target directory."
          ))
        },
        error = function(e){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        },
        warning = function(w){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        }
      )
    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please fill in all of the necessary fields and submit again!"
      ))
    }
  })





















  ################## GPS/Grouping Import File Function Buttons ###################

  # Show modal when button is clicked.
  shiny::observeEvent(input$ASVFile_button, {
    tryCatch(
      expr = {
        ASVFile$data <- file.choose()
        output$ASVFile_out <- shiny::renderText({as.character(ASVFile$data)})
      },
      error = function(e){
        print("Error - Data Input Submit choose file cancelled - 1")
      },
      warning = function(w){
        print("Warning - Data Input Submit choose file cancelled - 2")
      }
    )
  },ignoreInit = TRUE)

  # Show modal when button is clicked.
  shiny::observeEvent(input$MetadataFile_button, {

    tryCatch(
      expr = {
        metaDataFile$data <- file.choose()
        output$metaDataFile_out <- shiny::renderText({as.character(metaDataFile$data)})
      },
      error = function(e){
        print("Error - Data Input Submit choose file cancelled - 1")
      },
      warning = function(w){
        print("Warning - Data Input Submit choose file cancelled - 2")
      }
    )
  },ignoreInit = TRUE)

  ################## GPS/Grouping Reset Data Import Function #####################

  #Data reset button
  shiny::observeEvent(input$resetDataImport, {
    output$ASVFile_out <- shiny::renderText({as.character("No data file selected")})
    output$metaDataFile_out <- shiny::renderText({as.character("No data file selected")})

    ASVFileTable$data <- NA
    metaDataFileTable$data <- NA
    mergedTable$data <- NA
    markerColumns <- NULL
    firstLoadFlag$data = 0

    shiny::showModal(shiny::modalDialog(
      title = "Cleared Data",
      "The data selection has been cleared from this instance!"
    ))

  },ignoreInit = TRUE)

  ################## GPS/Grouping Data Processing Upon Hitting Submit ################

  shiny::observeEvent(input$submitDataImport, {

    #The below content is checking the files submitted are .tsv
    if(is.null(ASVFile$data)){
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Either the ASV file or the metadata file are missing. Please select both files and resubmit - 1"
      ))
    }else if(is.null(metaDataFile$data)){
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Either the ASV file or the metadata file are missing. Please select both files and resubmit - 2"
      ))
    }else if(is.na(ASVFile$data)){
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Either the ASV file or the metadata file are missing. Please select both files and resubmit - 3"
      ))
    }else if(is.na(metaDataFile$data)){
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Either the ASV file or the metadata file are missing. Please select both files and resubmit - 4"
      ))
    }else if(grepl("\\.[Tt][Ss][Vv]$", ASVFile$data)){

      if(grepl("\\.[Tt][Ss][Vv]$", metaDataFile$data)){

        #Loading in the data files
        tryCatch(
          expr = {

            shiny::showModal(shiny::modalDialog("Loading the data files, please wait...",  footer=NULL))

            #Load in the data to the formatted_metadata variable
            ASVFileTable$data<-read.table(ASVFile$data, header = TRUE, check.names=FALSE, sep="\t", dec=".")
print("Here loading data 1")
            metaDataFileTable$data<-read.table(metaDataFile$data, header = TRUE, check.names=FALSE, sep="\t", dec=".")
print("Here loading data 2")
            ########################Process the submitted data files###################

            #Get the columns with the '_MarkerResults' string which indicates that it was a combined marker file
            # from the combined reduced taxa assign function
            markerColumns <- grep("_MarkerResults", names(ASVFileTable$data), value = TRUE)
print(paste0("Here loading data 3 and the columns with _MarkerResults are...", markerColumns))

            if(length(markerColumns)>0){
print("Here loading data 4")
              #Get the last column number with the '_MarkerResults' in the title
              maxColNum <- max(which(names(ASVFileTable$data) %in% markerColumns))
print("Here loading data 5")
              #Get the names of the samples
              sampleNames <- names(ASVFileTable$data)[c((maxColNum+1):ncol(ASVFileTable$data))]
print("Here loading data 6")
              #Flatten the data.frame
              flatTable <- reshape(
                ASVFileTable$data,
                varying = list(sampleNames),
                v.names = "Abundance",
                direction = "long",
                times = sampleNames,
                timevar = "Sample",
                sep = ""
              )
print("Here loading data 7")
              #Remove all entries with 0 in the Abundance column
              flatTable <- flatTable[flatTable$Abundance != 0, ]
print("Here loading data 8")
              #Merge the flattened file with the GPS file.
              mergedTable$data <- merge(flatTable, metaDataFileTable$data, by = "Sample")
print("Here loading data 9")
              #Initialize the finalMergedTable with the first marker.
              finalMergedTable <- NULL
print("Here loading data 10")
              #Loop through the markers to add them on to the finalMergedTable
              for (numMarkerCol in 1:length(markerColumns)){

                #Initialize the tempMergedTable
                tempMergedTable <- mergedTable$data[, -which(names(mergedTable$data) %in% markerColumns)]

                #Add a column for the marker name
                tempMergedTable$Marker <- gsub("_MarkerResults", "", markerColumns[numMarkerCol])

                #Create a temporary mergedTable variable
                tempMergedTable <-cbind(tempMergedTable, mergedTable$data[,markerColumns[numMarkerCol]])

                #Rename the last column name to Final_Taxa
                names(tempMergedTable)[ncol(tempMergedTable)] <- "Final_Taxa"

                #Remove records based on NA in the Marker column
                tempMergedTable <- subset(tempMergedTable, !is.na(Final_Taxa))

                #Split the final column by ' - '
                tempMergedTable <- cbind(tempMergedTable[,-ncol(tempMergedTable)],do.call(rbind, strsplit(as.character(tempMergedTable$Final_Taxa), " - ", fixed = TRUE)))

                #Rename the last two columns to Final_Rank and Final_Taxa
                names(tempMergedTable)[(ncol(tempMergedTable)-1)] <- "Final_Rank"
                names(tempMergedTable)[ncol(tempMergedTable)] <- "Final_Taxa"

                #Add the data back onto a dataframe with each molecular Marker data in a marker column
                finalMergedTable <- rbind(finalMergedTable,tempMergedTable)

print(paste0("Here loading data 11 - ", numMarkerCol))
              }

              #Add a category column to place the data points on the map
              finalMergedTable$AbundanceCategory <- cut(
                finalMergedTable$Abundance,
                breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000),
                labels = c("10", "100", "1000", "10000", "100000", "1000000"),
                include.lowest = TRUE
              )
print("Here loading data 12")
              #Remove unnecessary variables
              mergedTable$data <- finalMergedTable
              remove(finalMergedTable)
              remove(tempMergedTable)
              remove(flatTable)

            } else{

              #Get the last column number with the 'Results' in the title
              maxColNum <- max(which(names(ASVFileTable$data) %in% "Results"))

              #Get the names of the samples
              sampleNames <- names(ASVFileTable$data)[c((maxColNum+1):ncol(ASVFileTable$data))]

              #Flatten the data.frame
              flatTable <- reshape(
                ASVFileTable$data,
                varying = list(sampleNames),
                v.names = "Abundance",
                direction = "long",
                times = sampleNames,
                timevar = "Sample",
                sep = ""
              )

              #Remove all entries with 0 in the Abundance column
              flatTable <- flatTable[flatTable$Abundance != 0, ]

              #Merge the flattened file with the GPS file.
              tempMergedTable <- merge(flatTable, metaDataFileTable$data, by = "Sample")

              # Remove brackets and contents from all values in the data frame
              tempMergedTable[,c(which(names(tempMergedTable) %in% "superkingdom"):which(names(tempMergedTable) %in% "species"))] <- lapply(tempMergedTable[,c(which(names(tempMergedTable) %in% "superkingdom"):which(names(tempMergedTable) %in% "species"))], function(x) gsub("\\(.*?\\)", "", x))

              #Add a column for the marker name
              tempMergedTable$Marker <- gsub("_MarkerResults", "", markerColumns[numMarkerCol])

              #Add a category column to place the data points on the map
              tempMergedTable$AbundanceCategory <- cut(
                tempMergedTable$Abundance,
                breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000),
                labels = c("10", "100", "1000", "10000", "100000", "1000000"),
                include.lowest = TRUE
              )

              #Remove unnecessary variables
              mergedTable$data <- tempMergedTable
              remove(flatTable)
              remove(tempMergedTable)

            }
print("Here loading data 13")

            #################### update the filters based on submitted ASV data ############

            shiny::updateSliderInput(session = session, inputId = "abundance",
                                     min=0,
                                     max = max(mergedTable$data$Abundance),
                                     value=c(0,max(mergedTable$data$Abundance))
            )

            #Dropdown menu for Final_Rank
            shinyWidgets::updatePickerInput(session, inputId = "finalRank",
                                            choices = sort(unique(mergedTable$data$Final_Rank), na.last = TRUE),
                                            selected = sort(unique(mergedTable$data$Final_Rank), na.last = TRUE))

            #Dropdown menu for superkingdom
            shinyWidgets::updatePickerInput(session, inputId = "kingdomFilterInput",
                                            choices = sort(unique(mergedTable$data$superkingdom), na.last = TRUE),
                                            selected = sort(unique(mergedTable$data$superkingdom), na.last = TRUE))

            #Dropdown menu for phylum
            shinyWidgets::updatePickerInput(session, inputId = "phylumFilterInput",
                                            choices = sort(unique(mergedTable$data$phylum), na.last = TRUE),
                                            selected = sort(unique(mergedTable$data$phylum), na.last = TRUE))

            #Dropdown menu for class
            shinyWidgets::updatePickerInput(session, inputId = "classFilterInput",
                                            choices = sort(unique(mergedTable$data$class), na.last = TRUE),
                                            selected = sort(unique(mergedTable$data$class), na.last = TRUE))

            #Dropdown menu for order
            shinyWidgets::updatePickerInput(session, inputId = "orderFilterInput",
                                            choices = sort(unique(mergedTable$data$order), na.last = TRUE),
                                            selected = sort(unique(mergedTable$data$order), na.last = TRUE))

            #Dropdown menu for family
            shinyWidgets::updatePickerInput(session, inputId = "familyFilterInput",
                                            choices = sort(unique(mergedTable$data$family), na.last = TRUE),
                                            selected = sort(unique(mergedTable$data$family), na.last = TRUE))

            #Dropdown menu for genus
            shinyWidgets::updatePickerInput(session, inputId = "genusFilterInput",
                                            choices = sort(unique(mergedTable$data$genus), na.last = TRUE),
                                            selected = sort(unique(mergedTable$data$genus), na.last = TRUE))

            #Dropdown menu for species
            shinyWidgets::updatePickerInput(session, inputId = "speciesFilterInput",
                                            choices = sort(unique(mergedTable$data$species), na.last = TRUE),
                                            selected = sort(unique(mergedTable$data$species), na.last = TRUE))

            #################### update the filters based on submitted Meta Data ###########

            #Dropdown menu for Sample
            shinyWidgets::updatePickerInput(session, inputId = "sampleFilterInput",
                                            choices = sort(unique(mergedTable$data$Sample), na.last = TRUE),
                                            selected = sort(unique(mergedTable$data$Sample), na.last = TRUE))

            #Dropdown menu for Run
            shinyWidgets::updatePickerInput(session, inputId = "runFilterInput",
                                            choices = sort(unique(mergedTable$data$Run), na.last = TRUE),
                                            selected = sort(unique(mergedTable$data$Run), na.last = TRUE))

            #Dropdown menu for Lab
            shinyWidgets::updatePickerInput(session, inputId = "labFilterInput",
                                            choices = sort(unique(mergedTable$data$Lab), na.last = TRUE),
                                            selected = sort(unique(mergedTable$data$Lab), na.last = TRUE))

            #Dropdown menu for Type
            shinyWidgets::updatePickerInput(session, inputId = "typeFilterInput",
                                            choices = sort(unique(mergedTable$data$Type), na.last = TRUE),
                                            selected = sort(unique(mergedTable$data$Type), na.last = TRUE))

            #Dropdown menu for Region
            shinyWidgets::updatePickerInput(session,inputId = "markerFilterInput",
                                            choices = sort(unique(mergedTable$data$Marker), na.last = TRUE),
                                            selected = sort(unique(mergedTable$data$Marker), na.last = TRUE))

            shiny::updateSliderInput(session = session, inputId = "dateInput",
                                     min = as.Date(min(as.Date(mergedTable$data$Date, "%Y-%m-%d")[!is.na(as.Date(mergedTable$data$Date, "%Y-%m-%d"))]),"%Y-%m-%d"),
                                     max = as.Date(max(as.Date(mergedTable$data$Date, "%Y-%m-%d")[!is.na(as.Date(mergedTable$data$Date, "%Y-%m-%d"))]),"%Y-%m-%d"),
                                     value=c(as.Date(min(as.Date(mergedTable$data$Date, "%Y-%m-%d")[!is.na(as.Date(mergedTable$data$Date, "%Y-%m-%d"))]),"%Y-%m-%d"),
                                              as.Date(max(as.Date(mergedTable$data$Date, "%Y-%m-%d")[!is.na(as.Date(mergedTable$data$Date, "%Y-%m-%d"))]),"%Y-%m-%d")),step = 1)

            #Remove the processing files modal
            removeModal()

          },
          error = function(e){
            shiny::showModal(shiny::modalDialog(
              title = "Incorrect File Type - 1",
              "Error - Loading one of the data files, please check the files and try again - 1"
            ))
          },
          warning = function(w){
            shiny::showModal(shiny::modalDialog(
              title = "Incorrect File Type - 2",
              "Error - Loading one of the data files, please check the files and try again - 2"
            ))
          }
        )#Closing trycatch
      }else{#Closing off the if grep metaDataFile
        shiny::showModal(shiny::modalDialog(
          title = "Incorrect File Type - 3",
          "Incorrect file type, please select a properly formatted '.tsv' file and resubmit."
        ))
      }
    }else{#The file extension was incorrect
      shiny::showModal(shiny::modalDialog(
        title = "Incorrect File Type - 4",
        "Incorrect file type, please select a properly formatted '.tsv' file and resubmit."
      ))
    }#End of if-else right extension

  },ignoreInit = TRUE)# end of the observeEvent Input Submit

  ############### Filtering Button ###########################################

  shiny::observeEvent(input$updateFilterMappingButton, {

    shiny::showModal(shiny::modalDialog("Updating, please standby...",  footer=NULL))

    workMergedTable <- mergedTable$data
    workMergedTable <- workMergedTable[workMergedTable$Abundance >= input$abundance[1] & workMergedTable$Abundance <= input$abundance[2],]
    workMergedTable <- workMergedTable[workMergedTable$Final_Rank %in% input$finalRank,]
    workMergedTable <- workMergedTable[workMergedTable$superkingdom %in% input$kingdomFilterInput,]
    workMergedTable <- workMergedTable[workMergedTable$phylum %in% input$phylumFilterInput,]
    workMergedTable <- workMergedTable[workMergedTable$class %in% input$classFilterInput,]
    workMergedTable <- workMergedTable[workMergedTable$order %in% input$orderFilterInput,]
    workMergedTable <- workMergedTable[workMergedTable$family %in% input$familyFilterInput,]
    workMergedTable <- workMergedTable[workMergedTable$genus %in% input$genusFilterInput,]
    workMergedTable <- workMergedTable[workMergedTable$species %in% input$speciesFilterInput,]
    workMergedTable <- workMergedTable[workMergedTable$Sample %in% input$sampleFilterInput,]
    workMergedTable <- workMergedTable[workMergedTable$Run %in% input$runFilterInput,]
    workMergedTable <- workMergedTable[workMergedTable$Lab %in% input$labFilterInput,]
    workMergedTable <- workMergedTable[workMergedTable$Type %in% input$typeFilterInput,]
    workMergedTable <- workMergedTable[workMergedTable$Marker %in% input$markerFilterInput,]
    workMergedTable <- workMergedTable[workMergedTable$Date >= input$dateInput[1] & workMergedTable$Date <= input$dateInput[2],]

    AVal <- input$abundance[1]
    BVal <- input$abundance[2]
    CVal <- input$finalRank
    DVal <- input$kingdomFilterInput
    EVal <- input$phylumFilterInput
    FVal <- input$classFilterInput
    GVal <- input$orderFilterInput
    HVal <- input$familyFilterInput
    IVal <- input$genusFilterInput
    JVal <- input$speciesFilterInput
    KVal <- input$sampleFilterInput
    LVal <- input$runFilterInput
    MVal <- input$labFilterInput
    NVal <- input$typeFilterInput
    OVal <- input$markerFilterInput
    PVal <- input$dateInput[1]
    QVal <- input$dateInput[2]

    shiny::updateSliderInput(session = session, inputId = "abundance",value=c(AVal,BVal),step = 1)
    shinyWidgets::updatePickerInput(session, "finalRank", choices = sort(unique(workMergedTable$Final_Rank), na.last = TRUE), selected = CVal)
    shinyWidgets::updatePickerInput(session, "kingdomFilterInput", choices = sort(unique(workMergedTable$superkingdom), na.last = TRUE), selected = DVal)
    shinyWidgets::updatePickerInput(session, "phylumFilterInput", choices = sort(unique(workMergedTable$phylum), na.last = TRUE), selected = EVal)
    shinyWidgets::updatePickerInput(session, "classFilterInput", choices = sort(unique(workMergedTable$class), na.last = TRUE), selected = FVal)
    shinyWidgets::updatePickerInput(session, "orderFilterInput", choices = sort(unique(workMergedTable$order), na.last = TRUE), selected = GVal)
    shinyWidgets::updatePickerInput(session, "familyFilterInput", choices = sort(unique(workMergedTable$family), na.last = TRUE), selected = HVal)
    shinyWidgets::updatePickerInput(session, "genusFilterInput", choices = sort(unique(workMergedTable$genus), na.last = TRUE), selected = IVal)
    shinyWidgets::updatePickerInput(session, "speciesFilterInput", choices = sort(unique(workMergedTable$species), na.last = TRUE), selected = JVal)
    shinyWidgets::updatePickerInput(session, "sampleFilterInput", choices = sort(unique(workMergedTable$Sample), na.last = TRUE), selected = KVal)
    shinyWidgets::updatePickerInput(session, "runFilterInput", choices = sort(unique(workMergedTable$Run), na.last = TRUE), selected = LVal)
    shinyWidgets::updatePickerInput(session, "labFilterInput", choices = sort(unique(workMergedTable$Lab), na.last = TRUE), selected = MVal)
    shinyWidgets::updatePickerInput(session, "typeFilterInput", choices = sort(unique(workMergedTable$Type), na.last = TRUE), selected = NVal)
    shinyWidgets::updatePickerInput(session, "markerFilterInput", choices = sort(unique(workMergedTable$Marker), na.last = TRUE), selected = OVal)
    shiny::updateSliderInput(session = session, inputId = "dateInput", value=c(as.Date(PVal,"%Y-%m-%d"),as.Date(QVal,"%Y-%m-%d")),step = 1)


    leaflet::leafletProxy("mymap", data = as.data.frame(workMergedTable)) %>%
      clearMarkers() %>%
      clearMarkerClusters() %>%
      clearPopups() %>%
      #Adding labels to markers on map
      addCircleMarkers(lng = ~North,
                       lat = ~West,
                       color = ~palette_map(workMergedTable$AbundanceCategory),
                       clusterOptions = markerClusterOptions(spiderfyDistanceMultiplier=1.5),
                       popup= paste("<strong><h5>Species:", workMergedTable$Final_Taxa, "</strong>",
                                    "<br><h6>Sample:", workMergedTable$Sample,
                                    "<br><h6>Marker:", workMergedTable$Marker,
                                    "<br><h6>Date:", workMergedTable$Date,
                                    "<br><h6>Run:", workMergedTable$Run,
                                    "<br><h6>Type:", workMergedTable$Type,
                                    "<br><h6>Lab:", workMergedTable$Lab,
                                    "<h6>Coord(Lat, Lon):", workMergedTable$North,",", workMergedTable$West))



    removeModal()

  },ignoreInit = TRUE)


  ############### Reset Filtering Button ###########################################

  shiny::observeEvent(input$resetFilterMappingButton, {

    shiny::showModal(shiny::modalDialog("Updating, please standby...",  footer=NULL))

    shiny::updateSliderInput(session = session, inputId = "abundance",
                             min=0,
                             max = max(mergedTable$data$Abundance),
                             value=c(0,max(mergedTable$data$Abundance))
    )

    #Dropdown menu for Final_Rank
    shinyWidgets::updatePickerInput(session, inputId = "finalRank",
                                    choices = sort(unique(mergedTable$data$Final_Rank), na.last = TRUE),
                                    selected = sort(unique(mergedTable$data$Final_Rank), na.last = TRUE))

    #Dropdown menu for superkingdom
    shinyWidgets::updatePickerInput(session, inputId = "kingdomFilterInput",
                                    choices = sort(unique(mergedTable$data$superkingdom), na.last = TRUE),
                                    selected = sort(unique(mergedTable$data$superkingdom), na.last = TRUE))

    #Dropdown menu for phylum
    shinyWidgets::updatePickerInput(session, inputId = "phylumFilterInput",
                                    choices = sort(unique(mergedTable$data$phylum), na.last = TRUE),
                                    selected = sort(unique(mergedTable$data$phylum), na.last = TRUE))

    #Dropdown menu for class
    shinyWidgets::updatePickerInput(session, inputId = "classFilterInput",
                                    choices = sort(unique(mergedTable$data$class), na.last = TRUE),
                                    selected = sort(unique(mergedTable$data$class), na.last = TRUE))

    #Dropdown menu for order
    shinyWidgets::updatePickerInput(session, inputId = "orderFilterInput",
                                    choices = sort(unique(mergedTable$data$order), na.last = TRUE),
                                    selected = sort(unique(mergedTable$data$order), na.last = TRUE))

    #Dropdown menu for family
    shinyWidgets::updatePickerInput(session, inputId = "familyFilterInput",
                                    choices = sort(unique(mergedTable$data$family), na.last = TRUE),
                                    selected = sort(unique(mergedTable$data$family), na.last = TRUE))

    #Dropdown menu for genus
    shinyWidgets::updatePickerInput(session, inputId = "genusFilterInput",
                                    choices = sort(unique(mergedTable$data$genus), na.last = TRUE),
                                    selected = sort(unique(mergedTable$data$genus), na.last = TRUE))

    #Dropdown menu for species
    shinyWidgets::updatePickerInput(session, inputId = "speciesFilterInput",
                                    choices = sort(unique(mergedTable$data$species), na.last = TRUE),
                                    selected = sort(unique(mergedTable$data$species), na.last = TRUE))

    #################### update the filters based on submitted Meta Data ###########

    #Dropdown menu for Sample
    shinyWidgets::updatePickerInput(session, inputId = "sampleFilterInput",
                                    choices = sort(unique(mergedTable$data$Sample), na.last = TRUE),
                                    selected = sort(unique(mergedTable$data$Sample), na.last = TRUE))

    #Dropdown menu for Run
    shinyWidgets::updatePickerInput(session, inputId = "runFilterInput",
                                    choices = sort(unique(mergedTable$data$Run), na.last = TRUE),
                                    selected = sort(unique(mergedTable$data$Run), na.last = TRUE))

    #Dropdown menu for Lab
    shinyWidgets::updatePickerInput(session, inputId = "labFilterInput",
                                    choices = sort(unique(mergedTable$data$Lab), na.last = TRUE),
                                    selected = sort(unique(mergedTable$data$Lab), na.last = TRUE))

    #Dropdown menu for Type
    shinyWidgets::updatePickerInput(session, inputId = "typeFilterInput",
                                    choices = sort(unique(mergedTable$data$Type), na.last = TRUE),
                                    selected = sort(unique(mergedTable$data$Type), na.last = TRUE))
    #Dropdown menu for Region
    shinyWidgets::updatePickerInput(session,inputId = "markerFilterInput",
                                    choices = sort(unique(mergedTable$data$Marker), na.last = TRUE),
                                    selected = sort(unique(mergedTable$data$Marker), na.last = TRUE))

    shiny::updateSliderInput(session = session, inputId = "dateInput",
                             min = as.Date(min(as.Date(mergedTable$data$Date, "%Y-%m-%d")[!is.na(as.Date(mergedTable$data$Date, "%Y-%m-%d"))]),"%Y-%m-%d"),
                             max = as.Date(max(as.Date(mergedTable$data$Date, "%Y-%m-%d")[!is.na(as.Date(mergedTable$data$Date, "%Y-%m-%d"))]),"%Y-%m-%d"),
                             value=c(as.Date(min(as.Date(mergedTable$data$Date, "%Y-%m-%d")[!is.na(as.Date(mergedTable$data$Date, "%Y-%m-%d"))]),"%Y-%m-%d"),
                                     as.Date(max(as.Date(mergedTable$data$Date, "%Y-%m-%d")[!is.na(as.Date(mergedTable$data$Date, "%Y-%m-%d"))]),"%Y-%m-%d")),step = 1)

    leaflet::leafletProxy("mymap", data = as.data.frame(mergedTable)) %>%
      clearMarkers() %>%
      clearMarkerClusters() %>%
      clearPopups() %>%
      #Adding labels to markers on map
      addCircleMarkers(lng = ~North,
                       lat = ~West,
                       color = ~palette_map(mergedTable$AbundanceCategory),
                       clusterOptions = markerClusterOptions(spiderfyDistanceMultiplier=1.5),
                       popup= paste("<strong><h5>Species:", mergedTable$Final_Taxa, "</strong>",
                                    "<br><h6>Sample:", mergedTable$Sample,
                                    "<br><h6>Marker:", mergedTable$Marker,
                                    "<br><h6>Date:", mergedTable$Date,
                                    "<br><h6>Run:", mergedTable$Run,
                                    "<br><h6>Type:", mergedTable$Type,
                                    "<br><h6>Lab:", mergedTable$Lab,
                                    "<h6>Coord(Lat, Lon):", mergedTable$North,",", mergedTable$West))

    removeModal()

  },ignoreInit = TRUE)

} # End of Server
