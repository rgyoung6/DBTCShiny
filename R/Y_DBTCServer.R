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

        mergedTableGlobal <<- mergedTable$data
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
  shiny::observeEvent(input$dadaDirectoryButton, {
    tryCatch(
      expr = {
        dadaLocation$data <- dirname(dirname(file.choose()))
        output$dadaDirectoryDisplay <- shiny::renderText({as.character(dadaLocation$data)})
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
  },ignoreInit = TRUE)

  # Get the file with the primer data for this analysis
  shiny::observeEvent(input$primerFileButton, {
    tryCatch(
      expr = {
        primerFile$data <- file.choose()
        output$primerFileDisplay <- shiny::renderText({as.character(primerFile$data)})
      },
      error = function(e){
        print("Error - Primer File Button choose file cancelled")
        primerFile$data <- NA
      },
      warning = function(w){
        print("Warning - Primer File Button choose file cancelled")
        primerFile$data <- NA
      }
    )
  },ignoreInit = TRUE)

  shiny::observeEvent(input$dadaSubmit, {
    suppressWarnings(if(!is.na(dadaLocation$data) && !is.na(primerFile$data)){

      # Create variables to call the dada_implement so that there are no conflicts
      # with the multithreading and the shiny app

       runFolderLoc <- force(as.character(dadaLocation$data))
       primerFile <- force(as.character(primerFile$data))
       fwdIdent <- force(input$fwdIdent)
       revIdent <- force(input$revIdent)
       nonMergeProcessing <- force(input$nonMergeProcessing)
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
                     nonMergeProcessing = nonMergeProcessing,
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
       }
       )
    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please select a primer file and try submitting again!"
      ))
    })
  },ignoreInit = TRUE)

  ################## Dada Combine Function ####################################
  #Get the location of the dada output files you would like to combine
  shiny::observeEvent(input$dadaCombineButton, {
    tryCatch(
      expr = {
        dadaCombineLoc$data <- dirname(file.choose())
        output$dadaCombineDisplay <- shiny::renderText({as.character(dadaCombineLoc$data)})
      },
      error = function(e){
        print("Error")
        dadaCombineLoc$data <- NA
      },
      warning = function(w){
        print("Warning")
        dadaCombineLoc$data <- NA
      }
    )
  },ignoreInit = TRUE)

  #Running the data combine
  shiny::observeEvent(input$dadaCombine, {

    suppressWarnings(if(!is.na(dadaCombineLoc$data)){
      # Create variables for the arguments to avoid conflicts between the multithreading
      # and the shiny
      fileLoc= force(dadaCombineLoc$data)
      minLen = force(input$dadaCombineMinLen)

      tryCatch(
        expr = {
          #Run the Dada function here.

          shiny::showModal(shiny::modalDialog(
            title = "Dada combine analysis results is underway.",
            "Processing, please stand by...", footer=""

          ))

        #Run the Dada combine function here.
        combine_dada_output(fileLoc = fileLoc,
                            minLen = minLen)

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
    })
  },ignoreInit = TRUE)

  ################## Make BLAST DB Function ###################################

  # Get the fasta file you want to use to build your db
  shiny::observeEvent(input$makeBlastDBFileLocButton, {
    tryCatch(
      expr = {
        makeBlastDBFileLoc$data <- file.choose()
        output$makeBlastDBFileLocDisplay <- shiny::renderText({as.character(makeBlastDBFileLoc$data)})
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


  # Select where the makeblastdb program is on your computer
  shiny::observeEvent(input$makeblastdbPathButton, {
    tryCatch(
      expr = {
        makeblastdbPath$data <- file.choose()
        output$makeblastdbPathDisplay <- shiny::renderText({as.character(makeblastdbPath$data)})
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

  # Point to the NCBI taxonomic data base
  shiny::observeEvent(input$makeBlastTaxaDBLocButton, {
    tryCatch(
      expr = {
        makeBlastTaxaDBLoc$data <- file.choose()
        output$makeBlastTaxaDBLocDisplay <- shiny::renderText({as.character(makeBlastTaxaDBLoc$data)})
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

    suppressWarnings(if(!is.na(makeBlastDBFileLoc$data) && !is.na(makeBlastTaxaDBLoc$data) && !is.na(input$inputFormat) && !is.na(input$dbName) && !is.na(input$makeBLASTDBMinLen)){
      if(is.na(makeblastdbPath$data)){
        makeblastdbPath$data <- "makeblastdb"
      }
      # Create local variables to avoid conflicts with shiny and multithread
      fileLoc = force(makeBlastDBFileLoc$data)
      makeblastdbPath = force(makeblastdbPath$data)
      taxaDBLoc = force(makeBlastTaxaDBLoc$data)
      inputFormat = force(input$inputFormat)
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
                        inputFormat = inputFormat,
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
            ""
          ))
        },
        warning = function(w){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            ""
          ))
        }
      )

    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please fill in all of the necessary fields and submit again!"
      ))
    })
  },ignoreInit = TRUE)

  ################## BLAST sequences Function #################################

  # Point to the NCBI taxonomic data base
  shiny::observeEvent(input$BLASTDatabasePathButton, {
    tryCatch(
      expr = {
        BLASTDatabasePath$data <- file.choose()
        output$BLASTDatabasePathDisplay <- shiny::renderText({as.character(BLASTDatabasePath$data)})
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

  # Get the path where all of the folders containing the fastq files are located
  shiny::observeEvent(input$blastnPathButton, {
    tryCatch(
      expr = {
        blastnPath$data <- file.choose()
        output$blastnPathDisplay <- shiny::renderText({as.character(blastnPath$data)})
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

  # Point to the location of the fasta files you want to BLAST
  shiny::observeEvent(input$querySeqPathButton, {
    tryCatch(
      expr = {
        querySeqPath$data <- file.choose()
        output$querySeqPathDisplay <- shiny::renderText({as.character(querySeqPath$data)})
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

    suppressWarnings(if(!is.na(BLASTDatabasePath$data) && !is.na(querySeqPath$data)&& !is.na(input$BLASTResults) && !is.na(input$blastSeqNumCores)){
      if(is.na(blastnPath$data)){
        blastnPath$data <- "blastn"
      }

      # Create local variables to avoid conflicts with shiny and multithread
        databasePath = force(BLASTDatabasePath$data)
        blastnPath = force(blastnPath$data)
        querySeqPath = force(querySeqPath$data)
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
    })
  },ignoreInit = TRUE)

  ################## Taxon Assign Function ####################################
  # Point to the BLAST output files to get taxonomic assignment
  shiny::observeEvent(input$taxaAssignFileLocButton, {
    tryCatch(
      expr = {
        taxaAssignFileLoc$data <- file.choose()
        output$taxaAssignFileLocDisplay <- shiny::renderText({as.character(taxaAssignFileLoc$data)})
      },
      error = function(e){
        print("Error")
        taxaAssignFileLoc$data <- NA
      },
      warning = function(w){
        print("Warning")
        taxaAssignFileLoc$data <- NA
      }
    )
  },ignoreInit = TRUE)

  # Point to the NCBI taxonomic data base
  shiny::observeEvent(input$taxaAssignDBLocButton, {
    tryCatch(
      expr = {
        taxaAssignDBLoc$data <- file.choose()
        output$taxaAssignDBLocDisplay <- shiny::renderText({as.character(taxaAssignDBLoc$data)})
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
    suppressWarnings(if(!is.na(taxaAssignFileLoc) && !is.na(taxaAssignDBLoc)){

print("Made it into the submission buttion past the no NA values if")

      # Create local variables to avoid conflicts with shiny and multithread
       fileLoc = force(taxaAssignFileLoc$data)
       taxaDBLoc = force(taxaAssignDBLoc$data)
       numCores = force(input$taxaAssignNumCores)
       coverage = force(input$coverage)
       ident = force(input$ident)
       propThres = force(input$propThres)
       coverReportThresh = force(input$coverReportThresh)
       identReportThresh = force(input$identReportThresh)
       includeAllDada = force(input$includeAllDada)

print("Made it into the sumbission buttion section after assiging values to local variables")

       tryCatch(
         expr = {
           #Run the function here.

print("Right before the modal processing message")

           shiny::showModal(shiny::modalDialog(
             title = "Taxonomic assingment is underway.",
             "Processing, please stand by...", footer=""

           ))
print("Right after the modal processing modal and right before running the function")
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

print("Right after running the function")

            removeModal()

print("Right after the removal of the modal and right before the next modal")

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
    })
  })

  ################## Combine Taxa Assign Function #############################

  # Point to the folder with the Taxon Assign records you would like to combine
  shiny::observeEvent(input$combineTaxaFileLocButton, {
    tryCatch(
      expr = {
        combineTaxaFileLoc$data <- file.choose()
        output$combineTaxaFileLocDisplay <- shiny::renderText({as.character(combineTaxaFileLoc$data)})
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
    suppressWarnings(if(!is.na(combineTaxaFileLoc$data)){
      # Create local variables to avoid conflicts with shiny and multithread
      fileLoc = force(combineTaxaFileLoc$data)
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
    })
  })

  ################## Reduce Taxa Assign Function ##############################

  # Point to the folder with the Taxon Assign records you would like to combine
  shiny::observeEvent(input$reduceTaxaFileLocButton, {
    tryCatch(
      expr = {
        reduceTaxaFileLoc$data <- file.choose()
        output$reduceTaxaFileLocDisplay <- shiny::renderText({as.character(reduceTaxaFileLoc$data)})
      },
      error = function(e){
        print("Error")
        reduceTaxaFileLoc$data <- NA
      },
      warning = function(w){
        print("Warning")
        reduceTaxaFileLoc$data <- NA
      }
    )
  },ignoreInit = TRUE)


  shiny::observeEvent(input$reduceTaxa, {
    suppressWarnings(if(!is.na(reduceTaxaFileLoc$data)){

      # Create local variables to avoid conflicts with shiny and multithread
      fileLoc = force(reduceTaxaFileLoc$data)
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
    })
  })

  ################## Combine Reduce Taxa Assign Function ######################

  # Point to the folder with the Taxon Assign records you would like to combine
  shiny::observeEvent(input$combineReducedTaxaFileLocButton, {
    tryCatch(
      expr = {
        combineReducedTaxaFileLoc$data <- file.choose()
        output$combineReducedTaxaFileLocDisplay <- shiny::renderText({as.character(combineReducedTaxaFileLoc$data)})
      },
      error = function(e){
        print("Error")
        combineReducedTaxaFileLoc$data <- NA
      },
      warning = function(w){
        print("Warning")
        combineReducedTaxaFileLoc$data <- NA
      }
    )
  },ignoreInit = TRUE)

  shiny::observeEvent(input$combineReduceTaxa, {
    suppressWarnings(if(!is.na(combineReducedTaxaFileLoc$data)){

      # Create local variables to avoid conflicts with shiny and multithread
      fileLoc = force(combineReducedTaxaFileLoc$data)
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
    })
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
mergedTableGlobal <<- mergedTable
print("Here loading data 14")
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

