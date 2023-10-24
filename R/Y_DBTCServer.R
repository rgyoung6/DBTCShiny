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
  maxReads <- reactiveValues(data = NA)

  ###################### Initialize the Map ######################################
  #Building the initial map and using the
  #Default static leaflet map before filtering parameters are applied.
  output$mymap <- leaflet::renderLeaflet({
    leaflet::leaflet() %>%
      leaflet::addTiles() %>%

      leaflet::addLegend(position = "topright",#"bottomright",
                pal = leaflet::colorFactor(palette = c("#fbd300",
                                                       "#ff8600",
                                                       "#ea5f94",
                                                       "#9d02d7",
                                                       "#0000ff"),
                                           levels = c("1",
                                                      "10",
                                                      "100",
                                                      "1000",
                                                      "10000")),
                values =  c("1",
                            "10",
                            "100",
                            "1000",
                            "10000"),
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
      }
    }
  })

  ####################### Data Export Download Button ###########################

  #Download button for downloading CSV of filtered mapping data
  #  output$downloadTable <- downloadHandler(
  #    filename = function() {
  #      paste0(format(Sys.time(), "%Y_%m_%d_%H%M"), "_MDMAPR_Data.tsv")
  #    },
  #    content = function(file) {
  #      write.csv(as.data.frame(filtered$value[2:ncol(uploaded_data$value)]), file,  row.names=FALSE)
  #    }
  #  )

  #  print("Here 10")


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

       tryCatch(
         expr = {
           #Run the function here.

           shiny::showModal(shiny::modalDialog(
             title = "Taxonomic assingment is underway.",
             "Processing, please stand by...", footer=""

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

      tryCatch(
        expr = {
          #Run the function here.

          shiny::showModal(shiny::modalDialog(
            title = "Combine reduced taxonomic results for multiple markers for the same samples is underway.",
            "Processing, please stand by...", footer=""

          ))

          #Run the function
          combine_reduced_output(fileLoc = fileLoc)

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

        #Loading in the MDMAPR data file
        tryCatch(
          expr = {

            shiny::showModal(shiny::modalDialog("Loading the data files, please wait...",  footer=NULL))

            #Load in the data to the formatted_metadata variable and check format
            ASVFileTable$data<-read.table(ASVFile$data, header = TRUE, check.names=FALSE, sep="\t", dec=".")
            metaDataFileTable$data<-read.table(metaDataFile$data, header = TRUE, check.names=FALSE, sep="\t", dec=".")

            removeModal()

            #################### update the filters based on submitted ASV data ############

            shiny::updateSliderInput(session = session, inputId = "numReads",
                                     max = max(ASVFileTable$data[,22:ncol(ASVFileTable$data)]),
                                     value=c(0,max(ASVFileTable$data[,22:ncol(ASVFileTable$data)]))
            )




            #################### update the filters based on submitted Meta Data ###########

            #Dropdown menu for Sample
            shinyWidgets::updatePickerInput(session, inputId = "sampleIDFilterInput",
                                            choices = unique(metaDataFileTable$data$SampleID),
                                            selected = "None")

            #Dropdown menu for Run
            shinyWidgets::updatePickerInput(session, inputId = "runIDFilterInput",
                                            choices = unique(metaDataFileTable$data$RunID),
                                            selected = "None")

            #Dropdown menu for Lab
            shinyWidgets::updatePickerInput(session, inputId = "labIDFilterInput",
                                            choices = unique(metaDataFileTable$data$LabID),
                                            selected = "None")

            #Dropdown menu for Site
            shinyWidgets::updatePickerInput(session, inputId = "siteIDFilterInput",
                                            choices = unique(metaDataFileTable$data$SiteID),
                                            selected = "None")

            #Dropdown menu for Type
            shinyWidgets::updatePickerInput(session, inputId = "typeFilterInput",
                                            choices = unique(metaDataFileTable$data$Type),
                                            selected = "None")

            #Dropdown menu for Date
            shinyWidgets::updatePickerInput(session,inputId = "dateFilterInput",
                                            choices = unique(metaDataFileTable$data$CollectDate),
                                            selected = "None")

            #Dropdown menu for Region
            shinyWidgets::updatePickerInput(session,inputId = "regionIDFilterInput",
                                            choices = unique(metaDataFileTable$data$RegionID),
                                            selected = "None")
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

  ################# Update Mapping Filters when selections are changed ###########

  #  numReads

  #  observe({
  #    updateSliderInput(session = session, inputId = "numReads",
  #                      choices = project_list(),
  #                      selected = project_list() )
  #  })


} # End of Server

