# Written by Rob Young at the University of Guelph in Ontario Canada, Sept, 2023
# ******************************************************************************
#Roxygen2 Documentation:

#' @export
#'
#' @title Combine Taxa Assignment for Same ASV Using Different Databases
#'
#' @author Robert G. Young
#'
#' @description
#' This function takes a file selection and then uses all 'taxaAssign' files in
#' that directory and combines them into a single output
#' 'taxaAssignCombined.tsv' file.
#'
#' @details
#' The User Input: This function requires a file in a directory where all 'taxaAssign'
#' files in that directory will be combined.
#'
#' @examples
#' \dontrun{
#' combine_assign_output()
#' combine_assign_output(fileLoc = NULL,   numCores = 1)
#' }
#'
#' @param fileLoc The location of a file in a directory where all of the 'taxaAssign'
#' files are located.
#' @param numCores The number of cores used to run the function (default = 1,
#' Windows systems can only use a single core).
#'
#' @returns
#' This function produces a '2023_08_03_0913_taxaAssignCombined.tsv' and
#' a '2023_08_03_0913_taxaAssignCombined.txt' file in the selected target
#' directory.
#'
#' @references
#' <https://github.com/rgyoung6/DBTC>
#' Young, R. G., Hanner, R. H. (Submitted October 2023). Title Here. Biodiversity Data Journal.
#'
#' @note
#' When running DBTCShiny functions the paths for the files selected cannot have
#' whitespace! File folder locations should be as short as possible (close to
#' the root directory) as some functions do not process long naming conventions.
#' Also, special characters should be avoided (including question mark, number
#' sign, exclamation mark). It is recommended that dashes be used for
#' separations in naming conventions while retaining underscores for use as
#' information delimiters (this is how DBTC functions use underscore). There
#' are several key character strings used in the DBTC pipeline, the presence of
#' these strings in file or folder names will cause errors when running DBTC
#' functions.
#'
#' The following strings are those used in DBTC and should not be used in file
#' or folder naming:
#' - _BLAST
#' - _taxaAssign
#' - _taxaCombined
#' - _taxaReduced
#'
#' @seealso
#' dada_implement()
#' combine_dada_output()
#' make_BLAST_DB()
#' seq_BLAST()
#' taxon_assign()
#' reduce_taxa()
#' combine_reduced_output()
#'
################################ COMBINE TWO REDUCED TAXA FILES INTO A SINGLE FILE ##################
combine_assign_output <- function(fileLoc = NULL,   numCores = 1){

  #Get the initial working directory
  start_wd <- getwd()
  on.exit(setwd(start_wd))

  #load in the files list
  if (is.null(fileLoc)){
    print("Select a file in the file folder with the taxa assigned files")
    print("you would like to combine (extension '_taxaAssign_YYYY_MM_DD_HHMM.tsv').")
    print("*** NOTE: all '_taxaAssign_' files in the folder location should")
    print("originate from the same dada output file but have outputs from different")
    print("BLAST sequence libraries and therefore conatin the same ASVs ***")
    fileLoc <- file.choose()
  }
  if (is.null(fileLoc)){

    print("********************************************************************************")
    print("The file location is required and needs to be submited as an argument 'fileLoc'")
    print("when calling the combine_ouput() function or when prompted to select the folder")
    print("through a popup window (where available).")
    print("Please rerun the function and provide a character string for the fileLoc")
    print("argument.")
    print(paste0("Current file location (fileLoc) is: ", fileLoc))
    print("********************************************************************************")

  }else{

    #Printing the start time
    print(paste0("Start time...", Sys.time()))
    startTime <- paste0("Start time...", Sys.time())
    dateStamp <- paste0(format(Sys.time(), "%Y_%m_%d_%H%M"))

    fileLoc <- dirname(fileLoc)
    #Set the working directory
    setwd(fileLoc)

    #Get the all files in the selected file folder with full paths of the files
    files <- as.data.frame(list.files(path = fileLoc, pattern = "*[.]*", full.names = TRUE))
    # Get the local paths
    files[,2] <- list.files(path = fileLoc, pattern = "*[.]*")
    #Get all files with the '_taxaAssign' string
    files <- files[grepl("_taxaAssign_.*", files[,2]),]
    #Remove files with _taxaAssignCombined_
    files <- files[!grepl("_taxaAssignCombined_.*", files[,2]),]
    # Get the names of the files
    files[,3] <- gsub("_taxaAssign_.*","",files[,2])
    # Get a unique number for each file to associate with the combined output results
    files[,4] <- c(1:nrow(files))

    #Add in headers for the file table
    colnames(files)<-c("File_Path", "File", "File_Name", "File_Unique_ID")

    #Write start time to the summary file
    suppressWarnings(write(paste0(dateStamp, " Start"), file = paste0(dateStamp, "_taxaAssignCombined.txt"), append = FALSE))

    #Write the files being combined
    suppressWarnings(write(paste0("The files being combined with an associated unique identifying number (File_Unique_ID column) used in the uniqueID column of the output file..."), file = paste0(dateStamp, "_taxaAssignCombined.txt"), append = TRUE))
    suppressWarnings(write.table(files, file = paste0(dateStamp, "_taxaAssignCombined.txt"), append = TRUE, na = "NA", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"))

    if(nrow(files)<2){

      print("********************************************************************************")
      print("There is only one file in the target directory for the file type and so ")
      print("the combination of files is not needed.")
      print("********************************************************************************")

    }else{

      #Initiate totalResults final data frame
      totalResults = NULL

      print(paste0("Reading in all of the submitted files at... ", Sys.time()))
      suppressWarnings(write(paste0("Reading in all of the submitted files at... ",Sys.time()), file = paste0(dateStamp, "_taxaAssignCombined.txt"), append = TRUE))

      #Load in all of the records
      for(records in 1:nrow(files)){

        #Read in the files for this loop
        loopResult <- as.data.frame(read.delim(files[records,1], header = TRUE, check.names=FALSE))

        #Add the file name to the front of the data frame
        loopResult <- cbind(File = files[records,"File_Name"], loopResult)

        if(is.null(totalResults)){
          totalResults <-loopResult
        }else{
          if(identical(colnames(loopResult),colnames(totalResults))){
            totalResults<-rbind(totalResults, loopResult)
          }else{
            stop("The files that you have submitted are not from the same root file as they do not have the same headers (samples). Please check the files and resubmit.")
          }
        }
      }#End of the loop through the files and loading the files

      print(paste0("Completed reading the submitted files at... ", Sys.time()))
      suppressWarnings(write(paste0("Completed reading the submitted files at... ",Sys.time()), file = paste0(dateStamp, "_taxaAssignCombined.txt"), append = TRUE))

      #Remove unnecessary loop variable
      rm(loopResult)

      #Get unique sequences
      uniqueRecords <- unique(totalResults$Sequence)

      #Initiate finalResults final data frame
      finalResults = NULL

      print(paste0("Begin looping through the unique reads at... ", Sys.time()))
      suppressWarnings(write(paste0("Begin looping through the unique reads at... ",Sys.time()), file = paste0(dateStamp, "_taxaAssignCombined.txt"), append = TRUE))


      taxaResults <- function(uniqueRecordsCount){
      #Loop through each of the unique reads and keep the best taxonomic result
#      for(uniqueRecordsCount in 1:length(uniqueRecords)){

        #subset the totalResults data frame to get the records for the loop sequence
        loopResult <- totalResults[totalResults$Sequence==uniqueRecords[uniqueRecordsCount],, drop = FALSE]

        #Check to see if there are records other than NA in the Final_Taxa column
        if(nrow(loopResult[!is.na(loopResult$Final_Taxa),, drop = FALSE]) == 0){


          #print(paste0("Record...", loopResult[1,"uniqueID"], " all NA with rows...", nrow(loopResult)))


          #If after the removal of duplicate records there is more than one record
          if(nrow(loopResult[!duplicated(loopResult[,3:ncol(loopResult)]), , drop=FALSE]) > 1){

            #print(paste0("Record...", loopResult[1,"uniqueID"], " more than 1 derep with rows...", nrow(loopResult)))


            #Choose the best record from this set where the taxa is NA
            loopResult<-loopResult[which.max(rowSums(!is.na(loopResult[,c(1:18)]))),, drop=FALSE]

            #If there are more than one with the same number of NA then take the first
            #and combine the file origins
            if(nrow(loopResult)>1){

              #Take the loopResultFiles and combine into a single variable and
              #Change the first row File variable to the combined value
              loopResult[1,"File"]<- paste0(unique(loopResult$File), collapse = ",")

              #Keep the first row of the results
              loopResult<-loopResult[1,, drop=FALSE]

            }

          }else{ #Else if there is only a single dereplicated records that is NA for taxa

            #print(paste0("Here is the loop record...", loopResult[1,"uniqueID"], " and the number of rows of loopResult...", nrow(loopResult)))

            #Take the loopResultFiles and combine into a single variable and
            #Change the first row File variable to the combined value
            loopResult[1,"File"]<- paste0(unique(loopResult$File), collapse = ",")

            #Keep the first row of the results
            loopResult<-loopResult[1,, drop=FALSE]

          }

        }else{ #Else if there are records that are not NA in the Final_Taxa

          #If there is only one record that is not NA then use this one
          if(nrow(loopResult[!is.na(loopResult$Final_Taxa),, drop = FALSE]) == 1){

            #Keep non NA Final_Taxa records
            loopResult <- loopResult[!is.na(loopResult$Final_Taxa),, drop = FALSE]

          } else if(nrow(loopResult[!is.na(loopResult$Final_Taxa),, drop = FALSE])> 1){

            #Choose best record where the final taxa is complete
            #If after the removal of duplicate records there is more than one record
            if(nrow(loopResult[!duplicated(loopResult[,3:ncol(loopResult)]), , drop=FALSE]) > 1){

              #Choose the best record where the taxa is NA by first getting the result with the fewest NA
              loopResult<-loopResult[which.max(rowSums(!is.na(loopResult[,c(1:18)]))),, drop=FALSE]

              #If there is still more than one record look at the contents of the Final_Taxa results
              if(nrow(loopResult)>1){

                #Split the column by ( and replace the close bracket )
                loopResultCompare <- as.data.frame(do.call('rbind', strsplit(as.character(loopResult[1,"Final_Taxa"]),'(',fixed = TRUE)))
                loopResultCompare[,2] <- gsub(")","", as.character(loopResultCompare[,2]))
                #Split the values column by ,
                loopResultCompare <- cbind(loopResultCompare[,1],data.frame(do.call('rbind', strsplit(as.character(loopResultCompare[,2]),',', fixed = TRUE))), row.names = NULL)
                #Name the total results temp columns
                colnames(loopResultCompare) <- c("Final_Taxa", "Num_Rec", "Coverage", "Identity", "Max_eVal")

                #Keep the row with the most coverage
                loopResult<-loopResult[which.max(loopResultCompare$Coverage),]

                #If there are more than one with the same number of NA then take the first
                #and combine the file origins
                if(nrow(loopResult)>1){

                  #Keep the row with the highest identity
                  loopResult<-loopResult[which.max(loopResultCompare$Identity),]

                  #If there are more than one with the same number of NA then take the first
                  #and combine the file origins
                  if(nrow(loopResult)>1){

                    #Take the loopResultFiles and combine into a single variable and
                    #Change the first row File variable to the combined value
                    loopResult[1,"File"]<- paste0(unique(loopResult$File), collapse = ",")

                    #Keep the first row of the results
                    loopResult<-loopResult[1,, drop=FALSE]
                  }
                }

              }

            }else{

              #Remove exact duplicate records.
              loopResult <- loopResult[!duplicated(loopResult), , drop=FALSE]

            }

          }#End of the if/else if where the records are not NA but there are records left
        }#End of the if/else if/else where we check for non NA values

        #Update the uniqueID to include the file name
        loopResult[1,"uniqueID"] <- paste0(loopResult[1,"uniqueID"],",",loopResult$File)

        #Remove the File column
        loopResult<-loopResult[,-1]

        return <- loopResult

      }#End of looping through the unique reads



      if(numCores==1){

        finalResults <- pbapply::pblapply(seq_len(length(uniqueRecords)), taxaResults)

      }else{

        finalResults <- pbapply::pblapply(seq_len(length(uniqueRecords)), taxaResults, cl = numCores)

      }

      finalResults <- do.call(rbind.data.frame, finalResults)

      print(paste0("The end looping through the unique reads at... ", Sys.time()))
      suppressWarnings(write(paste0("The end looping through the unique reads at... ",Sys.time()), file = paste0(dateStamp, "_taxaAssignCombined.txt"), append = TRUE))

      #Export the results
      suppressWarnings(write.table(finalResults, file = paste0(dateStamp, "_taxaAssignCombined.tsv"), append = FALSE, na = "NA", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"))

      print(paste0(startTime, " End time ", Sys.time()))
      suppressWarnings(write(paste0(startTime, " End time ", Sys.time()), file = paste0(dateStamp, "_taxaAssignCombined.txt"), append = TRUE))

    } #End of the if/else there are more than 2 files
  }# End of the if/else confirming file location
}# End of function




