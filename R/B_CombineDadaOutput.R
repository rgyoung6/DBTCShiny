# Written by Rob Young at the University of Guelph in Ontario Canada, April, 2024
# ******************************************************************************
# Roxygen2 Documentation:
#' @export
#'
#' @title Combine Dada Output
#'
#' @author Robert G. Young
#'
#' @description
#' This function uses DBTC dada_implement ASV output files
#' (YYYY_MM_DD_HH_MM_UserInputRunName_Merge,
#' YYYY_MM_DD_HH_MM_UserInputRunName_MergeFwdRev, and/or
#' YYYY_MM_DD_HH_MM_UserInputRunName_TotalTable) and combines them into a single
#' ASV table with accompanying fasta file. This function also produces a file
#' containing the processing information for the function. The main input
#' argument for this function is the location of a file in a folder containing
#' all ASV tables wanting to be combined. Output files are generated with the
#' naming convention YYYY_MM_DD_HH_MM_combinedDada.
#'
#' @details
#' Two or more files to be combined are required as input for this function.
#' These files need to be ASV files as outputted from the dada_implement() and
#' can include Merge, MergeFwdRev, or TotalTable.tsv files.
#' In addition, the user can input the desired minimum length of sequences that
#' are wanted in the output combined file.
#'
#' @examples
#' \dontrun{
#' combine_dada_output()
#' combine_dada_output(fileLoc = NULL, minLen = 100)
#' }
#'
#' @param fileLoc Select a file in the file folder with dada_implement() results you would
#' like to combine (YYYY_MM_DD_HHMM_FileName_MergeFwdRev OR
#' YYYY_MM_DD_HHMM_FileName_Merge both .tsv and .fas files (Default NULL).
#' @param minLen The minimum final desired length of the read (Default 100).

#' @returns
#' The output from this function includes three files.
#' 1. YYYY_MM_DD_HHMM_combinedDada.tsv - combined ASV table
#' 2. YYYY_MM_DD_HHMM_combinedDada.fas - combined fasta file
#' 3. YYYY_MM_DD_HHMM_combinedDada.tsv - Summary file from the combine_dada_output run
#'
#' @references
#' <https://github.com/rgyoung6/DBTC>
#' Young, R. G., Hanner, R. H. (Submitted October 2023). Dada-BLAST-Taxon Assign-Condense
#' Shiny Application (DBTCShiny). Biodiversity Data Journal.
#'
#' @note
#' When running DBTC functions the paths for the files selected cannot have
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
#' make_BLAST_DB()
#' seq_BLAST()
#' taxon_assign()
#' combine_assign_output()
#' reduce_taxa()
#' combine_reduced_output()

################################ COMBINE TWO REDUCED TAXA FILES INTO A SINGLE FILE ##################

combine_dada_output <- function(fileLoc = NULL, minLen = 100){

  #If there are issues and I need to audit the script make this 1
  auditScript=0

  #Get the initial working directory
  start_wd <- getwd()
  on.exit(setwd(start_wd))

  #load in the files list
  if (is.null(fileLoc)){
    print(paste0("Select a file in the file folder with dada files you would like to combine (extension '_Merge.tsv' OR '_MergeFwdRev.tsv' OR '_Forward.tsv')."))
    print("********** NOTE: all files being combined should have used the same protocols (molecular marker, dada arguments, etc.) **********")
    fileLoc <- file.choose()
  }

  #Get the directory
  fileLoc <- dirname(fileLoc)

  if (is.null(fileLoc)){

    print("********************************************************************************")
    print("The file location is required and needs to be submited as an argument 'fileLoc'")
    print("when calling the combine_ouput() function or when prompted to select a file in ")
    print("the folder of interest through a popup window.")
    print("Please rerun the function and provide a character string for the fileLoc")
    print("argument or select a file when prompted.")
    print(paste0("Current file location (fileLoc) is: ", fileLoc))
    print("********************************************************************************")

  }else{

    #Audit line
    if(auditScript>0){
      auditFile <- paste0(fileLoc,"/", format(Sys.time(), "%Y_%m_%d_%H%M"), "_audit.txt")
      print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 1"))
      suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 1"), file = auditFile, append = FALSE))
    }


    #Printing the start time
    print(paste0("Start time...", Sys.time()))
    startTime <- paste0("Start time...", Sys.time())
    dateStamp <- paste0(format(Sys.time(), "%Y_%m_%d_%H%M"))

    #Set the working directory
    setwd(fileLoc)

    #Write start time to the summary file
    suppressWarnings(write(paste0(dateStamp, " Start"), file = paste0(dateStamp, "_combinedDada.txt"), append = FALSE))

    #Get the all files in the selected file folder with full paths of the files
    files <- as.data.frame(list.files(path = fileLoc, pattern = "*[.]*", full.names = TRUE))
    # Get the local paths
    files[,2] <- list.files(path = fileLoc, pattern = "*[.]*")
    #Get all files with extension '_Merge.tsv' OR '_MergeFwdRev.tsv' OR '_Forward.tsv'
    files <- rbind(files[grepl("_Merge.tsv", files[,2]),],files[grepl("_MergeFwdRev.tsv", files[,2]),],files[grepl("_Forward.tsv", files[,2]),])

    # Get the names of the files
    files[,3] <- gsub("*.tsv","",files[,2])
    # Get a unique number for each file to associate with the combined output results
    files[,4] <- c(1:nrow(files))

    #Add in headers for the file table
    colnames(files)<-c("File_Path", "File", "File_Name", "File_Unique_ID")

    #Write the files being combined
    suppressWarnings(write(paste0("The files being combined with an associated unique identifying number (File_Unique_ID column) used in the uniqueID column of the output file..."), file = paste0(dateStamp, "_combinedDada.txt"), append = TRUE))
    suppressWarnings(write.table(files, file = paste0(dateStamp, "_combinedDada.txt"), append = TRUE, na = "NA", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"))

    #Audit line
    if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 2")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 2"), file = auditFile, append = TRUE))}

    if(nrow(files)<2){

      print("********************************************************************************")
      print("There is only one file in the target directory for the file type and so ")
      print("the combination of files is not needed.")
      print("********************************************************************************")

    }else{

      #Create a flag for the first loop
      flag = TRUE

      for(records in 1:nrow(files)){

        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 3")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 3"), file = auditFile, append = TRUE))}

        print("********************************************************************************")
        print(paste0("Starting file ", files[records,3], " at time...", Sys.time()))
        print("********************************************************************************")
        suppressWarnings(write(paste0("Starting file ", files[records,3], " at time...", Sys.time()), file = paste0(dateStamp, "_combinedDada.txt"), append = TRUE))

        #Read in the files for this loop
        loopResults <- read.delim(files[records,1], header = TRUE, check.names = FALSE)

        #subset the results to only include values with the included min length
        loopResults <- loopResults[loopResults$Length >= minLen, ]

        #Remove the unique ID variable
        loopResults<-loopResults[,-1]

        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 4")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 4"), file = auditFile, append = TRUE))}

        #Using the dataframe without quality metrics collapse all columns up to
        #the sample columns into a single variable
        mergeVariable<-paste0(loopResults$Sequence,
                              loopResults$Length,
                              loopResults$Results)

        #Create a loopResultsData table with sequence, length and results in columns
        loopResultsData<-data.frame(mergeVariable=mergeVariable,
                                    Sequence = loopResults$Sequence,
                                    Length = loopResults$Length,
                                    Results = loopResults$Results)

        #Remove the columns that made up the mergeVariable - Seq, Length, Results
        loopResults<-loopResults[,-c(1:3)]

        #Bind the combined column back onto the dataframe
        loopResults<-cbind(data.frame(mergeVariable=mergeVariable),loopResults)

        #update the unique sample names to add on the file name
        colnames(loopResults)[2:ncol(loopResults)] <- paste0(files[records, 3], "_SEP_", colnames(loopResults[2:ncol(loopResults)]))

        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 5")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 5"), file = auditFile, append = TRUE))}

        if (flag == TRUE){

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 6")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 6"), file = auditFile, append = TRUE))}

          #Make the total results equal the loop results
          totalResults <- loopResults

          #Make a data results data.frame
          totalResultsData<-loopResultsData

          #Reset the flag to add to the total results for the next loop
          flag = FALSE

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 7")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 7"), file = auditFile, append = TRUE))}

        }else{

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 8")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 8"), file = auditFile, append = TRUE))}

          #Merge the total and loop data
          totalResults <- merge(totalResults[,c(1:ncol(totalResults))], loopResults[,c(1:ncol(loopResults))], by = "mergeVariable", all=TRUE )

          #Createa a totalResultsData data frame
          totalResultsData<-rbind(totalResultsData, loopResultsData)

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 9")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 9"), file = auditFile, append = TRUE))}

        }#End of the if/else flag
      }#Closing off the loop going through the records loop
    }#End of the check if there are more than one file so combo makes sense
  }#End of the if checking that a file location value was submitted or chosen from the pop up

  #Audit line
  if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 10")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 10"), file = auditFile, append = TRUE))}

  #Change all NA to 0 in the sample section
  totalResults[,c(2:ncol(totalResults))][is.na(totalResults[,c(2:ncol(totalResults))])] <- 0

  #reduce all of the unique mergeVariables and sum the values in the sample columns
  totalResults <- aggregate(.~mergeVariable, totalResults, sum)

  #Reduce all of the totalResultsData to unique rows
  totalResultsData<-unique(totalResultsData)

  #Merge totalResults and totalResultsData
  totalResults <- merge(x = totalResultsData, y = totalResults, all=TRUE )

  #Remove the mergeVariable from the totalResults data frame
  totalResults<-totalResults[,-1]

  #Add uniqueID back onto the dataframe before printing to file
  totalResults <- data.frame(uniqueID = paste0("C",c(1:nrow(totalResults))), totalResults, check.names=FALSE)

  #Audit line
  if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 11")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 11"), file = auditFile, append = TRUE))}

  #Change the .(which were forced to change from the original naming convention of -) to - for output
  colnames(totalResults)<-gsub("\\.","-",colnames(totalResults))

  #print the results fo the combined datasets to the file
  suppressWarnings(write.table(totalResults, file = paste0(dateStamp, "_combinedDada.tsv"), append = FALSE, na = "NA", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"))

  #print to fasta file
  totalResults<-totalResults[,1:2]
  totalResults[,1]<-paste0(">",totalResults[,1])
  suppressWarnings(write.table(totalResults, file = paste0(dateStamp, "_combinedDada.fas"), append = FALSE, na = "", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n"))

  print(paste0(startTime, " end time ", Sys.time()))
  suppressWarnings(write(paste0(startTime, " end time ", Sys.time()), file = paste0(dateStamp, "_combinedDada.txt"), append = TRUE))

  #Audit line
  if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 12")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 12"), file = auditFile, append = TRUE))}

}#End of function
