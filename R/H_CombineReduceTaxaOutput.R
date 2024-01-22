# Written by Rob Young at the University of Guelph in Ontario Canada, Sept, 2023
# ******************************************************************************
#Roxygen2 Documentation:

#' @export
#'
#' @title Combine Reduce Taxa Files for the Same Biological Samples using Different Markers
#'
#' @author Robert G. Young
#'
#' @description
#' This function takes a file selection and then uses all 'taxaReduced' files
#' in that directory and combines them into a single taxa table file with presence
#' absence results._CombineTaxaReduced.tsv
#'
#' @details
#' The User Input: This function requires a file in a directory where all 'taxaReduced'
#' files in that directory will be combined. The output format will be a taxa table
#' with all taxa from all files combined into a single table with presence absence
#' (0 or 1) results. The value metrics for the identification of the taxa from each
#' combined file will remain in a column with the parenthetical results from the 'taxaReduced'
#' files ("Num_Rec", "Coverage", "Identity", "Max_eVal").
#'
#' @examples
#' \dontrun{
#' combine_reduced_output()
#' combine_reduced_output(fileLoc = NULL)
#' }
#'
#' @param fileLoc The location of a file in a directory where all of the 'taxa_assign'
#' and/or 'combined_taxa_assign' files are located.
#'
#' @returns
#' This function produces a single 'YYYY_MM_DD_HHMM_CombineTaxaReduced' file and
#' associated summary file in the target directory.
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
#' combine_assign_output()
#' reduce_taxa()
#'
################ COMBINE TWO REDUCED TAXA FILES INTO A SINGLE FILE #############
combine_reduced_output <- function(fileLoc = NULL, presenceAbsence = TRUE){

  #Get the initial working directory
  start_wd <- getwd()
  on.exit(setwd(start_wd))

  #load in the files list
  if (is.null(fileLoc)){
    print(paste0("Select a file in the file folder with the reduced taxa files files you would like to combine (extension '_taxaReduced_YYYY_MM_DD_HHMM.tsv')."))
    print("********** NOTE: all files being combined should have used the exact same samples **********")
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

    #get the directory of interest
    fileLoc <- dirname(fileLoc)

    #Set the working directory
    setwd(fileLoc)

    #Write start time to the summary file
    suppressWarnings(write(paste0(dateStamp, " Start"), file = paste0(dateStamp, "_CombineTaxaReduced.txt"), append = FALSE))

    #Get the all files in the selected file folder with full paths of the files
    files <- as.data.frame(list.files(path = fileLoc, pattern = "*[.]*", full.names = TRUE))
    # Get the local paths
    files[,2] <- list.files(path = fileLoc, pattern = "*[.]*")
    #Get all files with the '_taxaReduced' string
    files <- files[grepl("_taxaReduced_.*", files[,2]),]
    # Get the names of the files
    files[,3] <- gsub("_taxaReduced_.*","",files[,2])
    # Get a unique number for each file to associate with the combined output results
    files[,4] <- c(1:nrow(files))

    #Add in headers for the file table
    colnames(files)<-c("File_Path", "File", "File_Name", "File_Unique_ID")

    #Write the files being combined
    suppressWarnings(write(paste0("The files being combined with an associated unique identifying number (File_Unique_ID column) used in the uniqueID column of the output file..."), file = paste0(dateStamp, "_CombineTaxaReduced.txt"), append = TRUE))
    suppressWarnings(write.table(files, file = paste0(dateStamp, "_CombineTaxaReduced.txt"), append = TRUE, na = "NA", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"))

    if(nrow(files)<2){

      print("********************************************************************************")
      print("There are one or fewer files in the target directory for the file type and so ")
      print("the combination of files is not needed.")
      print("********************************************************************************")

    }else{

      #Create a flag for the first loop
      flag = TRUE

      for(records in 1:nrow(files)){

        print(paste0("Starting file ", files[records,3], " at time...", Sys.time()))
        suppressWarnings(write(paste0("Starting file ", files[records,3], " at time...", Sys.time()), file = paste0(dateStamp, "_CombineTaxaReduced.txt"), append = TRUE))

        #Read in the files for this loop
        loopResults <- read.delim(files[records,1], header = TRUE, check.names=FALSE)

        #Remove the bracket values
        mergeVariable <- data.frame(lapply(loopResults[,1:7], function(x) {gsub("\\s*\\([^\\)]+\\)", "", x)}), check.names=FALSE)

        #Change out the Greater than 20 to NA
        mergeVariable <- data.frame(lapply(mergeVariable, function(x) {gsub("Greater than 20", NA, x)}), check.names=FALSE)

        #Change all cells in the dataframe  with a comma to NA
        mergeVariable <- lapply(mergeVariable, function(x) {x[grepl(",", x, fixed = TRUE)] <- NA; x})

        #Create a merge variable to sort the data and combine for final output
        mergeVariable <- as.data.frame(paste0(mergeVariable$superkingdom, "|", mergeVariable$phylum, "|", mergeVariable$class, "|", mergeVariable$order, "|", mergeVariable$family, "|", mergeVariable$genus, "|", mergeVariable$species), drop=FALSE, check.names=FALSE)

        #Remove the ends of the strings after the first instance of NA
        mergeVariable <- as.data.frame(lapply(mergeVariable, function(x) {gsub("\\|NA.*","",x)}), drop=FALSE, check.names=FALSE)

        #Create two dataframes one with the taxonomic information and the mergeVariable
        #and one with the sample records and the mergeVariable
        loopResultsRecords<-as.data.frame(cbind(mergeVariable, loopResults[,20:ncol(loopResults)]), check.names=FALSE)
        loopResultsTaxa<-as.data.frame(cbind(mergeVariable, MarkerResults = paste0(loopResults[,10]," - ",loopResults[,11])), check.names=FALSE)

        #Change the name of the mergeVariable column
        colnames(loopResultsRecords)[1]<-"mergeVariable"
        colnames(loopResultsTaxa)[1]<-"mergeVariable"

        #Add the file name to the
        colnames(loopResultsTaxa)[2] <- paste0(files[records, 3], "_", colnames(loopResultsTaxa[2]))

        if (flag){

          #Make the total results equal the loop results
          totalLoopResultsRecords <- loopResultsRecords
          totalLoopResultsTaxa <- loopResultsTaxa

          #Reset the flag to add to the total results for the next loop
          flag = FALSE

        }else{

          #Make the total results equal the loop results
          totalLoopResultsRecords <- rbind.fill(totalLoopResultsRecords, loopResultsRecords)
          totalLoopResultsTaxa<-merge(totalLoopResultsTaxa, loopResultsTaxa, by = "mergeVariable", all=TRUE)

        }#end of if/else flag

        #Change all NA to 0 in the sample section
        totalLoopResultsRecords[,c(2:ncol(totalLoopResultsRecords))][is.na(totalLoopResultsRecords[,c(2:ncol(totalLoopResultsRecords))])] <- 0

        #Sum all values in the records data frame using the aggregate function
        totalLoopResultsRecords<-stats::aggregate(totalLoopResultsRecords[,2:ncol(totalLoopResultsRecords)], by=list(totalLoopResultsRecords$mergeVariable), FUN=sum, drop=TRUE)

        #We lose the column name in aggregate so we add it back in here...
        colnames(totalLoopResultsRecords)[1]<-"mergeVariable"

        #Merge the two dataframes with records and unique taxa into a single dataframe
        totalResults<-merge(totalLoopResultsTaxa,totalLoopResultsRecords, by = "mergeVariable", all = TRUE)
      }#Closing off the loop going through all the files
    }#End of the check if there are more than one file so combo makes sense
  }#End of the if checking that a file location value was submitted or chosen from the pop up

  #Finalizing the format for output
  #Add a dummy string at the top of mergeVariable and add headers
  #without file names to the top of the other columns
  totalResults<-rbind(c("King|Phy|Class|Order|Fam|Gen|Sp", colnames(totalResults)[2:ncol(totalResults)]),totalResults)

  # split the strings
  temp <- strsplit(totalResults$mergeVariable, split="[|]")
  # maximum length of the list items
  maxL <- max(sapply(temp, length))
  # construct data.frame with NAs as fills
  totalResults <- data.frame(cbind(do.call(rbind, lapply(temp, function(i) c(i, rep(NA, maxL-length(i)))))), totalResults[,2:ncol(totalResults)], check.names=FALSE)

  #Remove the first row which is the combining unique identifier
  totalResults<-totalResults[-1,]

  #Add column names to the newly created columns
  colnames(totalResults)[1:7]<-c("superkingdom", "phylum", "class", "order", "family", "genus", "species")

  #Change all values to binary results 0 and 1 to reflect presence absence because you can't look at num of reads after combining primers
  #Only do the replacement for the columns after the taxa names and the files so "superkingdom", "phylum", "class", "order",
  #"family", "genus", "species" is 7 plus one and then add the num of rows in the file variable

  # If the presenceAbsence is TRUE then change all of the read numbers to 1 or 0
  if (presenceAbsence == TRUE){

    #Subset and get the columns of interest
    ATempVariable<-totalResults[,as.numeric(8+nrow(files)):ncol(totalResults), drop = FALSE]

    #Change all positive values to 1
    ATempVariable[ATempVariable > 0] <- 1

    #Add the changed values back onto the final data frame to be printed to file
    totalResults<-cbind(totalResults[,1:as.numeric(7+nrow(files))],ATempVariable)

  }

  #Write the results to the file
  suppressWarnings(write.table(totalResults, file = paste0(dateStamp, "_CombineTaxaReduced.tsv"), append = FALSE, na = "NA", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"))

  print(paste0(startTime, " End time ", Sys.time()))
  suppressWarnings(write(paste0(startTime, " End time ", Sys.time()), file = paste0(dateStamp, "_CombineTaxaReduced.txt"), append = TRUE))

}#End of function
