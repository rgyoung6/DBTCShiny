# Written by Rob Young at the University of Guelph in Ontario Canada, Sept, 2023
# ******************************************************************************
#Roxygen2 Documentation:

#' @export
#'
#' @title BLAST Query File Against Local Database
#'
#' @author Robert G. Young
#'
#' @description
#' This function takes fasta files as input along with a user selected NCBI
#' formatted library to BLAST sequences against. The outcome of the function are
#' two files, a BLAST run file and a single file containing all of the BLAST
#' results in tab delimited format (Note: there are no headers but the columns
#' are, query sequence ID, search sequence ID, search taxonomic ID, query to
#' sequence coverage, percent identity, search scientific name, search common
#' name, query start, query end, search start, search end, e-value.
#'
#' @details
#' The user input provides a location for the BLAST database you would like to
#' use by selecting a file in the target directory. Then provide the location
#' of the query sequence files by indicating a file in a directory that contains
#' the fasta files. Provide the path for the blast+ blastn program. Finally, provide
#' the minimum query sequence length to BLAST (Default = 100), the depth of the BLAST
#' returned results (default = 250), and finally the number of cores to process
#' the function (default = 1, Windows implementation will only accept this value
#' as 1).
#'
#' @examples
#' \dontrun{
#' seq_BLAST()
#' seq_BLAST(databasePath = NULL, querySeqPath=NULL,  blastnPath="blastn",
#' minLen = 100, BLASTResults=200, numCores=1)
#' }
#'
#' @param databasePath The location of a file in a directory where the desired
#' BLAST database is located.
#' @param querySeqPath The local path for the directory containing all of the
#' fasta files wishing to be BLASTed
#' @param blastnPath The location of the NCBI blast+ blastn program (default = blastn).
#' @param minLen The minimum length of the sequences that will be BLASTed (default = 100).
#' @param BLASTResults The number of returned results, or the depth of the reported
#' results, saved from the BLAST (default = 250).
#' @param numCores The number of cores used to run the function (default = 1,
#' Windows systems can only use a single core).
#'
#' @returns
#' Two files are produced from this function, a BLAST run file and a BLAST results
#' file for each of the fasta files in the target directory.
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
#' taxon_assign()
#' combine_assign_output()
#' reduce_taxa()
#' combine_reduced_output()

##################################### BLAST FUNCTION ##############################################################
seq_BLAST <- function(databasePath = NULL, querySeqPath=NULL,  blastnPath="blastn", minLen = 100, BLASTResults=200, numCores=1){

  #If there are issues and I need to audit the script make this 1
  auditScript=0

  #Get the initial working directory
  start_wd <- getwd()
  on.exit(setwd(start_wd))

  #Printing the start time
  print(paste0("Start time...", Sys.time()))
  startTime <- paste0("Start time...", Sys.time())
  dateStamp <- paste0(format(Sys.time(), "%Y_%m_%d_%H%M"))

  if(is.null(databasePath)){
    #Get the directory with the BLAST data base and set the working directory to the location of the database
    print(paste0("Select a file in the folder with the NCBI database you would like to use."))
    databasePath <- file.choose()
  }

  if(is.null(querySeqPath)){
    print(paste0("Select a file in the folder with the fasta files you would like to BLAST."))
    querySeqPath<-file.choose()
  }


  if(is.null(blastnPath)){
    print(paste0("Select the path for the blastn command."))
    blastnPath<-file.choose()
  }

  if(is.null(databasePath)){
    print("********************************************************************************")
    print("A database location is required (databasePath).")
    print("Please rerun the function and provide a character string for the databasePath")
    print("argument or select a location when prompted.")
    print(paste0("Current database name (databasePath) is: ", databasePath))
    print("********************************************************************************")
  } else if(is.null(BLASTResults)){
    print("********************************************************************************")
    print("An integer value for the desired maximum number of BLAST returned results")
    print("(BLASTResults) is required.")
    print("Please rerun the function and provide a value for the BLASTResults argument.")
    print(paste0("Current database name (BLASTResults) is: ", BLASTResults))
    print("********************************************************************************")
  } else if (is.null(querySeqPath)){
    print("********************************************************************************")
    print("The query sequence path was not provided. Please select a file and run the function again.")
    print("********************************************************************************")
  } else if (grepl(" ",querySeqPath) | grepl(" ",databasePath) | grepl(" ",blastnPath)){
    print("********************************************************************************")
    print("Error: One or more of the file paths contains a space in the naming convention. Please change the naming and try again.")
    print("********************************************************************************")
  }else{

    querySeqPath<-dirname(querySeqPath)

    #Set the working directory to the location of the BLAST database
    setwd(dirname(databasePath))

    #Audit line
    if(auditScript>0){
      auditFile <- paste0(querySeqPath,"/", format(Sys.time(), "%Y_%m_%d_%H%M"), "_audit.txt")
      write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 1"), file = auditFile, append = FALSE)
      print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 1"))
    }

    #Setting up the database name
    databaseName <- list.files(path=dirname(databasePath), pattern = "*[.]nto$")

    if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 1A")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 1A"), file = auditFile, append = TRUE))}

    databaseName <- sub("\\..*","", databaseName)

    if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 1B")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 1B"), file = auditFile, append = TRUE))}

    if(length(databaseName) == 0 ){
      print("********************************************************************************")
      print("The database location (databasePath) does not appear to contain a BLAST formated database.")
      print("Please rerun the function and provide an appropriate character string for the databasePath")
      print("argument or select a proper location when prompted.")
      print(paste0("Current database name (databaseName) is: ", databaseName))
      print("********************************************************************************")

    }else {

      #Audit line
      if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 2")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 2"), file = auditFile, append = TRUE))}

      #Get the list of files in the target directory to loop through
      pathList <- list.files(path=querySeqPath, pattern = "*[.]([Ff][Aa][Ss]$)|([Ff][Aa][Ss][Tt][Aa]$)", full.names = TRUE)

      fileList <- list.files(path=querySeqPath, pattern = "*[.]([Ff][Aa][Ss]$)|([Ff][Aa][Ss][Tt][Aa]$)")
      fileNames <- sub("\\..*","", as.vector(fileList))

      if(length(pathList)==0){

        print("********************************************************************************")
        print("There are no fasta (*.FAS or *.FASTA) files in the target directory.")
        print("Please rerun the function and provide a folder with fasta files.")
        print(paste0("The current path is... ",querySeqPath ))
        print("********************************************************************************")

      }else{

        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 3")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 3"), file = auditFile, append = TRUE))}

        #Create a temp file to hold the reduced dataset after filtering
        temp_file <- tempfile(fileext = ".fas")

        #Build the run file.
        for(filesInFolder in 1:length(pathList)){

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 4")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 4"), file = auditFile, append = TRUE))}

          print(paste0("Beginning the BLASTing of file ",pathList[filesInFolder], " at - ", Sys.time()))

          #Read in the data from the target file
          seqTable <- read.delim(pathList[filesInFolder], header = FALSE)

          if(!grepl(">", seqTable[3,1], fixed = TRUE)){

            print("********************************************************************************")
            print("The submitted fasta file is not in the single line nucleotide format which is ")
            print("needed for this script. Please correct the format and rerun this script")
            print(paste0("for file ",pathList[filesInFolder]))
            print("********************************************************************************")

          }else{

            #Audit line
            if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 5")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 5"), file = auditFile, append = TRUE))}

            #taking the read in file and changing from Fasta to tab delimited
            seqTableTemp <- data.frame(Header = (seqTable[seq(from = 1, to = nrow(seqTable), by = 2), 1]))
            seqTableTemp["Sequence"] <- seqTable[seq(from = 2, to = nrow(seqTable), by = 2), 1]
            seqTable<-seqTableTemp
            #Remove the seqTableTemp to free memory
            rm(seqTableTemp)

            #Remove the records with fewer than submitted lower limit of nucleotides
            seqTable <- seqTable[nchar(as.character(seqTable[,2])) > minLen,]

            if(nrow(seqTable)==0){

              print("********************************************************************************")
              print("The submitted fasta file does not have any records to BLAST after applying the ")
              print("length filter. Please check the data and run it again.")
              print("********************************************************************************")

            }else{

              #Audit line
              if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 6")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 6"), file = auditFile, append = TRUE))}

              #save the reduced fasta to a temp file
              write.table(seqTable, file = temp_file, append = FALSE, na = "", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n")

              #Build the BLAST command
              BLASTCmdString<- paste0("\"",blastnPath, "\" -db ", databaseName, " -query \"", temp_file, "\" -max_target_seqs ", BLASTResults, " -outfmt \"6 qseqid sseqid staxid qcovs pident ssciname scomname qstart qend sstart send evalue\" -num_threads ", numCores, " -out \"", querySeqPath,"/", fileNames[filesInFolder], "_BLAST_", databaseName,"_",dateStamp, ".tsv\"")

              #Build the file to run based on the operating system
              if (.Platform$OS.type == "windows"){

                #Audit line
                if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 7")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 7"), file = auditFile, append = TRUE))}

                #Windows command file
                blastCommandFile <- paste0(fileNames[filesInFolder], "_BLAST_", databaseName,"_", dateStamp, ".bat")
                write(BLASTCmdString, file = blastCommandFile, append = FALSE)

                #Run the BLAST command in a system command
                system(blastCommandFile)

                #After running the command file move the file to the location of the fasta files being BLASTed
                file.rename(from = blastCommandFile,  to = paste0(querySeqPath, "/", blastCommandFile))

              } else{

                #Audit line
                if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 8")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 8"), file = auditFile, append = TRUE))}

                #linux and Mac OS command file
                blastCommandFile <- paste0(fileNames[filesInFolder], "_BLAST_", databaseName,"_", dateStamp, ".sh")

                #Initialize the file that will be run for the BLAST
                write("#!/biin/sh", file = blastCommandFile, append = FALSE)
                write("\n", file = blastCommandFile, append = TRUE)
                write(BLASTCmdString, file = blastCommandFile, append = TRUE)

                #Run the BLAST command in a system command
                BLASTOutput<-system(paste0("bash './", blastCommandFile, "'"))

                tryCatch(
                  expr = {
                    #After running the command file move the file to the location of the fasta files being BLASTed
                    file.rename(from = blastCommandFile,  to = paste0(querySeqPath, "/", blastCommandFile))
                  },
                  error = function(e){
                    print(paste0("Error - unable to move the BLAST run file (", blastCommandFile, " due to persmissions."))
                  },
                  warning = function(w){
                    print(paste0("Warning - unable to move the BLAST run file (", blastCommandFile, " due to persmissions."))
                  }
                )

                #Audit line
                if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 9")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 9"), file = auditFile, append = TRUE))}

              }#End of if checking for the operating system
            }#Close if where there are no fasta files in the selected folder
          }#Close the check if there are records in the file after length filter
        }#End of the checking the format of the file
      }#End of the loop going through the fasta files to be BLASTed
    }#End of the check for the length of the databaseName as obtained from the database location.
  }#Closing the if checking for the submitting arguments

  print(paste0("BLAST Complete - Started at ", startTime, " Ended at ", Sys.time()))

}#End of the function
