# Written by Rob Young at the University of Guelph in Ontario Canada, Sept, 2023
# ******************************************************************************
#Roxygen2 Documentation:

#' @export
#'
#' @title Dada Implement
#'
#' @author Robert G. Young
#'
#' @description
#' This function requires a main directory containing a folder(s) representing
#' sequencing runs which in-turn contain fastq files (the location of one of the
#' fastq files in one of the sequencing run folders is used as an input
#' argument). All sequencing folders in the main directory need to represent
#' data from sequencing runs that have used the same primers and protocols.
#' Output from this function includes all processing files and final main output
#' files in the form of fasta files and amplicon sequencing variant (ASV)
#' tables.
#'
#' @details
#' Two file types are required as input for the dada_implement() function.
#' The first are the fastq files in the appropriate folder structure (see below)
#' and the second is a file containing the primers used for the amplification of
#' the sequence reads.
#'
#' Fastq File Folder Structure
#'
#'          Parent Directory
#'                 |
#'                 |
#'          -----------------
#'          |               |
#'          |               |
#'    Run1 Directory     Run2 Directory
#'    -Fastq             -Fastq
#'    -Fastq             -Fastq
#'    ...                ...
#'
#' Format of the primer file
#'
#'  | Forward              | Reverse                          |
#'  | AGTGTGTAGTGATTG      | CGCATCGCTCAGACTGACTGC            |
#'  | GAGCCCTCGATCGCT      | GGTCGATAGCTACGCGCGCATACGACT      |
#'  |                      | GGTTCACATCGCATTCAT               |
#'
#'
#' @examples
#' \dontrun{
#' dada_implement()
#' dada_implement(runFolderLoc = NULL, primerFile = NULL,fwdIdent = "_R1_001",
#' revIdent = "_R2_001",unidirectional = FALSE, bidirectional = TRUE
#' maxPrimeMis = 2, fwdTrimLen = 0, revTrimLen = 0,maxEEVal=2, truncQValue = 2,
#' truncLenValueF = 0, truncLenValueR = 0,error = 0.1, nbases = 1e80,
#' maxMismatchValue = 0, minOverlapValue = 12,trimOverhang = FALSE,
#' minFinalSeqLen = 100)
#' }
#'
#' @param runFolderLoc Select a directory that contains the run folders with the
#' fastq files.
#' @param primerFile Select a file with the primers for this analysis.
#' @param fwdIdent Forward identifier naming substring.
#' @param revIdent Reverse identifier naming substring.
#' @param bidirectinoal Selection to process paired forward and reverse sequence
#' for analysis (Default TRUE).
#' @param unidirectional Seection to process files independently (Default FALSE).
#' @param maxPrimeMis Maximum number of mismatches allowed when pattern matching
#' trimming the primers from the ends of the reads for the ShortRead trimLRPatterns()
#' function (Default maxPrimeMis = 2).
#' @param fwdTrimLen Select a forward trim length for the Dada filterAndTrim()
#' function (Default fwdTrimLen = 0).
#' @param revTrimLen Select a reverse trim length for the Dada filterAndTrim()
#' function (Default revTrimLen = 0).
#' @param maxEEVal Maximum number of expected errors allowed in a read for the
#' Dada filterAndTrim() function (Default maxEEVal = 2).
#' @param truncQValue Truncation value use to trim ends of reads, nucleotides
#' with quality values less than this value will be used to trim the remainder
#' of the read for the Dada filterAndTrim() function (Default truncQValue = 2).
#' @param truncLenValueF Dada forward length trim value for the Dada filterAndTrim()
#' function. This function is set to 0 when the pattern matching trim function is
#' enabled (Default truncLenValueF = 0).
#' @param truncLenValueR Dada reverse length trim value for the Dada filterAndTrim()
#' function. This function is set to 0 when the pattern matching trim function is
#' enabled (Default truncLenValueR = 0).
#' @param error Percent of fastq files used to assess error rates for the Dada
#' learnErrors() function (Default error = 0.1).
#' @param nbases The total number of bases used to assess errors for the Dada
#' learnErrors() function (Default nbases = 1e80) NOTE: this value is set very
#' high to get all nucleotides in the error persent file subset. If the error
#' is to be assessed using total reads and not specific fastq files then set
#' the error to 1 and set this value to the desired number of reads.
#' @param maxMismatchValue Maximum number of mismatches allowed when merging
#' two reads for the Dada mergePairs() function (Default maxMismatchValue = 2).
#' @param minOverlapValue Minimum number of overlapping nucleotides for the
#' forward and reverse reads for the Dada mergePairs() function (Default
#' minOverlapValue = 12).
#' @param trimOverhang Trim merged reads past the start of the complimentary
#' primer regions for the Dada mergePairs() function (Default trimOverhang = FALSE).
#' @param minFinalSeqLen The minimum final desired length of the read (Default
#' minFinalSeqLen = 100).
#'
#' @returns
#' The output from this function includes four folders.
#' A_Qual - Contains quality pdf files for the input fastq files.
#' B_Filt - Contains dada filtered fastq files and a folder with the end trimmed
#' fastq files before quality filtering.
#' C_FiltQual - Contains quality pdf files for the filtered fastq files.
#' D_Output - This folder contains output files including and analysis summary,
#' an analysis summary table of processing values, foward and reverse error assessments,
#' and finally the output ASV and fasta files of obtained sequences.                                                                                                   -TotalTable.tsv
#'
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
#' combine_dada_output()
#' make_BLAST_DB()
#' seq_BLAST()
#' taxon_assign()
#' combine_assign_output()
#' reduce_taxa()
#' combine_reduced_output()

##################################### dada implementation FUNCTION ##############################################################
dada_implement <- function(runFolderLoc = NULL, primerFile = NULL,
                           fwdIdent = "_R1_001", revIdent = "_R2_001",
                           unidirectional = FALSE,
                           bidirectional = TRUE,
                           printQualityPdf = TRUE,
                           maxPrimeMis = 2, fwdTrimLen = 0, revTrimLen = 0,
                           maxEEVal=2, truncQValue = 2,
                           truncLenValueF = 0, truncLenValueR = 0,
                           error = 0.1, nbases = 1e80,
                           maxMismatchValue = 0, minOverlapValue = 12,
                           trimOverhang = FALSE, minFinalSeqLen = 100){

  #If there are issues and I need to audit the script make this 1
  auditScript=0

  if(unidirectional == FALSE & bidirectional == FALSE){
    print("One or both of the unidirectional or bidirectional selections need to be TRUE!")
  }else{

    #Get the initial working directory
    start_wd <- getwd()
    on.exit(setwd(start_wd))

    if(is.null(primerFile)){
      # prompting to choose the file of interest with the tab delimited primer info
      n <- substr(readline(prompt="Choose the file with the primers of interest. Hit enter key to continue..."),1,1)
      primerFile <- file.choose()
    }

    if(is.null(runFolderLoc)){
      # prompting to choose the file of interest with the tab delimited primer info
      n <- substr(readline(prompt="Select a fastq file in one of the run folders in the directory of interest (NOTE: all run folders with fastq data in the parent directory will be processed by DBTC. If this is not what you want please rearrange your folder structure)."),1,1)
      runFolderLoc <- dirname(dirname(file.choose()))
    }

    #Audit line
    if(auditScript>0){
      auditFile <- paste0(runFolderLoc,"/", format(Sys.time(), "%Y_%m_%d_%H%M"), "_audit.txt")
      print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 1"))
      suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 1"), file = auditFile, append = FALSE))
    }

    #Read in the primers from the primer file.
    primerFileData <- read.delim(primerFile, header = TRUE)

    #Set a list for the forward and reverse primers
    fwdPrimerList<-primerFileData[,1]
    fwdPrimerList<-fwdPrimerList[!(fwdPrimerList=="")]
    revPrimerList<-primerFileData[,2]
    revPrimerList<-revPrimerList[!(revPrimerList=="")]

    #Get the maximum length of the forward and reverse primers
    fwdPrimeMaxLen<-as.numeric(max(nchar(fwdPrimerList)))
    revPrimeMaxLen<-as.numeric(max(nchar(revPrimerList)))

    #Set up the Dada forward and reverse length trim to remove poor quality nucleotides
    truncLenValue = c(truncLenValueF, truncLenValueR)

    #Set the number of allowed unknown nucleotides in the analysis. Dada requires no
    #unknown nucleotides so it is set to 0
    maxNVal=0

    if(fwdPrimeMaxLen == 0){
      primerFileCheck = FALSE
    }else if(is.na(revPrimeMaxLen) & bidirectional==TRUE){
      primerFileCheck = FALSE
    }else{
      primerFileCheck = TRUE
    }

    if(is.null(runFolderLoc)){
      print("A target directory containing folders with runs and in these run folders the fastq files is necessary to run the program")
      print("These MiSeq runs need to be for a single set of primers. If another set of primers is desired then this script will need to be run multiple times.")
      print(paste0("Current file folder location: ", runFolderLoc))
    } else if(isFALSE(primerFileCheck) && (fwdTrimLen < 1 || revTrimLen < 1)){
      print("Primer sequences or primer end trimming lengths are required to run the program.")
      print(paste0("Please verify the current values included as arguments or in the Primer File: "))
      print(paste0("primerFile: ", primerFile))
      print(paste0("fwdTrimLen: ", fwdTrimLen))
      print(paste0("revTrimLen: ", revTrimLen))
    } else{

      #Audit line
      if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 2")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 2"), file = auditFile, append = TRUE))}

      #Get the list of files in the target directory to loop through
      runFolderLoc <- list.dirs(path = runFolderLoc, full.names = TRUE, recursive = FALSE)

      #Audit line
      if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 3")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 3"), file = auditFile, append = TRUE))}

      for(runCounter in 1:length(runFolderLoc)){

        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 4")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 4"), file = auditFile, append = TRUE))}

        workLoc <- runFolderLoc[runCounter]
        setwd(workLoc)

        #Extracting file folder name to append to the results
        dirName <- sub(".*/","", runFolderLoc[runCounter])

        #Printing the start time
        print("###########################################################################################")
        print(paste0(runFolderLoc[runCounter], " - Start time...", Sys.time()))
        startTime <- paste0("Start time...", Sys.time())
        dateStamp <- paste0(format(Sys.time(), "%Y_%m_%d_%H%M"), "_")

        #Initiating the list of files in the folders
        pathList <- list.files(path = workLoc, pattern = "*[.][Ff][Aa][Ss][Tt][Qq][.]*", full.names = TRUE)
        fileList <- list.files(path = workLoc, pattern = "*[.][Ff][Aa][Ss][Tt][Qq][.]*")
        fileNames <- sub("\\..*","", as.vector(fileList))
        fileNamesNoDir <- sub(fwdIdent,"", fileNames)
        fileNamesNoDir <- sub(revIdent,"", fileNamesNoDir)

        #Create a table with the files of interest
        filesTbl <- as.data.frame(cbind(pathList, fileList, fileNames, fileNamesNoDir))

        #Create a subfolder for filtered data files
        filteredFolder <-paste0(workLoc, "/", dateStamp, dirName,"_B_Filt")
        dir.create(filteredFolder)

        #Create a subfolder for filtered data files
        filteredSubFolder <-paste0(filteredFolder, "/Primer_Trim")
        dir.create(filteredSubFolder)

        #Create a subfolder for output files
        outFolder <- paste0(workLoc, "/", dateStamp, dirName,"_D_Output")
        dir.create(outFolder)

        if (printQualityPdf == TRUE){

          #Create a subfolder for pre-filter quality figures
          qualityFolder <- paste0(workLoc, "/", dateStamp, dirName,"_A_Qual")
          dir.create(qualityFolder)

          #Create a subfolder for filtered quality figures
          filterQualFolder <- paste0(workLoc, "/", dateStamp, dirName,"_C_FiltQual")
          dir.create(filterQualFolder)

        }

        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 5")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 5"), file = auditFile, append = TRUE))}










        #writing the start time to the log file here
        suppressWarnings(write(paste0(dateStamp, "Start"), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = FALSE))
        suppressWarnings(write(paste0(""), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0("Here are the variables used in the analyses..."), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0("********************************************************************************"), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0("Platform: ", .Platform$OS.type),file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0("The working directory: ", runFolderLoc[runCounter]), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0("Directional Processing: Uni - ", unidirectional, " ; bi - ", bidirectional),file=paste0(outFolder,"/", dateStamp,dirName, "_dadaSummary.txt"),append = TRUE))
        suppressWarnings(write(paste0("Parsing file directional identifiers: Forward(fwdIdent) - ", fwdIdent, "; Reverse(revIdent) - ",  revIdent), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0("Quality truncation value (truncQ): ", truncQValue), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0("Minimum read length for the Dada2 filterandtrim() function (truncLenValue): ", truncLenValue), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0("Forward and reverse trim length values used in the Dada2 filterandtrim() function: Forward (fwdTrimLen) = ", fwdTrimLen, "; Reverse (revTrimLen) = ", revTrimLen), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0("File used containing the forward and reverse trim sequences data: ", primerFile), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0("Forward and reverse trim sequences used in the trim() function from the Insect package: \nForward - "), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0(fwdPrimerList), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0("Reverse - "), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0(revPrimerList), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0("Minimum read length for the final results (minFinalSeqLen): ", minFinalSeqLen), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0("********************************************************************************"), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
        suppressWarnings(write(paste0(""), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))









        ############################ FILTER & TRIM ##########################################

        #Message starting the initial quality figures - Print to screen and log file
        print(paste0("Begin generating initial quality figures at time...", Sys.time()))
        suppressWarnings(write(paste0("Begin generating initial quality figures at time...", Sys.time()), file = paste0(outFolder, "/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))

        #initialize the variable to store the raw data quality values
        fwdInputTempPlotValuesOut <- NULL
        revInputTempPlotValuesOut <- NULL
        fwdFilteredTempPlotValuesOut <- NULL
        revFilteredTempPlotValuesOut <- NULL

        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 6")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 6"), file = auditFile, append = TRUE))}

        if(bidirectional==TRUE | (bidirectional==TRUE & unidirectional == TRUE)){

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 7")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 7"), file = auditFile, append = TRUE))}

          #Place files in the forward and reverse variables based on the user
          fwdFile <- filesTbl[grepl(fwdIdent, filesTbl[,2]),,drop = FALSE]
          revFile <- filesTbl[grepl(revIdent, filesTbl[,2]),,drop = FALSE]

          for (listCounter in 1:length(fwdFile$fileNamesNoDir)){
            fwdInputTempPlotValues<-suppressWarnings(dada2::plotQualityProfile(fwdFile[listCounter,2])) # Forward
            fwdInputTempPlotValues<-summary(fwdInputTempPlotValues$data$Score)
            fwdInputTempPlotValues<-c(fwdFile$fileNamesNoDir[listCounter], fwdInputTempPlotValues[c(1,3,4,6)])
            fwdInputTempPlotValuesOut<-cbind(fwdInputTempPlotValuesOut, fwdInputTempPlotValues)
            if (printQualityPdf == TRUE){
              ggplot2::ggsave(paste0(qualityFolder, "/", dateStamp, fwdFile[listCounter,4], "_fwdQual.pdf"))
            }
          }

          for (listCounter in 1:length(revFile$fileNamesNoDir)){
            revInputTempPlotValues<-suppressWarnings(dada2::plotQualityProfile(revFile[listCounter,1])) # Reverse
            revInputTempPlotValues<-summary(revInputTempPlotValues$data$Score)
            revInputTempPlotValues<-c(revFile$fileNamesNoDir[listCounter], revInputTempPlotValues[c(1,3,4,6)])
            revInputTempPlotValuesOut<-cbind(revInputTempPlotValuesOut, revInputTempPlotValues)
            if (printQualityPdf == TRUE){
              ggplot2::ggsave(paste0(qualityFolder, "/", dateStamp, revFile[listCounter,4], "_revQual.pdf"))
            }
          }

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 8")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 8"), file = auditFile, append = TRUE))}

        }else{

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 9")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 9"), file = auditFile, append = TRUE))}

          fwdFile <- filesTbl

          for (listCounter in 1:length(fwdFile$fileNamesNoDir)){
            fwdInputTempPlotValues<-suppressWarnings(dada2::plotQualityProfile(fwdFile[listCounter,1])) # Forward
            fwdInputTempPlotValues<-summary(fwdInputTempPlotValues$data$Score)
            fwdInputTempPlotValues<-c(fwdFile$fileNamesNoDir[listCounter], fwdInputTempPlotValues[c(1,3,4,6)])
            fwdInputTempPlotValuesOut<-cbind(fwdInputTempPlotValuesOut, fwdInputTempPlotValues)
            if (printQualityPdf == TRUE){
              ggplot2::ggsave(paste0(qualityFolder, "/", dateStamp, fwdFile[listCounter,4], "_fwdQual.pdf"))
            }
          }

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 10")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 10"), file = auditFile, append = TRUE))}

        }

        ########################### PATTERN TRIMMING SECTION ############################################

        #Message quality trimming - Print to screen and log file
        print(paste0("Begin pattern based trimming at time...", Sys.time()))
        suppressWarnings(write(paste0("Begin pattern based trimming at time...", Sys.time()), file = paste0(outFolder, "/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))

        if(bidirectional==TRUE | (bidirectional==TRUE & unidirectional == TRUE)){

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 11")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 11"), file = auditFile, append = TRUE))}

          #For loop for each file in the forward direction in the submitted set
          for (listCounter in 1:length(fwdFile$fileNamesNoDir)){
            #Read in the fastqfile in to a ShortreadQ format file
            readInFile<-ShortRead::readFastq(fwdFile[listCounter,1])

            for(i in 1:length(primerFileData[,1])){
              #Pattern trim the reads
              readInFile <- ShortRead::trimLRPatterns(Lpattern = as.character(primerFileData[i,1]), Rpattern = "", subject = readInFile, max.Lmismatch = maxPrimeMis, max.Rmismatch	 = maxPrimeMis)
            }
            #Write files into a folder in the A_folder
            ShortRead::writeFastq(readInFile, paste0(filteredSubFolder, "/", fwdFile[listCounter,3], "_primeTrim.fastq.gz"), compress=TRUE)
          }

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 12")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 12"), file = auditFile, append = TRUE))}

          #For loop for each file in the reverse direction in the submitted set
          for (listCounter in 1:length(revFile$fileNamesNoDir)){
            #Read in the fastqfile in to a ShortreadQ format file
            readInFile<-readFastq(revFile[listCounter,1])

            for(i in 1:length(primerFileData[,1])){
              #Pattern trim the reads
              readInFile<-ShortRead::trimLRPatterns(Lpattern = as.character(primerFileData[i,2]), Rpattern = "", subject = readInFile, max.Lmismatch = maxPrimeMis, max.Rmismatch = maxPrimeMis)
            }
            #Write files into a folder in the A_folder
            ShortRead::writeFastq(readInFile, paste0(filteredSubFolder, "/", revFile[listCounter,3], "_primeTrim.fastq.gz"), compress=TRUE)
          }

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 13")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 13"), file = auditFile, append = TRUE))}

        }else{

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 14")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 14"), file = auditFile, append = TRUE))}

          #For loop for each file in the forward direction in the submitted set
          for (listCounter in 1:length(fwdFile$fileNamesNoDir)){
            #Read in the fastqfile in to a ShortreadQ format file
            readInFile<-readFastq(fwdFile[listCounter,1])

            for(i in 1:length(primerFileData[,1])){
              #Pattern trim the reads
              readInFile <- ShortRead::trimLRPatterns(Lpattern = as.character(primerFileData[i,1]), Rpattern = "", subject = readInFile, max.Lmismatch = maxPrimeMis, max.Rmismatch	 = maxPrimeMis)
            }
            #Write files into a folder in the A_folder
            ShortRead::writeFastq(readInFile, paste0(filteredSubFolder, "/", fwdFile[listCounter,3], "_primeTrim.fastq.gz"), compress=TRUE)
          }

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 15")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 15"), file = auditFile, append = TRUE))}

        }

        ########################### QUALITY TRIMMING SECTION ############################################

        #Message quality trimming - Print to screen and log file
        print(paste0("Begin filtering and trimming based on quality arguments at time...", Sys.time()))
        suppressWarnings(write(paste0("Begin filtering and trimming based on quality arguments at time...", Sys.time()), file = paste0(outFolder, "/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))

        if(bidirectional==TRUE | (bidirectional==TRUE & unidirectional == TRUE)){

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 16")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 16"), file = auditFile, append = TRUE))}

          #Variable names for input prime trimmed files
          inputFwdFilt <- paste0(filteredSubFolder, "/", fwdFile[,3], "_primeTrim.fastq.gz")
          inputRevFilt <- paste0(filteredSubFolder, "/", revFile[,3], "_primeTrim.fastq.gz")

          #Variable names for filtered & trimmed fastq files
          fwdFilt <- paste0(filteredFolder, "/", fwdFile[,4], "_fwdFilt.fastq.gz")
          revFilt <- paste0(filteredFolder, "/", revFile[,4], "_revFilt.fastq.gz")

          #Filter and trim
          if (.Platform$OS.type == "windows"){
            filterOut <- filterAndTrim(fwd = inputFwdFilt, filt = fwdFilt, rev = inputRevFilt, filt.rev = revFilt, compress = TRUE, truncQ = truncQValue, truncLen = truncLenValue, trimLeft = c(fwdTrimLen,revTrimLen), maxN = maxNVal, maxEE = maxEEVal, verbose = TRUE, multithread = FALSE)
          } else{
            filterOut <- filterAndTrim(fwd = inputFwdFilt, filt = fwdFilt, rev = inputRevFilt, filt.rev = revFilt, compress = TRUE, truncQ = truncQValue, truncLen = truncLenValue, trimLeft = c(fwdTrimLen,revTrimLen), maxN = maxNVal, maxEE = maxEEVal, verbose = TRUE, multithread = TRUE)
          }

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 17")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 17"), file = auditFile, append = TRUE))}

          #Get the list of files in the filtered folder
          fwdFiltOutFiles<-list.files(path = filteredFolder, pattern = "*fwdFilt[.]*", full.names = TRUE)
          revFiltOutFiles<-list.files(path = filteredFolder, pattern = "*revFilt[.]*", full.names = TRUE)

          #remove the ending of the read in files to create a list to use in the naming below
          #Create dataframe with filtered file paths and names
          fwdFiltOutFilesName <- sub(".*/","", fwdFiltOutFiles)
          fwdFiltOutFilesName <- sub("_fwdFilt.*","", fwdFiltOutFilesName)

          revFiltOutFilesName <- sub(".*/","", revFiltOutFiles)
          revFiltOutFilesName <- sub("_revFilt.*","", revFiltOutFilesName)

          #get the records that are in common
          totalFiltOutFiles<-intersect(fwdFiltOutFilesName, revFiltOutFilesName)

          #Report the files that were not able pass the initial quality filtering
          poorQual<-filesTbl[!(filesTbl[,4] %in% totalFiltOutFiles) ,2]

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 18")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 18"), file = auditFile, append = TRUE))}

          if(length(poorQual)>0){
            print(paste0("There were low quality samples not able to pass the quality filtering at time ", Sys.time(), "..."))
            print(poorQual)
            suppressWarnings(write(paste0("There were low quality samples not able to pass the quality filtering at time ", Sys.time(), "..."), file = paste0(outFolder, "/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
            suppressWarnings(write(poorQual, file = paste0(outFolder, "/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
          }

          #Remove from the list of files if they don't have both F and R
          fwdFiltOutFiles<-fwdFiltOutFiles[grepl(paste(totalFiltOutFiles, collapse="|"), fwdFiltOutFiles)]
          revFiltOutFiles<-revFiltOutFiles[grepl(paste(totalFiltOutFiles, collapse="|"), revFiltOutFiles)]

          #Create dataframe with filtered file paths and names
          fwdFiltTbl <- as.data.frame(cbind(fwdFiltOutFiles, totalFiltOutFiles))
          revFiltTbl <- as.data.frame(cbind(revFiltOutFiles, totalFiltOutFiles))

          if (printQualityPdf == TRUE){
            #Message starting the initial quality figures - Print to screen and log file
            print(paste0("Begin generating quality filtered figures at time...", Sys.time()))
            suppressWarnings(write(paste0("Begin generating quality filtered figures at time...", Sys.time()), file = paste0(outFolder, "/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
          }

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 19")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 19"), file = auditFile, append = TRUE))}

          #Output the quality for the files post filter and trim in a figure
          for (listCounter in 1:length(totalFiltOutFiles)){
            fwdFilteredTempPlotValues<-suppressWarnings(plotQualityProfile(fwdFiltTbl[listCounter,1])) # Forward
            fwdFilteredTempPlotValues<-summary(fwdFilteredTempPlotValues$data$Score)
            fwdFilteredTempPlotValues<-c(totalFiltOutFiles[listCounter], fwdFilteredTempPlotValues[c(1,3,4,6)])
            fwdFilteredTempPlotValuesOut<-cbind(fwdFilteredTempPlotValuesOut, fwdFilteredTempPlotValues)
            if (printQualityPdf == TRUE){
              ggsave(paste0(filterQualFolder, "/", dateStamp,dirName,fwdFiltTbl[listCounter,2], "_filtFwdQual.pdf"))
            }
          }

          for (listCounter in 1:nrow(revFiltTbl)){
            revFilteredTempPlotValues<-suppressWarnings(plotQualityProfile(revFiltTbl[listCounter,1])) # Reverse
            revFilteredTempPlotValues<-summary(revFilteredTempPlotValues$data$Score)
            revFilteredTempPlotValues<-c(totalFiltOutFiles[listCounter], revFilteredTempPlotValues[c(1,3,4,6)])
            revFilteredTempPlotValuesOut<-cbind(revFilteredTempPlotValuesOut, revFilteredTempPlotValues)
            if (printQualityPdf == TRUE){
              ggsave(paste0(filterQualFolder, "/", dateStamp,dirName, revFiltTbl[listCounter,2], "_filtRevQual.pdf"))
            }
          }

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 20")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 20"), file = auditFile, append = TRUE))}

          #Getting values for reporting
          #Get the column names
          fwdInputTempPlotValuesOutColNames<-fwdInputTempPlotValuesOut[1,]
          revInputTempPlotValuesOutColNames<-revInputTempPlotValuesOut[1,]
          fwdFilteredTempPlotValuesOutColNames<-fwdFilteredTempPlotValuesOut[1,]
          revFilteredTempPlotValuesOutColNames<-revFilteredTempPlotValuesOut[1,]

          #remove the first row from the qual data frame
          fwdInputTempPlotValuesOut<-as.data.frame(fwdInputTempPlotValuesOut[-1,])
          revInputTempPlotValuesOut<-as.data.frame(revInputTempPlotValuesOut[-1,])
          fwdFilteredTempPlotValuesOut<-as.data.frame(fwdFilteredTempPlotValuesOut[-1,])
          revFilteredTempPlotValuesOut<-as.data.frame(revFilteredTempPlotValuesOut[-1,])

          #add in the col names
          colnames(fwdInputTempPlotValuesOut)<-fwdInputTempPlotValuesOutColNames
          colnames(revInputTempPlotValuesOut)<-revInputTempPlotValuesOutColNames
          colnames(fwdFilteredTempPlotValuesOut)<-fwdFilteredTempPlotValuesOutColNames
          colnames(revFilteredTempPlotValuesOut)<-revFilteredTempPlotValuesOutColNames

          #Add in the row names
          row.names(fwdInputTempPlotValuesOut)<-c("fwdInputQualMin", "fwdInputQualMedian", "fwdInputQualMean", "fwdInputQualMax")
          row.names(revInputTempPlotValuesOut)<-c("revInputQualMin", "revInputQualMedian", "revInputQualMean", "revInputQualMax")
          row.names(fwdFilteredTempPlotValuesOut)<-c("fwdFilteredQualMin", "fwdFilteredQualMedian", "fwdFilteredQualMean", "fwdFilteredQualMax")
          row.names(revFilteredTempPlotValuesOut)<-c("revFilteredQualMin", "revFilteredQualMedian", "revFilteredQualMean", "revFilteredQualMax")

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 21")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 21"), file = auditFile, append = TRUE))}

        }else{

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 22")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 22"), file = auditFile, append = TRUE))}

          #Variable names for input prime trimmed files
          inputFwdFilt <- paste0(filteredSubFolder, "/", fwdFile[,3], "_primeTrim.fastq.gz")

          #Variable names for filtered & trimmed fastq files
          fwdFilt <- paste0(filteredFolder, "/", fwdFile[,4], "_fwdFilt.fastq.gz")

          #Filter and trim
          if (.Platform$OS.type == "windows"){
            filterOut <- filterAndTrim(fwd = inputFwdFilt, filt = fwdFilt, compress = TRUE, truncQ = truncQValue, truncLen = truncLenValueF, trimLeft = fwdTrimLen, maxN = maxNVal, maxEE = maxEEVal, verbose = TRUE, multithread = FALSE)
          } else{
            filterOut <- filterAndTrim(fwd = inputFwdFilt, filt = fwdFilt, compress = TRUE, truncQ = truncQValue, truncLen = truncLenValueF, trimLeft = fwdTrimLen, maxN = maxNVal, maxEE = maxEEVal, verbose = TRUE, multithread = TRUE)
          }

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 23")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 23"), file = auditFile, append = TRUE))}

          #Get the list of files in the filtered folder
          fwdFiltOutFiles<-list.files(path = filteredFolder, pattern = "*fwdFilt[.]*", full.names = TRUE)

          #remove the ending of the read in files to create a list to use in the naming below
          #Create dataframe with filtered file paths and names
          fwdFiltOutFilesName <- sub(".*/","", fwdFiltOutFiles)
          fwdFiltOutFilesName <- sub("_fwdFilt.*","", fwdFiltOutFilesName)

          #get the records that are in common
          totalFiltOutFiles<-fwdFiltOutFilesName

          #Report the files that were not able pass the initial quality filtering
          poorQual<-filesTbl[!(filesTbl[,4] %in% totalFiltOutFiles) ,2]

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 24")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 24"), file = auditFile, append = TRUE))}

          if(length(poorQual)>0){
            print(paste0("There were low quality samples not able to pass the quality filtering at time ", Sys.time(), "..."))
            print(poorQual)
            suppressWarnings(write(paste0("There were low quality samples not able to pass the quality filtering at time ", Sys.time(), "..."), file = paste0(outFolder, "/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
            suppressWarnings(write(poorQual, file = paste0(outFolder, "/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
          }

          #Remove from the list of files if they don't have both F and R
          fwdFiltOutFiles<-fwdFiltOutFiles[grepl(paste(totalFiltOutFiles, collapse="|"), fwdFiltOutFiles)]

          #Create dataframe with filtered file paths and names
          fwdFiltTbl <- as.data.frame(cbind(fwdFiltOutFiles, totalFiltOutFiles))

          #Message starting the initial quality figures - Print to screen and log file
          print(paste0("Begin generating quality filtered figures at time...", Sys.time()))
          suppressWarnings(write(paste0("Begin generating quality filtered figures at time...", Sys.time()), file = paste0(outFolder, "/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))

          #Output the quality for the files post filter and trim in a figure
          for (listCounter in 1:length(totalFiltOutFiles)){
            fwdFilteredTempPlotValues<-suppressWarnings(plotQualityProfile(fwdFiltTbl[listCounter,1])) # Forward
            fwdFilteredTempPlotValues<-summary(fwdFilteredTempPlotValues$data$Score)
            fwdFilteredTempPlotValues<-c(totalFiltOutFiles[listCounter], fwdFilteredTempPlotValues[c(1,3,4,6)])
            fwdFilteredTempPlotValuesOut<-cbind(fwdFilteredTempPlotValuesOut, fwdFilteredTempPlotValues)
            if (printQualityPdf == TRUE){
              ggsave(paste0(filterQualFolder, "/", dateStamp,dirName,fwdFiltTbl[listCounter,2], "_filtFwdQual.pdf"))
            }
          }

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 25")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 25"), file = auditFile, append = TRUE))}

          #Getting values for reporting
          #Get the column names
          fwdInputTempPlotValuesOutColNames<-fwdInputTempPlotValuesOut[1,]
          fwdFilteredTempPlotValuesOutColNames<-fwdFilteredTempPlotValuesOut[1,]

          #remove the first row from the qual data frame
          fwdInputTempPlotValuesOut<-as.data.frame(fwdInputTempPlotValuesOut[-1,])
          fwdFilteredTempPlotValuesOut<-as.data.frame(fwdFilteredTempPlotValuesOut[-1,])

          #add in the col names
          colnames(fwdInputTempPlotValuesOut)<-fwdInputTempPlotValuesOutColNames
          colnames(fwdFilteredTempPlotValuesOut)<-fwdFilteredTempPlotValuesOutColNames

          #Add in the row names
          row.names(fwdInputTempPlotValuesOut)<-c("fwdInputQualMin", "fwdInputQualMedian", "fwdInputQualMean", "fwdInputQualMax")
          row.names(fwdFilteredTempPlotValuesOut)<-c("fwdFilteredQualMin", "fwdFilteredQualMedian", "fwdFilteredQualMean", "fwdFilteredQualMax")

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 26")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 26"), file = auditFile, append = TRUE))}

        }

        #Section to get estimated error rates by looking at user specified percent of the files for this run
        ############################ DETERMINING EXPECTED ERROR RATES FOR THE RUN ##########################################

        #Print to screen and log file
        print(paste0("Begin getting error rates at time...", Sys.time()))
        suppressWarnings(write(paste0("Begin getting error rates at time...", Sys.time()), file = paste0(outFolder, "/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))

        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 27")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 27"), file = auditFile, append = TRUE))}

        #If there are 3 or fewer (times the error) use the number of samples that passed the quality filter,
        #anything more than 3 use the floor of the percentage of records as supplied by user
        if(length(totalFiltOutFiles)* error < 3){
          if(length(totalFiltOutFiles)>2){
            numErrFiles = 3
          }else{
            numErrFiles = length(totalFiltOutFiles)
          }
        } else {
          numErrFiles = floor(length(totalFiltOutFiles) * error)
        }

        #Random subset of the samples for the user specified percent
        errCheckFiles <- sample(totalFiltOutFiles, numErrFiles)

        #Print the error checking files to the screen and log file
        print(paste0("Here are the error checking files used..."))
        print(filesTbl[filesTbl$fileNamesNoDir %in% errCheckFiles,]$fileList)
        errCheckNote <- c(paste0("Here are the error checking files used..."), filesTbl[filesTbl$fileNamesNoDir %in% errCheckFiles,]$fileList)
        suppressWarnings(write(errCheckNote, file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))

        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 28")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 28"), file = auditFile, append = TRUE))}

        if(bidirectional==TRUE | (bidirectional==TRUE & unidirectional == TRUE)){

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 29")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 29"), file = auditFile, append = TRUE))}

          #Extract the filtered error check files
          fwdFiltError <- fwdFiltTbl$fwdFiltOutFiles[fwdFiltTbl$totalFiltOutFiles %in% errCheckFiles]
          revFiltError <- revFiltTbl$revFiltOutFiles[revFiltTbl$totalFiltOutFiles %in% errCheckFiles]

          #Estimating the error present in the reads and print to screen and log file
          print(paste0("At the beginning of the forward read error rate calculation at ", Sys.time()))
          write.table(paste0("At the beginning of the forward read error rate calculation at ", Sys.time()), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE, na = "", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n")
          errF <- learnErrors(fwdFiltError, nbases = nbases, multithread = TRUE)

          print(paste0("At the beginning of the reverse read error rate calculation at ", Sys.time()))
          write.table(paste0("At the beginning of the reverse read error rate calculation at ", Sys.time()), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE, na = "", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n")
          errR <- learnErrors(revFiltError, nbases = nbases, multithread = TRUE)

          invisible(plotErrors(errF, nominalQ = TRUE))
          suppressWarnings(ggsave(paste0(outFolder, "/", dateStamp, dirName, "_ErrorForward.pdf")))

          invisible(plotErrors(errR, nominalQ = TRUE))
          suppressWarnings(ggsave(paste0(outFolder, "/", dateStamp, dirName, "_ErrorReverse.pdf")))

          #Print to screen and log file
          print(paste0("At the end of the error rate calculations at ", Sys.time()))
          suppressWarnings(write(paste0("At the end of the error rate calculations at ", Sys.time()), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))

          ########## DADA FUNCTIONS SECTION (QUALITY FILTER, DEREP, ERROR RATE FILTER, MERGE) ##################

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 30")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 30"), file = auditFile, append = TRUE))}

          #Start chimera section
          logStr <- paste0("Start dereplicating at ", Sys.time())
          print(logStr)
          suppressWarnings(write(logStr, file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))

          #De-replicating the filtered and trimmed sequences
          fwdDerep <- derepFastq(fwdFiltOutFiles, verbose = TRUE)
          revDerep <- derepFastq(revFiltOutFiles, verbose = TRUE)

          ################## PREPARING THE FINAL DATA FILE #########################################

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 31")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 31"), file = auditFile, append = TRUE))}

          #Infer sample composition using the estimated error rates from the previous step
          logStr <- paste0("Start forward dada at ", Sys.time())
          print(logStr)
          suppressWarnings(write(logStr, file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
          dadaF1 <- dada(fwdDerep, err = errF, multithread = TRUE)

          logStr <- paste0("Start reverse dada at ", Sys.time())
          print(logStr)
          suppressWarnings(write(logStr, file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
          dadaR1 <- dada(revDerep, err = errR, multithread = TRUE)

          #Merge forward and reverse reads
          logStr <- paste0("Start merged at ", Sys.time())
          print(logStr)
          suppressWarnings(write(logStr, file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
          merged <- mergePairs(dadaF1, fwdDerep, dadaR1, revDerep, trimOverhang = TRUE, returnRejects = FALSE, maxMismatch = maxMismatchValue, verbose = TRUE, minOverlap = minOverlapValue)

          #Prepare tables for export
          finalTable <- t(makeSequenceTable(merged))

          #Create a final table dataframe for reporting general numbers at the end.
          mergedTable<-data.frame(Sequence = row.names(finalTable), finalTable, check.names=FALSE)

          if(nrow(finalTable)==0){

            logStr <- paste0("No sequences were successfully merged - ", Sys.time())
            print(logStr)
            suppressWarnings(write(logStr, file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
            finalTable <-cbind(finalTable, rep(0,length(finalTable)))
            finalTable<-as.data.frame(cbind(Results="Merged", finalTable), drop=FALSE, check.names = FALSE)

          }else{
            #Add in a column to contain analysis information
            finalTable<-as.data.frame(cbind(Results="Merged", finalTable), drop=FALSE, check.names = FALSE)
          }

          ####################### CHIMERA SECTION ###################################################

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 32")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 32"), file = auditFile, append = TRUE))}

          #Start chimera section
          logStr <- paste0("Start chimera checking at ", Sys.time())
          print(logStr)
          suppressWarnings(write(logStr, file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))

          #Remove chimeric sequences
          if (.Platform$OS.type == "windows"){
            noChim <- removeBimeraDenovo(merged, multithread = FALSE, verbose = TRUE)
          } else{
            noChim <- removeBimeraDenovo(merged, multithread = TRUE, verbose = TRUE)
          }

          #Prepare tables for export
          nochimTblTrimmed <- t(makeSequenceTable(noChim))

          #Create a final table dataframe for reporting general numbers at the end.
          nochimTrimmedOut<- data.frame(Sequence = rownames(nochimTblTrimmed), nochimTblTrimmed, check.names=FALSE)

          #Get the reads that were removed as chimeric
          chim<-finalTable[!rownames(finalTable) %in% rownames(nochimTblTrimmed),]

          if(nrow(chim)>0){
            #Set the final Table Results column to indicate merged reads that were chimeric
            finalTable[rownames(finalTable) %in% rownames(chim), "Results"] <- "ChimRemoved"
          }

          #Add a sequence column to the final table
          finalTable<-data.frame(Sequence = row.names(finalTable),Length = nchar(as.character(row.names(finalTable))), finalTable, check.names=FALSE)

          ###################### Output the data files with the merged data ############
          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 33")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 33"), file = auditFile, append = TRUE))}

          #Get a subset with records that are not too short
          notShort<-finalTable[finalTable$Length > minFinalSeqLen,, drop=FALSE]

          #Remove the file extension from the header names and change the . back to - if needed
          colnames(notShort)<-gsub("_fwdFilt.fastq.gz","",colnames(notShort))
          colnames(notShort)<-gsub("\\.","-",colnames(notShort))

          #Remove NoChim records
          notShort<-notShort[notShort$Results != "ChimRemoved",, drop=FALSE]

          #add unique ids to the records
          notShort <- data.frame(uniqueID = paste0("MRG",c(1:nrow(notShort))), notShort, check.names=FALSE)

          #Print to file the final dataset with short sequences removed and in the format for use in the taxonomic assignment step
          write.table(notShort, file = paste0(outFolder, "/", dateStamp, dirName, "_Merge.tsv"), append = FALSE, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

          #add unique ids to the records
          seqOut <- data.frame(uniqueID = paste0(">MRG",c(1:nrow(notShort))), notShort$Sequence, check.names=FALSE)

          #print to file
          write.table(seqOut, file = paste0(outFolder, "/", dateStamp, dirName, "_Merge.fas"), append = FALSE, na = "", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n")

          ####################### COMPLETE UNMERGED ANALYSIS ###################################################

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 34")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 34"), file = auditFile, append = TRUE))}

          if(unidirectional == TRUE){

            #Audit line
            if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 35")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 35"), file = auditFile, append = TRUE))}

            #Start trimming section
            logStr <- paste0("Completing the unmerged analyses at... ", Sys.time())
            print(logStr)
            suppressWarnings(write(logStr, file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))

            #Get the Forward table
            fwdTblUnpairedTmp <- t(makeSequenceTable(dadaF1))

            #Preparing the forward data file to be combine with the final data file
            #Add in a column to contain analysis information
            fwdTblUnpairedTmp<-cbind(Results="Forward", fwdTblUnpairedTmp)

            #Add a column for the sequence lengths
            fwdTblUnpairedTmp <- data.frame(Length = nchar(as.character(row.names(fwdTblUnpairedTmp))), fwdTblUnpairedTmp, check.names=FALSE)

            #Make the row names the first column
            fwdTblUnpairedTmp <- data.frame(Sequence = as.character(row.names(fwdTblUnpairedTmp)), fwdTblUnpairedTmp, check.names=FALSE)

            #Create a data frame with a single variable to aggregate
            mergeVariable<-paste0(fwdTblUnpairedTmp$Sequence,"|",fwdTblUnpairedTmp$Length,"|",fwdTblUnpairedTmp$Results)

            #build a temp dataframe
            fwdTblUnpairedTmpTemp<-as.data.frame(cbind(mergeVariable,fwdTblUnpairedTmp[,4:ncol(fwdTblUnpairedTmp)]))
            colnames(fwdTblUnpairedTmpTemp)[1]<-"mergeVariable"

            #make the columns numeric
            fwdTblUnpairedTmpTemp <- cbind(fwdTblUnpairedTmpTemp[,1, drop=FALSE], apply(fwdTblUnpairedTmpTemp[,2:ncol(fwdTblUnpairedTmpTemp), drop=FALSE] , 2, function(x) as.numeric(as.character(x))))

            #reduce all of the unique mergeVariables and sum the values in the sample columns
            fwdTblUnpairedTmpTemp <- stats::aggregate(. ~ mergeVariable,fwdTblUnpairedTmpTemp, sum)

            #Take the sequence_id column and split into columns
            fwdTblUnpairedTmpTemp <- data.frame(cbind(do.call('rbind', strsplit(as.character(fwdTblUnpairedTmpTemp[,1]),'|', fixed = TRUE)),fwdTblUnpairedTmpTemp[,2:ncol(fwdTblUnpairedTmpTemp)]), check.names=FALSE)

            #add in row names
            colnames(fwdTblUnpairedTmpTemp)<-colnames(fwdTblUnpairedTmp)

            #make the temp equal to the finaTable
            fwdTblUnpairedTmp<-fwdTblUnpairedTmpTemp

            #Audit line
            if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 36")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 36"), file = auditFile, append = TRUE))}

            #Get the Reverse table
            revTblUnpairedTmp <- t(makeSequenceTable(dadaR1))

            #Preparing the forward data file to be combine with the final data file
            #Add in a column to contain analysis information
            revTblUnpairedTmp<-cbind(Results="Reverse", revTblUnpairedTmp)

            #Add a column for the sequence lengths
            revTblUnpairedTmp <- data.frame(Length = nchar(as.character(row.names(revTblUnpairedTmp))), revTblUnpairedTmp, check.names=FALSE)

            #Make the row names the first column
            revTblUnpairedTmp <- data.frame(Sequence = as.character(row.names(revTblUnpairedTmp)), revTblUnpairedTmp, check.names=FALSE)

            #Create a data frame with a single variable to aggregate
            mergeVariable<-paste0(revTblUnpairedTmp$Sequence,"|",revTblUnpairedTmp$Length,"|",revTblUnpairedTmp$Results)

            #build a temp dataframe
            revTblUnpairedTmpTemp<-as.data.frame(cbind(mergeVariable,revTblUnpairedTmp[,4:ncol(revTblUnpairedTmp)]), check.names=FALSE)
            colnames(revTblUnpairedTmpTemp)[1]<-"mergeVariable"

            #make the columns numeric
            revTblUnpairedTmpTemp <- cbind(revTblUnpairedTmpTemp[,1, drop=FALSE], apply(revTblUnpairedTmpTemp[,2:ncol(revTblUnpairedTmpTemp), drop=FALSE] , 2, function(x) as.numeric(as.character(x))))

            #reduce all of the unique mergeVariables and sum the values in the sample columns
            revTblUnpairedTmpTemp <- stats::aggregate(. ~ mergeVariable,revTblUnpairedTmpTemp, sum)

            #Take the sequence_id column and split into columns
            revTblUnpairedTmpTemp <- data.frame(cbind(do.call('rbind', strsplit(as.character(revTblUnpairedTmpTemp[,1]),'|', fixed = TRUE)),revTblUnpairedTmpTemp[,2:ncol(revTblUnpairedTmpTemp)]), check.names=FALSE)

            #add in row names
            colnames(revTblUnpairedTmpTemp)<-colnames(revTblUnpairedTmp)

            #make the temp equal to the finaTable
            revTblUnpairedTmp<-revTblUnpairedTmpTemp

            #Audit line
            if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 37")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 37"), file = auditFile, append = TRUE))}

            ########################## ADD IN THE FORWARD AND REVERSE ANALYSIS TABLES IN TO THE FINAL TABLE #############
            #Remove exact sequences or sub sequences from the forward and reverse tables if the exist in the merged table
            fwdTblUnpairedTmp<-fwdTblUnpairedTmp[!(fwdTblUnpairedTmp$Sequence %in% finalTable$Sequence),, drop=FALSE ]
            revTblUnpairedTmp<-revTblUnpairedTmp[!(revTblUnpairedTmp$Sequence %in% finalTable$Sequence),, drop=FALSE ]

            #Get the maximum length of merged sequences so that we only check forward and reverse subsequences if they are less than the max
            maxMergedSeq<-as.numeric(max(nchar(finalTable$Sequence)))
            fwdOutput<-NULL
            revOutput<-NULL

            #These could be apply
            if (nrow(fwdTblUnpairedTmp)>0){

              for (fwdFinalTableReduceCount in 1:nrow(fwdTblUnpairedTmp)){

                if(as.numeric(nchar(fwdTblUnpairedTmp$Sequence[fwdFinalTableReduceCount]))>maxMergedSeq){

                  fwdOutput<-rbind(fwdOutput, fwdTblUnpairedTmp[fwdFinalTableReduceCount,])

                } else if (nrow(fwdTblUnpairedTmp[grepl(fwdTblUnpairedTmp$Sequence[fwdFinalTableReduceCount], finalTable$Sequence),])==0){

                  fwdOutput<-rbind(fwdOutput, fwdTblUnpairedTmp[fwdFinalTableReduceCount,])

                }
              }
            }
            fwdTblUnpairedTmp <- fwdOutput

            #Audit line
            if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 38")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 38"), file = auditFile, append = TRUE))}

            #These could be apply
            if (nrow(revTblUnpairedTmp)>0){

              for (revFinalTableReduceCount in 1:nrow(revTblUnpairedTmp)){

                if(as.numeric(nchar(revTblUnpairedTmp$Sequence[revFinalTableReduceCount]))>maxMergedSeq){

                  revOutput<-rbind(revOutput, revTblUnpairedTmp[revFinalTableReduceCount,])

                } else if (nrow(revTblUnpairedTmp[grepl(revTblUnpairedTmp$Sequence[revFinalTableReduceCount], finalTable$Sequence),])==0){

                  revOutput<-rbind(revOutput, revTblUnpairedTmp[revFinalTableReduceCount,])

                }
              }
            }
            revTblUnpairedTmp <- revOutput

            #Audit line
            if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 39")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 39"), file = auditFile, append = TRUE))}

            ##reverse complement to correct orientation
            revTblUnpairedTmp$Sequence <- dada2::rc(revTblUnpairedTmp$Sequence) #Dada2 package function

            #Make all column headers equal
            colnames(revTblUnpairedTmp)<-colnames(finalTable)

            #Add the forward and reverse results to the finalTable
            finalTable<-rbind(finalTable, fwdTblUnpairedTmp, revTblUnpairedTmp)

            #Audit line
            if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 40")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 40"), file = auditFile, append = TRUE))}

            #Get a subset with records that are not too short
            notShort<-finalTable[finalTable$Length > minFinalSeqLen,, drop=FALSE]

            #Remove the file extension from the header names and change the . back to - if needed
            colnames(notShort)<-gsub("_fwdFilt.fastq.gz","",colnames(notShort))
            colnames(notShort)<-gsub("\\.","-",colnames(notShort))

            #Remove NoChim records
            notShort<-notShort[notShort$Results != "ChimRemoved",, drop=FALSE]

            #add unique ids to the records
            notShort <- data.frame(uniqueID = paste0("MFR",c(1:nrow(notShort))), notShort, check.names=FALSE)

            #Print to file the final dataset with short sequences removed and in the format for use in the taxonomic assignment step
            write.table(notShort, file = paste0(outFolder, "/", dateStamp, dirName, "_MergeFwdRev.tsv"), append = FALSE, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

            #add unique ids to the records
            seqOut <- data.frame(uniqueID = paste0(">MFR",c(1:nrow(notShort))), notShort$Sequence, check.names=FALSE)

            #print to file
            write.table(seqOut, file = paste0(outFolder, "/", dateStamp, dirName, "_MergeFwdRev.fas"), append = FALSE, na = "", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n")

            #Audit line
            if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 41")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 41"), file = auditFile, append = TRUE))}


          }#end of if for including non_paired merged analysis

        }else{ #Else for Bidirectional or Bi and unidirectional if

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 41")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 41"), file = auditFile, append = TRUE))}

          #Extract the filtered error check files
          fwdFiltError <- fwdFiltTbl$fwdFiltOutFiles[fwdFiltTbl$totalFiltOutFiles %in% errCheckFiles]

          #Estimating the error present in the reads and print to screen and log file
          print(paste0("At the beginning of the forward read error rate calculation at ", Sys.time()))
          write.table(paste0("At the beginning of the forward read error rate calculation at ", Sys.time()), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE, na = "", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n")
          errF <- learnErrors(fwdFiltError, nbases = nbases, multithread = TRUE)

          invisible(plotErrors(errF, nominalQ = TRUE))
          suppressWarnings(ggsave(paste0(outFolder, "/", dateStamp, dirName, "_ErrorForward.pdf")))

          #Print to screen and log file
          print(paste0("At the end of the error rate calculations at ", Sys.time()))
          suppressWarnings(write(paste0("At the end of the error rate calculations at ", Sys.time()), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))

          ########## DADA FUNCTIONS SECTION (QUALITY FILTER, DEREP, ERROR RATE FILTER, MERGE) ##################

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 42")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 42"), file = auditFile, append = TRUE))}

          #Start chimera section
          logStr <- paste0("Start dereplicating at ", Sys.time())
          print(logStr)
          suppressWarnings(write(logStr, file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))

          #De-replicating the filtered and trimmed sequences
          fwdDerep <- derepFastq(fwdFiltOutFiles, verbose = TRUE)

          ################## PREPARING THE FINAL DATA FILE #########################################

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 43")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 43"), file = auditFile, append = TRUE))}

          #Infer sample composition using the estimated error rates from the previous step
          logStr <- paste0("Start forward dada at ", Sys.time())
          print(logStr)
          suppressWarnings(write(logStr, file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))
          dadaF1 <- dada(fwdDerep, err = errF, multithread = TRUE)

          #Get the Forward table
          fwdTblUnpairedTmp <- t(makeSequenceTable(dadaF1))

          #Preparing the forward data file to be combine with the final data file
          #Add in a column to contain analysis information
          fwdTblUnpairedTmp<-cbind(Results="Forward", fwdTblUnpairedTmp)

          #Add a column for the sequence lengths
          fwdTblUnpairedTmp <- data.frame(Length = nchar(as.character(row.names(fwdTblUnpairedTmp))), fwdTblUnpairedTmp, check.names=FALSE)

          #Make the row names the first column
          fwdTblUnpairedTmp <- data.frame(Sequence = as.character(row.names(fwdTblUnpairedTmp)), fwdTblUnpairedTmp, check.names=FALSE)

          #Create a data frame with a single variable to aggregate
          mergeVariable<-paste0(fwdTblUnpairedTmp$Sequence,"|",fwdTblUnpairedTmp$Length,"|",fwdTblUnpairedTmp$Results)

          #build a temp dataframe
          fwdTblUnpairedTmpTemp<-as.data.frame(cbind(mergeVariable,fwdTblUnpairedTmp[,4:ncol(fwdTblUnpairedTmp)]))
          colnames(fwdTblUnpairedTmpTemp)[1]<-"mergeVariable"

          #make the columns numeric
          fwdTblUnpairedTmpTemp <- cbind(fwdTblUnpairedTmpTemp[,1, drop=FALSE], apply(fwdTblUnpairedTmpTemp[,2:ncol(fwdTblUnpairedTmpTemp), drop=FALSE] , 2, function(x) as.numeric(as.character(x))))

          #reduce all of the unique mergeVariables and sum the values in the sample columns
          fwdTblUnpairedTmpTemp <- stats::aggregate(. ~ mergeVariable,fwdTblUnpairedTmpTemp, sum)

          #Take the sequence_id column and split into columns
          fwdTblUnpairedTmpTemp <- data.frame(cbind(do.call('rbind', strsplit(as.character(fwdTblUnpairedTmpTemp[,1]),'|', fixed = TRUE)),fwdTblUnpairedTmpTemp[,2:ncol(fwdTblUnpairedTmpTemp)]), check.names=FALSE)

          #add in row names
          colnames(fwdTblUnpairedTmpTemp)<-colnames(fwdTblUnpairedTmp)

          #make the temp equal to the finaTable
          finalTable<-fwdTblUnpairedTmpTemp

          #Get a subset with records that are not too short
          notShort<-finalTable[finalTable$Length > minFinalSeqLen,, drop=FALSE]

          #Remove the file extension from the header names and change the . back to - if needed
          colnames(notShort)<-gsub("_fwdFilt.fastq.gz","",colnames(notShort))
          colnames(notShort)<-gsub("\\.","-",colnames(notShort))

          #add unique ids to the records
          notShort <- data.frame(uniqueID = paste0("UNI",c(1:nrow(notShort))), notShort, check.names=FALSE)

          #Print to file the final dataset with short sequences removed and in the format for use in the taxonomic assignment step
          write.table(notShort, file = paste0(outFolder, "/", dateStamp, dirName, "_Foward.tsv"), append = FALSE, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

          #add unique ids to the records
          seqOut <- data.frame(uniqueID = paste0(">UNI",c(1:nrow(notShort))), notShort$Sequence, check.names=FALSE)

          #print to file
          write.table(seqOut, file = paste0(outFolder, "/", dateStamp, dirName, "_Forward.fas"), append = FALSE, na = "", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n")

          #Audit line
          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 44")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 44"), file = auditFile, append = TRUE))}

        }#End of bi or uni and bi directional if

        #output the final table with all data.

        #add unique ids to the records
        finalTableOut <- data.frame(uniqueID = paste0("T",c(1:nrow(finalTable))), finalTable, check.names=FALSE)

        #Remove the file extension from the header names
        colnames(finalTableOut)<-gsub("_fwdFilt.fastq.gz","",colnames(finalTableOut))
        colnames(finalTableOut)<-gsub("\\.","-",colnames(finalTableOut))

        #print to file the total data set table
        write.table(finalTableOut, file = paste0(outFolder, "/", dateStamp, dirName, "_TotalTable.tsv"), append = FALSE, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 45")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 45"), file = auditFile, append = TRUE))}

        ########################################## RUN REPORTING #############################
        #Start reporting section
        logStr <- paste0("Running the reporting for the analyses... ", Sys.time())
        print(logStr)
        suppressWarnings(write(logStr, file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))

        ## Make functions to calculate the other reporting values - Stats functions
        getN <- function(x) sum(getUniques(x))
        getSeqQual <- function(seqRow) mean(na.omit(as.numeric(seqRow$quality)))
        getMeanSeqLen <- function(seqRow) mean(nchar(seqRow$sequence))
        getSeqLength <- function(seqRow) length(seqRow$sequence)

        #Get the filterOut values ready for reporting
        row.names(filterOut)<-gsub("_primeTrim.fastq.gz","",row.names(filterOut))
        row.names(filterOut)<-gsub(fwdIdent,"",row.names(filterOut))
        row.names(filterOut)<-gsub(revIdent,"",row.names(filterOut))
        filterOut<-t(filterOut)

        #Add in the initial quality and filtered quality numbers for reporting
        track<-rbind.fill(as.data.frame(filterOut[1,,drop=FALSE], check.names=FALSE),as.data.frame(fwdInputTempPlotValuesOut, check.names=FALSE))

        #Row names vector
        rowNamesTemp <- c("inputReads", gsub("_primeTrim","",rownames(fwdInputTempPlotValuesOut)))

        if(bidirectional==TRUE){

          #Add in the rev input values
          track<-rbind.fill(as.data.frame(track, check.names=FALSE),as.data.frame(revInputTempPlotValuesOut, check.names=FALSE))

          #Row names vector
          rowNamesTemp <- c(rowNamesTemp, rownames(revInputTempPlotValuesOut))

        }

        #Add in the filtered read numbers
        track<-rbind.fill(as.data.frame(track, check.names=FALSE),as.data.frame(filterOut[2,,drop=FALSE], check.names=FALSE))

        #Row names vector
        rowNamesTemp <- c(rowNamesTemp, "filteredReads")

        #Check to see if there is more than one sample being analysed.
        if (ncol(track)>1){
          #Get the fwd denoised values
          fwdDenoisedReads = as.data.frame(sapply(dadaF1, getN), check.names=FALSE)
          row.names(fwdDenoisedReads)<-gsub("_fwdFilt.fastq.gz","",row.names(fwdDenoisedReads))
          fwdDenoisedReads<-as.data.frame(t(fwdDenoisedReads), check.names=FALSE)
          track<-rbind.fill(track, fwdDenoisedReads)
          rowNamesTemp<-c(rowNamesTemp,"fwdDenoisedReads")

          #get the denoised quality values
          fwdDenoisedAvgQual <- round(as.data.frame(mapply(getSeqQual, dadaF1)),2)
          row.names(fwdDenoisedAvgQual)<-gsub("_fwdFilt.fastq.gz","",row.names(fwdDenoisedAvgQual))
          fwdDenoisedAvgQual<-as.data.frame(t(fwdDenoisedAvgQual), check.names=FALSE)
          track<-rbind.fill(track, fwdDenoisedAvgQual)
          rowNamesTemp<-c(rowNamesTemp,"fwdDenoisedAvgQual")

          #Get the denoised average length
          fwdDenoisedAvgLen <- round(as.data.frame(mapply(getMeanSeqLen, dadaF1)),2)
          row.names(fwdDenoisedAvgLen)<-gsub("_fwdFilt.fastq.gz","",row.names(fwdDenoisedAvgLen))
          fwdDenoisedAvgLen<-as.data.frame(t(fwdDenoisedAvgLen), check.names=FALSE)
          track<-rbind.fill(track, fwdDenoisedAvgLen)
          rowNamesTemp<-c(rowNamesTemp,"fwdDenoisedAvgLen")

          #Get the denoised unique sequences
          fwdDenoisedUniqueReads <- as.data.frame(mapply(getSeqLength, dadaF1), check.names=FALSE)
          row.names(fwdDenoisedUniqueReads)<-gsub("_fwdFilt.fastq.gz","",row.names(fwdDenoisedUniqueReads))
          fwdDenoisedUniqueReads<-as.data.frame(t(fwdDenoisedUniqueReads), check.names=FALSE)
          track<-rbind.fill(track, fwdDenoisedUniqueReads)
          rowNamesTemp<-c(rowNamesTemp,"fwdDenoisedUniqueReads")

          if(bidirectional==TRUE){

            #Get the rev denoised values
            revDenoisedReads = as.data.frame(sapply(dadaR1, getN), check.names=FALSE)
            row.names(revDenoisedReads)<-gsub("_revFilt.fastq.gz","",row.names(revDenoisedReads))
            revDenoisedReads<-as.data.frame(t(revDenoisedReads), check.names=FALSE)
            track<-rbind.fill(track, revDenoisedReads)
            rowNamesTemp<-c(rowNamesTemp,"revDenoisedReads")

            #get the denoised quality values
            revDenoisedAvgQual <- round(as.data.frame(mapply(getSeqQual, dadaR1), check.names=FALSE),2)
            row.names(revDenoisedAvgQual)<-gsub("_revFilt.fastq.gz","",row.names(revDenoisedAvgQual))
            revDenoisedAvgQual<-as.data.frame(t(revDenoisedAvgQual), check.names=FALSE)
            track<-rbind.fill(track, revDenoisedAvgQual)
            rowNamesTemp<-c(rowNamesTemp,"revDenoisedAvgQual")

            #Get the denoised average length
            revDenoisedAvgLen <- round(as.data.frame(mapply(getMeanSeqLen, dadaR1), check.names=FALSE),2)
            row.names(revDenoisedAvgLen)<-gsub("_revFilt.fastq.gz","",row.names(revDenoisedAvgLen))
            revDenoisedAvgLen<-as.data.frame(t(revDenoisedAvgLen), check.names=FALSE)
            track<-rbind.fill(track, revDenoisedAvgLen)
            rowNamesTemp<-c(rowNamesTemp,"revDenoisedAvgLen")

            #Get the denoised unique sequences
            revDenoisedUniqueReads <- as.data.frame(mapply(getSeqLength, dadaR1), check.names=FALSE)
            row.names(revDenoisedUniqueReads)<-gsub("_revFilt.fastq.gz","",row.names(revDenoisedUniqueReads))
            revDenoisedUniqueReads<-as.data.frame(t(revDenoisedUniqueReads), check.names=FALSE)
            track<-rbind.fill(track, revDenoisedUniqueReads)
            rowNamesTemp<-c(rowNamesTemp,"revDenoisedUniqueReads")

          }

        }else{

          if(bidirectional==TRUE){

            #Add the values to the reporting track dataframe for a single sample, first get the values
            trackTemp <- as.data.frame(c(getN(dadaF1), getSeqQual(dadaF1), getMeanSeqLen(dadaF1), getSeqLength(dadaF1), getN(dadaR1), getSeqQual(dadaR1), getMeanSeqLen(dadaR1), getSeqLength(dadaR1)))

            #Add rownames to the temp data frame
            rowNamesTemp <- c(rowNamesTemp, "fwdDenoisedReads", "fwdDenoisedAvgQual", "fwdDenoisedAvgLen", "fwdDenoisedUniqueReads", "revDenoisedReads", "revDenoisedAvgQual", "revDenoisedAvgLen", "revDenoisedUniqueReads")

          }else{

            #Add the values to the reporting track dataframe for a single sample, first get the values
            trackTemp <- as.data.frame(c(getN(dadaF1), getSeqQual(dadaF1), getMeanSeqLen(dadaF1), getSeqLength(dadaF1)))

            #Add rownames to the temp data frame
            rowNamesTemp <- c(rowNamesTemp, "fwdDenoisedReads", "fwdDenoisedAvgQual", "fwdDenoisedAvgLen", "fwdDenoisedUniqueReads")

          }

          #make the trackTemp colname the same as track
          colnames(trackTemp)<-colnames(track)
          #Add the values onto the track data frame
          track <- rbind(track, trackTemp)

        }

  ################################################################################
        if(bidirectional==TRUE){
          #Get the values for the merged reads coming out of the analysis
          if(nrow(mergedTable)==0){
            mergedReads=as.data.frame(rep(0,ncol(mergedTable)-1), check.names=FALSE)
            mergedReadsAvgLen=as.data.frame(rep(0,ncol(mergedTable)-1), check.names=FALSE)
          }else{
            #get the col name
            mergedReads = as.data.frame(colSums(sapply(as.data.frame(mergedTable[,-c(1)]), as.numeric)), check.names=FALSE)
            #If there is only one record then the name of the sample is dropped in the colSums so we need to add it back
            if(length(fwdFile$fileNamesNoDir)==1){
              row.names(mergedReads)<-fwdFile$fileNamesNoDir
            }

            #mergedReadsAvgLen <- lapply(mergedTable[,c(2:ncol(mergedTable))], function(x) rownames(mergedTable)[x>0])
            #Get the sample columns in to a variable
            mergedReadsAvgLen <- mergedTable[,2:ncol(mergedTable),drop=FALSE]
            #Make the mergedReadsAvgLen values numeric
            mergedReadsAvgLen <- data.frame(apply(mergedReadsAvgLen, 2, function(x) as.numeric(as.character(x))), check.names=FALSE)
            row.names(mergedReadsAvgLen)<-row.names(mergedTable)

            #Place the length of the read in the dataframe where the number of reads go first get the indices
            k <- which(mergedReadsAvgLen>0, arr.ind=TRUE)
            #Then replace all positive values with the number of characters in row names
            mergedReadsAvgLen[k] <- nchar(row.names(k))
            #Now set all 0 values to NA
            mergedReadsAvgLen[mergedReadsAvgLen==0] <-NA

            #Get the mean length of the columns which are the samples
            mergedReadsAvgLen <- as.data.frame(colMeans(mergedReadsAvgLen, na.rm = TRUE), check.names=FALSE)

            #If there is only one record then the name of the sample is dropped in the colMeans so we need to add it back
            if(length(fwdFile$fileNamesNoDir)==1){
              row.names(mergedReadsAvgLen)<-fwdFile$fileNamesNoDir
            }
        }

          row.names(mergedReads)<-gsub("_fwdFilt.fastq.gz","",row.names(mergedReads))
          row.names(mergedReadsAvgLen)<-gsub("_fwdFilt.fastq.gz","",row.names(mergedReadsAvgLen))
          row.names(mergedReads)<-gsub("\\.","-",row.names(mergedReads))
          row.names(mergedReadsAvgLen)<-gsub("\\.","-",row.names(mergedReadsAvgLen))
          mergedReads<-as.data.frame(t(mergedReads), check.names=FALSE)
          mergedReadsAvgLen<-as.data.frame(t(mergedReadsAvgLen), check.names=FALSE)
          track<-rbind.fill(track, mergedReads)
          track<-rbind.fill(track, mergedReadsAvgLen)
          rowNamesTemp<-c(rowNamesTemp,"mergedReads")
          rowNamesTemp<-c(rowNamesTemp,"mergedReadsAvgLen")

    ################################################################################

          #Get the reads that were trimmed using the pattern matching
          noChimReadsdf<-as.data.frame(finalTable[finalTable$Results =="Merged",,drop=F], check.names=FALSE)
          noChimReadsdf<-as.data.frame(noChimReadsdf[noChimReadsdf$Results != "ChimRemoved",,drop=F], check.names=FALSE)

          #Get the values for the no chimeric reads coming out of the analysis
          if(nrow(noChimReadsdf)==0){
            noChimReads=as.data.frame(rep(0,ncol(noChimReadsdf)-3), check.names=FALSE)
          }else{
             noChimReads = as.data.frame(colSums(sapply(noChimReadsdf[,4:ncol(noChimReadsdf), drop=FALSE], as.numeric)), check.names=FALSE)
            if(length(fwdFile$fileNamesNoDir)==1){
              row.names(noChimReads)<-fwdFile$fileNamesNoDir
            }
          }

          row.names(noChimReads)<-gsub("_fwdFilt.fastq.gz","",row.names(noChimReads))
          row.names(noChimReads)<-gsub("_revFilt.fastq.gz","",row.names(noChimReads))
          row.names(noChimReads)<-gsub("\\.","-",row.names(noChimReads))
          noChimReads<-as.data.frame(t(noChimReads), check.names=FALSE)
          track<-rbind.fill(track, noChimReads)
          rowNamesTemp<-c(rowNamesTemp,"noChimReads")

          #No Chimeric above min length reporting
          noChimTrimReadsOverMinLendf<-data.frame(noChimReadsdf[as.numeric(noChimReadsdf$Length) > minFinalSeqLen,,drop=F], check.names=FALSE)
          if(nrow(noChimTrimReadsOverMinLendf)==0){

            noChimTrimReadsOverMinLen=as.data.frame(rep(0,ncol(noChimTrimReadsOverMinLendf)-3), check.names=FALSE)
            noChimTrimUniqueReadsOverMinLen=as.data.frame(rep(0,ncol(noChimTrimReadsOverMinLendf)-3), check.names=FALSE)
            noChimTrimAvgLenOverMinLen=as.data.frame(rep(0,ncol(noChimTrimReadsOverMinLendf)-3), check.names=FALSE)
          }else if(nrow(noChimTrimReadsOverMinLendf)==1){

            noChimTrimReadsOverMinLen = t(as.data.frame(noChimTrimReadsOverMinLendf[,4:ncol(noChimTrimReadsOverMinLendf)], drop=FALSE, check.names=FALSE))
            noChimTrimUniqueReadsOverMinLen = t(as.data.frame(noChimTrimReadsOverMinLendf[,4:ncol(noChimTrimReadsOverMinLendf)], drop=FALSE, check.names=FALSE))
            noChimTrimUniqueReadsOverMinLen[noChimTrimUniqueReadsOverMinLen>0]<-1
            noChimTrimAvgLenOverMinLen = t(as.data.frame(noChimTrimReadsOverMinLendf[,4:ncol(noChimTrimReadsOverMinLendf)], drop=FALSE, check.names=FALSE))
            noChimTrimAvgLenOverMinLen[noChimTrimAvgLenOverMinLen>0]<-unique(noChimTrimReadsOverMinLendf$Length)
          }else{

            noChimTrimReadsOverMinLen = as.data.frame(colSums(sapply(noChimTrimReadsOverMinLendf[,4:ncol(noChimTrimReadsOverMinLendf), drop=FALSE], as.numeric)), check.names=FALSE)
            noChimTrimUniqueReadsOverMinLen = as.data.frame(colSums(noChimTrimReadsOverMinLendf[,4:ncol(noChimTrimReadsOverMinLendf), drop=FALSE] != 0), check.names=FALSE)

            #Get the sample columns in to a variable
            noChimTrimAvgLenOverMinLen <- noChimTrimReadsOverMinLendf[,-c(1,2,3), drop=FALSE]

            #Place the length of the read in the dataframe where the number of reads go first get the indices
            k <- which(noChimTrimAvgLenOverMinLen>0, arr.ind=TRUE)
            #Then for the indices place the length of the row name
            noChimTrimAvgLenOverMinLen[k] <- nchar(row.names(k))
            #Now set all 0 values to NA
            noChimTrimAvgLenOverMinLen[noChimTrimAvgLenOverMinLen==0] <-NA

            #Make the noChimTrimAvgLenOverMinLen numeric
            noChimTrimAvgLenOverMinLen <- data.frame(apply(noChimTrimAvgLenOverMinLen, 2, function(x) as.numeric(as.character(x))), check.names=FALSE)
            #Get the mean length of the columns which are the samples
            noChimTrimAvgLenOverMinLen <- data.frame(colMeans(noChimTrimAvgLenOverMinLen, na.rm = TRUE), check.names=FALSE)

            #If there is only one record then the name of the sample is dropped in the colMeans so we need to add it back
            if(length(fwdFile$fileNamesNoDir)==1){
              print("Here 30E")
              row.names(noChimTrimReadsOverMinLen)<-fwdFile$fileNamesNoDir
              row.names(noChimTrimUniqueReadsOverMinLen)<-fwdFile$fileNamesNoDir
              row.names(noChimTrimAvgLenOverMinLen)<-fwdFile$fileNamesNoDir
            }
          }

          row.names(noChimTrimReadsOverMinLen)<-gsub("_fwdFilt.fastq.gz","",row.names(noChimTrimReadsOverMinLen))
          row.names(noChimTrimReadsOverMinLen)<-gsub("\\.","-",row.names(noChimTrimReadsOverMinLen))
          row.names(noChimTrimUniqueReadsOverMinLen)<-gsub("_fwdFilt.fastq.gz","",row.names(noChimTrimUniqueReadsOverMinLen))
          row.names(noChimTrimUniqueReadsOverMinLen)<-gsub("\\.","-",row.names(noChimTrimUniqueReadsOverMinLen))
          row.names(noChimTrimAvgLenOverMinLen)<-gsub("_fwdFilt.fastq.gz","",row.names(noChimTrimAvgLenOverMinLen))
          row.names(noChimTrimAvgLenOverMinLen)<-gsub("\\.","-",row.names(noChimTrimAvgLenOverMinLen))

    ################################################################################

          noChimTrimReadsOverMinLen<-as.data.frame(t(noChimTrimReadsOverMinLen), check.names=FALSE)
          noChimTrimUniqueReadsOverMinLen<-as.data.frame(t(noChimTrimUniqueReadsOverMinLen), check.names=FALSE)
          noChimTrimAvgLenOverMinLen<-as.data.frame(t(noChimTrimAvgLenOverMinLen), check.names=FALSE)
          track<-rbind.fill(track, noChimTrimReadsOverMinLen)
          track<-rbind.fill(track, noChimTrimUniqueReadsOverMinLen)
          track<-rbind.fill(track, noChimTrimAvgLenOverMinLen)
          rowNamesTemp<-c(rowNamesTemp,"noChimTrimReadsOverMinLen")
          rowNamesTemp<-c(rowNamesTemp,"noChimTrimUniqueReadsOverMinLen")
          rowNamesTemp<-c(rowNamesTemp,"noChimTrimAvgLenOverMinLen")

          #Add on the row names
          row.names(track)<-rowNamesTemp

          #Get the final percent of reads from the initial number and then add this to the output tracking data frame
          finalCleanedPerReads <- (round((as.numeric(track[rownames(track) %in% "noChimTrimReadsOverMinLen",])/as.numeric(track[rownames(track) %in% "inputReads",])), 4))

        }else{#End of the reporting if bidirectional is TRUE
          #Add on the row names
          row.names(track)<-rowNamesTemp
          finalCleanedPerReads <- (round((as.numeric(track[rownames(track) %in% "fwdDenoisedReads",])/as.numeric(track[rownames(track) %in% "inputReads",])), 4))
        }

        finalCleanedPerReads <- t(data.frame(colnames(track), finalCleanedPerReads = finalCleanedPerReads, check.names=FALSE) )
        colnames(finalCleanedPerReads)<-finalCleanedPerReads[1,]
        finalCleanedPerReads<-as.data.frame(finalCleanedPerReads[-1,,drop=FALSE], check.names=FALSE)

        #Add the percent reads to the tracking table for reporting
        track<-rbind.fill(track,finalCleanedPerReads)
        #Add on the row names
        rowNamesTemp<-c(rowNamesTemp,"finalCleanedPerReads")
        row.names(track)<-rowNamesTemp

        #Print to log file
        suppressWarnings(write.table(t(track), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummaryTable.tsv"), append = FALSE, na = "NA", row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t"))
        print(paste0(startTime, " - End time...", Sys.time()))
        suppressWarnings(write(paste0("At the end...", Sys.time()), file = paste0(outFolder,"/", dateStamp, dirName, "_dadaSummary.txt"), append = TRUE))


      }# End of loop going through the different files

       #Clean up the variables and memory from this loop
      suppressWarnings(rm(chim,
          dadaF1,
          dadaR1,
          dateStamp,
          dirName,
          errCheckFiles,
          errF,
          errR,
          fileList,
          fileNames,
          fileNamesNoDir,
          filesTbl,
          filteredFolder,
          filterOut,
          filterQualFolder,
          finalCleanedPerReads,
          finalTable,
          finalTableOut,
          fwdDenoisedAvgLen,
          fwdDenoisedAvgQual,
          fwdDenoisedReads,
          fwdDenoisedUniqueReads,
          fwdDerep,
          fwdFile,
          fwdFilt,
          fwdFilteredTempPlotValues,
          fwdFilteredTempPlotValuesOut,
          fwdFilteredTempPlotValuesOutColNames,
          fwdFiltError,
          fwdFiltOutFiles,
          fwdFiltOutFilesName,
          fwdFiltTbl,
          fwdInputTempPlotValues,
          fwdInputTempPlotValuesOut,
          fwdInputTempPlotValuesOutColNames,
          merged,
          mergedReads,
          mergedReadsAvgLen,
          mergedTable,
          noChim,
          noChimReads,
          noChimReadsdf,
          noChimTrimReadsOverMinLendf,
          nochimTblTrimmed,
          nochimTrimmedOut,
          noChimTrimReadsOverMinLen,
          noChimTrimUniqueReadsOverMinLen,
          notShort,
          numErrFiles,
          outFolder,
          pathList,
          poorQual,
          qualityFolder,
          revDenoisedAvgLen,
          revDenoisedAvgQual,
          revDenoisedReads,
          revDenoisedUniqueReads,
          revDerep,
          revFile,
          revFilteredTempPlotValues,
          revFilteredTempPlotValuesOut,
          revFilteredTempPlotValuesOutColNames,
          revFiltError,
          revFiltOutFiles,
          revFiltOutFilesName,
          revFiltTbl,
          revInputTempPlotValues,
          revInputTempPlotValuesOut,
          revInputTempPlotValuesOutColNames,
          rowNamesTemp,
          seqOut,
          startTime,
          totalFiltOutFiles,
          track,
          workLoc
       ))

      #Clean memory
      invisible(gc())

    }#End of the main if statements that checks that the proper data are submitted to run the script
  }#End of the if else where the unidirectional and bidirectional are checked.
} #End of the Function
