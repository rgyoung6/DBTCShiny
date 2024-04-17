# Written by Rob Young at the University of Guelph in Ontario Canada, April, 2024
# ******************************************************************************
# Roxygen2 Documentation:
#' @export
#'
#' @title Assign Taxa using BLAST Results
#'
#' @author Robert G. Young
#'
#' @description
#' This function takes a BLAST result file and associated fasta files (either on
#' their own or with accompanying ASV files generated from the dada_implement
#' function) and collapses the multiple BLAST results into as single result for
#' each query sequence. When an ASV table is present the taxonomic results will
#' be combined with the ASV table.
#'
#' @details
#' This function requires a BLAST output file and an associated fasta file. In
#' addition, if present an ASV file will also be used and combined with the
#' taxonomic results when present. The BLAST results are reduced to a single
#' result for each read. At each taxonomic level there may be one or more
#' taxonomic assignments. Each assignment has quality metrics in parentheses after
#' the name. These values ("Num_Rec", "Coverage", "Identity", "Max_eVal") represent
#' the number of records with this taxonomic placement, the minimum coverage and
#' identity, and the maximum eValue for the reported taxa.
#'
#' @examples
#' \dontrun{
#' taxon_assign()
#' taxon_assign(fileLoc = NULL, taxaDBLoc = NULL, numCores = 1, coverage = 95,
#' ident = 95, propThres = 0.95, coverReportThresh=0, identReportThresh=0, includeAllDada=TRUE)
#' }
#'
#' @param fileLoc The location of a file in a directory where all of the paired
#' fasta and BLAST (and potentially ASV) files are located (Default NULL).
#' @param taxaDBLoc The location of the NCBI taxonomic data base (Default NULL;
#' for accessionTaxa.sql see the main DBTC page for details).
#' @param numCores The number of cores used to run the function (Default 1,
#' Windows systems can only use a single core)
#' @param coverage The percent coverage used for taxonomic assignment for the
#' above threshold results (Default 95).
#' @param ident The percent identity used for the taxonomic assignment for above
#' threshold results (Default 95).
#' @param propThres The proportional threshold flags the final result based on
#' the preponderance of the data. So if the threshold is set to 0.95, results
#' will be flagged if the taxa directly below the assigned taxa has fewer than
#' 0.95 percent of the records causing the upward taxonomic placement (Default 0.95).
#' @param coverReportThresh The percent coverage threshold used for reporting
#' flags below this threshold (Default 95).
#' @param identReportThresh The percent identity threshold used for reporting
#' flags below this threshold (Default 95).
#' @param includeAllDada When paired Dada ASV tables are present, when set to
#' FALSE, this will exclude records without taxonomic assignment (Default TRUE).
#'
#' @returns
#' This function produces a taxa_reduced file for each submitted BLAST-fasta submission.
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
#' combine_dada_output()
#' make_BLAST_DB()
#' seq_BLAST()
#' combine_assign_output()
#' reduce_taxa()
#' combine_reduced_output()

##################################### taxon_assign FUNCTION ##############################################################
taxon_assign<- function(fileLoc = NULL, taxaDBLoc = NULL, numCores = 1, coverage = 95, ident = 95, propThres = 0.95, coverReportThresh=0, identReportThresh=0, includeAllDada=TRUE){

  #If there are issues and I need to audit the script make this 1
  auditScript=0

  #Get the initial working directory
  start_wd <- getwd()
  on.exit(setwd(start_wd))

  if (is.null(fileLoc)){
    print(paste0("Select a file in the file folder with the BLAST output file(s), associated fas file, and associated Dada ASV output files (if present)."))
    fileLoc <- file.choose()
  }

  #Set the location of the BLAST taxonomic database
  if(is.null(taxaDBLoc)){
    print(paste0("Select the data base file with the NCBI taxonomic database (i.e. accessionTaxa.sql)"))
    taxaDBLoc <- file.choose()
  }

  #Audit line
  if(auditScript>0){
    auditFile <- paste0(dirname(fileLoc),"/", format(Sys.time(), "%Y_%m_%d_%H%M"), "_audit.txt")
    write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 1"), file = auditFile, append = FALSE)
    print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 1"))
  }

  #Get the directory of interest
  fileLoc <- dirname(fileLoc)
  #Set the current working directory to be the location where the data base file is located.
  setwd(dirname(taxaDBLoc))

  #Printing the start time
  print(paste0("Start time...", Sys.time()))
  startTime <- paste0("Start time...", Sys.time())
  dateStamp <- paste0(format(Sys.time(), "%Y_%m_%d_%H%M"))

  #Setting up the table of BLAST results to work through (with at times associated asv tables)
  files <- list.files(path = fileLoc, pattern = "*[.]*")
  fullPathFiles<-list.files(path = fileLoc, pattern = "*[.]*", full.names = TRUE)
  files <- as.data.frame(cbind(fullPathFiles, files))

  if(max(lengths(regmatches(fullPathFiles, gregexpr("_BLAST", fullPathFiles))))<=1){

    #Get the BLAST files and generate the dada output files from the BLAST naming conventions
    filesList <- files[grepl("_BLAST", files[,2]),c(1,2)]

    filesList <- filesList[!grepl(".sh$", filesList[,2]),c(1,2)]
    filesList <- filesList[!grepl(".bat$", filesList[,2]),c(1,2)]

    #Add the fileNames to the front
    filesList <- cbind(gsub("_BLAST.*","",filesList[,2]), filesList)

    #Add the asvFilesFullPath and asvFiles
    filesList <- cbind(filesList, gsub("_BLAST.*",".tsv",filesList[,2]), gsub("_BLAST.*",".tsv",filesList[,3]))

    #Add in a column to verify that the ASV file exists
    filesList <- cbind(filesList, lapply(filesList[4], file.exists))

    #Add the fasFilesFullPath and fasFiles
    filesList <- cbind(filesList, gsub("_BLAST.*",".fas",filesList[,2]), gsub("_BLAST.*",".fas",filesList[,3]))

    #Add in a column to verify that the FAS file exists
    filesList <- cbind(filesList, lapply(filesList[7], file.exists))
    filesList <- cbind(filesList, gsub("_.*", "", gsub(".*(_BLAST_)", "", filesList[,3])))
    colnames(filesList) <- c("fileNames", "blastFilesFullPath", "blastFiles", "asvFilesFullPath", "asvFiles", "asvExists", "fasFilesFullPath", "fasFiles","fasExists", "databaseName" )

    #Remove the rows with FALSE in the fasExists column
    filesList <- filesList[filesList$fasExists, ]

    if(nrow(filesList)>0){

      #Check to see if there are files where there are both no fasta or asv file associated with the BLAST results
      positions <- which(filesList$asvExists == FALSE & filesList$fasExists == FALSE)

      if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 2")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 2"), file = auditFile, append = TRUE))}

      if(length(positions)==0){

        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 3")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 3"), file = auditFile, append = TRUE))}

        #initialize the output variable
        finalCondensedOut<-NULL

        #Set up looping through the files in the submitted folders by pairing the fasta to the Dada ASV tables where present
        for(fileCounter in 1:nrow(filesList)){

          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 4")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 4"), file = auditFile, append = TRUE))}

          #Read in the BLAST results for the file in this loop
          blastResults <- read.delim(filesList[fileCounter,2], header = FALSE)

          #Add in the headers
          colnames(blastResults) <- c("uniqueID","sequence_id","taxa_id","query_coverage","percent_identity","s_name","c_name","q_start","q_end","s_start","s_end","evalue")

          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 5")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 5"), file = auditFile, append = TRUE))}

          # This if checks to see if the output is a MACER formatted database or a NCBI
          # database by checking the number of |. There are 3 in Custom and 4 in NCBI
          if(sum(gregexpr("|", blastResults[1,2], fixed=TRUE)[[1]] > 0)<4){

            if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 6")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 6"), file = auditFile, append = TRUE))}

            #Take the sequence_id column and split into columns
            blastResultsTemp <- data.frame(do.call('rbind', strsplit(as.character(blastResults[,2]),'|', fixed = TRUE)))

            #Take the second column from the sequence_id variables that holds the taxa_id's and place in the taxa id column
            blastResults$taxa_id <- blastResultsTemp[,2]

            #Take columns 3 and 4 and place them in the s_name column
            blastResults$s_name <- paste0(blastResultsTemp[,3]," ",blastResultsTemp[,4])

            #Remove the blastResultsTemp to clear memory
            rm(blastResultsTemp)

          }

          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 7")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 7"), file = auditFile, append = TRUE))}

          # get all of the unique uniqueID's
          blastResultsuniqueID <- unique(blastResults$uniqueID)

          #Get the start time to use in the reporting of the progress throughout
          startTime <- Sys.time()

          print("********************************************************************************")
          print(paste0("Starting analysis ", filesList[fileCounter,3], ": ", fileCounter, " of ", nrow(filesList)," at ", Sys.time()))
          print("********************************************************************************")

          #Create the output dataframe
          condensedOut <- data.frame(uniqueID = character(), superkingdom = character(), phylum = character(), class = character(), order = character(), family = character(), genus = character(), species = character(), Top_BLAST = character(), Lowest_Single_Rank = character(), Lowest_Single_Taxa = character(), Lowest_Single_Rank_Above_Thres = character(), Lowest_Single_Taxa_Above_Thres = character(), Final_Common_Names = character(), Final_Rank = character(), Final_Taxa = character(), Final_Rank_Taxa_Thres = character(), Result_Code = character())
          write.table(as.data.frame(condensedOut), file = paste0(fileLoc, "/", as.vector(filesList[fileCounter,1]), "_", filesList[fileCounter,10], "_taxaAssign_", dateStamp, ".tsv"), row.names = FALSE, col.names = TRUE, append = FALSE, sep = "\t", quote = FALSE)

          #Make a list of the ranks
          taxaRank <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")

          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 8")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 8"), file = auditFile, append = TRUE))}

#          for(records in 1:length(blastResultsuniqueID)){
          taxaResults <- function(records){
            #subset the data for the target unique ID for the loop and identity and coverage values above define threshold
            blastResultsTarget <- blastResults[blastResults$uniqueID==blastResultsuniqueID[records],, drop=FALSE]

            if(nrow(blastResultsTarget) > 0){

              Top_BLAST <- paste0(blastResultsTarget[1,c("s_name")], "(", blastResultsTarget[1,c("query_coverage")], ",", blastResultsTarget[1,c("percent_identity")],",", blastResultsTarget[1,c("evalue")],")")

              #Create the output dataframe
              condensedOut <- data.frame(uniqueID = character(), superkingdom = character(), phylum = character(), class = character(), order = character(), family = character(), genus = character(), species = character(), Top_BLAST = character(), Lowest_Single_Rank = character(), Lowest_Single_Taxa = character(), Lowest_Single_Rank_Above_Thres = character(), Lowest_Single_Taxa_Above_Thres = character(), Final_Common_Names = character(), Final_Rank = character(), Final_Taxa = character(), Result_Code = character())
              condensedOutThres <- data.frame(uniqueID = character(), superkingdom = character(), phylum = character(), class = character(), order = character(), family = character(), genus = character(), species = character(), Top_BLAST = character(), Lowest_Single_Rank = character(), Lowest_Single_Taxa = character(), Lowest_Single_Rank_Above_Thres = character(), Lowest_Single_Taxa_Above_Thres = character(), Final_Common_Names = character(), Final_Rank = character(), Final_Taxa = character(), Result_Code = character())

              #Add in the unique identifier for this read to the condensedOut data frame
              condensedOut[1,"uniqueID"]<-blastResultsuniqueID[records]
              condensedOutThres[1,"uniqueID"]<-blastResultsuniqueID[records]

              #Adding the first BLAST hit to the output data
              condensedOut[1,"Top_BLAST"] <- Top_BLAST
              condensedOutThres[1,"Top_BLAST"] <- Top_BLAST

              #Get the total number of returned results for the current unique query ID
              numBLASTResults<-nrow(blastResultsTarget)

              #Run the taxonomizr taxonomy function
              suppressWarnings(taxa <- taxonomizr::getTaxonomy(blastResultsTarget$taxa_id, taxaDBLoc))
              taxa <- as.data.frame(taxa)

              if(nrow(taxa) == 1){
                #get the row name for use later
                oneRowID = row.names(taxa)
                #Setting any result with a non alpha character to NA
                taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".*[[:punct:]].*", NA, y)))
                taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".*[[:digit:]].*", NA, y)))
                #with one row the data get transformed so this is just placing it back in the proper dataframe format
                taxa <- as.data.frame(t(taxa))
                row.names(taxa) <- oneRowID
              }else{
                #Setting any result with a non alpha character to NA
                taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".*[[:punct:]].*", NA, y)))
                taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".*[[:digit:]].*", NA, y)))
              }

              #create a function to return FALSE if there is more than two spaces in the species name
              countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) } != 1
              #using the countSpaces function change all species records with two spaces to NA
              taxa$species[countSpaces(taxa$species)] = NA

              #create a function to return FALSE if there is more than one space in the higher level taxa
              countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) } != 0
              #using the countSpaces function change all higher taxa records with one spaces to NA
              taxa$genus[countSpaces(taxa$genus)] = NA
              taxa$family[countSpaces(taxa$family)] = NA
              taxa$order[countSpaces(taxa$order)] = NA
              taxa$class[countSpaces(taxa$class)] = NA
              taxa$phylum[countSpaces(taxa$phylum)] = NA
              taxa$superkingdom[countSpaces(taxa$superkingdom)] = NA

              #Change all taxonomic ranks below the highest level of NA to NA
              taxa[t(apply(is.na(taxa), 1, cumsum)) > 0 ] <- NA

              #Add on the raw BLAST results to the taxonomy output table
              taxa <- cbind(taxa, blastResultsTarget)

              #Subset the taxa dataframe to only have records with greater than the defined query coverage and identity
              taxaAboveFilter <- taxa[taxa$query_coverage > coverage,, drop=FALSE]
              taxaAboveFilter <- taxaAboveFilter[taxaAboveFilter$percent_identity > ident,, drop=FALSE]

              #Loop through each of the ranks
              for(taxaRankCount in 1:length(taxaRank)){
                #Get the unique taxa for this rank from the main dataset.
                uniqueTaxa <- unique(taxa[[taxaRank[taxaRankCount]]])[!unique(taxa[[taxaRank[taxaRankCount]]]) %in% NA]

                #Main results using the main dataset
                if(length(uniqueTaxa) > 0){
                  #subset taxa for the target rank and name to get the working dataframe and the number of records at the rank(s)
                  taxaWork <- taxa[taxa[[taxaRank[taxaRankCount]]] %in% uniqueTaxa,c(taxaRank[taxaRankCount], "query_coverage", "percent_identity", "evalue")]
                  taxaNumRecords <- as.data.frame(table(taxaWork[[taxaRank[taxaRankCount]]]))
                  colnames(taxaNumRecords) <- c(taxaRank[taxaRankCount], "Num_records")

                  #Get the min query coverage value for the taxa specific dataset
                  taxaQuery <- as.data.frame(tapply(as.vector(as.numeric(taxaWork$query_coverage)), taxaWork[[taxaRank[taxaRankCount]]], min))
                  #make the row name the first column
                  taxaQuery <- data.frame(uniqueID = row.names(taxaQuery), taxaQuery)
                  colnames(taxaQuery) <- c(taxaRank[taxaRankCount],"taxaQuery")
                  taxaTotal <- merge(taxaNumRecords, taxaQuery, by = taxaRank[taxaRankCount])

                  #Get the min percent identity for the target taxa specific dataset
                  taxaPerIdent <- as.data.frame(tapply(as.vector(as.numeric(taxaWork$percent_identity)), taxaWork[[taxaRank[taxaRankCount]]], min))
                  #make the row name the first column
                  taxaPerIdent <- data.frame(uniqueID = row.names(taxaPerIdent), taxaPerIdent)
                  colnames(taxaPerIdent) <- c(taxaRank[taxaRankCount],"taxaPerIdent")
                  taxaTotal <- merge(taxaTotal, taxaPerIdent, by = taxaRank[taxaRankCount])

                  #Get the max e value for the target taxa specific dataset
                  taxaEval <- as.data.frame(tapply(as.vector(as.numeric(taxaWork$evalue)), taxaWork[[taxaRank[taxaRankCount]]], max))
                  #make the row name the first column
                  taxaEval <- data.frame(uniqueID = row.names(taxaEval), taxaEval)
                  colnames(taxaEval) <- c(taxaRank[taxaRankCount],"taxaEval")
                  taxaTotal <- merge(taxaTotal, taxaEval, by = taxaRank[taxaRankCount])

                  #add the number of records, query coverage, percent identity, and e values on to the taxa specific dataset
                  taxaTotal <- as.data.frame(cbind(taxaTotal, paste0("(",round(taxaTotal$Num_records,2),",",round(taxaTotal$taxaQuery,2),",",round(taxaTotal$taxaPerIdent,2),",",format(taxaTotal$taxaEval, digits = 2),")")))

                  # Remove all white spaces from the column
                  taxaTotal[,ncol(taxaTotal)]<- gsub("\\s+", "", taxaTotal[,ncol(taxaTotal)])

                  #adding the values in brackets back onto the taxa in a dataframe for merging below
                  taxaTotal <- taxaTotal[,c(1, ncol(taxaTotal))]
                  condenseWork <- paste0(taxaTotal[,1],taxaTotal[,2])

                  if(length(na.omit(uniqueTaxa)) > 1){
                    #Place the results in the correct rank column when multiple results are present
                    condensedOut[1,taxaRank[taxaRankCount]] <- paste0(condenseWork, collapse = ", ")
                  }else if (length(na.omit(uniqueTaxa))== 1){
                    #Place the results in the correct rank column when a single result is present
                    condensedOut[1,c(taxaRank[taxaRankCount])] <- condenseWork
                    #If the ranks above the current loop rank do not contain NA or multiple results then assign the lowest values
                    valueCheck<-condensedOut[1,c(which(colnames(condensedOut) == "superkingdom"):which(colnames(condensedOut)== taxaRank[taxaRankCount])),drop=FALSE]
                    if(sum(is.na(valueCheck[1,]))==0){
                      if(sum(grep("),",valueCheck[1,]))==0){
                        condensedOut[1,"Lowest_Single_Rank"] <- taxaRank[taxaRankCount]
                        condensedOut[1,"Lowest_Single_Taxa"] <- paste0(condenseWork, collapse = ", ")
                      }
                    }
                  }

                }#end of the if uniqueTaxa> 0

            ############################ Repeat the above with filters applied (where data is available) #############

                #Get the unique taxa for this rank for the data above the query and identity thresholds
                uniqueTaxaAboveFilter <- unique(taxaAboveFilter[[taxaRank[taxaRankCount]]])[!unique(taxaAboveFilter[[taxaRank[taxaRankCount]]]) %in% NA]

                #Main results using the main dataset with filters
                if(length(uniqueTaxaAboveFilter) > 0){

                  #subset taxa for the target rank and name to get the working dataframe and the number of records at the rank(s)
                  taxaWork <- taxaAboveFilter[taxaAboveFilter[[taxaRank[taxaRankCount]]] %in% uniqueTaxaAboveFilter, c(taxaRank[taxaRankCount], "query_coverage", "percent_identity", "evalue")]
                  taxaNumRecords <- as.data.frame(table(taxaWork[[taxaRank[taxaRankCount]]]))
                  colnames(taxaNumRecords) <- c(taxaRank[taxaRankCount],"Num_records")

                  #Get the min query coverage value for the taxa specific dataset
                  taxaQuery <- as.data.frame(tapply(as.vector(as.numeric(taxaWork$query_coverage)), taxaWork[[taxaRank[taxaRankCount]]], min))
                  #make the row name the first column
                  taxaQuery <- data.frame(uniqueID = row.names(taxaQuery), taxaQuery)
                  colnames(taxaQuery) <- c(taxaRank[taxaRankCount],"taxaQuery")
                  taxaTotal <- merge(taxaNumRecords, taxaQuery, by = taxaRank[taxaRankCount])

                  #Get the min percent identity for the target taxa specific dataset
                  taxaPerIdent <- as.data.frame(tapply(as.vector(as.numeric(taxaWork$percent_identity)), taxaWork[[taxaRank[taxaRankCount]]], min))
                  #make the row name the first column
                  taxaPerIdent <- data.frame(uniqueID = row.names(taxaPerIdent), taxaPerIdent)
                  colnames(taxaPerIdent) <- c(taxaRank[taxaRankCount],"taxaPerIdent")
                  taxaTotal <- merge(taxaTotal, taxaPerIdent, by = taxaRank[taxaRankCount])

                  #Get the max e value for the target taxa specific dataset
                  taxaEval <- as.data.frame(tapply(as.vector(as.numeric(taxaWork$evalue)),taxaWork[[taxaRank[taxaRankCount]]], max))
                  #make the row name the first column
                  taxaEval <- data.frame(uniqueID = row.names(taxaEval), taxaEval)
                  colnames(taxaEval) <- c(taxaRank[taxaRankCount],"taxaEval")
                  taxaTotal <- merge(taxaTotal, taxaEval, by = taxaRank[taxaRankCount])

                  #add the number of records, query coverage, percent identity, and e values on to the taxa specific dataset
                  taxaTotal <- as.data.frame(cbind(taxaTotal, paste0("(",round(taxaTotal$Num_records,2),",", round(taxaTotal$taxaQuery,2),",", round(taxaTotal$taxaPerIdent,2),",", format(taxaTotal$taxaEval, digits = 2),")")))

                  # Remove all white spaces from the column
                  taxaTotal[,ncol(taxaTotal)]<- gsub("\\s+", "", taxaTotal[,ncol(taxaTotal)])

                  #adding the values in brackets back onto the species in a dataframe for merging below
                  taxaTotal <- taxaTotal[,c(1,ncol(taxaTotal))]
                  condenseWork <- paste0(taxaTotal[,1],taxaTotal[,2])

                  if(length(na.omit(uniqueTaxaAboveFilter)) > 1){
                    #Place the results in the correct rank column when multiple results are present
                    condensedOutThres[1,taxaRank[taxaRankCount]] <- paste0(condenseWork, collapse = ", ")
                  }else if (length(na.omit(uniqueTaxaAboveFilter))== 1){
                    #Place the results in the correct rank column when a single result is present
                    condensedOutThres[1,c(taxaRank[taxaRankCount])] <- condenseWork
                    #If the ranks above the current loop rank do not contain NA or multiple results then assign the lowest values
                    valueCheck<-condensedOutThres[1,c(which(colnames(condensedOutThres) == "superkingdom"):which(colnames(condensedOutThres)== taxaRank[taxaRankCount])),drop=FALSE]
                    if(sum(is.na(valueCheck[1,]))==0){
                      if(sum(grep("),",valueCheck[1,]))==0){
                        condensedOut[1,"Lowest_Single_Rank_Above_Thres"] <- taxaRank[taxaRankCount]
                        condensedOut[1,"Lowest_Single_Taxa_Above_Thres"] <- paste0(condenseWork, collapse = ", ")
                      }
                    }
                  }

                }# end of the if uniqueTaxaAboveFilter >0

              } #end of the loop going through each of the ranks

              #add two columns with the lowest final identification as taken from total lowest common taxa
              #and lowest common taxa with records above the per coverage and per ident final identification
              if (is.na(condensedOut[1,"Lowest_Single_Rank"])){
                if (is.na(condensedOut[1,"Lowest_Single_Rank_Above_Thres"])){
                  condensedOut[1,"Final_Rank"]<-NA
                  condensedOut[1,"Final_Taxa"]<-NA
                  condensedOut[1,"Final_Common_Names"]<-NA
                  #Add a column with the coverage and identity thresholds of 0
                  condensedOut<-data.frame(append(condensedOut, c("Final_Rank_Taxa_Thres"=paste0("0,0")), after=which( colnames(condensedOut)=="Final_Taxa" )))
                }else{
                  condensedOut[1,"Final_Rank"]<-condensedOut[1,"Lowest_Single_Rank_Above_Thres"]
                  condensedOut[1,"Final_Taxa"]<-condensedOut[1,"Lowest_Single_Taxa_Above_Thres"]
                  condensedOut[1,"Final_Common_Names"]<-paste0(head(unique(na.omit(taxaAboveFilter[taxaAboveFilter$uniqueID==blastResultsuniqueID[records],]$c_name)),20), collapse = ", ")
                  #update all taxonomic results to reflect those of the above threshold final rank and taxa
                  condensedOut[1, c("superkingdom","phylum","class","order","family","genus","species")]<-condensedOutThres[1, c("superkingdom","phylum","class","order","family","genus","species")]
                  #Add a column with the coverage and identity thresholds
                  condensedOut<-data.frame(append(condensedOut, c("Final_Rank_Taxa_Thres"=paste0(coverage,",",ident)), after=which( colnames(condensedOut)=="Final_Taxa" )))
                }
              }else {
                if (is.na(condensedOut[1,"Lowest_Single_Rank_Above_Thres"])){
                  condensedOut[1,"Final_Rank"]<-condensedOut[1,"Lowest_Single_Rank"]
                  condensedOut[1,"Final_Taxa"]<-condensedOut[1,"Lowest_Single_Taxa"]
                  condensedOut[1,"Final_Common_Names"]<-paste0(head(unique(na.omit(taxa[taxa$uniqueID==blastResultsuniqueID[records],]$c_name)),20), collapse = ", ")
                  #Add a column with the coverage and identity thresholds of 0
                  condensedOut<-data.frame(append(condensedOut, c(Final_Rank_Taxa_Thres=paste0("0,0")), after=which( colnames(condensedOut)=="Final_Taxa" )))

                }else{#So both lowest single rank and lowest single rank above threshold have values
                  if (which(colnames(condensedOut)==condensedOut[1,"Lowest_Single_Rank"]) > which(colnames(condensedOut)==condensedOut[1,"Lowest_Single_Rank_Above_Thres"])){
                    condensedOut[1,"Final_Rank"]<-condensedOut[1,"Lowest_Single_Rank"]
                    condensedOut[1,"Final_Taxa"]<-condensedOut[1,"Lowest_Single_Taxa"]
                    condensedOut[1,"Final_Common_Names"]<-paste0(head(unique(na.omit(taxa[taxa$uniqueID==blastResultsuniqueID[records],]$c_name)),20), collapse = ", ")
                    #Add a column with the coverage and identity thresholds of 0
                    condensedOut<-data.frame(append(condensedOut, c("Final_Rank_Taxa_Thres"=paste0("0,0")), after=which( colnames(condensedOut)=="Final_Taxa" )))

                  }else{
                    condensedOut[1,"Final_Rank"]<-condensedOut[1,"Lowest_Single_Rank_Above_Thres"]
                    condensedOut[1,"Final_Taxa"]<-condensedOut[1,"Lowest_Single_Taxa_Above_Thres"]
                    condensedOut[1,"Final_Common_Names"]<-paste0(head(unique(na.omit(taxaAboveFilter[taxaAboveFilter$uniqueID==blastResultsuniqueID[records],]$c_name)),20), collapse = ", ")
                    #update all taxonomic results to reflect those of the above threshold final rank and taxa
                    condensedOut[1, c("superkingdom","phylum","class","order","family","genus","species")]<-condensedOutThres[1, c("superkingdom","phylum","class","order","family","genus","species")]
                    #Add a column with the coverage and identity thresholds
                    condensedOut<-data.frame(append(condensedOut, c("Final_Rank_Taxa_Thres"=paste0(coverage,",",ident)), after=which( colnames(condensedOut)=="Final_Taxa" )))
                  }
                }
              }


              #Get the final taxonomic rank and use the uniqueTaxaAboveFilter dataframe to evaluate the number of different values in the rank right below.
              #If one of the returned values has greater than 95% of the returned results update the final result to be this value
              # TBAT: Taxa Below Assigned Taxa
              if(!is.na(condensedOut[1,"Final_Rank"])){

                #if the rank is species then there isn't anything lower to do. Add dash to the proportional columns
                if(condensedOut[1,"Final_Rank"] != "species"){
                  #Get the results for the rank right below the final rank
                  propResults <-condensedOut[1,which(colnames(condensedOut) == condensedOut[1,"Final_Rank"])+1]

                  # Split the string into separate elements
                  propResults <- unlist(strsplit(unlist(propResults), ", "))

                  # Define a function to remove substring until the first occurrence of "(" and the "(" itself,
                  # and then remove all characters after the first comma ","
                  remove_prefix <- function(string) {
                    sub("^[^(]*\\(([^,]*).*", "\\1", string)
                  }

                  # Apply the function to each element of the list
                  propResults <- unlist(lapply(propResults, remove_prefix))
                  if(length(propResults)>1){
                    if (max(as.numeric(propResults)/sum(as.numeric(propResults))) > propThres){

                      #Add in code here to flag and replace or just to flag.
                      #Perhaps have a check value to indicate that for all records over the
                      #Thresholds as shown in the 'Final_Rabk_Taxa_Thres' column
                      #That have a TBAT flag then automatically place the assignment
                      #To the rank below at the most occurance taxa


                      condensedOut[1,"Result_Code"] <- paste0("TBAT(",propThres,")")
                    }else{condensedOut[1,"Result_Code"] <- "-"}
                  }else{condensedOut[1,"Result_Code"] <- "-"}
                }else{condensedOut[1,"Result_Code"] <- "-"}
              }else{condensedOut[1,"Result_Code"] <- "-"}

              #Where the user has defined the reporting threshold values coverReportThresh and identReportThresh
              if(!is.na(condensedOut[1,"Final_Taxa"])){
                #Subset the final dataframe to only have records with greater than the defined query coverage and identity
                #Split the column by ( and replace the close bracket )
                entries <- as.data.frame(do.call('rbind', strsplit(as.character(condensedOut[1,"Final_Taxa"]),'(',fixed = TRUE)))
                entries[,2] <- gsub(")","", as.character(entries[,2]))
                #Split the values column by ,
                entries <- cbind(entries[,1],data.frame(do.call('rbind', strsplit(as.character(entries[,2]),',', fixed = TRUE))), row.names = NULL)
                #Name the total results temp columns
                colnames(entries) <- c("Final_Taxa", "Num_Rec", "Coverage", "Identity", "Max_eVal")

                #Add a result code to all records that had coverage or identity below the submitted thresholds
                if(as.numeric(entries[1,3]) < coverReportThresh){
                  #BCRT: Final taxa below the nucleotide coverage reporting threshold
                  if(condensedOut[1,"Result_Code"] == "-"){
                     condensedOut[1,"Result_Code"] <- paste0("BCRT(",coverReportThresh,")")
                  }else{
                    condensedOut[1,"Result_Code"] <- paste0(condensedOut[1,"Result_Code"], ", BCRT(",coverReportThresh,")")
                  }
                }

                if(as.numeric(entries[1,4]) < identReportThresh){
                  #BIRT: Final taxa result is below the identity reporting threshold
                  if(condensedOut[1,"Result_Code"] == "-"){
                    condensedOut[1,"Result_Code"] <- paste0("BIRT(",identReportThresh,")")
                  }else{
                    condensedOut[1,"Result_Code"] <- paste0(condensedOut[1,"Result_Code"], ", BIRT(",identReportThresh,")")
                  }
                }

                if(nrow(condensedOut)!=0){
                  if(!is.na(condensedOut[1,"Lowest_Single_Taxa"])){
                    numHitsNF <- unique(taxa[,which( colnames(taxa)==condensedOut[1,"Final_Rank"] )])
                    if(length(numHitsNF) == 1){
                      #SANF: Saturated non-filtered (same taxa for all hits in non-filtered Lowest_Single_Taxa)
                      if(condensedOut[1,"Result_Code"] == "-"){
                        condensedOut[1,"Result_Code"] <- paste0("SANF(",coverage, ",",ident,")")
                      }else{
                        condensedOut[1,"Result_Code"] <- paste0(condensedOut[1,"Result_Code"], ", SANF(",coverage, ",",ident,")")
                      }
                    }
                  }
                }

                if(!is.na(condensedOut[1,"Lowest_Single_Taxa_Above_Thres"])){
                  numHitsNF <- unique(taxaAboveFilter[,which( colnames(taxaAboveFilter)==condensedOut[1,"Final_Rank"] )])
                  if(length(numHitsNF) == 1){
                    #SFAT: Saturated filtered taxa above threshold
                    if(condensedOut[1,"Result_Code"] == "-"){
                      condensedOut[1,"Result_Code"] <- paste0("SFAT(",coverage, ",",ident,")")
                    }else{
                      condensedOut[1,"Result_Code"] <- paste0(condensedOut[1,"Result_Code"], ", SFAT(",coverage, ",",ident,")")
                    }
                  }
                }

                #Where the final result has a blank in the Final_Common_Names column place an NA
                condensedOut$Final_Common_Names[condensedOut$Final_Common_Names==""] <- NA

            } #Closing the if nrow condensedOut !=0

              write.table(as.data.frame(condensedOut), file = paste0(fileLoc, "/", as.vector(filesList[fileCounter,1]), "_", filesList[fileCounter,10], "_taxaAssign_", dateStamp, ".tsv"), row.names = FALSE, col.names = FALSE, append = TRUE, sep = "\t", quote = FALSE)

              return(as.data.frame(condensedOut))

            }#End of if there is a returned result for the BLAST

          }#End of the taxaResults function

          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 9")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 9"), file = auditFile, append = TRUE))}
          if(numCores==1){
            if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 10")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 10"), file = auditFile, append = TRUE))}
            finalCondensedOut <- pbapply::pblapply(seq_len(length(blastResultsuniqueID)), taxaResults)
            if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 10A")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 10A"), file = auditFile, append = TRUE))}
          }else{
            if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 11")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 11"), file = auditFile, append = TRUE))}
            finalCondensedOut <- pbapply::pblapply(seq_len(length(blastResultsuniqueID)), taxaResults, cl = numCores)
          }

          #Take the list from the pblapply and place it in to a data frame
          finalCondensedOut <- as.data.frame(do.call(rbind,finalCondensedOut))

          #Print the results fo the condensed taxa assign to a file
          colnames(finalCondensedOut)<-gsub("\\.","-",colnames(finalCondensedOut))
          write.table(as.data.frame(finalCondensedOut, check.names=FALSE), file = paste0(fileLoc, "/", as.vector(filesList[fileCounter,1]), "_", filesList[fileCounter,10], "_taxaAssign_", dateStamp, ".tsv"), row.names = FALSE, col.names=TRUE, append = FALSE, sep="\t", quote = FALSE)

          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 12")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 12"), file = auditFile, append = TRUE))}

          #Check to see if the BLAST results has an accompanying dada output file
          if(file.exists(filesList[fileCounter,4])){

            if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 13")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 13"), file = auditFile, append = TRUE))}

            dadaOutputTable <- read.delim(filesList[fileCounter,4], header = TRUE, check.names=FALSE)

            #Read in the taxa condense output file.
            taxaTotalFile <- read.delim(file = paste0(fileLoc, "/", as.vector(filesList[fileCounter,1]), "_", filesList[fileCounter,10], "_taxaAssign_", dateStamp, ".tsv"), header = TRUE, check.names=FALSE)

            #Merge the taxa with the records per sample
            if(includeAllDada==TRUE){
              totalDataset <- merge(taxaTotalFile, dadaOutputTable, "uniqueID", all = TRUE)
            }else{
                totalDataset <- merge(taxaTotalFile, dadaOutputTable, "uniqueID", all.x = TRUE)
            }
            #Change the .(which were forced to change from the original naming convention of -) to - for output
            colnames(totalDataset)<-gsub("\\.","-",colnames(totalDataset))
            write.table(as.data.frame(totalDataset, check.names=FALSE), file = paste0(fileLoc, "/", as.vector(filesList[fileCounter,1]), "_", filesList[fileCounter,10], "_taxaAssign_", dateStamp, ".tsv"), row.names = FALSE, col.names=TRUE, append = FALSE, sep="\t", quote = FALSE)

            if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 14")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 14"), file = auditFile, append = TRUE))}

          }else if(file.exists(filesList[fileCounter,7])){ #if no dada asv file check for a fasta file

            if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 15")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 15"), file = auditFile, append = TRUE))}

            #Read in the data from the target file
            seqTable <- read.delim(filesList[fileCounter,7], header = FALSE)

            if(!grepl(">", seqTable[3,1], fixed = TRUE)){

              if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 16")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 16"), file = auditFile, append = TRUE))}

              print("*****ERROR**********************************************************************")
              print("The submitted fasta file is not in the single line nucleotide format which is ")
              print("needed for this script. Please correct the format and rerun this script.")
              print(" Sequence data was not added to this taxaAssign file. Please correct and try again!")
              print("********************************************************************************")

            }else{

              if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 17")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 17"), file = auditFile, append = TRUE))}

              #taking the read in file and changing from Fasta to tab delimited
              seqTableTemp <- data.frame(uniqueID = (seqTable[seq(from = 1, to = nrow(seqTable), by = 2), 1]))
              seqTableTemp["Sequence"] <- seqTable[seq(from = 2, to = nrow(seqTable), by = 2), 1]
              seqTable<-seqTableTemp
              #Remove the seqTableTemp to free memory
              rm(seqTableTemp)

              #Remove the > from the uniqueID
              seqTable$uniqueID <- gsub(">", "", seqTable$uniqueID)

              #Read in the taxa condense output file.
              taxaTotalFile <- read.delim(file = paste0(fileLoc, "/", as.vector(filesList[fileCounter,1]), "_", filesList[fileCounter,10], "_taxaAssign_", dateStamp, ".tsv"), header = TRUE, check.names=FALSE)

              #Merge the taxa with the records per sample
              if(includeAllDada==TRUE){
                totalDataset <- merge(taxaTotalFile, seqTable, "uniqueID", all = TRUE)
              }else{
                totalDataset <- merge(taxaTotalFile, seqTable, "uniqueID", all.x = TRUE)
              }
              #Change the .(which were forced to change from the original naming convention of -) to - for output
              colnames(totalDataset)<-gsub("\\.","-",colnames(totalDataset))
              write.table(as.data.frame(totalDataset, check.names=FALSE), file = paste0(fileLoc, "/", as.vector(filesList[fileCounter,1]), "_", filesList[fileCounter,10], "_taxaAssign_", dateStamp, ".tsv"), row.names = FALSE, col.names=TRUE, append = FALSE, sep="\t", quote = FALSE)

              if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 18")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 18"), file = auditFile, append = TRUE))}

            }

          }

          if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 19")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 19"), file = auditFile, append = TRUE))}
          print(paste0("End of this loop ", fileCounter, " of ", nrow(filesList), " at ", Sys.time()))

        }#end of the loop
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 20")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 20"), file = auditFile, append = TRUE))}
        print(paste0("Start at ", startTime, " end at ", Sys.time()))
      }else{ #End of the checking to see if there is a fas or asv associated file
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 21")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 21"), file = auditFile, append = TRUE))}
        print("ERROR - There is no associated fasta or ASV file with one or more of the BLAST files in the target folder. Please try again.")
      }
      if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 22")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 22"), file = auditFile, append = TRUE))}

    }else{# End of if checking to see if there is a fasta to go with the _BLAST file

      if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 23")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 23"), file = auditFile, append = TRUE))}
      print("ERROR - There is no associated .fas file please ensure both the _BLAST file and the .fas file are located in the selected location and try again.")

    }
  }else{#End of the if checking to see if there is more than one _BLAST in the file path
    if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 24")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 24"), file = auditFile, append = TRUE))}
    print("ERROR - remove the '_BLAST' from the naming convention in the file structure and try again.")

  }
}#end of the function
