# Written by Rob Young at the University of Guelph in Ontario Canada, Sept, 2023
# ******************************************************************************
#Roxygen2 Documentation:

#' @export
#'
#' @title Reduce Taxa Assignment
#'
#' @author Robert G. Young
#'
#' @description
#' This function takes a file selection and then uses all 'taxaAssign' and/or
#' 'taxaAssignCombine' files in that directory and reduces all ASV with the same
#' taxonomic assignment into a single result and places these results in a
#' 'taxaReduced' file for each of the target files in the directory.
#'
#' @details
#' This function requires a file in a directory where all 'taxaAssign'
#' and/or 'taxaAssignCombine' files in that directory will be combined. All records
#' with the same taxonomic result will be combined. The BLAST values in parentheses
#' ("Num_Rec", "Coverage", "Identity", "Max_eVal") are combine by the mean number of records,
#' the minimum coverage and identity, and the maximum eValue.
#'
#' @examples
#' \dontrun{
#' reduce_taxa()
#' reduce_taxa(fileLoc = NULL,   numCores = 1)
#' }
#'
#' @param fileLoc The location of a file in a directory where all of the 'taxaAssign'
#' and/or 'taxaAssignCombine' files are located.
#' @param numCores The number of cores used to run the function (default = 1,
#' Windows systems can only use a single core).
#'
#' @returns
#' This function produces a 'taxa_reduced' file for every 'taxaAssign' or
#' 'taxaAssignCombine' present in the target directory.
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
#' combine_reduced_output()
#'
################# reduce_taxa FUNCTION #########################################
reduce_taxa<- function(fileLoc = NULL,   numCores = 1){

  #Get the initial working directory
  start_wd <- getwd()
  on.exit(setwd(start_wd))

  #load in the files list
  if (is.null(fileLoc)){
    print(paste0("Select a file in the file folder with taxon assign files or the taxon assign combined files you would like to reduce (extension '_taxaAssign_YYYY_MM_DD_HHMM.tsv' or '_YYYY_MM_DD_HHMM_taxaAssignCombine.tsv' )."))
    fileLoc <- file.choose()
  }

  #Printing the start time
  print(paste0("Start time...", Sys.time()))
  startTime <- paste0("Start time...", Sys.time())
  dateStamp <- paste0(format(Sys.time(), "%Y_%m_%d_%H%M"))

  #Get the directory for the location
  fileLoc <- dirname(fileLoc)
  #Set the working directory
  setwd(fileLoc)

  #Get the all files in the selected file folder with full paths of the files
  files <- as.data.frame(list.files(path = fileLoc, pattern = "*[.]*", full.names = TRUE), check.names=FALSE)
  # Get the local paths
  files[,2] <- list.files(path = fileLoc, pattern = "*[.]*")
  #Get all files with the '_taxaAssign' string
  files <- files[grepl("_taxaAssign.*", files[,2]),]
  #Remove files with _taxaAssignCombine.txt
  files <- files[!grepl("_taxaAssignCombined.txt", files[,2]),]
  # Get the names of the files
  files[,3] <- gsub("_taxaAssign.*","",files[,2])

  for(records in 1:nrow(files)){

    flag = 0

    print(paste0("Started file ", files[records,2], " at ",  format(Sys.time(), "%Y_%m_%d_%H%M")))
    print(paste0("Reading in file...", files[records,2], " at ",  format(Sys.time(), "%Y_%m_%d_%H%M")))

    #Read in the files for this loop
    totalResults <- read.delim(files[records,1], header = TRUE, check.names=FALSE)

    print(paste0("Finished reading in file ", records, " of ", nrow(files), " - ", files[records,2], " at ",  format(Sys.time(), "%Y_%m_%d_%H%M")))

    #Get the column indices
    colIndex <- which(names(totalResults) %in% c("uniqueID", "Lowest_Single_Rank", "Lowest_Single_Taxa", "Lowest_Single_Rank_Above_Thres", "Lowest_Single_Taxa_Above_Thres"))
    #Drop unnecessary columns
    totalResults <- totalResults[ , -c(colIndex)]

    #add in a column for the taxa without the values
    totalResults <- cbind(taxa = sub("\\(.*", "", totalResults$Final_Taxa), totalResults, row.names = NULL)

    #Remove records with NA or two results in the kingdom and/or NA in the taxa or sequence
    #two results in kingdom
    totalResults <- totalResults[!grepl("),", totalResults$superkingdom),, drop = FALSE]
    #NA in kingdom
    totalResults <- totalResults[!is.na(totalResults$superkingdom),, drop = FALSE ]
    #NA in taxa
    totalResults <- totalResults[!is.na(totalResults$taxa), , drop = FALSE]
    #NA in sequence
    totalResults <- totalResults[!is.na(totalResults$Sequence),, drop = FALSE ]

    #Get a list of unique taxa
    uniqueTaxa <- unique(totalResults$taxa)

    taxaResults <- function(numUniqueTaxa){
    #loop through the unique taxa
#    for(numUniqueTaxa in 1:length(uniqueTaxa)){

      #subset totalResults for the target unique taxa
      uniqueTotalResults <- totalResults[totalResults$taxa == uniqueTaxa[numUniqueTaxa],]

      #Remove NA if present in the taxa column (not sure why I am picking these up?)
      uniqueTotalResults <- uniqueTotalResults[!is.na(uniqueTotalResults$taxa),]

      #Remove the Final_Rank_Taxa_Thres column for this set
      colIndex <- which(names(uniqueTotalResults) %in% c("Final_Rank_Taxa_Thres"))
      uniqueTotalResults <- uniqueTotalResults[ , -c(colIndex)]

      #kingdom
      uniqueTotalResultsTemp<-data.frame(superkingdom = combine_row(uniqueTotalResults$superkingdom), check.names=FALSE)
      #phylum
      uniqueTotalResultsTemp["phylum"] <- combine_row(uniqueTotalResults$phylum)
      #class
      uniqueTotalResultsTemp["class"] <- combine_row(uniqueTotalResults$class)
      #Order
      uniqueTotalResultsTemp["order"] <- combine_row(uniqueTotalResults$order)
      #Family
      uniqueTotalResultsTemp["family"] <- combine_row(uniqueTotalResults$family)
      #Genus
      uniqueTotalResultsTemp["genus"] <- combine_row(uniqueTotalResults$genus)
      #species
      uniqueTotalResultsTemp["species"] <- combine_row(uniqueTotalResults$species)
      #top BLAST hit
      if(length(unique(uniqueTotalResults$Top_BLAST)) < 21){
        uniqueTotalResultsTemp["Top_BLAST"] <- paste0(unique(uniqueTotalResults$Top_BLAST), collapse = ", ")
      }else{
        uniqueTotalResultsTemp["Top_BLAST"] <- "Greater than 20"
      }

      #Final_Common_Names
      if(length(unique(uniqueTotalResults$Final_Common_Names)) == 1){
        uniqueTotalResultsTemp["Final_Common_Names"] <-uniqueTotalResults[1,"Final_Common_Names"]
      }else if(length(unique(uniqueTotalResults$Final_Common_Names)) < 21){
        uniqueTotalResultsTemp["Final_Common_Names"] <- paste0(unique(uniqueTotalResults$Final_Common_Names), collapse = ", ")
      } else{
        uniqueTotalResultsTemp["Final_Common_Names"] <- "Greater than 20"
      }

      #Final_Rank
      if(unique(uniqueTotalResults$Final_Rank)==1){
      uniqueTotalResultsTemp["Final_Rank"] <- unique(uniqueTotalResults$Final_Rank)
      }else{
        uniqueTotalResultsTemp["Final_Rank"] <- paste0(unique(uniqueTotalResults$Final_Rank), collapse = ", ")
      }

      #Final_Taxa
      uniqueTotalResultsTemp["Final_Taxa"] <- combine_row(uniqueTotalResults$Final_Taxa)

      #Get all results in the Result_Code column and place into a vector
      result_code_list <- unique(unlist(strsplit(uniqueTotalResults$Result_Code, split = ", ")))
#        result_code_list <- gsub(" ","",result_code_list)
#        result_code_list <- sort(unique(unlist(strsplit(result_code_list, split = ","))))
#        result_code_list<- paste(result_code_list, collapse=",")
#        result_code_list <- gsub("-,","",result_code_list)
#        uniqueTotalResultsTemp["Result_Code"] <- result_code_list

      #order the Result_Code vector so that all of the output are in the same format
      result_code_list <- sort(result_code_list, decreasing = TRUE)

      #Get all results in the Result_Code column and place into a vector
      result_code_list <- unique(unlist(strsplit(result_code_list, split = ", ")))
      # Use grep to identify elements containing a dash
      elements_with_dash <- grep("-", result_code_list, value = TRUE)
      # Filter out elements with dashes
      result_code_list <- result_code_list[!result_code_list %in% elements_with_dash]
      # Combine the elements into a comma-separated single variable
      result_code_list <- paste(result_code_list, collapse = ", ")

      #Add - back into the empty rows.
      if(result_code_list==""){
        uniqueTotalResultsTemp["Result_Code"] <- "-"
      }else{
        uniqueTotalResultsTemp["Result_Code"] <- result_code_list
      }

      #Use the dereplicate_sequences function to get the longest dereplicated sequence
      reducedSeq<-dereplicate_sequences(uniqueTotalResults$Sequence[!is.na(uniqueTotalResults$Sequence)])
      uniqueTotalResultsTemp["RepSequence"]<-reducedSeq[1]
      #Get the number of unique sequences for the taxa
      uniqueTotalResultsTemp["Number_ASV"]<-length(unique(uniqueTotalResults$Sequence[!is.na(uniqueTotalResults$Sequence)]))
      #Get the average length of the sequence
      uniqueTotalResultsTemp["Average_ASV_Length"]<- mean(nchar(uniqueTotalResults$Sequence))
      # Get the number of occurrences for the taxa
      uniqueTotalResultsTemp["Number_Occurrences"]<-sum(colSums(uniqueTotalResults[,c(as.numeric(which(colnames(uniqueTotalResults) == "Results")+1):ncol(uniqueTotalResults))]) != 0)
      # Get the average number of reads
      per_sample_sums<-colSums(uniqueTotalResults[,c(as.numeric(which(colnames(uniqueTotalResults) == "Results")+1):ncol(uniqueTotalResults))])
      uniqueTotalResultsTemp["Average_ASV_Per_Sample"]<-mean(per_sample_sums[per_sample_sums!=0])
      #Get the median number of reads
      uniqueTotalResultsTemp["Median_ASV_Per_Sample"]<-stats::median(per_sample_sums[per_sample_sums!=0])

      #Results
      if(length(unique(uniqueTotalResults$Results))==1){
        uniqueTotalResultsTemp["Results"] <- unique(uniqueTotalResults$Results)
      }else{
        uniqueTotalResultsTemp["Results"] <- paste0(unique(uniqueTotalResults$Results), collapse = ", ")
      }

      #Add on the sequence reads per sample on to the end of the reporting table
      sampleSums<-as.data.frame(colSums(uniqueTotalResults[,c(as.numeric(which(colnames(uniqueTotalResults) == "Results")+1):ncol(uniqueTotalResults))]), check.names=FALSE)
      sampleSumNames<-rownames(sampleSums)
      sampleSums<-t(sampleSums)
      colnames(sampleSums)<-sampleSumNames
      uniqueTotalResultsTemp<-cbind(uniqueTotalResultsTemp, sampleSums)

      return(uniqueTotalResultsTemp)

    } #Closing off the loop through unique taxa

    if(numCores==1){
      finalResults <- pbapply::pblapply(seq_len(length(uniqueTaxa)), taxaResults)
    }else{
      finalResults <- pbapply::pblapply(seq_len(length(uniqueTaxa)), taxaResults, cl = numCores)
    }
    finalResults <- do.call(rbind, finalResults)

    #Change the .(which were forced to change from the original naming convention of -) to - for output
    colnames(finalResults)<-gsub("\\.","-",colnames(finalResults))
    #print the output variable to file for the first loop including the headers
    write.table(as.data.frame(finalResults, check.names=FALSE), file = paste0(files[records,3], "_taxaReduced_",dateStamp,".tsv"), row.names = FALSE, col.names = TRUE, append = FALSE, sep = "\t", quote = FALSE)

  }#Closing off looping through each file in the target directory

  print(paste0(startTime, " end at ", Sys.time()))

}#Closing off the function reduce_taxa


################################ ACCEPT A NUCLEOTIDE STRING AND REVERSE COMPLIMENT AND RETURN A NUCLEOTIDE STRING ###########

rev_comp <- function(sequence){

  #breaking the sequence string into columns
  seq <- unlist(strsplit(toupper(sequence), NULL))

  #Reversing the order of the columns
  seq <- rev(seq)

  #The paste collapses the columns to a single string before returning and the ifs swap the complimentary nucleotides
  paste(unlist(lapply(seq, function(swap){
    if(swap=="A") compNucleotide <- "T"
    if(swap=="C") compNucleotide <- "G"
    if(swap=="G") compNucleotide <- "C"
    if(swap=="T") compNucleotide <- "A"
    if(swap=="U") compNucleotide <- "A"
    if(swap=="Y") compNucleotide <- "R"
    if(swap=="R") compNucleotide <- "Y"
    if(swap=="K") compNucleotide <- "M"
    if(swap=="M") compNucleotide <- "K"
    if(swap=="B") compNucleotide <- "V"
    if(swap=="D") compNucleotide <- "H"
    if(swap=="H") compNucleotide <- "D"
    if(swap=="V") compNucleotide <- "B"
    if(swap=="N") compNucleotide <- "N"
    return(compNucleotide)
  })), collapse = "")
}

##################################### Dereplicate sequences FUNCTION ##############################################################
dereplicate_sequences <- function(sequences){

  #First do a simple dereplication (this should do nothing as the table already represents ASVs)
  sequences <- unique(sequences)

  if(length(sequences) > 1){

    #generate a column with the number of characters
    sequences <- cbind(sequences, nchar(sequences))

    #Order the remaining sequences
    sequences <- sequences[order(as.numeric(sequences[,2]), decreasing = FALSE),]

    #Starting at the first sequences loop through all remaining sequences and test against the first.
    #set the working variable to the sequences variable.
    workSeq <- sequences

    for (numOfSeq in 1:nrow(sequences)){
      #Check if the counter element sequence in the sequences dataframe is present in the workSeq
      #if yes then remove the sequence from the workSeq dataframe. The if is greater than 1 (so 2 or more) because
      #it would find itself and another sequence where it was a sub sequence
      if(length(grep(sequences[numOfSeq,1], workSeq[,1], fixed = TRUE)) > 1){
        #remove the substring if present in the other strings
        workSeq <- as.data.frame(workSeq[-numOfSeq,, drop = FALSE], check.names=FALSE)
      }#End of the if seeing if the current sequence is a substring in another sequence

      #Check if the rev_comp counter element sequence in the sequences dataframe is present in the workSeq
      #if yes then remove the sequence from the workSeq dataframe. The if is only greater than 0 (so 1 or more) because when checking
      #for the reverse compliment sequence it will not find itself.
      if(length(grep(rev_comp(sequences[numOfSeq,1]), workSeq[,1], fixed = TRUE))>0){
        #remove the substring if present in the other strings
        workSeq <- as.data.frame(workSeq[-numOfSeq,, drop = FALSE], check.names=FALSE)
      }#End of the if seeing if the current sequence is a substring in another sequence
    }#end of the for loop going through each sequence element in the sequences dataframe
    sequences <- workSeq[,1]
  }#End of if checking to see if there is more than one sequence in the initial submission
  return(sequences)
}#end of dereplicate_sequences function

############################################ combine_row FUNCTION #############################
combine_row <- function(entries){

  #Remove any that are NA
  entries <- entries[!is.na(entries)]

  if(length(entries) == 0){
    entries="NA"
  }else {

    if(length(grep("),", entries)) == 0){
      #Split the column by ( and replace the close bracket )
      entries <- as.data.frame(do.call('rbind', strsplit(as.character(entries),'(',fixed = TRUE)), check.names=FALSE)
      entries[,2] <- gsub(")","", as.character(entries[,2]))

      #Split the values column by ,
      entries <- cbind(entries[,1],data.frame(do.call('rbind', strsplit(as.character(entries[,2]),',', fixed = TRUE)), check.names=FALSE), row.names = NULL)

      #Name the total results temp columns
      colnames(entries) <- c("Final_Taxa", "Num_Rec", "Coverage", "Identity", "Max_eVal")

      if(length(unique(entries[,1])) == 1){
        if(length(entries[,1]) == 0){
          entries="NA"
        }else if(length(entries[,1]) == 1){
          #Create the final value to return
          entries<-paste0(unique(entries$Final_Taxa), "(", entries$Num_Rec, ",", entries$Coverage, ",", entries$Identity, ",", entries$Max_eVal,")")
        }else {
          entries<-paste0(unique(entries$Final_Taxa), "(", round(mean(as.numeric(entries$Num_Rec)),1), ",", round(min(as.numeric(entries$Coverage)),1), ",", round(min(as.numeric(entries$Identity)),1), ",", format(max(as.numeric(entries$Max_eVal)), digits = 2, scientific = TRUE),")")
        }#Closing off the check if there are more than one entry after building data frame
      }else{
        entries<-NA
      }

    }else{ # else to the if where there are more than one taxa per row for one of the multiple rows
      entries <- NA
    }
  }#Closing off the if the length of the total dataset is 0 meaning after removal of NA values there were no results.
  return(entries)
}
