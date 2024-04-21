# Dada-BLAST-Taxon Assign-Condense Shiny Application (DBTCShiny)
DBTCShiny is an R shiny implementation of the DBTC Metabarcode analysis pipeline.

# Description

This repository contains the DBTCShiny package located at [rgyoung6/DBTCShiny](https://github.com/rgyoung6/DBTCShiny) . The Dada-BLAST-Taxon Assign-Condense Shiny package contains the foundational [DBTC](https://github.com/rgyoung6/DBTC) functions (able to be run through the command line) which have been wrapped in a Shiny application for easy user interface. The DBTC functions have four main outcomes...

  - [Fastq](https://en.wikipedia.org/wiki/FASTQ_format) file processing using Dada in R
  - Using the Basic Local Alignment Search Tool ([BLAST](https://en.wikipedia.org/wiki/BLAST_(biotechnology))) amplicon sequence variants ([ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant)) can be searched against local NCBI or custom libraries
  - Assign taxa to the unique reads using NCBI taxon database through taxa names and/or taxaID's
  - Condense the resulting [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) taxonomic assignment tables to unique taxa with the ability to combine datasets (using different sequence libraries for the same reads, or results from the same samples for different molecular regions) into a combined results table

Before using these files and working through this tutorial please install the [DBTCShiny](https://github.com/rgyoung6/DBTCShiny) package.

**NOTE:** While the DBTCShiny package has been built for the analysis of high-thoughput sequencing results, the BLAST and taxonomic assignment, taxonomic condense, and mapping functions can be utilized with single specimen Sanger sequencing data.

# Table of Contents  
- [Installation](#installation)
- [Package Dependencies](#package-dependencies)
  * [External R Element Dependencies](#external-r-element-dependencies)
    * [NCBI BLAST+ local program to run BLAST on local databases](#ncbi-blast+-local-program-to-run-blast-on-local-databases)
    * [R package taxonomizr to establish a local NCBI taxonomy database](#r-package-taxonomizr-to-establish-a-local-ncbi-taxonomy-database)
    * [Establish a local NCBI prepared sequence database](#establish-a-local-ncbi-prepared-sequence-database)
    * [Create a custom sequence database to BLAST against](#create-a-custom-sequence-database-to-blast-against)
    * [Create a local NCBI taxonomy database to assign taxonomic identifications to BLAST results](#create-a-local-ncbi-taxonomy-database-to-assign-taxonomic-identifications-to-blast-results)
  * [R Packages Dependencies](#r-packages-dependencies)
     * [Bioconductor ShortRead and Dada2 packages](#bioconductor-shortread-and-dada2-packages)
     * [CRAN](#cran)
- [Function Descriptions](#function-descriptions)
- [Naming Convention Rules](#naming-convention-rules)
- [Citation](#citation)
- [Package Function Details](#package-function-details)
  * [Dada Implement](#dada-implement)
  * [Combine Dada Output](#combine-dada-output)
  * [Make BLAST DB](#make-blast-db)
  * [Sequence BLAST](#sequence-blast)
  * [Taxon Assignment](#taxon-assignment)
  * [Combine Assignment Output](#combine-assignment-output)
  * [Reduce Taxa](#reduce-taxa)
  * [Combine Reduced Output](#combine-reduced-output)
- [Import GPS and grouping data](#import-gps-and-grouping-data)
- [Mapping Dashboard](#mapping-dashboard)
  * [Mapping](#mapping)
  * [Data Filtering](#data-filtering)
  * [Mapped Data Table](#mapped-data-table) 

# Installation 

DBTCShiny can be installed three ways.

## 1. Install from CRAN

Coming Soon.... install.packages('DBTCShiny’)

## 2. Install via GitHub
Run the following commands in your R terminal...<br/>
```
if(!require(devtools)) install.packages('devtools')
library('devtools')
devtools::install_github('rgyoung6/DBTCShiny')
library('DBTCShiny')
```
**Note:** the first command to install the "devtools" may not be necessary if already installed.<br/>

## 3. Install through download from GitHub
Navigate to the [DBTCShiny](https://github.com/rgyoung6/DBTCShiny) GitHub page. Download the files associated with this page to your local computer and place them somewhere in the main file folder named DBTCShiny. Then run the following command pointing to that location on your local computer by replacing the HERE with the path in the below command...<br/>
```
library("DBTCShiny", lib.loc="HERE")
```

([Back to Top](#table-of-contents))
***
</details>

# Package Dependencies
 
## External R Element Dependencies

### NCBI BLAST+ local program to run BLAST on local databases

Follow the instructions on the NCBI [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#blast-executables) executables page to obtain a local version of the BLAST tools. The list of the latest installation files can be found [here](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

Note: It is best to download and install the most recent versions of the blast+ suit to your computer and place the programs in your computers path so you can access the program from any folder. However, the program files for both blastn and makeblastdb have been included in the DBTCShinyTutorial GitHub page for ease of use (please note that these may not be the most recent versions).

### R package taxonomizr to establish a local NCBI taxonomy database
The R package taxonomizr is used to establish a NCBI taxaID database (NOTE: this package is also required when using the taxon assignment elements in the DBTC pipeline).
```
install.packages('taxonomizr')
library('taxonomizr')
```

### Establish a local NCBI prepared sequence database
NCBI preformatted databases can be established through two methods.

1. Download your desired preformatted NCBI database by using the 'update_blastdb.pl' (found in the NCBI BLAST+ local install folder). NOTE: Perl programming langugage needs to be installed on your local machine. Instructions can be found at [Get NCBI BLAST databases](https://www.ncbi.nlm.nih.gov/books/NBK569850/).

2. You can download your desired preformatted NCBI database manually with instructions at [BLAST FTP Site](https://www.ncbi.nlm.nih.gov/books/NBK62345/#blast_ftp_site.The_blastdb_subdirectory) and a list of the available databases at [Index of /blast/db](https://ftp.ncbi.nlm.nih.gov/blast/db/). 

### Create a custom sequence database to BLAST against
In addition to the NCBI resources, DBTC can also use custom databases. To establish these databases you will requre a [Fasta](https://en.wikipedia.org/wiki/FASTA_format) file with the desired records with MACER formatted headers. The MACER R package and instructions can be found at either of the two locations:

[MACER CRAN](https://cran.r-project.org/web/packages/MACER/index.html)

[MACER GitHub](https://github.com/rgyoung6/MACER) (will have the most recent version and development versions)

### Create a local NCBI taxonomy database to assign taxonomic identifications to BLAST results
In the 'Preparation' section of the [taxonomizr website](https://cran.r-project.org/web/packages/taxonomizr/vignettes/usage.html), use the instructions and the prepareDatabase('accessionTaxa.sql', getAccessions = FALSE) taxonomizr command to establish a local taxonomy database.
```
prepareDatabase('accessionTaxa.sql', getAccessions = FALSE)
```
Note: Along with the accessionTaxa.sql two other files nodes.dmp and names.dmp files are downloaded. These two files are not necessary for downstream analyses and can be deleted.

([Back to Top](#table-of-contents))
## R Packages Dependencies

### Bioconductor ShortRead and Dada2 packages

The [ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html) package is required to run elements of the DBTC pipeline and can be obtained through Bioconductor.

The [dada2](https://www.bioconductor.org/packages/release/bioc/html/dada2.html) package is main package to process the raw fastq fules and can be obtained from Bioconductor. There is also a [dada2](https://benjjneb.github.io/dada2/) GitHub resource. 

```
if (!require('BiocManager', quietly = TRUE))
install.packages('BiocManager')
BiocManager::install('ShortRead')
BiocManager::install('dada2')
library(dada2)
```

### CRAN
Each of below CRAN packages and their dependencies are required for the DBTCShiny package.
```
install.packages(c('ggplot2',
                   'leaflet',
                   'leaflet.extras',
                   'magrittr',
                   'pbapply',
                   'plyr',
                   'shiny',
                   'shinycssloaders',
                   'shinydashboard',
                   'shinyWidgets',
                   'taxonomizr'))
library(c('ggplot2',
                   'leaflet',
                   'leaflet.extras',
                   'magrittr',
                   'pbapply',
                   'plyr',
                   'shiny',
                   'shinycssloaders',
                   'shinydashboard',
                   'shinyWidgets',
                   'taxonomizr'))
```

([Back to Top](#table-of-contents))
***
</details>

# Run DBTCShiny

After installation of the DBTCShiny and all of its dependencies you need to load the packge and then run the Shiny Graphical User Interface (GUI) using the following commands...

```
library(DBTCShiny)
launchDBTCShiny()
```

***
# Function Descriptions

## dada_implement()
This function requires a main directory containing a folder(s) representing sequencing runs which in-turn contains fastq files (the location of one of the fastq files in one of the sequencing run folders is used as an input argument). A run is a group of results processed at the same time on the same machine representing the same molecular methods. All sequencing folders in the main directory need to represent data from sequencing runs that have used the same primers and protocols. Output from this function includes all processing files and final main output files in the form of [Fasta](https://en.wikipedia.org/wiki/FASTA_format) files and amplicon sequencing variant ([ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant)) tables. 
    
## combine_dada_output()
This function uses DBTC Dada [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) output files (YYYY_MM_DD_HH_MM_UserInputRunName_Merge, YYYY_MM_DD_HH_MM_UserInputRunName_MergeFwdRev, and/or YYYY_MM_DD_HH_MM_UserInputRunName_TotalTable) and combines them into a single [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table with accompanying [Fasta](https://en.wikipedia.org/wiki/FASTA_format) file. This function also produces a file containing the processing information for the function. The main input argument for this function is the location of a file in a folder containing all [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) tables wanting to be combined. Output files are generated with the naming convention YYYY_MM_DD_HH_MM_combinedDada.

## make_BLAST_DB()
This function takes a [Fasta](https://en.wikipedia.org/wiki/FASTA_format) file with headers in the [MACER](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8443542/) format and establishes a database upon which a BLAST search can be completed. There are also NCBI preformatted databases available where the make_BLAST_DB() function can then be skipped ([NCBI BLAST databases](https://www.ncbi.nlm.nih.gov/books/NBK62345/#blast_ftp_site.The_blastdb_subdirectory)). The outcome of the function is a folder with a BLASTable NCBI formatted sequence database. Note: The there are three main required elements for MACER formatted records which include a the Unique Identifier 
- The MACER [Fasta](https://en.wikipedia.org/wiki/FASTA_format) header format - ```>UniqueID|OtherInformation|Genus|species|OtherInformation|Marker```


## seq_BLAST()
This function takes [Fasta](https://en.wikipedia.org/wiki/FASTA_format) files as input along with a user selected NCBI formatted database to BLAST sequences against. The outcome of the function are two files, a BLAST run file and a single file containing all of the BLAST results in tab delimited format. There are no headers in the BLAST results file but the columns are: query sequence ID, search sequence ID, search taxonomic ID, query to sequence coverage, percent identity, search scientific name, search common name, query start, query end, search start, search end, e-value.

## taxon_assign()
This function takes a BLAST result file and associated [Fasta](https://en.wikipedia.org/wiki/FASTA_format) files (either on their own or with accompanying [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) files generated from the dada_implement function) and collapses the multiple BLAST results into as single result for each query sequence. When an [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table is present the taxonomic results will be combined with the [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table.

## combine_assign_output()
This function takes a file selection and then uses all 'taxaAssign' files in that directory and combines them into a single output 'taxaAssignCombined.tsv' file.

## reduce_taxa()
This function takes a file selection and then uses all 'taxaAssign' and/or 'taxaAssignCombine' files in that directory and reduces all [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) with the same taxonomic assignment into a single result and places these results in a 'taxaReduced' file for each of the target files in the directory.

## combine_reduced_output()
This function takes a file selection and then uses all 'taxaReduced' files in that directory and combines them into a single 'CombineTaxaReduced' taxa table file with presence absence results.

([Back to Top](#table-of-contents))
***

# Naming convention rules

WARNING - NO WHITESPACE!

When running DBTCShiny functions the paths for the files selected cannot have whitespace! File folder locations should be as short as possible (close to the root directory) as some functions do not process long naming conventions. 

Also, special characters should be avoided (including question mark, number sign, exclamation mark). It is reccommended that dashes be used for separations in naming conventions while retaining underscores for use as information deliminters (this is how DBTC functions use underscore). 

There are several key character strings used in the DBTC pipeline, the presence of these strings in file or folder names will cause errors when running DBTC functions. 
The following strings are those used in DBTC and should not be used in file or folder naming:
  - _BLAST
  - _taxaAssign
  - _taxaCombined
  - _taxaReduced

([Back to Top](#table-of-contents))
***

# Citation
Young RG, et al., Hanner RH (2023) A Scalable, Open Source, Cross Platform, MetaBarcode Analysis Method using Dada2 and BLAST. Biodiversity Data Journal (In progress)

([Back to Top](#table-of-contents))

# Package Function Details

## Dada Implement
dada_implement() - Process metabarcode raw fastq files by run using Dada2 (Note: molecular markers are independently analysed and combined at the end of the analysis pipeline using the reduce_taxa() function).

### Input 
Two file types are required as input for the dada_implement() function. The first are the fastq files in the appropriate folder structure (see below) and the second is a file containing the primers used for the amplification of the sequence reads (tab separated file). To select all of the fastq files simply select one file in one of the Run directories to point to the desired data files, as long as the file and folder structure is correct.

**Fastq File Folder Structure**

```
            Parent Directory
                  |
                  |
          -----------------
          |               |
          |               |
    Run1 Directory     Run2 Directory
    -Fastq             -Fastq
    -Fastq             -Fastq
    ...                ...
```

**Format of the primer file**

| Forward        | Reverse           | 
| :------------- |:-------------| 
| AGTGTGTAGTGATTG      | CGCATCGCTCAGACTGACTGC | 
| GAGCCCTCGATCGCT      | GGTCGATAGCTACGCGCGCATACGACT      |  
|  | GGTTCACATCGCATTCAT      |   


### Arguments
- <strong>runFolderLoc -</strong> Select a directory that contains the run folders with the fastq files.
- <strong>primerFile -</strong> Select a file with the primers for this analysis.
- <strong>fwdIdent </strong> Foward identifier naming substring.
- <strong>revIdent </strong> Reverse identifier naming substring.
- <strong>nonMergeProcessing </strong> Non-Merge Processing (Default TRUE).
- <strong>maxPrimeMis</strong> Maximum number of mismatches allowed when pattern matching trimming the primers from the ends of the reads for the ShortRead trimLRPatterns() function (Default maxPrimeMis = 2).
- <strong>fwdTrimLen</strong> Select a forward trim length for the Dada filterAndTrim() function (Default fwdTrimLen = 0).
- <strong>revTrimLen</strong> Select a reverse trim length for the Dada filterAndTrim() function (Default revTrimLen = 0).
- <strong>maxEEVal</strong> Maximum number of expected errors allowed in a read for the Dada filterAndTrim() function (Default maxEEVal = 2).
- <strong>truncQValue</strong> Truncation value use to trim ends of reads, nucleotide with quality values less than this value will be used to trim the remainder of the read for the Dada filterAndTrim() function (Default truncQValue = 2).
- <strong>truncLenValueF</strong> Dada forward length trim value for the Dada filterAndTrim() function. This function is set to 0 when the pattern matching trim function is enabled (Default truncLenValueF = 0).
- <strong>truncLenValueR</strong> Dada reverse length trim value for the Dada filterAndTrim() function. This function is set to 0 when the pattern matching trim function is enabled (Default truncLenValueR = 0).
- <strong>error</strong> Percent of fastq files used to assess error rates for the Dada learnErrors() function (Default error = 0.1).
- <strong>nbases</strong> The total number of bases used to assess errors for the Dada learnErrors() function (Default nbases = 1e80) NOTE: this value is set very high to get all nucleotides in the error persent file subset. If the error is to be assessed using total reads and not specific fastq files then set the error to 1 and set this value to the desired number of reads.
- <strong>maxMismatchValue</strong> Maximum number of mismatches allowed when merging two reads for the Dada mergePairs() function (Default maxMismatchValue = 2).
- <strong>minOverlapValue</strong> Minimum number of overlapping nucleotides for the forward and reverse reads for the Dada mergePairs() function (Default minOverlapValue = 12).
- <strong>trimOverhang</strong> Trim merged reads past the start of the complimentary primer regions for the Dada mergePairs() function (Default trimOverhang = FALSE).
- <strong>minFinalSeqLen</strong> The minimum final desired length of the read (Default minFinalSeqLen = 100).

### Output

The output files from this function appear in four folders. See the below diagram for the saved file structure.
```            
                              Parent Directory
                                    |
                                    |
                      -----------------------------
                      |                           |
                      |                           |
                Run1 Directory                 Run2 Directory
                - Fastq                        - Fastq
                - Fastq                        - Fastq
                ...                            ...
                      |
                      |
            -------------------------------------------------------------------------------------------------------   
            |                                 |                                 |                                 |  
            |                                 |                                 |                                 |   
2023_10_16_1511_Run1_A_Qual       2023_10_16_1511_Run1_B_Filt    2023_10_16_1511_Run1_C_FiltQual    2023_10_16_1511_Run1_D_Output
  -fwdQual.pdf                      -fwdFilt.fastq                 -filtFwdQual.pdf                   -dadaSummary.txt
  -revQual.pdf                      -revFilt.fastq                 -filtRevQual.pdf                   -dadaSummaryTable.tsv
  ...                               ...                            ...                                -ErrorForward.pdf
                                    Primer_Trim                                                       -ErrorReverse.pdf
                                      -primeTrim.fastq                                                -MergeFwdRev.tsv
                                      ...                                                             -MergeFwdRev.fas
                                                                                                      -Merge.tsv
                                                                                                      -Merge.fas
                                                                                                      -TotalTable.tsv
```

### Intrepretation
Quality pdf's in the A_Qual folder represent the quality metrics for the raw fastq files
Files in the B_Filt folder represent the trimming (in the Primer_Trim folder) and trimmed and cleaned in the larger folder.
Quality pdf's in the C_FiltQual folder represent the quality metrics for the trimmed and cleaned fastq files.
There are numerous files in the D_Output folder. These include:
- dadaSummary.txt file which provides all of the information on the running of the dada_implement() function.
- dadaSummaryTable.tsv contains a table with summary information for the processing of the samples in the run.
- ErrorForward.pdf and ErrorReverse.pdf provide visualizations on the assessed sequencing error for the sequencing run.
- MergeFwdRev.tsv is the [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table with the data from the sequencing run and the MergeFwdRev.fas is a [Fasta](https://en.wikipedia.org/wiki/FASTA_format) file with the reads from the samples. The MergeFwdRev files include reads that were able to be merged, as well as reads that were not able to be merged.
  NOTE: The merged, forward, and reverse reads are obtained in parallel analyses and combined into a single file so MergeFwdRev files will represent triplicate molecular processing results. These files are present to see if there are reads that are not represented or poorly represented across merged and unidirectional results, perhaps indicating issues with one of the primers.
- Merge.tsv is the [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table with the read data able to be merged from the sequencing run. The companion Merge.fas is a [Fasta](https://en.wikipedia.org/wiki/FASTA_format) file with the reads from the samples in [Fasta](https://en.wikipedia.org/wiki/FASTA_format) format.
- TotalTable.tsv is an [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table with all of the merged, forward, and reverse results as well as the retained results for the merged reads removed due to being suspected chimeric combinations.

### Dependencies
- dada2
- ShortRead readFastq()
- ShortRead writeFastq()
  
([Back to Top](#table-of-contents))
***

## Combine Dada Output
combine_dada_output() - Combine Dada2 [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) tables (YYYY_MM_DD_HHMM_FileName_MergeFwdRev.tsv OR YYYY_MM_DD_HHMM_FileName_Merge.tsv) into a single [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) output file.

### Input 
Two or more files to be combined are required as input for this function. These files need to be [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) files as outputted from the dada_inplement() and can include Merge, MergeFwdRev, or TotalTable.tsv files.

### Arguments
- <strong> fileLoc -</strong> Select a file in the file folder with dada results you would like to combine (YYYY_MM_DD_HHMM_FileName_MergeFwdRev.tsv OR YYYY_MM_DD_HHMM_FileName_Merge.tsv OR 2023_10_16_1511_Run1_TotalTable.tsv files.
- <strong> minLen -</strong> The minimum final desired length of the read (Default minLen = 100).

### Output
The output from this function includes three files.
  1. YYYY_MM_DD_HHMM_combinedDada.tsv - combined [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table
  2. YYYY_MM_DD_HHMM_combinedDada.fas - combined [Fasta](https://en.wikipedia.org/wiki/FASTA_format) file
  3. YYYY_MM_DD_HHMM_combinedDada.tsv - Summary file from the combine_dada_output run

### Intrepretation
Outputted data files will come in the same [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table format as the output dada_implement() [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) files. 

### Dependencies
- Base R
  
([Back to Top](#table-of-contents))
***

## Make BLAST DB
make_BLAST_DB() - Create a local database to BLAST against.

### Input 
This function takes a [Fasta](https://en.wikipedia.org/wiki/FASTA_format) file (in MACER format) and establishes a database upon which a BLAST search can be completed. The outcome of the function is a folder with an NCBI database.
- The MACER [Fasta](https://en.wikipedia.org/wiki/FASTA_format) header format - ```>UniqueID|OtherInformation|Genus|species|OtherInformation|Marker```
- An example of the header format output from the MACER program is ```>GenBankAccessionOrBOLDID|GenBankAccession|Genus|species|UniqueID|Marker```
  
### Arguments
- <strong>fileLoc -</strong> The location of a file in a directory where all [Fasta](https://en.wikipedia.org/wiki/FASTA_format) files will be used to construct a BLASTable database. Default = NULL makeblastdbPath The local path for the blast+ 
- <strong>makeblastdbPath -</strong> program taxaDBLoc The location of the NCBI taxonomic data base (accessionTaxa.sql see the main DBTCShiny page for details).
- <strong>inputFormat -</strong> This will either be NCBI formatted [Fasta](https://en.wikipedia.org/wiki/FASTA_format) file (header example ```>uniqueID|NCBItaxaID|Genus|species```) or a MACER (header example ```>uniqueID|other_ID|Genus|species|Other_info|markerOrDatabase```) formatted file.
- <strong>dbName -</strong> A short 6-8 alpha character name used when building a database.
- <strong>minLen -</strong> The minimum sequence length used to construct the BLAST database.

### Output
The output from this function includes a folder with the BLAST database named according to the submitted dbName.

### Intrepretation
The constructed database can then be used with the seq_BLAST() function.

### Dependencies
- taxonomizr()
  
([Back to Top](#table-of-contents))
***

## Sequence BLAST
seq_BLAST() - Search [Fasta](https://en.wikipedia.org/wiki/FASTA_format) files of unknown sequences against a BLAST formatted database.

### Input 
Provide a location for the BLAST database you would like to use by selecting a file in the target directory (This can be a built database using the make_BLAST_DB() function or a preformatted [NCBI BLAST database](https://www.ncbi.nlm.nih.gov/books/NBK62345/#blast_ftp_site.The_blastdb_subdirectory)). Then provide the location of the query sequence files by indicating a file in a directory that contains the [Fasta](https://en.wikipedia.org/wiki/FASTA_format) files. Provide the path for the blast+ blastn program. Finally, provide the minimum query sequence length to BLAST (Default = 100), the maximum depth of the BLAST returned results (default = 200), and finally the number of cores to process the function (default = 1, Windows implementation can only use a single core and will default to this value when running on Windows).

### Arguments
- <strong>databasePath -</strong> The location of a file in a directory where the desired BLAST database is located.
- <strong>querySeqPath -</strong> The local path for the directory containing all of the [Fasta](https://en.wikipedia.org/wiki/FASTA_format) files wishing to be BLASTed
- <strong>blastnPath -</strong> The location of the NCBI blast+ blastn program (default = blastn).
- <strong>minLen -</strong> The minimum length of the sequences that will be BLASTed (default = 100).
- <strong>BLASTResults -</strong> The maximum number of returned results, or the depth of the reported results, saved from the BLAST (default = 250).
- <strong>numCores -</strong> The number of cores used to run the function (default = 1, Windows systems can only use a single core).

### Output
Two files are produced from this function, a BLAST run file and a BLAST results file for each of the [Fasta](https://en.wikipedia.org/wiki/FASTA_format) files in the target directory.

### Intrepretation
The BLAST run file contains the command used to run the BLAST search. The BLAST results file includes all results in a tab delimited .tsv file format with the columns qseqid, sseqid, staxid, qcovs, pident, ssciname, scomname, qstart, qend, sstart, send, evalue.

### Dependencies
- Base R
  
([Back to Top](#table-of-contents))
***

## Taxon Assignment
taxon_assign() - Using BLAST results to construct a table with taxonomic assignments for each unique sequence.

### Input 
This function requires a BLAST output file and an associated [Fasta](https://en.wikipedia.org/wiki/FASTA_format) file. In addition, if present an [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) file will also be used to combine the taxonomic results when present. The BLAST results are reduced to a single result for each read. 

### Arguments
- <strong>fileLoc -</strong> The location of a file in a directory where all of the paired [Fasta](https://en.wikipedia.org/wiki/FASTA_format) and BLAST (and potentially [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant)) files are located.
- <strong>taxaDBLoc -</strong> The location of the NCBI taxonomic data base (accessionTaxa.sql, see the main DBTCShiny page for details). The local path for the directory containing all of the [Fasta](https://en.wikipedia.org/wiki/FASTA_format) files wishing to be BLASTed.
- <strong>numCores -</strong> The number of cores used to run the function (default = 1, Windows systems can only use a single core).
- <strong>coverage -</strong> The percent coverage used for taxonomic assignment for the above threshold results (Default coverage = 95).
- <strong>ident -</strong> The percent identity used for the taxonomic assignment for above threshold results (Default ident = 95).
- <strong>propThres -</strong> The proportional threshold flags the final result based on the preponderance of the data. So if the threshold is set to 0.95, results will be flagged if the taxa directly below the assigned taxa has fewer than 0.95 percent of the records causing the upward taxonomic placement (Default ident = 0.95).
- <strong>coverReportThresh -</strong> The percent coverage threshold used for reporting flags below this threshold (Default coverReportThresh = 95).
- <strong>identReportThresh -</strong> The percent identity threshold used for reporting flags below this threshold (Default identReportThresh = 95).
- <strong>includeAllDada  -</strong> When paired Dada [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) tables are present, when set to FALSE, this will exclude records without taxonomic assignment (Default includeAllDada = TRUE).

### Output
A single taxonomic assignment file is created contains the string 'taxaAssign_YYYY_MM_DD_HHMM'.
<a id="taxon-assignment-intrepretation"></a>
### Intrepretation
The number of returned BLAST results is dictated by the seq_BLAST() BLASTResults argument. The taxon_assign() function takes into account all returned BLAST results for each read. At each taxonomic level assignmenta have quality metrics in parentheses after the name. These values ("Num_Rec", "Coverage", "Identity", "Max_eVal") represent the number of records with this taxonomic placement, the minimum coverage and identity, and the maximum eValue for the reported taxa.

Column headers for the resulting taxonomic assignments include...

uniqueID, superkingdom, phylum, class, order, family, genus, species, Top_BLAST, Lowest_Single_Ran, Lowest_Single_Taxa, Lowest_Single_Rank_Above_Thres, Lowest_Single_Taxa_Above_Thres, Final_Common_Names, Final_Rank, Final_Taxa, Final_Rank_Taxa_Thres, Result_Code, Sequence, Length, Results, Followed by columns of samples containing [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) values

There are three columns that deserve special explaination. 

The Final_Rank_Taxa_Thres column contains the threshold values (Coverage, Identitiy) applied to the final rank and taxonomic values for the associated records.  
The Results column contains the reference to the record if it was able to be merged, or if it is representing a forward or reverse unidirectional read.
The Result_Code column contains flags placed on the results to better understand the quality of the resulting taxonomic assignments. Below is a list of codes. T

  - SFAT(coverage, ident): Saturated filtered taxa above threshold
  - SANF(coverage, ident): Saturated non-filtered
  - BIRT(identReportThresh): Final taxa result is below the identity reporting threshold
  - BCRT(coverReportThresh): Final taxa below the nucleotide coverage reporting threshold
  - TBAT(propThres): Taxa Below Assigned Taxa Threshold

Note: Records with BIRT, BCRT, and TBAT flags should be highly scruntized. SANF results should be explored and the size of the database, the trust placed in the records in the database, and the depth of the BLAST results should be considered when assessing records with this flag. Records with SFAT are among the least concerning as the BLAST results were saturated but this taxonomic assignment saturation occurred above your set quality coverage and identity threshold. Concerns with records with this result could be that the depth of the BLAST analysis was not low enough for very large databases, or that the database is not complete (taxonomic breadth) enough for smaller databases.

### Dependencies
- taxonomizr()
- pbapply()
  
([Back to Top](#table-of-contents))
***

## Combine Assignment Output
combine_assign_output() - Using results from the taxon_assign() function, combine all files with the string 'taxaAssign_YYYY_MM_DD_HHMM' in to a single .tsv file.

### Input 
Select a file in a folder with the taxa assigned files you would like to combine (extension '_taxaAssign_YYYY_MM_DD_HHMM.tsv'). NOTE: all '_taxaAssign_' files in the folder location should originate from the same dada output file but have outputs from different BLAST sequence libraries and therefore contain the same [ASV's](https://en.wikipedia.org/wiki/Amplicon_sequence_variant).

### Arguments
- <strong>fileLoc -</strong> The location of a file in a directory where all of the 'taxaAssign' files are located.
- <strong>numCores -</strong> The number of cores used to run the function (default = 1, Windows systems can only use a single core).

### Output
This function produces a '2023_08_03_0913_taxaAssignCombined.tsv' and a '2023_08_03_0913_taxaAssignCombined.txt' file in the selected target directory.

### Intrepretation
The intrepretation of the output file for the combine_assign_output() 'taxaAssignCombined' files is the same as the is the same as the taxon_assign() 'taxaAssign' files (see above description).

### Dependencies
- pbapply()
  
([Back to Top](#table-of-contents))
***

## Reduce Taxa
reduce_taxa() - Using results from taxon_assign() and/or combine_assign_output() this function combines all reads with the same taxonomic assignment into a single result.

### Input 
This function requires a file in a directory where all 'taxaAssign' and/or 'taxaAssignCombine' files in that directory will be combined. All records with the same taxonomic result will be combined. The BLAST values in parentheses ("Num_Rec", "Coverage", "Identity", "Max_eVal") are combine by the mean number of records, the minimum coverage and identity, and the maximum eValue.

### Arguments
- <strong>fileLoc -</strong>  The location of a file in a directory where all of the 'taxaAssign' and/or 'taxaAssignCombine' files are located.
- <strong>numCores -</strong>  The number of cores used to run the function (default = 1, Windows systems can only use a single core).

### Output
This function produces a '_CombineTaxaReduced.tsv' file for every 'taxaAssign' or 'taxaAssignCombine' present in the target directory.

### Intrepretation
Reduced taxonomic assignment files have fewer columns in the main taxa_reduced.tsv file than the taxaAssign files as columns are collapsed. In addition, the values in the taxonomic columns in parentheses represent the average values across all of the results with the same taxonomic assignment (see taxon_assign() intrepretation above).

The columns include, superkingdom, phylum, class, order, family, genus, species, Top_BLAST, Final_Common_Names, Final_Rank, Final_Taxa, Result_Code, RepSequence, Number_ASV, Average_ASV_Length, Number_Occurrences, Average_ASV_Per_Sample, Median_ASV_Per_Sample, Results.


Add in a note about the representative sequence for the taxa and how I got it.



### Dependencies
- pbapply()

([Back to Top](#table-of-contents))
***

## Combine Reduced Output
combine_reduced_output() - This function takes 'taxaReduced' files generated from the same biological samples but representing different amplified molecular markers and collapses these data into a single file. The outcome of this process results in a presence absence matrix for all taxa and markers.

### Input 
Select a file in a folder with 'taxaReduced' files representing data for the same biological samples but representing different amplified molecular markers.

### Arguments
There are only two arguments necessary for this function. The first is the location of a file in the folder with all of the files that are wanting to be combined. The second is the combineReduceTaxa function which will reduce all read count values to binary 1 or 0 entries. 
- <strong>fileLoc -</strong>  The location of a file in a directory where all of the 'taxaReducedAssign' files are located.
- <strong>presenceAbsence -</strong>  A TRUE or FALSE value used to indicate if the read values should be replaced with presence/absence (1/0) data. This change is necessary when combining taxa for the same named samples across molecular markers (TRUE) but is not necessary when combining results for taxa with all unique sample names (FALSE). Note: For visulization on the DBTCShiny mapping component displaying read intensity all samples should have unique names and this value should be set to FALSE (default = TRUE to replace values with presence/absence 1/0 values).

### Output
Two files, a CombineTaxaReduced.tsv result file and a CombineTaxaReduced.txt run summary file are generated from this function. The result file contains presence/absence data in a matrix that associates the data with samples, taxa, and molecular marker. The column headers in the results file includes the following, superkingdom, phylum, class, order, family, genus, species, markers(n number of columns), samples (n number).

### Intrepretation
There is no specific unique intrepretation for this file.

### Dependencies
- plyr() rbind.fill

([Back to Top](#table-of-contents))
***

# Mapping Dashboard

## Mapping
The interactive map where loaded data can be viewed.

## Data Import
Data import buttons to load [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) data files generated by the DBTCShiny pipeline along with provenance data to visualize on the map.

## Data Filtering
Filtering options to selectively view different data on the map.

## Data Table
A tabular display of the data loaded and filtered and appearing on the map. 

([Back to Top](#table-of-contents))
