# Written by Rob Young at the University of Guelph in Ontario Canada, Sept, 2023
#********************************************Main program section***********************************************
#Roxygen2 Documentation:
#' @export
#'
#' @title Launch DBTCShiny
#'
#' @author Robert G. Young
#'
#' @description
#' This function launches the DBTCShiny Application
#'
#' @details
#' This function launches a DBTCShiny Application which allows a user to run the DBTC
#' functions and process high throughput sequencing data.
#'
#' @examples
#' \dontrun{
#' launchDBTCShiny()
#' }
#'
#' @returns
#' There are no values or files returned from this function
#'
#' @references
#' <https://github.com/rgyoung6/DBTC>
#' Young, R. G., Hanner, R. H. (Submitted October 2023). Metabarcoding analysis
#' using Dada-BLAST-Taxon Assign-Condense Shiny Application (DBTCShiny). Biodiversity Data Journal.
#'
#' @note
#' This is a wrapper function which launches the DBTCShiny package as a shiny
#' application in the systems default browser program
#'
#' @seealso
#' dada_implement()
#' combine_dada_output()
#' make_BLAST_DB()
#' seq_BLAST()
#' taxon_assign()
#' combine_assign_output()
#' reduce_taxa()
#' combine_reduced_output()


# wrapper for shiny::shinyApp()
launchDBTCShiny <- function() {
  shiny::shinyApp(ui = shinyAppUI, server = shinyAppServer, options = list(launch.browser = TRUE))
}
