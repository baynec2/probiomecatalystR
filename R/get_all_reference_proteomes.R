#' get_all_reference_proteomes
#'
#' This function returns all current reference proteomes from uniprot. This
#' gets updated every 8 weeks per their documentation. See
#' https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README
#' for more info.
#'
#' Contains the following columns:
#' Proteome_ID,Tax_ID,OSCODE,SUPERREGNUM,#(1),#(2),#(3),Species Name
#' #(1),#(2),#(3) are statistics as follows:
# (1) Number of entries in main fasta (canonical)
# (2) Number of entries in additional fasta (isoforms)
# (3) Number of entries in gene2acc mapping file
#'
#' @returns a tibble
#' @export
#'
#' @examples
#' rps <- get_all_reference_proteomes()
get_all_reference_proteomes <- function() {
  current_reference_proteomes <- readr::read_tsv("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README",
    skip = 153,skip_empty_rows = TRUE
  )
  if (sum(names(current_reference_proteomes) == c(
    "Proteome_ID", "Tax_ID", "OSCODE", "SUPERREGNUM", "#(1)", "#(2)", "#(3)",
    "Species Name"
  )) == 8) {
    return(current_reference_proteomes)
  } else {
    stop("Unexpected column names are returned from the file")
  }
}
