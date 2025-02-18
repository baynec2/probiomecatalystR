#' extract fasta info
#'
#' extracts organism name, organism id and sequence from fasta file
#'
#' @param fasta_file
#'
#' @returns a data frame with organism name, organism id and sequence
#' @export
#'
#' @examples
extract_fasta_info <- function(fasta_file) {
  # Read the FASTA file
  fasta_data <- Biostrings::readAAStringSet(fasta_file)
  sequences <- as.character(fasta_data)  # Convert the sequences to character strings

  # Extract headers and sequences
  headers <- names(fasta_data)

  # Extract organism name (OS=...)
  organism_name <- ifelse(
    grepl("OS=", headers),
    sub(".*OS=([^=]+?)\\s*OX=.*", "\\1", headers),
    NA
  )
  organism_name <- trimws(organism_name)  # Remove extra spaces

  # Extract organism ID (OX=...)
  organism_id <- ifelse(
    grepl("OX=", headers),
    sub(".*OX=([0-9]+).*", "\\1", headers),
    NA
  )
  organism_id <- as.integer(organism_id)  # Convert to integer

  # Extract protein ID (first field in UniProt-style headers)
  protein_id <- ifelse(
    grepl("\\|", headers),
    sub("^.*\\|([^|]+)\\|.*$", "\\1", headers),
    NA
  )


  # Create a DataFrame
  df <- data.frame(
    protein_id = protein_id,
    organism_name = organism_name,
    organism_id = organism_id,
    sequence = sequences,
    stringsAsFactors = FALSE
  )

  return(df)
}
