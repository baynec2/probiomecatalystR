#' create_fasta_from_organism_ids
#'
#' Given a vector of organism ids, pull down a list of proteomes from uniprot
#' for each organism. The first one from the list is then downloaded as a fasta
#' file. This is then repeated for each organism id and concatintated into a
#' combined fasta file.
#'
#' This is particularly usefull for generating a metaproteomics .fasta database
#'
#' @param organism_ids vector of organism ids
#' @param destination_fp path to the destination fasta file. Must be .fasta
#' @param additonal_organism_ids vector of any additional organism ids to include
#'
#' @returns a downloaded .fasta file to user specified file path.
#' @export
#'
#' @examples
#' # Making a database with Human and E.coli
#'create_fasta_from_organism_ids(organism_ids = c(9606, 562))
get_fasta_from_organism_ids <- function(organism_ids,
                                           destination_fp = paste0(getwd(),"/",
                                                                   Sys.Date(),
                                                                   ".fasta"),
                                           additonal_organism_ids = NULL) {
  if (file.exists(destination_fp)) {
    # Delete any FASTA file at the output location if it already exists
    unlink(destination_fp)
    cat("File deleted at", destination_fp, "\n")
  }
  # making sure organism_ids are unique
  organism_ids <- unique(organism_ids)
  # Adding any additional specified organism ids to the list
  if (!is.null(additonal_organism_ids)) {
    organism_ids <- c(organism_ids, additonal_organism_ids)
  }
  # Dowloading list of reference proteomes. This will be used to
  rp = get_all_reference_proteomes() |>
    dplyr::pull(.data$Proteome_ID)
  # Create temp directory at the destination of database
  temp_dir <- paste0(dirname(destination_fp), "/temp/")
  dir.create(temp_dir)

  ##############################################################################
  # Download a fasta file for each organism id
  ##############################################################################
  for (i in seq_along(organism_ids)) {
    organism_id <- organism_ids[i]
    tryCatch(
      {
        # UniProt API request to get the Proteome ID
        request_url <- paste0(
          "https://rest.uniprot.org/proteomes/search?query=organism_id:",
          organism_id, "&format=tsv"
        )
        org_to_rp <- read.csv(request_url, header = TRUE, sep = "\t")
        # Taking first proteome id
        if(org_to_rp$Proteome.Id[1] %in% rp){
          first_rp <- org_to_rp |>
            dplyr::filter(.data$Proteome.Id %in% rp) |>
            dplyr::pull(.data$Proteome.Id)
        } else {
          # If we cant find a reference, take the first one returned.
          first_rp <- org_to_rp$Proteome.Id[1]
        }
        if (!is.na(first_rp)) {
          # UniProt API request to get the Proteome ID
          request_url <- paste0(
            "https://rest.uniprot.org/uniprotkb/stream?query=proteome:",
            first_rp, "&format=fasta"
          )
          # Downloading reference proteome
          reference_proteome <- download.file(request_url,
            destfile = paste0(temp_dir, first_rp, ".fasta")
          )
        }
      },
      error = function(e) {
        message(paste("Error processing organism ID", organism_id, ":",
                      e$message))
      }
    )
  }
  ##############################################################################
  # Concatenate all FASTA files in the temp directory into one file
  ##############################################################################
  # List all FASTA files in the temp directory
  fasta_files <- list.files(temp_dir, pattern = "\\.fasta$", full.names = TRUE)
  # Check if there are FASTA files to concatenate
  if (length(fasta_files) > 0) {
    # Loop through each FASTA file and append its content
    for (i in fasta_files) {
      # Read the content of each file
      fasta_content <- readLines(i)
      if (file.exists(destination_fp)) {
        readr::write_lines(fasta_content, destination_fp, append = TRUE)
      } else {
        readr::write_lines(fasta_content, destination_fp)
      }
    }
    # Print informative message
    cat("FASTA files concatenated into:", destination_fp, "\n")

    # Delete the temp directory and its contents
    unlink(temp_dir, recursive = TRUE)
    cat("Temporary directory deleted.\n")
  } else {
    cat("No FASTA files found in the temporary directory.\n")
  }
}
