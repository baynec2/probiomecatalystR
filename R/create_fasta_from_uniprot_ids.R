#' create_fasta_from_uniprot_ids
#' Creates a .fasta database from a list of uniprot ids.
#' file from the first proteome returned via uniprots API.
#'
#' @param uniprot_ids = vector of uniprot ids
#' @param fasta_out_fp = file path to write the .fasta database too
#'
#' @returns a downloaded concatinated
#' @export
#'
#' @examples
#' # Human and mouse insulin, will pull down the human and mouse proteome.
#' create_fasta_from_uniprot_ids(c("P01325", "P01308"),fasta_out_fp = "test.fasta")
create_fasta_from_uniprot_ids <- function(uniprot_ids,
                                          fasta_out_fp = paste0(getwd(),
                                                                "/",
                                                                Sys.Date(),
                                                                ".fasta")) {
  ## Pull the required data for each uniprot id from Uniprot.
  annotated_up <- annotate_uniprot_ids(uniprot_ids,
    batch_size = 100,
    columns = "accession,id,protein_name,organism_name,organism_id,gene_primary,protein_existence,sequence_version,sequence"
  )
  # Filter out missing sequences
  na_omit <- annotated_up %>% dplyr::filter(sequence != "")

  # Warning message about # of IDs with no sequence in current release.
  message(paste0("Note that ", nrow(annotated_up) - nrow(na_omit), " uniprot ids in your sequence set have been deleted from the most recent uniprot release"),
    appendLF = TRUE
  )

  # Format data frame appropriately
  df <- na_omit %>%
    dplyr::mutate(header = paste0(
      ">sp|", accession, "|", id, " ", protein_name,
      " OS=", organism_name, " OX=", organism_id, " GN=",
      gene_primary, " PE=", protein_existence, " SV=", sequence_version
    ), ) %>%
    dplyr::select(header, sequence)

  # writing to file
  fasta_out_dir <- dirname(fasta_out_fp)
  file_conn <- file(paste0(fasta_out_dir, "/", fasta_out_fp), open = "w")

  # Loop through each row in the tibble
  for (i in 1:nrow(df)) {
    # Write the header line to the file
    writeLines(df$header[i], file_conn)
    # Write the sequence to the file
    writeLines(df$sequence[i], file_conn)
  }
  # Close the file connection
  close(file_conn)
  # Print where the file was written to
  message(paste0("Fasta file written to: ", fasta_out_dir, "/", fasta_out_fp))
}
