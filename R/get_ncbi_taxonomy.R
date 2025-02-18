#' get_ncbi_taxonomy
#'
#' given a list of ncbi taxonomy, fetch complete taxonomic information.
#'
#' @param ncbi_ids a vector of ncbi taxonomy ids
#'
#' @returns
#' @export
#'
#' @examples
get_ncbi_taxonomy <- function(ncbi_ids) {
  taxonomy_list <- list()
  for (id in ncbi_ids) {
    tryCatch({
      # Query NCBI Taxonomy
      xml_data <- rentrez::entrez_fetch(db = "taxonomy", id = id, rettype = "xml")

      # Parse the XML data using XML package
      xml_parsed <- XML::xmlParse(xml_data)

      # Extract taxonomic lineage
      ranks <- XML::xpathSApply(xml_parsed,
                                "//LineageEx/Taxon/Rank",
                                XML::xmlValue)
      names <- XML::xpathSApply(xml_parsed,
                                "//LineageEx/Taxon/ScientificName",
                                XML::xmlValue)

      # Extract scientific name of the queried organism
      species <- XML::xpathSApply(xml_parsed,
                                        "//TaxaSet/Taxon/ScientificName",
                                        XML::xmlValue)

      # Extract tax ID
      organism_id <- XML::xpathSApply(xml_parsed,
                                 "//TaxaSet/Taxon/TaxId",
                                 XML::xmlValue)

      # Combine all results into one dataframe
      # Define expected taxonomy ranks
      expected_ranks <- c("superkingdom", "kingdom", "phylum", "class", "order",
                          "family", "genus")

      taxonomy_df <- data.frame(
        organism_id = organism_id,
        species = species,
        rank = ranks,
        name = names,
        stringsAsFactors = FALSE
      ) |>
        # Ensure ranks are only those in expected_ranks
        dplyr::mutate(rank = ifelse(rank %in% expected_ranks, rank, NA)) |>
        tidyr::drop_na(rank) |>  # Remove any entries with NA rank
        # Ensure all expected ranks are present, filling missing ones with NA
        tidyr::complete(rank = expected_ranks, fill = list(organism_id = organism_id,
                                                           species = species,
                                                           name = NA_character_))

      taxonomy_list[[id]] <- taxonomy_df
    }, error = function(e) {
      message(paste("Error fetching taxonomy for ID:", id))
    })

}

  final_df <- dplyr::bind_rows(taxonomy_list) |>
    tidyr::pivot_wider(names_from = rank,values_from = name) |>
    dplyr::select(c("organism_id","superkingdom", "kingdom", "phylum", "class",
                    "order", "family", "genus","species"))

  return(final_df)
}
