#' get_kegg_in_batches
#'
#' Given a vector of kegg ids, get all of the functional information using the
#' KEGG API.
#'
#' @param kegg_ids = kegg ides in the format "hsa:1234"
#'
#' @returns
#' @export
#'
#' @examples
#'
#' insulin = get_kegg_in_batches("hsa:3630")
#'
get_kegg_in_batches <- function(kegg_ids, batch_size = 10) {
  all_results <- list()

  # Split IDs into batches
  batches <- split(kegg_ids, ceiling(seq_along(kegg_ids) / batch_size))

  for (i in seq_along(batches)) {
    message(sprintf("Processing batch %d of %d", i, length(batches)))

    batch <- batches[[i]]

    # Fetch KEGG entries
    entries <- lapply(batch, function(id) {
      tryCatch(KEGGREST::keggGet(id), error = function(e) NULL)
    })

    # Flatten list and remove NULLs
    entries <- unlist(entries, recursive = FALSE)
    entries <- entries[!sapply(entries, is.null)]

    # Process each entry and store results
    if (length(entries) > 0) {
      results <- purrr::map_dfr(entries, function(entry) {
        if (is.null(entry)) return(NULL)  # Handle missing entries gracefully

        gene_num <- entry$ENTRY  # Gene ID
        org_id <- names(entry$ORGANISM)  # Organism ID
        gene_id <- paste0(org_id, ":", gene_num)  # Full gene ID

        # Extract KO ID safely
        ko <- if (!is.null(entry$ORTHOLOGY)) names(entry$ORTHOLOGY) else NA_character_

        # Extract KEGG Pathway
        kegg_pathway <- if (!is.null(entry$PATHWAY)) {
          paste0(paste0(entry$PATHWAY, " [", names(entry$PATHWAY), "]"), collapse = "; ")
        } else {
          NA_character_
        }

        # # Extract BRITE hierarchy
        # brite_info <- if (!is.null(entry$BRITE)) {
        #   paste0(unlist(entry$BRITE), collapse = "; ")
        # } else {
        #   NA_character_
        # }

        # Return extracted info as a tibble
        tibble::tibble(
          kegg_id = gene_id,
          ko = ko,
          kegg_pathway = kegg_pathway#,
          # brite_info = brite_info
        )
      })

      all_results[[i]] <- results
    }

    # Avoid hitting KEGG rate limits (3x per second)
    Sys.sleep(0.34)
  }

  # Combine all batch results into a single dataframe
  out = dplyr::bind_rows(all_results) |>
   # Formatting into a useful format
    dplyr::mutate(kegg_pathway_list= stringr::str_split(kegg_pathway, "; ")) |>
    tidyr::unnest(kegg_pathway_list) |>
    dplyr::mutate(
      kegg_pathway = stringr::str_extract(kegg_pathway_list, "^[^\\[]+"),  # Extract pathway name (before bracket)
      code = stringr::str_extract(kegg_pathway_list, "\\[([^\\]]+)\\]")  # Extract code (inside brackets)
    ) |>
    dplyr::mutate(
      pathway = stringr::str_trim(kegg_pathway),  # Clean up any leading/trailing spaces
      code = stringr::str_remove_all(code, "[\\[\\]]")  # Remove brackets from code
    ) |>
    dplyr::select(kegg_pathway,kegg_id,ko,code)

}

