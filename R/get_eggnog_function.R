#' get_eggnog_function
#'
#' Get Eggnog functional go terms given a vector of EGGNOG ids
#'
#' @param eggnog_ids
#'
#' @returns
#' @export
#'
#' @examples
#' insulin = get_eggnog_function("ENOG502S5P5")
get_eggnog_function <- function(eggnog_ids) {

  # Initialize an empty results tibble
  results <- tibble::tibble(
    eggnog_id = character(), go_type = character(), go_id = character(),
    go_description = character(), evidence_codes = character(),
    seq_count = numeric(), frequency = numeric(), annotation_count = numeric()
  )

  # Loop through each EggNOG ID
  for (id in eggnog_ids) {
    url <- paste0("http://eggnogapi5.embl.de/nog_data/json/go_terms/", id)
    response <- httr::GET(url)

    if (httr::status_code(response) == 200) {
      data <- jsonlite::fromJSON(httr::content(response, "text",encoding = "UTF-8"))

      if ("go_terms" %in% names(data)) {
        col_names <- c("go_id", "go_description", "evidence_codes",
                       "seq_count", "frequency", "annotation_count")

        # Extract GO terms per category
        go_data <- list(
          MF = data$go_terms$`Molecular Function`,
          BP = data$go_terms$`Biological Process`,
          CC = data$go_terms$`Cellular Component`
        )

        go_df <- dplyr::bind_rows(lapply(names(go_data), function(go_type) {
          if (!is.null(go_data[[go_type]])) {
            as.data.frame(go_data[[go_type]]) |>
              setNames(col_names) |>
              dplyr::mutate(
                eggnog_id = id, go_type = go_type, .before = 1,
                seq_count = as.numeric(seq_count),
                frequency = as.numeric(frequency),
                annotation_count = as.numeric(annotation_count)
              )
          } else {
            tibble::tibble(eggnog_id = id, go_type = go_type,
                           go_id = NA_character_,go_description = NA_character_,
                           evidence_codes = NA_character_,seq_count = NA_real_,
                           frequency = NA_real_, annotation_count = NA_real_)
          }
        }))

        results <- dplyr::bind_rows(results, go_df)
        cat("GO terms downloaded for Eggnog ID:", id, "\n")
      } else {
        # No GO terms found, add empty row
        results <- dplyr::bind_rows(results, tibble(
          eggnog_id = id, go_type = NA_character_, go_id = NA_character_,
          go_description = NA_character_, evidence_codes = NA_character_,
          seq_count = NA_real_, frequency = NA_real_, annotation_count = NA_real_
        ))
        cat("Warning: No GO terms found for Eggnog ID:", id, "\n")
      }
    } else {
      # Failed API request, add empty row
      results <- bind_rows(results, tibble(
        eggnog_id = id, go_type = NA_character_, go_id = NA_character_,
        go_description = NA_character_, evidence_codes = NA_character_,
        seq_count = NA_real_, frequency = NA_real_, annotation_count = NA_real_
      ))
      cat("Error: Request failed for Eggnog ID:", id, "\n")
    }
  }
  out = tibble::as_tibble(results)
}
