#' annotate_uniprot_ids
#' Function to interact with the uniprot api and get data corresponding to each
#' uniprot id.
#' @param uniprot_ids uniprot acession numbers
#' @param columns This is a string  of the column names that you would like to include followed by commas.
#' Column names accepted by the API can be found here: https://www.uniprot.org/help/return_fields.
#' If not specified, the default columns are returned.
#' @param batch_size this is the size to include in a batch search
#'
#' @returns a tibble containing requested information.
#' @export
#'
#' @examples
#' annotate_uniprot_ids("P01308")
annotate_uniprot_ids <- function(uniprot_ids,
                                 columns = NULL,
                                 batch_size = 150) {
  # Define the regex pattern for valid UniProt accessions
  uniprot_regex <- "^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"

  # Function to validate UniProt accessions
  is_valid_uniprot_accession <- function(uniprot_ids) {
    valid_acessions <- uniprot_ids[grepl(uniprot_regex, uniprot_ids)]
    n_invalid <- length(uniprot_ids) - length(valid_acessions)
    message(paste0(n_invalid, " of the ids in list are not valid uniprot numbers"))
    return(valid_acessions)
  }

  uniprot_ids_filtered <- is_valid_uniprot_accession(uniprot_ids)

  # Base URL for request
  baseUrl <- "https://rest.uniprot.org/uniprotkb/search?query=accession:"

  # Check for long lists and break them into batches if necessary
  if (length(uniprot_ids_filtered) > batch_size) {
    num_iterations <- ceiling(length(uniprot_ids_filtered) / batch_size)
    # Calculate the start and end indices for each batch
    start <- seq(1, length(uniprot_ids_filtered), by = batch_size)
    end <- c(
      seq(batch_size, length(uniprot_ids_filtered), by = batch_size),
      length(uniprot_ids_filtered)
    )
  } else {
    num_iterations <- 1
    start <- 1
    end <- length(uniprot_ids_filtered)
  }

  # Set up the progress bar
  pb <- txtProgressBar(
    min = 0, max = num_iterations, style = 3,
    width = 50, char = "="
  )

  # Initialize an empty data frame to hold results
  output <- data.frame()

  # Iterate through the batches
  for (i in 1:num_iterations) {
    # Extract the current batch of protein IDs
    temp_list <- uniprot_ids_filtered[start[i]:end[i]]
    # Join them into a single query string
    temp_string_list <- paste(temp_list, collapse = "+OR+")

    # Construct the request URL
    if (is.null(columns)) {
      request <- paste0(
        baseUrl, temp_string_list, "&format=tsv",
        "&size=", length(temp_list)
      )
    } else {
      request <- paste0(
        baseUrl, temp_string_list, "&format=tsv",
        "&fields=", columns, "&size=", length(temp_list)
      )
    }

    # Fetch the data, handling errors for each individual protein
    returned_data <- tryCatch(
      {
        # Attempt to read the data
        read.csv(request, header = TRUE, sep = "\t")
      },
      error = function(cond) {
        # Print the error and return NULL for the batch
        message(paste("Request failed for batch:", i))
        return(NULL)
      }
    )

    # If the batch was successful, append the data to the output
    if (!is.null(returned_data)) {
      output <- dplyr::bind_rows(output, returned_data)
    } else {
      # If the batch fails, search each protein individually
      message("Batch failed, querying proteins individually...")
      # Set up the progress bar
      pb2 <- txtProgressBar(
        min = 0, max = length(temp_list), style = 3,
        width = 50, char = "="
      )

      for (protein_id in temp_list) {
        individual_request <- paste0(
          baseUrl, protein_id, "&format=tsv",
          "&size=1"
        )

        individual_data <- tryCatch(
          {
            # Try fetching data for the individual protein
            read.csv(individual_request, header = TRUE, sep = "\t")
          },
          error = function(cond) {
            # In case of an error, print a message
            message(paste("Request failed for protein:", protein_id))
            setTxtProgressBar(pb2, i)
            return(NULL)
          }
        )

        # Append individual protein data to output
        if (!is.null(individual_data)) {
          output <- dplyr::bind_rows(output, individual_data)
        }
      }
    }

    # Update the progress bar
    setTxtProgressBar(pb, i)
  }

  # Convert the result to a tibble and return
  output <- tibble::as_tibble(output)

  # fix names
  # Split the string by commas
  if (!is.null(columns)) {
    new_names <- unlist(strsplit(columns, split = ","))
    names(output) <- new_names
  }
  return(output)
}
