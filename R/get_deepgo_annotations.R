library(httr)
library(jsonlite)
library(tibble)
library(dplyr)

get_deepgo_annotations <- function(protein_sequence, release = "33", data_format = "fasta", threshold = 0.3) {
  # DeepGO API Endpoint
  url <- "https://deepgo.cbrc.kaust.edu.sa/deepgo/api/create"

  # Construct request payload
  payload <- list(
    release = release,
    data_format = data_format,
    threshold = threshold,
    data = protein_sequence
  )

  # Send POST request
  response <- POST(url, body = payload, encode = "json", content_type_json())

  # Check if request was successful
  if (status_code(response) == 200) {
    # Parse JSON response
    result <- fromJSON(content(response, "text", encoding = "UTF-8"))

    # Extract predictions
    if (!is.null(result$predictions$functions)) {
      predictions <- result$predictions$functions

      # Convert to tibble
      prediction_tibble <- as_tibble(predictions)
      return(prediction_tibble)
    } else {
      cat("No GO annotations found.\n")
      return(NULL)
    }
  } else {
    stop("Error: Failed to retrieve DeepGO annotations.")
  }
}

# # Example usage
# protein_seq <- ">sp|P80386|AAKB1_RAT 5'-AMP-activated protein kinase subunit beta-1 OS=Rattus norvegicus GN=Prkab1 PE=1 SV=4\nMGNTSSERAALERQAGHKTPRRDSSGGTKDGDRPKILMDSPEDADIFHTEEMKAPEKEEF\nLAWQHDLEVNEKAPAQARPTVFRWTGGGKEVYLSGSFNNWSKLPLTRSQNNFVAILDLPE\nGEHQYKFFVDGQWTHDPSEPIVTSQLGTVNNIIQVKKTDFEVFDALMVDSQKCSDVSELS\nSSPPGPYHQEPYISKPEERFKAPPILPPHLLQVILNKDTGISCDPALLPEPNHVMLNHLY\nALSIKDGVMVLSATHRYKKKYVTTLLYKPI"
#
# deepgo_results <- get_deepgo_annotations(protein_seq)
# print(deepgo_results)
#
# # Example usage
# protein_sequence <- ">sp|P80386|AAKB1_RAT 5'-AMP-activated protein kinase subunit beta-1 OS=Rattus norvegicus GN=Prkab1 PE=1 SV=4
# MGNTSSERAALERQAGHKTPRRDSSGGTKDGDRPKILMDSPEDADIFHTEEMKAPEKEEF
# LAWQHDLEVNEKAPAQARPTVFRWTGGGKEVYLSGSFNNWSKLPLTRSQNNFVAILDLP
# EGEHQYKFFVDGQWTHDPSEPIVTSQLGTVNNIIQVKKTDFEVFDALMVDSQKCSDVSE
# LSSSPPGPYHQEPYISKPEERFKAPPILPPHLLQVILNKDTGISCDPALLPEPNHVMLN
# HLYALSIKDGVMVLSATHRYKKKYVTTLLYKPI"
#
# annotations <- get_deepgo_annotations(protein_sequence)

