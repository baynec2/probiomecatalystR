#' perform_limma_analysis
#'
#' Allows user to perfom limma analysis on any "assay" in the QFeatures object
#'
#' @param qf_object qfeatures object
#' @param assay_name assay name
#' @param design_df design_data frame
#' @param formula formula
#' @param contrast contrasts to look for
#' @param normalize whether to log2 normalise
#'
#' @returns
#' @export
#'
#' @examples
perform_limma_analysis <- function(qf_object,
                                   assay_name,
                                   design_df,
                                   formula,
                                   contrast,
                                   normalize = TRUE) {

  # Extract the selected assay from QFeatures
  if (!(assay_name %in% names(qf_object))) {
    stop("Assay name not found in QFeatures object.")
  }

  # Extract a specific assay (e.g., "species")
  assay <- SummarizedExperiment::assay(qf_object[[assay_name]])

  rd <- SummarizedExperiment::rowData(qf_object[[assay_name]]) %>%
    as.data.frame() %>%
    dplyr::mutate(across(1, ~ ifelse(is.na(.), "NA", .)))

  # Combine row data and assay values
  exprs_matrix <- assay

  rownames(exprs_matrix) <- rd[,assay_name]

  if(assay_name == "go_description"){
    rownames(exprs_matrix) <- paste0(rd[,"organism_type"],"_", rd[,assay_name])
  } else{
    rownames(exprs_matrix) <- rd[,assay_name]
  }

  # Check row names (should be protein or taxon identifiers)
  if (is.null(rownames(exprs_matrix))) {
    stop("Row names (identifiers) are missing in the selected assay.")
  }

  # Log transformation if needed
  if (normalize) {
    exprs_matrix <- log2(exprs_matrix + 1)  # Avoid log(0) issues
  }

  # Ensure design_df has row names matching column names of exprs_matrix
  if (!all(colnames(exprs_matrix) %in% rownames(design_df))) {
    stop("Column names of expression matrix must match row names of design matrix.")
  }

  # Convert formula to a model matrix
  design_matrix <- model.matrix(formula, data = design_df)

  # Fit the limma model
  fit <- limma::lmFit(exprs_matrix, design_matrix)

  # Apply empirical Bayes smoothing
  fit <- limma::eBayes(fit)

  # Extract contrast results
  contrast_matrix <- limma::makeContrasts(contrasts = contrast,
                                          levels = design_matrix)
  fit_contrasted <- limma::contrasts.fit(fit, contrast_matrix)
  fit_contrasted <- limma::eBayes(fit_contrasted)

  # Get results
  results <- limma::topTable(fit_contrasted, coef = 1, number = Inf,
                             adjust.method = "BH") %>%
    tibble::rownames_to_column(var = "ID") %>%
    tibble::as_tibble()

  return(results)
}
