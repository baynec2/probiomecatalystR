#' prepare_taxonomy_matricies
#'
#' @param full_taxonomy_fp file path to taxonomy file
#' @param pr_group_matrix_fp file path to pr_group matrix from DIAnn
#' @param fasta_fp file path to fasta used for experiment.
#'
#' @returns taxonomy_matrix.tsv for each level of taxonomy
#' @export
#'
#' @examples
prepare_taxonomy_matricies = function(full_taxonomy_fp,
                                      pr_group_matrix_fp ,
                                      fasta_fp,
                                      output_dir = getwd()){
  # Reading in the files
  full_taxonomy = readr::read_delim(full_taxonomy_fp)
  peptide_groups = readr::read_tsv(pr_group_matrix_fp)
  fasta_info = GLabR::extract_fasta_info(fasta_fp)

  # Joining data
  peptide_and_org = dplyr::inner_join(peptide_groups,fasta_info,
                                      by = c("Protein.Ids" = "protein_id"))

  peptide_and_taxa = dplyr::inner_join(peptide_and_org,full_taxonomy,
                                       by = c("organism_id" = "tax_id"))


  taxa_units = c("superkingdom","kingdom","phylum","class","order","family",
  "genus","species")

  # Initialize an empty list to store results

  for (i in taxa_units) {
    df <- peptide_and_taxa %>%
      tidyr::pivot_longer(cols = contains("raw"), names_to = "colData",
                          values_to = "intensity") %>%
      dplyr::select(i,Protein.Ids,colData,intensity) %>%
      dplyr::group_by(colData, !!dplyr::sym(i)) %>%
      dplyr::summarise(intensity = sum(intensity, na.rm = TRUE) + 1,
                       .groups = "drop") %>%
      dplyr::mutate(log2 = log2(intensity)) %>%
      dplyr::select(-intensity) %>%
      tidyr::pivot_wider(names_from = colData, values_from = log2)

    readr::write_tsv(df, paste0(output_dir,"/",i,"_matrix.tsv"))
  }
}
