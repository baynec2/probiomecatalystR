#' get_go_kegg_eggnogg_annotation
#'
#' given the pg_matrix output of diann, assign go terms to each protein group
#' and summarise.
#'
#' @param pg_matrix_fp pg_matrix output of diann
#' @param outdir out directory
#'
#' @returns tsv at output directory location
#' @export
#'
#' @examples
#' insulin = get_annotations_from_uniprot("P01308")
#' t = get_annotations_from_uniprot("P47340")
get_annotations_from_uniprot = function(uniprot_ids){

  columns = "organism_id,accession,go,xref_kegg,xref_eggnog,cc_subcellular_location"

  # We can get both GO and Kegg data from the uniprot API
  go_kegg_data = annotate_uniprot_ids(uniprot_ids,
                                 columns = columns,
                                 batch_size = 100) |>
    dplyr::mutate(
      go = dplyr::na_if(go, " "),
      go = dplyr::na_if(go, ""),
      xref_kegg = dplyr::na_if(xref_kegg, " "),
      xref_kegg = dplyr::na_if(xref_kegg, ""),
      xref_eggnog = dplyr::na_if(xref_eggnog, " "),
      xref_eggnog = dplyr::na_if(xref_eggnog, "")
    )

  return(go_kegg_data)
}


