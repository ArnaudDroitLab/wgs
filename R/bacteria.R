#' Genbank's bacterias metadata
#'
#' Fetch the metadata for all the bacterias in Genbank
#'
#' Data is obtained from the \code{assembly_summary.txt} file:
#' ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
#'
#' @return A \code{data.frame} with the following columns:
#' \itemize{
#'   \item assembly_accession
#'   \item asm_name
#'   \item refseq_category
#'   \item organism_name
#'   \item infraspecific_name
#'   \item assembly_level
#'   \item version_status
#'   \item ftp_path
#' }
#'
#' @examples
#' \dontrun{
#'   bacteria <- fetch_bacteria_metadata()
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
#' @importFrom dplyr select
#'
#' @export
fetch_bacteria_metadata <- function() {
  col_types <- rep("c", 22) %>% paste(collapse = "")

  base_url <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria"
  url <- paste0(base_url, "/assembly_summary.txt")

  read_tsv(url, skip = 1, col_types = col_types, quote = "") %>%
      dplyr::select(assembly_accession = `# assembly_accession`,
                    asm_name,
                    refseq_category,
                    organism_name,
                    infraspecific_name,
                    assembly_level,
                    version_status,
                    ftp_path)
}
