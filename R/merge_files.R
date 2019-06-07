#' Merge multiple assemblies
#'
#' @param n The names of the bacteria to download, correcponds to the
#' \code{organism_name} column from the \code{metadata} parameter.
#' @param output The filename to save the merged fasta.
#' @param dir The download directory. Default: \code{.}.
#' @param metadata The metadata file downloaded with the
#' \code{fetch_bacteria_metadata} function. If \code{NULL}, the pre-computed
#' metadatas will be used. Default: \code{NULL}
#' @param strict If \code{TRUE}, the names must be identical to the
#' \code{organism_name} column of the metadata object (by default, it's the
#' \code{bacteria} object). Otherwise, partial match will be allowed and the
#' matchning will not be case sensitive. Default: \code{TRUE}.
#' @param force If \code{TRUE}, remove file if it already exists. Default:
#' \code{FALSE}.
#'
#' @return Invisibly return the subset of metadata that was selected to
#' download.
#'
#' @examples
#' \dontrun{
#'   merge("Campylobacter avium", output = "c_avium.fna")
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom purrr walk
#' @importFrom stringr str_extract
#' @importFrom stringr str_detect
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings writeXStringSet
#'
#' @export

merge_files <- function(n, output, dir = ".", metadata = NULL, strict = FALSE,
                  force = FALSE) {
    stopifnot(is(output, "character"))
    stopifnot(dir.exists(dirname(output)))

    current_metadata <- filter_metadata(metadata, n, strict)

    fna_files <- paste0(dir, "/",
                        stringr::str_extract(current_metadata$ftp_path, "[^\\/]*$"),
                        "_genomic.fna.gz")
    stopifnot(all(purrr::map_lgl(fna_files, file.exists)))

    if (file.exists(output)) {
        if (!force) {
            msg <- paste0("output file (", output, ") exists.\n")
            msg <- paste0(msg, "Please change output value, remove file")
            msg <- paste0(msg, " or use force = TRUE.")
            stop(msg)
        } else {
            file.remove(output)
        }
    }

    # Right now it's very simple, but this implementation will make it easier
    # to modify the files when merging in the future
    merge_file <- function(x) {
        readDNAStringSet(x) %>%
            writeXStringSet(output, append = TRUE)
    }
    purrr::walk(fna_files, merge_file)

    invisible(current_metadata)
}
