#' Download fasta files
#'
#' Download the \code{_cds_from_genomic.fna.gz} files for a list of bacteria
#'
#' This script will download all the \code{_cds_from_genomic.fna.gz} based on a
#' list of bacteria names.
#'
#' The results will be saved in the `code{dir/fna/} directory.
#'
#' To avoid redownloading the same files over and over, this function will
#' check by default if the file is already present with the correct md5sum
#' before downloading.
#'
#' @param n The names of the bacteria to download, correcponds to the
#' \code{organism_name} column from the \code{metadata} parameter.
#' @param metadata The metadata file downloaded with the
#' \code{fetch_bacteria_metadata} function. If \code{NULL}, the pre-computed
#' metadatas will be used. Default: \code{NULL}
#' @param dir The output directory. Default: \code{.}.
#' @param strict If \code{TRUE}, the names must be identical to the
#' \code{organism_name} column of the metadata object (by default, it's the
#' \code{bacteria} object). Otherwise, partial match will be allowed and the
#' matchning will not be case sensitive. Default: \code{TRUE}.
#' @param force If \code{TRUE}, the file will be downloaded even if it is
#' already present in the \code{dir/fna} directory. Default: \code{FALSE}.
#' @param verbose Print details of the downloading process. Default:
#' \code{FALSE}.
#' @param cores Number of cores to use for the download. Default: 1.
#'
#' @return Invisibly return the subset of metadata that was selected to
#' download.
#'
#' @examples
#' \dontrun{
#'   download_files("Campylobacter_avium")
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom tools md5sum
#' @importFrom readr read_delim
#' @importFrom stringr str_detect
#' @importFrom stringr str_extract
#' @importFrom stringr regex
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom parallel mclapply
#'
#' @export
download_files <- function(n, metadata = NULL, dir = ".", merge = NULL,
                           strict = TRUE, force = FALSE, verbose = FALSE,
                           cores = 1) {
    stopifnot(is(n, "character"))
    stopifnot(length(n) >= 1)
    stopifnot(dir.exists(dir))
    stopifnot(is(strict, "logical"))
    stopifnot(is(force, "logical"))
    stopifnot(is(verbose, "logical"))

    current_metadata <- filter_metadata(metadata, n, strict)
    stopifnot(nrow(current_metadata) >= 1)

    for (i in 1:nrow(current_metadata)) {
        download_file(current_metadata[i,], dir, force, verbose)
    }
    tmp <- split(current_metadata, 1:nrow(current_metadata)) %>%
        parallel::mclapply(function(x) download_file(x, dir, force, verbose),
                 mc.cores = cores)

    invisible(metadata)
}

download_file <- function(metadata, dir, force, verbose) {
    stopifnot(is(metadata, "data.frame"))
    stopifnot(nrow(metadata) == 1)

    cds_file <- paste0(metadata$assembly_accession, "_",
                       metadata$asm_name,
                       "_genomic.fna.gz")
    current_filename <- paste0(dir, "/", cds_file)
    current_url <- paste0(metadata$ftp_path, "/", cds_file)

    print_verbose(paste0("\nCurrent file: ", current_filename), verbose)

    md5_url <- paste0(metadata$ftp_path, "/md5checksums.txt")
    expected_md5 <- readr::read_delim(md5_url,
                                      col_names = c("md5", "file"),
                                      delim = "  ",
                                      col_types = "cc",
                                      progress = FALSE) %>%
        dplyr::filter(file == paste0(" ./", cds_file)) %>%
        dplyr::pull(md5)
    stopifnot(length(expected_md5) == 1)

    invalid_md5 <- TRUE
    if (!force) {
        if (file.exists(current_filename)) {
            print_verbose(paste0("    ", current_filename, " exists"), verbose)

            observed_md5 <- tools::md5sum(current_filename)

            invalid_md5 <- observed_md5 != expected_md5

            if (invalid_md5) {
                print_verbose("    md5sum is invalid", verbose)
            } else {
                print_verbose("    md5sum is valid, won't download.", verbose)
            }
        } else {
            to_print <- paste0("    ", current_filename, " does not exists")
            print_verbose(to_print, verbose)
        }
    } else {
        print_verbose("    force == TRUE", verbose)
    }

    if (!file.exists(current_filename) | force | invalid_md5) {
        print_verbose(paste0("    Downloading: ", current_filename), verbose)
        download.file(current_url,
                      quiet = TRUE,
                      destfile = current_filename,
                      method = "curl",
                      extra = "-L")

        observed_md5 <- tools::md5sum(current_filename)
        if (observed_md5 != expected_md5) {
            stop("Invalid md5")
        }
    }
    invisible(current_filename)
}

print_verbose <- function(msg, verbose) {
    if (verbose) message(msg)
}

filter_metadata <- function(metadata, n, strict) {
    metadata <- validate_metadata(metadata)
    if (strict) {
        stopifnot(all(n %in% metadata$organism_name))
        current_metadata <- dplyr::filter(metadata, organism_name %in% n)
    } else {
        organism_name <- tolower(metadata$organism_name) %>%
            stringr::str_extract("^[^ ]* [^ ]*")
        reg <- paste(tolower(n), collapse = "|")
        i <- stringr::str_detect(organism_name, reg)
        i[is.na(i)] <- FALSE
        current_metadata <- metadata[i,]
    }
    current_metadata
}
