#' Extract fasta sequences
#'
#' Extract the sequences obtained after a blast analysis.
#'
#' The blast results must contain an header with the following columns:
#'
#' \enumerate{
#'   \item qseqid
#'   \item sseqid
#'   \item evalue
#'   \item qlen
#'   \item qlen
#'   \item pident
#'   \item qcovs
#'   \item length
#'   \item sstart
#'   \item send
#' }
#'
#' @param blast_out The filename of the blast results
#' @param metadata The subset of the metadata obtained with the
#' \code{merge_files} or the \code{download_files} functions.
#' @param merged_fasta The fasta obtained with the \code{merge_files} function
#'
#' @return TODO
#'
#' @examples
#' \dontrun{
#'   metadata <- get_demo_metadata()
#'   extract_fasta("test.out", metadata, "merged.fna")
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
#' @importFrom dplyr left_join
#' @importFrom dplyr if_else
#' @importFrom dplyr mutate
#' @importFrom tibble tibble
#' @importFrom GenomicRanges GRanges 
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom Biostrings readDNAStringSet
#'
#' @export
extract_fasta <- function(blast_out, metadata, merged_fasta) {
    stopifnot(file.exists(blast_out))
    blast_res <- readr::read_tsv(blast_out, col_types = "ccddddddd")
    validate_blast(blast_res)

    stopifnot(!is.null(metadata))
    validate_metadata(metadata)

    fna <- Biostrings::readDNAStringSet(merged_fasta)
    fna_desc <- tibble::tibble(id = stringr::str_extract(names(fna), "^[^ ]*"),
                               full_name = names(fna))
    blast <- dplyr::left_join(blast_res, fna_desc, by = c("sseqid" = "id")) %>%
        dplyr::mutate(tmp = sstart,
                      sstart = dplyr::if_else(sstart < send, sstart, send),
                      send = dplyr::if_else(tmp < send, send, tmp))

    gr <- GenomicRanges::GRanges(blast$full_name,
                                 IRanges::IRanges(blast$sstart, blast$send))
    subset_fna <- fna[gr]
    names(subset_fna) <- paste0(names(subset_fna), " ", BiocGenerics::start(gr),
                                "-", BiocGenerics::end(gr))
    subset_fna
}

#' Annotate blast results
#'
#' This function will add metadata to the blast results.
#'
#' The blast results must contain an header with the following columns:
#'
#' \enumerate{
#'   \item qseqid
#'   \item sseqid
#'   \item evalue
#'   \item qlen
#'   \item qlen
#'   \item pident
#'   \item qcovs
#'   \item length
#'   \item sstart
#'   \item send
#' }
#'
#' @param blast_out The filename of the blast results
#' @param metadata The subset of the metadata obtained with the
#' \code{merge_files} or the \code{download_files} functions.
#' @param dir The directory where the fasta files are located.
#'
#' @return A \code{data.frame} with the annotated results
#'
#' @examples
#' \dontrun{
#'   metadata <- get_demo_metadata()
#'   annotate_blast("test.out", metadata, dir = ".")
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
#' @importFrom purrr map_dfr
#' @importFrom dplyr left_join
#' @importFrom tibble tibble
#' @importFrom Biostrings readDNAStringSet
#'
#' @export
annotate_blast <- function(blast_out, metadata, dir = ".") {
    stopifnot(file.exists(blast_out))
    blast_res <- readr::read_tsv(blast_out, col_types = "ccddddddd")
    validate_blast(blast_res)

    stopifnot(!is.null(metadata))
    validate_metadata(metadata)

    anno <- split(metadata, 1:nrow(metadata)) %>%
        purrr::map_dfr(extract_anno, dir)
    dplyr::left_join(blast_res, anno, by = c("sseqid" = "id"))
}

extract_anno <- function(metadata, dir) {
    stopifnot(nrow(metadata) == 1)

    fna <- paste0(dir, "/", basename(metadata$ftp_path),
                  "_genomic.fna.gz")
                 
    n <- names(Biostrings::readDNAStringSet(fna))

    tibble::tibble(id = stringr::str_extract(n ,"^[^ ]*"),
                   assembly_accession = assembly_accession,
                   asm_name = asm_name,
                   description = stringr::str_replace(n , "^[^ ]* ", ""))
}
