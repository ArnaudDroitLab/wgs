validate_metadata <- function(metadata) {
    stopifnot(is.null(metadata) | is(metadata, "data.frame"))
    if (is.null(metadata)) {
        metadata <- bacteria
    }
    stopifnot(nrow(metadata) >= 1)
    stopifnot("ftp_path" %in% colnames(metadata))
    stopifnot("assembly_accession" %in% colnames(metadata))
    stopifnot("asm_name" %in% colnames(metadata))
    invisible(TRUE)
}

validate_blast <- function(blast_res) {
    expected_cols <- c("qseqid", "sseqid", "evalue", "qlen", "pident", "qcovs",
                       "length", "sstart", "send")
    stopifnot(nrow(blast_res) >= 1)
    stopifnot(all(expected_cols %in% colnames(blast_res)))
    invisible(TRUE)

}
