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
