require(tidyverse)
require(rvest)
require(usethis)

base_url <- "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria"
bacteria_html <- read_html(base_url)

bacteria <- bacteria_html %>%
    html_nodes("a") %>%
    html_text %>%
    str_replace("\\/", "")

bacteria[!str_detect(bacteria, "\\.txt$")][-1] %>%
    usethis::use_data
