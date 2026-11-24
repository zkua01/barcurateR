library(tidyverse)
library(seqinr)
library(stringr)

local_raw <- read.fasta("data-raw/slovenian_fish/raw/ribe_sava_mlinarec_20261103.fasta", as.string = TRUE)

local_csv <- data.frame(species_name = sapply(local_raw, attr, "Annot"), sequence = sapply(local_raw, paste, collapse = ""), stringsAsFactors = FALSE) |>
  mutate(species_name = str_trim(str_replace_all(species_name, ">", "")))

write.csv(local_csv, "data-raw/slovenian_fish/raw/ribe_sava_mlinarec_20261103.csv", row.names = FALSE)
