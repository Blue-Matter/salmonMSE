
glossary <- read.csv("vignettes/tables/glossary.csv")[, -1]
usethis::use_data(glossary, overwrite = TRUE)
