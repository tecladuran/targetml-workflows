# Analyte names

normalize_name <- function(x) {
  x <- tolower(x)                 # lowercase
  x <- gsub("[0-9]+", "", x)      # remove numbers
  x <- gsub("\\(.*?\\)", "", x)   # remove (...) content
  x <- gsub("\\[.*?\\]", "", x)   # remove [...] content
  x <- gsub("[^a-z]+", "_", x)    # replace non-letters by _
  x <- gsub("^_|_$", "", x)       # trim leading/trailing _
  x <- gsub("_+", "_", x)         # collapse multiple _
  x
}



analyte_names <- c(
  anisole   = "Anisole",
  heptanone = "2-Heptanone",
  pentanone = "2-Pentanone"
)

display_name <- function(x) {
  if (x %in% names(analyte_names)) {
    analyte_names[[x]]
  } else {
    x
  }
}