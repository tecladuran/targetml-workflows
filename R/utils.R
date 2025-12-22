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

get_x_cols <- function(df) {
  x_cols <- df %>%
    dplyr::select(where(is.numeric)) %>%
    dplyr::select(-any_of(c("patient_id"))) %>%
    colnames()
  
  x_cols
}

# Function to check normality (Shapiro)
check_normal <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(NA)
  shapiro.test(x)$p.value
} 
