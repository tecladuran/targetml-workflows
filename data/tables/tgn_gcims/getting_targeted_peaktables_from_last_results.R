library(dplyr)
library(readr)
library(rprojroot)

# Find repository root (git)
repo_root <- rprojroot::find_root(rprojroot::is_git_root)

# Define root folder for results
root_folder <- file.path(
  repo_root,
  "data",
  "raw",
  "tgn_results_gcims",
  "urine"
)

# List folders starting with "results_"
folders <- list.dirs(
  root_folder,
  full.names = TRUE,
  recursive = FALSE
)

valid_folders <- folders[
  grepl("results_\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}", basename(folders)) &
    file.exists(file.path(folders, "relevant_results"))
]

if (length(valid_folders) == 0) {
  stop("No valid results folder with relevant_results found.")
}

latest_folder <- sort(valid_folders, decreasing = TRUE)[1]

cat("Using results from:", latest_folder, "\n")

# Read peak table

peak_table_relevant <- read_csv(
  file.path(
    latest_folder,
    "relevant_results",
    "peak_table_tgn_relevant.csv"
  )
)

peak_table_relevant <- peak_table_relevant %>%
  dplyr::select(
    SampleID,
    class,
    matrix,
    patient_id,
    patient_condition,
    control_level,
    `Anisole (Monomer)`,
    `2-Heptanone (Dimer)`
  ) %>%
  dplyr::rename(
    !!normalize_name("Anisole (Monomer)")   := `Anisole (Monomer)`,
    !!normalize_name("2-Heptanone (Dimer)") := `2-Heptanone (Dimer)`
  )

peak_table_patients <- peak_table_relevant %>%
  dplyr::filter(matrix == "urine") %>% 
  dplyr::select(-class, -matrix, -control_level)
  


output_dir <- file.path(repo_root, "data", "tables", "tgn_gcims")

write_csv(
  peak_table_relevant,
  file.path(output_dir, "tgn_gcims_targeted_all.csv")
)

write_csv(
  peak_table_patients,
  file.path(output_dir, "tgn_gcims_targeted_patients.csv")
)

