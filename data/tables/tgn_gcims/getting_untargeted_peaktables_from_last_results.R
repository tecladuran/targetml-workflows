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

results_folders <- folders[
  grepl("results_\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}", folders)
]

# Sort by timestamp (descending) and select latest
latest_folder <- sort(results_folders, decreasing = TRUE)[1]
cat("Results:", latest_folder, "\n")

# Read peak table
peak_table_corrected <- read.csv(
  file.path(
    latest_folder,
    "all_results",
    "peak_table_corrected_tgn.csv"
  )
)

# QC pool data
qcs <- peak_table_corrected %>%
  filter(matrix == "pool") %>%
  dplyr::select(-matrix, -patient_id, -patient_condition)

# Patient data
patient_data <- peak_table_corrected %>%
  filter(matrix == "urine") %>%
  dplyr::select(-matrix, -control_level, -day)

# Output directory
output_dir <- file.path(
  repo_root,
  "data",
  "tables",
  "tgn_gcims"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Write outputs
write_csv(
  peak_table_corrected,
  file.path(output_dir, "tgn_peak_table.csv")
)

write_csv(
  qcs,
  file.path(output_dir, "tgn_qcs_pool.csv")
)

write_csv(
  patient_data,
  file.path(output_dir, "tgn_peak_table_patient.csv")
)
