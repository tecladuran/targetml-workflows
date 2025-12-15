# To run after getting_peaktables_from_last_results.R and cluster_filtering_summary

peak_table<-read_csv("/storage/users/tduran/Projects/targetml-gcims-tools/data/peak_tables_targetml/tgn_peak_table.csv")
patients<-read_csv("/storage/users/tduran/Projects/targetml-gcims-tools/data/peak_tables_targetml/tgn_peak_table_patient.csv")

excluded_clusters <- read_lines("/storage/users/tduran/Projects/targetml-gcims-tools/data/peak_tables_targetml/excluded_clusters.txt")

peak_table_filtered <- peak_table %>%
  select(-any_of(excluded_clusters))

patients_filtered <- patients %>%
  select(-any_of(excluded_clusters))

out_path_1 <- "/storage/users/tduran/Projects/targetml-gcims-tools/data/peak_tables_targetml/tgn_peak_table_patient_filtered.csv"
out_path_2 <- "/storage/users/tduran/Projects/targetml-gcims-tools/data/peak_tables_targetml/tgn_peak_table_filtered.csv"

write_csv(patients_filtered, out_path_1)
write_csv(data_filtered, out_path_2)

