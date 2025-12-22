Cluster Filtering in Peak Table
================
Tecla Duran Fort
2025-12-22

- <a href="#0-load-data" id="toc-0-load-data">0. Load Data</a>
- <a href="#1-filter-by-representation-across-samples"
  id="toc-1-filter-by-representation-across-samples">1. Filter By
  Representation Across Samples</a>
  - <a href="#defining-a-threshold-20"
    id="toc-defining-a-threshold-20">Defining a threshold (20%)</a>
  - <a href="#exclude-them-drom-the-datasets"
    id="toc-exclude-them-drom-the-datasets">Exclude them drom the
    datasets</a>
- <a href="#save" id="toc-save">Save</a>
- <a href="#2-filter-by-qc-stability" id="toc-2-filter-by-qc-stability">2.
  Filter By QC Stability</a>

# 0. Load Data

``` r
data <- read.csv(file.path(latest_folder, "all_results/peak_table_corrected_tgn_old_new.csv"))

peak_table_duplicity <- read_csv(file.path(latest_folder,"peak_table_duplicity.csv"))[ , -1]

patient_data<-data%>%filter(matrix=="urine")%>%dplyr::select(-matrix,-control_level,-day)
```

    ## Data from (Latest Results Folder): data/raw/campaign_22//2025-12-22_14-39

    ## data (peak_table_corrected_tgn_old_new.csv)
    ##  Rows: 62  | Cols: 517 
    ##  Description: Corrected peak table for all TGN samples.

    ## peak_table_duplicity (peak_table_duplicity.csv)
    ##  Rows: 508  | Cols: 63 
    ##  Description: Duplicated peak information.

    ## patient_data
    ##  Rows: 48  | Cols: 514 
    ##  Description: Only samples corresponding to patients.

# 1. Filter By Representation Across Samples

    ## Data Frame of Representation

    ## # A tibble: 6 × 4
    ##   cluster    n_present total_samples representation_ratio
    ##   <chr>          <dbl>         <int>                <dbl>
    ## 1 Cluster1           3            62               0.0484
    ## 2 Cluster10          2            62               0.0323
    ## 3 Cluster100        15            62               0.242 
    ## 4 Cluster101         9            62               0.145 
    ## 5 Cluster102         2            62               0.0323
    ## 6 Cluster103        61            62               0.984

## Defining a threshold (20%)

``` r
rare_clusters <- representation_df %>%
  filter(representation_ratio < 0.2) %>%
  pull(cluster)
```

    ## Number of rare clusters (Under 20% representation): 293

## Exclude them drom the datasets

``` r
data_filtered <- data %>%
  dplyr::select(-any_of(rare_clusters))

peak_table_duplicity_filtered <- peak_table_duplicity %>%
  dplyr::select(-any_of(rare_clusters))

patient_data_filtered <- patient_data %>%
  dplyr::select(-any_of(rare_clusters))
```

# Save

# 2. Filter By QC Stability

(no esta fet perquè no hi ha QCs en pool a la campanya 22)
