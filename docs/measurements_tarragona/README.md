# Tarragona Measurements

This directory contains the documentation and workflows related to the Tarragona GC-IMS measurement campaign.  
It covers the full pipeline from raw data preprocessing to peak table preparation, which serves as the input for all subsequent analyses.

---

## [preprocessing](https://github.com/tecladuran/targetml-workflows/tree/main/docs/measurements_tarragona/preprocessing)

Preprocessing workflows starting from raw GC-IMS measurements.  
These documents describe the generation and preprocessing of datasets used in the Tarragona campaign.

The preprocessing scripts rely on raw data files that are not included in this repository.  
As a result, they can be executed as-is by group members with access to the cluster environment, without requiring any path modification.  
If a different dataset is to be used, the corresponding data paths must be adapted accordingly.

---

## [identification](https://github.com/tecladuran/targetml-workflows/tree/main/docs/measurements_tarragona/identification)

[Workflow for analyte identification](https://github.com/tecladuran/targetml-workflows/blob/main/docs/measurements_tarragona/identification/identification_tgn_urine.md) based on drift time and retention time.  
Analytes can be easily added or modified, and the workflow automatically identifies the corresponding clusters and generates a table containing the intensities with the names of the analytes as column names.

---

## [cluster_filtering_summary.md](https://github.com/tecladuran/targetml-workflows/blob/main/docs/measurements_tarragona/cluster_filtering_summary.md)

Summary of the cluster filtering strategy applied to the peak table.  
This document describes how clusters are filtered based on their representation across samples and their stability in quality control (QC) measurements.  
Underrepresented clusters are excluded using a minimum representation threshold, and additional filtering is applied based on QC variability using the relative standard deviation (RSD).  
The outcome of this process is a filtered peak table, together with a record of the excluded clusters, which is used in subsequent analyses.

---

## [data_transformation](https://github.com/tecladuran/targetml-workflows/tree/main/docs/measurements_tarragona/data_transformation)

Evaluation of different data transformations and normalization strategies.  
This section analyzes how individual transformations and their combinations affect the data, using multiple figures of merit their impact.

---

## Notes

The peak tables generated in this directory are used as input for all downstream targeted and untargeted analyses.
